/**************************************************************************
 * 
 * Copyright 2003 Tungsten Graphics, Inc., Cedar Park, Texas.
 * Copyright 2009 Intel Corporation.
 * All Rights Reserved.
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sub license, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 * 
 * The above copyright notice and this permission notice (including the
 * next paragraph) shall be included in all copies or substantial portions
 * of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NON-INFRINGEMENT.
 * IN NO EVENT SHALL TUNGSTEN GRAPHICS AND/OR ITS SUPPLIERS BE LIABLE FOR
 * ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 * 
 **************************************************************************/

#include "main/glheader.h"
#include "main/enums.h"
#include "main/image.h"
#include "main/mtypes.h"
#include "main/arrayobj.h"
#include "main/attrib.h"
#include "main/blend.h"
#include "main/bufferobj.h"
#include "main/buffers.h"
#include "main/depth.h"
#include "main/enable.h"
#include "main/macros.h"
#include "main/matrix.h"
#include "main/polygon.h"
#include "main/texstate.h"
#include "main/shaders.h"
#include "main/stencil.h"
#include "main/varray.h"
#include "glapi/dispatch.h"
#include "swrast/swrast.h"

#include "intel_context.h"
#include "intel_blit.h"
#include "intel_chipset.h"
#include "intel_clear.h"
#include "intel_fbo.h"
#include "intel_pixel.h"

#define FILE_DEBUG_FLAG DEBUG_BLIT

#define TRI_CLEAR_COLOR_BITS (BUFFER_BIT_BACK_LEFT |			\
			      BUFFER_BIT_FRONT_LEFT |			\
			      BUFFER_BIT_COLOR0 |			\
			      BUFFER_BIT_COLOR1 |			\
			      BUFFER_BIT_COLOR2 |			\
			      BUFFER_BIT_COLOR3 |			\
			      BUFFER_BIT_COLOR4 |			\
			      BUFFER_BIT_COLOR5 |			\
			      BUFFER_BIT_COLOR6 |			\
			      BUFFER_BIT_COLOR7)


/**
 * Per-context one-time init of things for intl_clear_tris().
 * Basically set up a private array object for vertex/color arrays.
 */
static void
init_clear(GLcontext *ctx)
{
   struct intel_context *intel = intel_context(ctx);
   struct gl_array_object *arraySave = NULL;
   const GLuint arrayBuffer = ctx->Array.ArrayBufferObj->Name;
   const GLuint elementBuffer = ctx->Array.ElementArrayBufferObj->Name;

   /* create new array object */
   intel->clear.arrayObj = _mesa_new_array_object(ctx, ~0);

   /* save current array object, bind new one */
   _mesa_reference_array_object(ctx, &arraySave, ctx->Array.ArrayObj);
   _mesa_reference_array_object(ctx, &ctx->Array.ArrayObj, intel->clear.arrayObj);

   /* one-time setup of vertex arrays (pos, color) */
   _mesa_BindBufferARB(GL_ARRAY_BUFFER_ARB, 0);
   _mesa_BindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB, 0);
   _mesa_ColorPointer(4, GL_FLOAT, 4 * sizeof(GLfloat), intel->clear.color);
   _mesa_VertexPointer(3, GL_FLOAT, 3 * sizeof(GLfloat), intel->clear.vertices);
   _mesa_Enable(GL_COLOR_ARRAY);
   _mesa_Enable(GL_VERTEX_ARRAY);

   /* restore original array object */
   _mesa_reference_array_object(ctx, &ctx->Array.ArrayObj, arraySave);
   _mesa_reference_array_object(ctx, &arraySave, NULL);

   /* restore original buffer objects */
   _mesa_BindBufferARB(GL_ARRAY_BUFFER_ARB, arrayBuffer);
   _mesa_BindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB, elementBuffer);
}



/**
 * Perform glClear where mask contains only color, depth, and/or stencil.
 *
 * The implementation is based on calling into Mesa to set GL state and
 * performing normal triangle rendering.  The intent of this path is to
 * have as generic a path as possible, so that any driver could make use of
 * it.
 */
void
intel_clear_tris(GLcontext *ctx, GLbitfield mask)
{
   struct intel_context *intel = intel_context(ctx);
   GLfloat dst_z;
   struct gl_framebuffer *fb = ctx->DrawBuffer;
   int i;
   GLboolean saved_fp_enable = GL_FALSE, saved_vp_enable = GL_FALSE;
   GLuint saved_shader_program = 0;
   unsigned int saved_active_texture;
   struct gl_array_object *arraySave = NULL;

   if (!intel->clear.arrayObj)
      init_clear(ctx);

   assert((mask & ~(TRI_CLEAR_COLOR_BITS | BUFFER_BIT_DEPTH |
		    BUFFER_BIT_STENCIL)) == 0);

   _mesa_PushAttrib(GL_COLOR_BUFFER_BIT |
		    GL_CURRENT_BIT |
		    GL_DEPTH_BUFFER_BIT |
		    GL_ENABLE_BIT |
		    GL_POLYGON_BIT |
		    GL_STENCIL_BUFFER_BIT |
		    GL_TRANSFORM_BIT |
		    GL_CURRENT_BIT);
   saved_active_texture = ctx->Texture.CurrentUnit;

   /* Disable existing GL state we don't want to apply to a clear. */
   _mesa_Disable(GL_ALPHA_TEST);
   _mesa_Disable(GL_BLEND);
   _mesa_Disable(GL_CULL_FACE);
   _mesa_Disable(GL_FOG);
   _mesa_Disable(GL_POLYGON_SMOOTH);
   _mesa_Disable(GL_POLYGON_STIPPLE);
   _mesa_Disable(GL_POLYGON_OFFSET_FILL);
   _mesa_Disable(GL_LIGHTING);
   _mesa_Disable(GL_CLIP_PLANE0);
   _mesa_Disable(GL_CLIP_PLANE1);
   _mesa_Disable(GL_CLIP_PLANE2);
   _mesa_Disable(GL_CLIP_PLANE3);
   _mesa_Disable(GL_CLIP_PLANE4);
   _mesa_Disable(GL_CLIP_PLANE5);
   _mesa_PolygonMode(GL_FRONT_AND_BACK, GL_FILL);
   if (ctx->Extensions.ARB_fragment_program && ctx->FragmentProgram.Enabled) {
      saved_fp_enable = GL_TRUE;
      _mesa_Disable(GL_FRAGMENT_PROGRAM_ARB);
   }
   if (ctx->Extensions.ARB_vertex_program && ctx->VertexProgram.Enabled) {
      saved_vp_enable = GL_TRUE;
      _mesa_Disable(GL_VERTEX_PROGRAM_ARB);
   }
   if (ctx->Extensions.ARB_shader_objects && ctx->Shader.CurrentProgram) {
      saved_shader_program = ctx->Shader.CurrentProgram->Name;
      _mesa_UseProgramObjectARB(0);
   }

   if (ctx->Texture._EnabledUnits != 0) {
      int i;

      for (i = 0; i < ctx->Const.MaxTextureUnits; i++) {
	 _mesa_ActiveTextureARB(GL_TEXTURE0 + i);
	 _mesa_Disable(GL_TEXTURE_1D);
	 _mesa_Disable(GL_TEXTURE_2D);
	 _mesa_Disable(GL_TEXTURE_3D);
	 if (ctx->Extensions.ARB_texture_cube_map)
	    _mesa_Disable(GL_TEXTURE_CUBE_MAP_ARB);
	 if (ctx->Extensions.NV_texture_rectangle)
	    _mesa_Disable(GL_TEXTURE_RECTANGLE_NV);
	 if (ctx->Extensions.MESA_texture_array) {
	    _mesa_Disable(GL_TEXTURE_1D_ARRAY_EXT);
	    _mesa_Disable(GL_TEXTURE_2D_ARRAY_EXT);
	 }
      }
   }

   /* save current array object, bind our private one */
   _mesa_reference_array_object(ctx, &arraySave, ctx->Array.ArrayObj);
   _mesa_reference_array_object(ctx, &ctx->Array.ArrayObj, intel->clear.arrayObj);

   intel_meta_set_passthrough_transform(intel);

   for (i = 0; i < 4; i++) {
      COPY_4FV(intel->clear.color[i], ctx->Color.ClearColor);
   }

   /* convert clear Z from [0,1] to NDC coord in [-1,1] */
   dst_z = -1.0 + 2.0 * ctx->Depth.Clear;

   /* Prepare the vertices, which are the same regardless of which buffer we're
    * drawing to.
    */
   intel->clear.vertices[0][0] = fb->_Xmin;
   intel->clear.vertices[0][1] = fb->_Ymin;
   intel->clear.vertices[0][2] = dst_z;
   intel->clear.vertices[1][0] = fb->_Xmax;
   intel->clear.vertices[1][1] = fb->_Ymin;
   intel->clear.vertices[1][2] = dst_z;
   intel->clear.vertices[2][0] = fb->_Xmax;
   intel->clear.vertices[2][1] = fb->_Ymax;
   intel->clear.vertices[2][2] = dst_z;
   intel->clear.vertices[3][0] = fb->_Xmin;
   intel->clear.vertices[3][1] = fb->_Ymax;
   intel->clear.vertices[3][2] = dst_z;

   while (mask != 0) {
      GLuint this_mask = 0;
      GLuint color_bit;

      color_bit = _mesa_ffs(mask & TRI_CLEAR_COLOR_BITS);
      if (color_bit != 0)
	 this_mask |= (1 << (color_bit - 1));

      /* Clear depth/stencil in the same pass as color. */
      this_mask |= (mask & (BUFFER_BIT_DEPTH | BUFFER_BIT_STENCIL));

      /* Select the current color buffer and use the color write mask if
       * we have one, otherwise don't write any color channels.
       */
      if (this_mask & BUFFER_BIT_FRONT_LEFT)
	 _mesa_DrawBuffer(GL_FRONT_LEFT);
      else if (this_mask & BUFFER_BIT_BACK_LEFT)
	 _mesa_DrawBuffer(GL_BACK_LEFT);
      else if (color_bit != 0)
	 _mesa_DrawBuffer(GL_COLOR_ATTACHMENT0 +
			  (color_bit - BUFFER_COLOR0 - 1));
      else
	 _mesa_ColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);

      /* Control writing of the depth clear value to depth. */
      if (this_mask & BUFFER_BIT_DEPTH) {
	 _mesa_DepthFunc(GL_ALWAYS);
	 _mesa_Enable(GL_DEPTH_TEST);
      } else {
	 _mesa_Disable(GL_DEPTH_TEST);
	 _mesa_DepthMask(GL_FALSE);
      }

      /* Control writing of the stencil clear value to stencil. */
      if (this_mask & BUFFER_BIT_STENCIL) {
	 _mesa_Enable(GL_STENCIL_TEST);
	 _mesa_StencilOpSeparate(GL_FRONT_AND_BACK,
				 GL_REPLACE, GL_REPLACE, GL_REPLACE);
	 _mesa_StencilFuncSeparate(GL_FRONT_AND_BACK, GL_ALWAYS,
				   ctx->Stencil.Clear,
				   ctx->Stencil.WriteMask[0]);
      } else {
	 _mesa_Disable(GL_STENCIL_TEST);
      }

      _mesa_DrawArrays(GL_TRIANGLE_FAN, 0, 4);

      mask &= ~this_mask;
   }

   intel_meta_restore_transform(intel);

   _mesa_ActiveTextureARB(GL_TEXTURE0 + saved_active_texture);
   if (saved_fp_enable)
      _mesa_Enable(GL_FRAGMENT_PROGRAM_ARB);
   if (saved_vp_enable)
      _mesa_Enable(GL_VERTEX_PROGRAM_ARB);

   if (saved_shader_program)
      _mesa_UseProgramObjectARB(saved_shader_program);

   _mesa_PopAttrib();

   /* restore current array object */
   _mesa_reference_array_object(ctx, &ctx->Array.ArrayObj, arraySave);
   _mesa_reference_array_object(ctx, &arraySave, NULL);
}

static const char *buffer_names[] = {
   [BUFFER_FRONT_LEFT] = "front",
   [BUFFER_BACK_LEFT] = "back",
   [BUFFER_FRONT_RIGHT] = "front right",
   [BUFFER_BACK_RIGHT] = "back right",
   [BUFFER_DEPTH] = "depth",
   [BUFFER_STENCIL] = "stencil",
   [BUFFER_ACCUM] = "accum",
   [BUFFER_AUX0] = "aux0",
   [BUFFER_COLOR0] = "color0",
   [BUFFER_COLOR1] = "color1",
   [BUFFER_COLOR2] = "color2",
   [BUFFER_COLOR3] = "color3",
   [BUFFER_COLOR4] = "color4",
   [BUFFER_COLOR5] = "color5",
   [BUFFER_COLOR6] = "color6",
   [BUFFER_COLOR7] = "color7",
};

/**
 * Called by ctx->Driver.Clear.
 */
static void
intelClear(GLcontext *ctx, GLbitfield mask)
{
   struct intel_context *intel = intel_context(ctx);
   const GLuint colorMask = *((GLuint *) & ctx->Color.ColorMask);
   GLbitfield tri_mask = 0;
   GLbitfield blit_mask = 0;
   GLbitfield swrast_mask = 0;
   struct gl_framebuffer *fb = ctx->DrawBuffer;
   GLuint i;

   if (0)
      fprintf(stderr, "%s\n", __FUNCTION__);

   /* HW color buffers (front, back, aux, generic FBO, etc) */
   if (colorMask == ~0) {
      /* clear all R,G,B,A */
      /* XXX FBO: need to check if colorbuffers are software RBOs! */
      blit_mask |= (mask & BUFFER_BITS_COLOR);
   }
   else {
      /* glColorMask in effect */
      tri_mask |= (mask & (BUFFER_BIT_FRONT_LEFT | BUFFER_BIT_BACK_LEFT));
   }

   /* HW stencil */
   if (mask & BUFFER_BIT_STENCIL) {
      const struct intel_region *stencilRegion
         = intel_get_rb_region(fb, BUFFER_STENCIL);
      if (stencilRegion) {
         /* have hw stencil */
         if (IS_965(intel->intelScreen->deviceID) ||
	     (ctx->Stencil.WriteMask[0] & 0xff) != 0xff) {
	    /* We have to use the 3D engine if we're clearing a partial mask
	     * of the stencil buffer, or if we're on a 965 which has a tiled
	     * depth/stencil buffer in a layout we can't blit to.
	     */
            tri_mask |= BUFFER_BIT_STENCIL;
         }
         else {
            /* clearing all stencil bits, use blitting */
            blit_mask |= BUFFER_BIT_STENCIL;
         }
      }
   }

   /* HW depth */
   if (mask & BUFFER_BIT_DEPTH) {
      /* clear depth with whatever method is used for stencil (see above) */
      if (IS_965(intel->intelScreen->deviceID) ||
	  tri_mask & BUFFER_BIT_STENCIL)
         tri_mask |= BUFFER_BIT_DEPTH;
      else
         blit_mask |= BUFFER_BIT_DEPTH;
   }

   /* If we're doing a tri pass for depth/stencil, include a likely color
    * buffer with it.
    */
   if (mask & (BUFFER_BIT_DEPTH | BUFFER_BIT_STENCIL)) {
      int color_bit = _mesa_ffs(mask & TRI_CLEAR_COLOR_BITS);
      if (color_bit != 0) {
	 tri_mask |= blit_mask & (1 << (color_bit - 1));
	 blit_mask &= ~(1 << (color_bit - 1));
      }
   }

   /* SW fallback clearing */
   swrast_mask = mask & ~tri_mask & ~blit_mask;

   for (i = 0; i < BUFFER_COUNT; i++) {
      GLuint bufBit = 1 << i;
      if ((blit_mask | tri_mask) & bufBit) {
         if (!fb->Attachment[i].Renderbuffer->ClassID) {
            blit_mask &= ~bufBit;
            tri_mask &= ~bufBit;
            swrast_mask |= bufBit;
         }
      }
   }

   if (blit_mask) {
      if (INTEL_DEBUG & DEBUG_BLIT) {
	 DBG("blit clear:");
	 for (i = 0; i < BUFFER_COUNT; i++) {
	    if (blit_mask & (1 << i))
	       DBG(" %s", buffer_names[i]);
	 }
	 DBG("\n");
      }
      intelClearWithBlit(ctx, blit_mask);
   }

   if (tri_mask) {
      if (INTEL_DEBUG & DEBUG_BLIT) {
	 DBG("tri clear:");
	 for (i = 0; i < BUFFER_COUNT; i++) {
	    if (tri_mask & (1 << i))
	       DBG(" %s", buffer_names[i]);
	 }
	 DBG("\n");
      }
      intel_clear_tris(ctx, tri_mask);
   }

   if (swrast_mask) {
      if (INTEL_DEBUG & DEBUG_BLIT) {
	 DBG("swrast clear:");
	 for (i = 0; i < BUFFER_COUNT; i++) {
	    if (swrast_mask & (1 << i))
	       DBG(" %s", buffer_names[i]);
	 }
	 DBG("\n");
      }
      _swrast_Clear(ctx, swrast_mask);
   }
}


void
intelInitClearFuncs(struct dd_function_table *functions)
{
   functions->Clear = intelClear;
}
