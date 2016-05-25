/**************************************************************************
 * 
 * Copyright 2003 Tungsten Graphics, Inc., Cedar Park, Texas.
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


#include "main/imports.h"
#include "main/buffers.h"
#include "main/context.h"
#include "main/framebuffer.h"
#include "main/matrix.h"
#include "main/renderbuffer.h"
#include "main/scissor.h"
#include "main/viewport.h"
#include "st_context.h"
#include "st_cb_fbo.h"
#include "st_public.h"
#include "pipe/p_defines.h"
#include "pipe/p_context.h"


struct st_framebuffer *
st_create_framebuffer( const __GLcontextModes *visual,
                       enum pipe_format colorFormat,
                       enum pipe_format depthFormat,
                       enum pipe_format stencilFormat,
                       uint width, uint height,
                       void *private)
{
   struct st_framebuffer *stfb = ST_CALLOC_STRUCT(st_framebuffer);
   if (stfb) {
      int samples = st_get_msaa();

      if (visual->sampleBuffers)
         samples = visual->samples;

      _mesa_initialize_framebuffer(&stfb->Base, visual);

      {
         /* fake frontbuffer */
         /* XXX allocation should only happen in the unusual case
            it's actually needed */
         struct gl_renderbuffer *rb
            = st_new_renderbuffer_fb(colorFormat, samples, FALSE);
         _mesa_add_renderbuffer(&stfb->Base, BUFFER_FRONT_LEFT, rb);
      }

      if (visual->doubleBufferMode) {
         struct gl_renderbuffer *rb
            = st_new_renderbuffer_fb(colorFormat, samples, FALSE);
         _mesa_add_renderbuffer(&stfb->Base, BUFFER_BACK_LEFT, rb);
      }

      if (depthFormat == stencilFormat && depthFormat != PIPE_FORMAT_NONE) {
         /* combined depth/stencil buffer */
         struct gl_renderbuffer *depthStencilRb
            = st_new_renderbuffer_fb(depthFormat, samples, FALSE);
         /* note: bind RB to two attachment points */
         _mesa_add_renderbuffer(&stfb->Base, BUFFER_DEPTH, depthStencilRb);
         _mesa_add_renderbuffer(&stfb->Base, BUFFER_STENCIL, depthStencilRb);
      }
      else {
         /* separate depth and/or stencil */

         if (visual->depthBits == 32) {
            /* 32-bit depth buffer */
            struct gl_renderbuffer *depthRb
               = st_new_renderbuffer_fb(depthFormat, samples, FALSE);
            _mesa_add_renderbuffer(&stfb->Base, BUFFER_DEPTH, depthRb);
         }
         else if (visual->depthBits == 24) {
            /* 24-bit depth buffer, ignore stencil bits */
            struct gl_renderbuffer *depthRb
               = st_new_renderbuffer_fb(depthFormat, samples, FALSE);
            _mesa_add_renderbuffer(&stfb->Base, BUFFER_DEPTH, depthRb);
         }
         else if (visual->depthBits > 0) {
            /* 16-bit depth buffer */
            struct gl_renderbuffer *depthRb
               = st_new_renderbuffer_fb(depthFormat, samples, FALSE);
            _mesa_add_renderbuffer(&stfb->Base, BUFFER_DEPTH, depthRb);
         }

         if (visual->stencilBits > 0) {
            /* 8-bit stencil */
            struct gl_renderbuffer *stencilRb
               = st_new_renderbuffer_fb(stencilFormat, samples, FALSE);
            _mesa_add_renderbuffer(&stfb->Base, BUFFER_STENCIL, stencilRb);
         }
      }

      if (visual->accumRedBits > 0) {
         /* 16-bit/channel accum */
         /* TODO: query the pipe screen for accumulation buffer format support */
         struct gl_renderbuffer *accumRb
            = st_new_renderbuffer_fb(PIPE_FORMAT_R16G16B16A16_SNORM, 0, TRUE);
         _mesa_add_renderbuffer(&stfb->Base, BUFFER_ACCUM, accumRb);
      }

      stfb->Base.Initialized = GL_TRUE;
      stfb->InitWidth = width;
      stfb->InitHeight = height;
      stfb->Private = private;
   }
   return stfb;
}


void st_resize_framebuffer( struct st_framebuffer *stfb,
                            uint width, uint height )
{
   if (stfb->Base.Width != width || stfb->Base.Height != height) {
      GET_CURRENT_CONTEXT(ctx);
      if (ctx) {
         _mesa_check_init_viewport(ctx, width, height);

         _mesa_resize_framebuffer(ctx, &stfb->Base, width, height);

         assert(stfb->Base.Width == width);
         assert(stfb->Base.Height == height);
      }
   }
}


void st_unreference_framebuffer( struct st_framebuffer *stfb )
{
   _mesa_reference_framebuffer((struct gl_framebuffer **) &stfb, NULL);
}



/**
 * Set/replace a framebuffer surface.
 * The user of the state tracker can use this instead of
 * st_resize_framebuffer() to provide new surfaces when a window is resized.
 */
void
st_set_framebuffer_surface(struct st_framebuffer *stfb,
                           uint surfIndex, struct pipe_surface *surf)
{
   GET_CURRENT_CONTEXT(ctx);
   static const GLuint invalid_size = 9999999;
   struct st_renderbuffer *strb;
   GLuint width, height, i;

   assert(surfIndex < BUFFER_COUNT);

   strb = st_renderbuffer(stfb->Base.Attachment[surfIndex].Renderbuffer);

   /* fail */
   if (!strb) return;

   /* replace the renderbuffer's surface/texture pointers */
   pipe_surface_reference( &strb->surface, surf );
   pipe_texture_reference( &strb->texture, surf->texture );

   if (ctx) {
      /* If ctx isn't set, we've likely not made current yet.
       * But when we do, we need to start setting this dirty bit
       * to ensure the renderbuffer attachements are up-to-date
       * via update_framebuffer.
       */
      ctx->st->dirty.st |= ST_NEW_FRAMEBUFFER;
   }

   /* update renderbuffer's width/height */
   strb->Base.Width = surf->width;
   strb->Base.Height = surf->height;

   /* Try to update the framebuffer's width/height from the renderbuffer
    * sizes.  Before we start drawing, all the rbs _should_ be the same size.
    */
   width = height = invalid_size;
   for (i = 0; i < BUFFER_COUNT; i++) {
      if (stfb->Base.Attachment[i].Renderbuffer) {
         if (width == invalid_size) {
            width = stfb->Base.Attachment[i].Renderbuffer->Width;
            height = stfb->Base.Attachment[i].Renderbuffer->Height;
         }
         else if (width != stfb->Base.Attachment[i].Renderbuffer->Width ||
                  height != stfb->Base.Attachment[i].Renderbuffer->Height) {
            /* inconsistant renderbuffer sizes, bail out */
            return;
         }
      }
   }

   if (width != invalid_size) {
      /* OK, the renderbuffers are of a consistant size, so update the
       * parent framebuffer's size.
       */
      stfb->Base.Width = width;
      stfb->Base.Height = height;
   }
}



/**
 * Return the pipe_surface for the given renderbuffer.
 */
int
st_get_framebuffer_surface(struct st_framebuffer *stfb, uint surfIndex, struct pipe_surface **surface)
{
   struct st_renderbuffer *strb;

   assert(surfIndex <= ST_SURFACE_DEPTH);

   /* sanity checks, ST tokens should match Mesa tokens */
   assert(ST_SURFACE_FRONT_LEFT == BUFFER_FRONT_LEFT);
   assert(ST_SURFACE_BACK_RIGHT == BUFFER_BACK_RIGHT);

   strb = st_renderbuffer(stfb->Base.Attachment[surfIndex].Renderbuffer);
   if (strb) {
      *surface = strb->surface;
      return GL_TRUE;
   }

   *surface = NULL;
   return GL_FALSE;
}

int
st_get_framebuffer_texture(struct st_framebuffer *stfb, uint surfIndex, struct pipe_texture **texture)
{
   struct st_renderbuffer *strb;

   assert(surfIndex <= ST_SURFACE_DEPTH);

   /* sanity checks, ST tokens should match Mesa tokens */
   assert(ST_SURFACE_FRONT_LEFT == BUFFER_FRONT_LEFT);
   assert(ST_SURFACE_BACK_RIGHT == BUFFER_BACK_RIGHT);

   strb = st_renderbuffer(stfb->Base.Attachment[surfIndex].Renderbuffer);
   if (strb) {
      *texture = strb->texture;
      return GL_TRUE;
   }

   *texture = NULL;
   return GL_FALSE;
}

/**
 * This function is to be called prior to SwapBuffers on the given
 * framebuffer.  It checks if the current context is bound to the framebuffer
 * and flushes rendering if needed.
 */
void
st_notify_swapbuffers(struct st_framebuffer *stfb)
{
   GET_CURRENT_CONTEXT(ctx);

   if (ctx && ctx->DrawBuffer == &stfb->Base) {
      st_flush( ctx->st, 
		PIPE_FLUSH_RENDER_CACHE | 
		PIPE_FLUSH_SWAPBUFFERS |
		PIPE_FLUSH_FRAME,
                NULL );
      if (st_renderbuffer(stfb->Base.Attachment[BUFFER_BACK_LEFT].Renderbuffer))
         ctx->st->frontbuffer_status = FRONT_STATUS_COPY_OF_BACK;
   }
}


void *st_framebuffer_private( struct st_framebuffer *stfb )
{
   return stfb->Private;
}

void st_get_framebuffer_dimensions( struct st_framebuffer *stfb,
				    uint *width,
				    uint *height)
{
   *width = stfb->Base.Width;
   *height = stfb->Base.Height;
}
