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

#include "i915_context.h"
#include "i915_blit.h"
#include "i915_state.h"
#include "pipe/p_defines.h"
#include "pipe/p_inlines.h"
#include "pipe/p_inlines.h"
#include "pipe/internal/p_winsys_screen.h"
#include "util/u_tile.h"
#include "util/u_rect.h"


/* Assumes all values are within bounds -- no checking at this level -
 * do it higher up if required.
 */
static void
i915_surface_copy(struct pipe_context *pipe,
		  struct pipe_surface *dst,
		  unsigned dstx, unsigned dsty,
		  struct pipe_surface *src,
		  unsigned srcx, unsigned srcy, unsigned width, unsigned height)
{
   struct i915_texture *dst_tex = (struct i915_texture *)dst->texture;
   struct i915_texture *src_tex = (struct i915_texture *)src->texture;

   assert( dst != src );
   assert( dst_tex->base.block.size == src_tex->base.block.size );
   assert( dst_tex->base.block.width == src_tex->base.block.height );
   assert( dst_tex->base.block.height == src_tex->base.block.height );
   assert( dst_tex->base.block.width == 1 );
   assert( dst_tex->base.block.height == 1 );

   i915_copy_blit( i915_context(pipe),
                   FALSE,
                   dst_tex->base.block.size,
                   (unsigned short) src_tex->stride, src_tex->buffer, src->offset,
                   (unsigned short) dst_tex->stride, dst_tex->buffer, dst->offset,
                   (short) srcx, (short) srcy, (short) dstx, (short) dsty, (short) width, (short) height );
}


static void
i915_surface_fill(struct pipe_context *pipe,
		  struct pipe_surface *dst,
		  unsigned dstx, unsigned dsty,
		  unsigned width, unsigned height, unsigned value)
{
   struct i915_texture *tex = (struct i915_texture *)dst->texture;

   assert(tex->base.block.width == 1);
   assert(tex->base.block.height == 1);

   i915_fill_blit( i915_context(pipe),
                   tex->base.block.size,
                   (unsigned short) tex->stride,
                   tex->buffer, dst->offset,
                   (short) dstx, (short) dsty,
                   (short) width, (short) height,
                   value );
}


void
i915_init_surface_functions(struct i915_context *i915)
{
   i915->pipe.surface_copy = i915_surface_copy;
   i915->pipe.surface_fill = i915_surface_fill;
}
