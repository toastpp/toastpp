/**************************************************************************
 * 
 * Copyright 2007 Tungsten Graphics, Inc., Cedar Park, Texas.
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

/* Authors:  Keith Whitwell <keith@tungstengraphics.com>
 */

#include "sp_context.h"
#include "sp_state.h"
#include "sp_surface.h"
#include "sp_tile_cache.h"

#include "draw/draw_context.h"


/**
 * XXX this might get moved someday
 * Set the framebuffer surface info: color buffers, zbuffer, stencil buffer.
 * Here, we flush the old surfaces and update the tile cache to point to the new
 * surfaces.
 */
void
softpipe_set_framebuffer_state(struct pipe_context *pipe,
                               const struct pipe_framebuffer_state *fb)
{
   struct softpipe_context *sp = softpipe_context(pipe);
   uint i;

   for (i = 0; i < PIPE_MAX_COLOR_BUFS; i++) {
      /* check if changing cbuf */
      if (sp->framebuffer.cbufs[i] != fb->cbufs[i]) {
         /* flush old */
         sp_flush_tile_cache(sp, sp->cbuf_cache[i]);

         /* assign new */
         sp->framebuffer.cbufs[i] = fb->cbufs[i];

         /* update cache */
         sp_tile_cache_set_surface(sp->cbuf_cache[i], fb->cbufs[i]);
      }
   }

   sp->framebuffer.nr_cbufs = fb->nr_cbufs;

   /* zbuf changing? */
   if (sp->framebuffer.zsbuf != fb->zsbuf) {
      /* flush old */
      sp_flush_tile_cache(sp, sp->zsbuf_cache);

      /* assign new */
      sp->framebuffer.zsbuf = fb->zsbuf;

      /* update cache */
      sp_tile_cache_set_surface(sp->zsbuf_cache, fb->zsbuf);
   }

#if 0
   /* XXX combined depth/stencil here */

   /* sbuf changing? */
   if (sp->framebuffer.sbuf != fb->sbuf) {
      /* flush old */
      sp_flush_tile_cache(sp, sp->sbuf_cache_sep);

      /* assign new */
      sp->framebuffer.sbuf = fb->sbuf;

      /* update cache */
      if (fb->sbuf != fb->zbuf) {
         /* separate stencil buf */
         sp->sbuf_cache = sp->sbuf_cache_sep;
         sp_tile_cache_set_surface(sp->sbuf_cache, fb->sbuf);
      }
      else {
         /* combined depth/stencil */
         sp->sbuf_cache = sp->zbuf_cache;
         sp_tile_cache_set_surface(sp->sbuf_cache, fb->sbuf);
      }
   }
#endif

   /* Tell draw module how deep the Z/depth buffer is */
   {
      int depth_bits;
      double mrd;
      if (sp->framebuffer.zsbuf) {
         depth_bits = pf_get_component_bits(sp->framebuffer.zsbuf->format,
                                            PIPE_FORMAT_COMP_Z);
      }
      else {
         depth_bits = 0;
      }
      if (depth_bits > 16) {
         mrd = 0.0000001;
      }
      else {
         mrd = 0.00002;
      }
      draw_set_mrd(sp->draw, mrd);
   }

   sp->framebuffer.width = fb->width;
   sp->framebuffer.height = fb->height;

   sp->dirty |= SP_NEW_FRAMEBUFFER;
}
