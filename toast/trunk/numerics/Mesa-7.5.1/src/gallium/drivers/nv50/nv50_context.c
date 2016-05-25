/*
 * Copyright 2008 Ben Skeggs
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
 * THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF
 * OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include "draw/draw_context.h"
#include "pipe/p_defines.h"
#include "pipe/internal/p_winsys_screen.h"

#include "nv50_context.h"
#include "nv50_screen.h"

static void
nv50_flush(struct pipe_context *pipe, unsigned flags,
	   struct pipe_fence_handle **fence)
{
	struct nv50_context *nv50 = (struct nv50_context *)pipe;
	
	FIRE_RING(nv50->screen->nvws->channel);
}

static void
nv50_destroy(struct pipe_context *pipe)
{
	struct nv50_context *nv50 = (struct nv50_context *)pipe;

	draw_destroy(nv50->draw);
	FREE(nv50);
}


static void
nv50_set_edgeflags(struct pipe_context *pipe, const unsigned *bitfield)
{
}

static unsigned int
nv50_is_texture_referenced( struct pipe_context *pipe,
			    struct pipe_texture *texture,
			    unsigned face, unsigned level)
{
   /**
    * FIXME: Optimize.
    */

   return PIPE_REFERENCED_FOR_READ | PIPE_REFERENCED_FOR_WRITE;
}

static unsigned int
nv50_is_buffer_referenced( struct pipe_context *pipe,
			   struct pipe_buffer *buf)
{
   /**
    * FIXME: Optimize.
    */

   return PIPE_REFERENCED_FOR_READ | PIPE_REFERENCED_FOR_WRITE;
}

struct pipe_context *
nv50_create(struct pipe_screen *pscreen, unsigned pctx_id)
{
	struct pipe_winsys *pipe_winsys = pscreen->winsys;
	struct nv50_screen *screen = nv50_screen(pscreen);
	struct nv50_context *nv50;

	nv50 = CALLOC_STRUCT(nv50_context);
	if (!nv50)
		return NULL;
	nv50->screen = screen;
	nv50->pctx_id = pctx_id;

	nv50->pipe.winsys = pipe_winsys;
	nv50->pipe.screen = pscreen;

	nv50->pipe.destroy = nv50_destroy;

	nv50->pipe.set_edgeflags = nv50_set_edgeflags;
	nv50->pipe.draw_arrays = nv50_draw_arrays;
	nv50->pipe.draw_elements = nv50_draw_elements;
	nv50->pipe.clear = nv50_clear;

	nv50->pipe.flush = nv50_flush;

	nv50->pipe.is_texture_referenced = nv50_is_texture_referenced;
	nv50->pipe.is_buffer_referenced = nv50_is_buffer_referenced;

	nv50_init_surface_functions(nv50);
	nv50_init_state_functions(nv50);
	nv50_init_query_functions(nv50);

	nv50->draw = draw_create();
	assert(nv50->draw);
	draw_set_rasterize_stage(nv50->draw, nv50_draw_render_stage(nv50));

	return &nv50->pipe;
}

		
