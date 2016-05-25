#include "draw/draw_context.h"
#include "pipe/p_defines.h"
#include "pipe/internal/p_winsys_screen.h"

#include "nv30_context.h"
#include "nv30_screen.h"

static void
nv30_flush(struct pipe_context *pipe, unsigned flags,
	   struct pipe_fence_handle **fence)
{
	struct nv30_context *nv30 = nv30_context(pipe);
	
	if (flags & PIPE_FLUSH_TEXTURE_CACHE) {
		BEGIN_RING(rankine, 0x1fd8, 1);
		OUT_RING  (2);
		BEGIN_RING(rankine, 0x1fd8, 1);
		OUT_RING  (1);
	}

	FIRE_RING(fence);
}

static void
nv30_destroy(struct pipe_context *pipe)
{
	struct nv30_context *nv30 = nv30_context(pipe);

	if (nv30->draw)
		draw_destroy(nv30->draw);
	FREE(nv30);
}

static unsigned int
nv30_is_texture_referenced( struct pipe_context *pipe,
			    struct pipe_texture *texture,
			    unsigned face, unsigned level)
{
   /**
    * FIXME: Optimize.
    */

   return PIPE_REFERENCED_FOR_READ | PIPE_REFERENCED_FOR_WRITE;
}

static unsigned int
nv30_is_buffer_referenced( struct pipe_context *pipe,
			   struct pipe_buffer *buf)
{
   /**
    * FIXME: Optimize.
    */

   return PIPE_REFERENCED_FOR_READ | PIPE_REFERENCED_FOR_WRITE;
}

struct pipe_context *
nv30_create(struct pipe_screen *pscreen, unsigned pctx_id)
{
	struct nv30_screen *screen = nv30_screen(pscreen);
	struct pipe_winsys *ws = pscreen->winsys;
	struct nv30_context *nv30;
	struct nouveau_winsys *nvws = screen->nvws;

	nv30 = CALLOC(1, sizeof(struct nv30_context));
	if (!nv30)
		return NULL;
	nv30->screen = screen;
	nv30->pctx_id = pctx_id;

	nv30->nvws = nvws;

	nv30->pipe.winsys = ws;
	nv30->pipe.screen = pscreen;
	nv30->pipe.destroy = nv30_destroy;
	nv30->pipe.draw_arrays = nv30_draw_arrays;
	nv30->pipe.draw_elements = nv30_draw_elements;
	nv30->pipe.clear = nv30_clear;
	nv30->pipe.flush = nv30_flush;

	nv30->pipe.is_texture_referenced = nv30_is_texture_referenced;
	nv30->pipe.is_buffer_referenced = nv30_is_buffer_referenced;

	nv30_init_query_functions(nv30);
	nv30_init_surface_functions(nv30);
	nv30_init_state_functions(nv30);

	/* Create, configure, and install fallback swtnl path */
	nv30->draw = draw_create();
	draw_wide_point_threshold(nv30->draw, 9999999.0);
	draw_wide_line_threshold(nv30->draw, 9999999.0);
	draw_enable_line_stipple(nv30->draw, FALSE);
	draw_enable_point_sprites(nv30->draw, FALSE);
	draw_set_rasterize_stage(nv30->draw, nv30_draw_render_stage(nv30));

	return &nv30->pipe;
}
	
