#include "draw/draw_context.h"
#include "pipe/p_context.h"
#include "pipe/p_state.h"

#include "nv20_context.h"
#include "nv20_state.h"

#include "nouveau/nouveau_channel.h"
#include "nouveau/nouveau_pushbuf.h"

boolean nv20_draw_elements( struct pipe_context *pipe,
                    struct pipe_buffer *indexBuffer,
                    unsigned indexSize,
                    unsigned prim, unsigned start, unsigned count)
{
	struct nv20_context *nv20 = nv20_context( pipe );
	struct draw_context *draw = nv20->draw;
	unsigned i;

	nv20_emit_hw_state(nv20);

	/*
	 * Map vertex buffers
	 */
	for (i = 0; i < PIPE_MAX_ATTRIBS; i++) {
		if (nv20->vtxbuf[i].buffer) {
			void *buf
				= pipe->winsys->buffer_map(pipe->winsys,
						nv20->vtxbuf[i].buffer,
						PIPE_BUFFER_USAGE_CPU_READ);
			draw_set_mapped_vertex_buffer(draw, i, buf);
		}
	}
	/* Map index buffer, if present */
	if (indexBuffer) {
		void *mapped_indexes
			= pipe->winsys->buffer_map(pipe->winsys, indexBuffer,
					PIPE_BUFFER_USAGE_CPU_READ);
		draw_set_mapped_element_buffer(draw, indexSize, mapped_indexes);
	}
	else {
		/* no index/element buffer */
		draw_set_mapped_element_buffer(draw, 0, NULL);
	}

	draw_set_mapped_constant_buffer(draw,
					nv20->constbuf[PIPE_SHADER_VERTEX],
					nv20->constbuf_nr[PIPE_SHADER_VERTEX]);

	/* draw! */
	draw_arrays(nv20->draw, prim, start, count);

	/*
	 * unmap vertex/index buffers
	 */
	for (i = 0; i < PIPE_MAX_ATTRIBS; i++) {
		if (nv20->vtxbuf[i].buffer) {
			pipe->winsys->buffer_unmap(pipe->winsys, nv20->vtxbuf[i].buffer);
			draw_set_mapped_vertex_buffer(draw, i, NULL);
		}
	}
	if (indexBuffer) {
		pipe->winsys->buffer_unmap(pipe->winsys, indexBuffer);
		draw_set_mapped_element_buffer(draw, 0, NULL);
	}

	draw_flush(nv20->draw);
	return TRUE;
}

boolean nv20_draw_arrays( struct pipe_context *pipe,
				 unsigned prim, unsigned start, unsigned count)
{
	return nv20_draw_elements(pipe, NULL, 0, prim, start, count);
}



