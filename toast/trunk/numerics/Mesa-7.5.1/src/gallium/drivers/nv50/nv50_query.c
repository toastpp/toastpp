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

#include "pipe/p_context.h"
#include "pipe/p_inlines.h"

#include "nv50_context.h"

struct nv50_query {
	struct pipe_buffer *buffer;
	unsigned type;
	boolean ready;
	uint64_t result;
};

static INLINE struct nv50_query *
nv50_query(struct pipe_query *pipe)
{
	return (struct nv50_query *)pipe;
}

static struct pipe_query *
nv50_query_create(struct pipe_context *pipe, unsigned type)
{
	struct pipe_screen *screen = pipe->screen;
	struct nv50_query *q = CALLOC_STRUCT(nv50_query);

	assert (q->type == PIPE_QUERY_OCCLUSION_COUNTER);
	q->type = type;

	q->buffer = screen->buffer_create(screen, 256, 0, 16);
	if (!q->buffer) {
		FREE(q);
		return NULL;
	}

	return (struct pipe_query *)q;
}

static void
nv50_query_destroy(struct pipe_context *pipe, struct pipe_query *pq)
{
	struct nv50_query *q = nv50_query(pq);

	if (q) {
		pipe_buffer_reference(&q->buffer, NULL);
		FREE(q);
	}
}

static void
nv50_query_begin(struct pipe_context *pipe, struct pipe_query *pq)
{
	struct nv50_context *nv50 = nv50_context(pipe);
	struct nouveau_channel *chan = nv50->screen->nvws->channel;
	struct nouveau_grobj *tesla = nv50->screen->tesla;
	struct nv50_query *q = nv50_query(pq);

	BEGIN_RING(chan, tesla, 0x1530, 1);
	OUT_RING  (chan, 1);
	BEGIN_RING(chan, tesla, 0x1514, 1);
	OUT_RING  (chan, 1);

	q->ready = FALSE;
}

static void
nv50_query_end(struct pipe_context *pipe, struct pipe_query *pq)
{
	struct nv50_context *nv50 = nv50_context(pipe);
	struct nouveau_channel *chan = nv50->screen->nvws->channel;
	struct nouveau_grobj *tesla = nv50->screen->tesla;
	struct nv50_query *q = nv50_query(pq);
	struct nouveau_bo *bo = nv50->screen->nvws->get_bo(q->buffer);

	WAIT_RING (chan, 5);
	BEGIN_RING(chan, tesla, 0x1b00, 4);
	OUT_RELOCh(chan, bo, 0, NOUVEAU_BO_VRAM | NOUVEAU_BO_WR);
	OUT_RELOCl(chan, bo, 0, NOUVEAU_BO_VRAM | NOUVEAU_BO_WR);
	OUT_RING  (chan, 0x00000000);
	OUT_RING  (chan, 0x0100f002);
	FIRE_RING (chan);
}

static boolean
nv50_query_result(struct pipe_context *pipe, struct pipe_query *pq,
		  boolean wait, uint64_t *result)
{
	struct pipe_winsys *ws = pipe->winsys;
	struct nv50_query *q = nv50_query(pq);

	/*XXX: Want to be able to return FALSE here instead of blocking
	 *     until the result is available..
	 */

	if (!q->ready) {
		uint32_t *map = ws->buffer_map(ws, q->buffer,
					       PIPE_BUFFER_USAGE_CPU_READ);
		q->result = map[1];
		q->ready = TRUE;
		ws->buffer_unmap(ws, q->buffer);
	}

	*result = q->result;
	return q->ready;
}

void
nv50_init_query_functions(struct nv50_context *nv50)
{
	nv50->pipe.create_query = nv50_query_create;
	nv50->pipe.destroy_query = nv50_query_destroy;
	nv50->pipe.begin_query = nv50_query_begin;
	nv50->pipe.end_query = nv50_query_end;
	nv50->pipe.get_query_result = nv50_query_result;
}
