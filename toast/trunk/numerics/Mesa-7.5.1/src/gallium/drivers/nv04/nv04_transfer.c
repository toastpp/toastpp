#include <pipe/p_state.h>
#include <pipe/p_defines.h>
#include <pipe/p_inlines.h>
#include <util/u_memory.h>
#include <nouveau/nouveau_winsys.h>
#include "nv04_context.h"
#include "nv04_screen.h"
#include "nv04_state.h"

struct nv04_transfer {
	struct pipe_transfer base;
	struct pipe_surface *surface;
	bool direct;
};

static unsigned nv04_usage_tx_to_buf(unsigned tx_usage)
{
	switch (tx_usage) {
		case PIPE_TRANSFER_READ:
			return PIPE_BUFFER_USAGE_CPU_READ;
		case PIPE_TRANSFER_WRITE:
			return PIPE_BUFFER_USAGE_CPU_WRITE;
		case PIPE_TRANSFER_READ_WRITE:
			return PIPE_BUFFER_USAGE_CPU_READ_WRITE;
		default:
			assert(0);
	}

	return -1;
}

static void
nv04_compatible_transfer_tex(struct pipe_texture *pt, unsigned level,
                             struct pipe_texture *template)
{
	memset(template, 0, sizeof(struct pipe_texture));
	template->target = pt->target;
	template->format = pt->format;
	template->width[0] = pt->width[level];
	template->height[0] = pt->height[level];
	template->depth[0] = 1;
	template->block = pt->block;
	template->nblocksx[0] = pt->nblocksx[level];
	template->nblocksy[0] = pt->nblocksx[level];
	template->last_level = 0;
	template->nr_samples = pt->nr_samples;

	template->tex_usage = PIPE_TEXTURE_USAGE_DYNAMIC |
	                      NOUVEAU_TEXTURE_USAGE_LINEAR;
}

static struct pipe_transfer *
nv04_transfer_new(struct pipe_screen *pscreen, struct pipe_texture *pt,
		  unsigned face, unsigned level, unsigned zslice,
		  enum pipe_transfer_usage usage,
		  unsigned x, unsigned y, unsigned w, unsigned h)
{
	struct nv04_miptree *mt = (struct nv04_miptree *)pt;
	struct nv04_transfer *tx;
	struct pipe_texture tx_tex_template, *tx_tex;

	tx = CALLOC_STRUCT(nv04_transfer);
	if (!tx)
		return NULL;

	pipe_texture_reference(&tx->base.texture, pt);
	tx->base.format = pt->format;
	tx->base.x = x;
	tx->base.y = y;
	tx->base.width = w;
	tx->base.height = h;
	tx->base.block = pt->block;
	tx->base.nblocksx = pt->nblocksx[level];
	tx->base.nblocksy = pt->nblocksy[level];
	tx->base.stride = mt->level[level].pitch;
	tx->base.usage = usage;
	tx->base.face = face;
	tx->base.level = level;
	tx->base.zslice = zslice;

	/* Direct access to texture */
	if ((pt->tex_usage & PIPE_TEXTURE_USAGE_DYNAMIC ||
	     debug_get_bool_option("NOUVEAU_NO_TRANSFER", TRUE/*XXX:FALSE*/)) &&
	    pt->tex_usage & NOUVEAU_TEXTURE_USAGE_LINEAR)
	{
		tx->direct = true;
		tx->surface = pscreen->get_tex_surface(pscreen, pt,
	                                               0, 0, 0,
	                                               nv04_usage_tx_to_buf(usage));
		return &tx->base;
	}

	tx->direct = false;

	nv04_compatible_transfer_tex(pt, level, &tx_tex_template);

	tx_tex = pscreen->texture_create(pscreen, &tx_tex_template);
	if (!tx_tex)
	{
		FREE(tx);
		return NULL;
	}

	tx->surface = pscreen->get_tex_surface(pscreen, tx_tex,
	                                       face, level, zslice,
	                                       nv04_usage_tx_to_buf(usage));

	pipe_texture_reference(&tx_tex, NULL);

	if (!tx->surface)
	{
		pipe_surface_reference(&tx->surface, NULL);
		FREE(tx);
		return NULL;
	}

	if (usage != PIPE_TRANSFER_WRITE) {
		struct nv04_screen *nvscreen = nv04_screen(pscreen);
		struct pipe_surface *src;

		src = pscreen->get_tex_surface(pscreen, pt,
	                                       face, level, zslice,
	                                       PIPE_BUFFER_USAGE_GPU_READ);

		/* TODO: Check if SIFM can deal with x,y,w,h when swizzling */
		/* TODO: Check if SIFM can un-swizzle */
		nvscreen->eng2d->copy(nvscreen->eng2d,
		                      tx->surface, 0, 0,
		                      src, 0, 0,
		                      src->width, src->height);

		pipe_surface_reference(&src, NULL);
	}

	return &tx->base;
}

static void
nv04_transfer_del(struct pipe_transfer *ptx)
{
	struct nv04_transfer *tx = (struct nv04_transfer *)ptx;

	if (!tx->direct && ptx->usage != PIPE_TRANSFER_READ) {
		struct pipe_screen *pscreen = ptx->texture->screen;
		struct nv04_screen *nvscreen = nv04_screen(pscreen);
		struct pipe_surface *dst;

		dst = pscreen->get_tex_surface(pscreen, ptx->texture,
	                                       ptx->face, ptx->level, ptx->zslice,
	                                       PIPE_BUFFER_USAGE_GPU_WRITE);

		/* TODO: Check if SIFM can deal with x,y,w,h when swizzling */
		nvscreen->eng2d->copy(nvscreen->eng2d,
		                      dst, 0, 0,
		                      tx->surface, 0, 0,
		                      dst->width, dst->height);

		pipe_surface_reference(&dst, NULL);
	}

	pipe_surface_reference(&tx->surface, NULL);
	pipe_texture_reference(&ptx->texture, NULL);
	FREE(ptx);
}

static void *
nv04_transfer_map(struct pipe_screen *pscreen, struct pipe_transfer *ptx)
{
	struct nv04_transfer *tx = (struct nv04_transfer *)ptx;
	struct nv04_surface *ns = (struct nv04_surface *)tx->surface;
	struct nv04_miptree *mt = (struct nv04_miptree *)tx->surface->texture;
	void *map = pipe_buffer_map(pscreen, mt->buffer,
	                            nv04_usage_tx_to_buf(ptx->usage));

	return map + ns->base.offset +
	       ptx->y * ns->pitch + ptx->x * ptx->block.size;
}

static void
nv04_transfer_unmap(struct pipe_screen *pscreen, struct pipe_transfer *ptx)
{
	struct nv04_transfer *tx = (struct nv04_transfer *)ptx;
	struct nv04_miptree *mt = (struct nv04_miptree *)tx->surface->texture;

	pipe_buffer_unmap(pscreen, mt->buffer);
}

void
nv04_screen_init_transfer_functions(struct pipe_screen *pscreen)
{
	pscreen->get_tex_transfer = nv04_transfer_new;
	pscreen->tex_transfer_destroy = nv04_transfer_del;
	pscreen->transfer_map = nv04_transfer_map;
	pscreen->transfer_unmap = nv04_transfer_unmap;
}
