#include "pipe/p_context.h"
#include "pipe/p_format.h"
#include "util/u_memory.h"

#include "nouveau/nouveau_winsys.h"
#include "nouveau/nouveau_util.h"
#include "nv04_surface_2d.h"

static INLINE int
nv04_surface_format(enum pipe_format format)
{
	switch (format) {
	case PIPE_FORMAT_A8_UNORM:
		return NV04_CONTEXT_SURFACES_2D_FORMAT_Y8;
	case PIPE_FORMAT_R16_SNORM:
	case PIPE_FORMAT_R5G6B5_UNORM:
		return NV04_CONTEXT_SURFACES_2D_FORMAT_R5G6B5;
	case PIPE_FORMAT_X8R8G8B8_UNORM:
	case PIPE_FORMAT_A8R8G8B8_UNORM:
		return NV04_CONTEXT_SURFACES_2D_FORMAT_A8R8G8B8;
	case PIPE_FORMAT_Z24S8_UNORM:
		return NV04_CONTEXT_SURFACES_2D_FORMAT_Y32;
	default:
		return -1;
	}
}

static INLINE int
nv04_rect_format(enum pipe_format format)
{
	switch (format) {
	case PIPE_FORMAT_A8_UNORM:
		return NV04_GDI_RECTANGLE_TEXT_COLOR_FORMAT_A8R8G8B8;
	case PIPE_FORMAT_R5G6B5_UNORM:
		return NV04_GDI_RECTANGLE_TEXT_COLOR_FORMAT_A16R5G6B5;
	case PIPE_FORMAT_A8R8G8B8_UNORM:
	case PIPE_FORMAT_Z24S8_UNORM:
		return NV04_GDI_RECTANGLE_TEXT_COLOR_FORMAT_A8R8G8B8;
	default:
		return -1;
	}
}

static INLINE int
nv04_scaled_image_format(enum pipe_format format)
{
	switch (format) {
	case PIPE_FORMAT_A1R5G5B5_UNORM:
		return NV04_SCALED_IMAGE_FROM_MEMORY_COLOR_FORMAT_A1R5G5B5;
	case PIPE_FORMAT_A8R8G8B8_UNORM:
		return NV04_SCALED_IMAGE_FROM_MEMORY_COLOR_FORMAT_A8R8G8B8;
	case PIPE_FORMAT_X8R8G8B8_UNORM:
		return NV04_SCALED_IMAGE_FROM_MEMORY_COLOR_FORMAT_X8R8G8B8;
	case PIPE_FORMAT_R5G6B5_UNORM:
	case PIPE_FORMAT_R16_SNORM:
		return NV04_SCALED_IMAGE_FROM_MEMORY_COLOR_FORMAT_R5G6B5;
	default:
		return -1;
	}
}

static INLINE unsigned
nv04_swizzle_bits(unsigned x, unsigned y)
{
	unsigned u = (x & 0x001) << 0 |
	             (x & 0x002) << 1 |
	             (x & 0x004) << 2 |
	             (x & 0x008) << 3 |
	             (x & 0x010) << 4 |
	             (x & 0x020) << 5 |
	             (x & 0x040) << 6 |
	             (x & 0x080) << 7 |
	             (x & 0x100) << 8 |
	             (x & 0x200) << 9 |
	             (x & 0x400) << 10 |
	             (x & 0x800) << 11;

	unsigned v = (y & 0x001) << 1 |
	             (y & 0x002) << 2 |
	             (y & 0x004) << 3 |
	             (y & 0x008) << 4 |
	             (y & 0x010) << 5 |
	             (y & 0x020) << 6 |
	             (y & 0x040) << 7 |
	             (y & 0x080) << 8 |
	             (y & 0x100) << 9 |
	             (y & 0x200) << 10 |
	             (y & 0x400) << 11 |
	             (y & 0x800) << 12;
	return v | u;
}

static int
nv04_surface_copy_swizzle(struct nv04_surface_2d *ctx,
			  struct pipe_surface *dst, int dx, int dy,
			  struct pipe_surface *src, int sx, int sy,
			  int w, int h)
{
	struct nouveau_channel *chan = ctx->nvws->channel;
	struct nouveau_grobj *swzsurf = ctx->swzsurf;
	struct nouveau_grobj *sifm = ctx->sifm;
	struct nouveau_bo *src_bo = ctx->nvws->get_bo(ctx->buf(src));
	struct nouveau_bo *dst_bo = ctx->nvws->get_bo(ctx->buf(dst));
	const unsigned src_pitch = ((struct nv04_surface *)src)->pitch;
	const unsigned max_w = 1024;
	const unsigned max_h = 1024;
	const unsigned sub_w = w > max_w ? max_w : w;
	const unsigned sub_h = h > max_h ? max_h : h;
	unsigned cx;
	unsigned cy;

	/* POT or GTFO */
	assert(!(w & (w - 1)) && !(h & (h - 1)));
	/* That's the way she likes it */
	assert(src_pitch == ((struct nv04_surface *)dst)->pitch);

	BEGIN_RING(chan, swzsurf, NV04_SWIZZLED_SURFACE_DMA_IMAGE, 1);
	OUT_RELOCo(chan, dst_bo,
	                 NOUVEAU_BO_GART | NOUVEAU_BO_VRAM | NOUVEAU_BO_WR);

	BEGIN_RING(chan, swzsurf, NV04_SWIZZLED_SURFACE_FORMAT, 1);
	OUT_RING  (chan, nv04_surface_format(dst->format) |
	                 log2i(w) << NV04_SWIZZLED_SURFACE_FORMAT_BASE_SIZE_U_SHIFT |
	                 log2i(h) << NV04_SWIZZLED_SURFACE_FORMAT_BASE_SIZE_V_SHIFT);
 
	BEGIN_RING(chan, sifm, NV04_SCALED_IMAGE_FROM_MEMORY_DMA_IMAGE, 1);
	OUT_RELOCo(chan, src_bo,
	                 NOUVEAU_BO_GART | NOUVEAU_BO_VRAM | NOUVEAU_BO_RD);
	BEGIN_RING(chan, sifm, NV04_SCALED_IMAGE_FROM_MEMORY_SURFACE, 1);
	OUT_RING  (chan, swzsurf->handle);

	for (cy = 0; cy < h; cy += sub_h) {
	  for (cx = 0; cx < w; cx += sub_w) {
	    BEGIN_RING(chan, swzsurf, NV04_SWIZZLED_SURFACE_OFFSET, 1);
	    OUT_RELOCl(chan, dst_bo, dst->offset + nv04_swizzle_bits(cx, cy) *
			     dst->texture->block.size, NOUVEAU_BO_GART |
			     NOUVEAU_BO_VRAM | NOUVEAU_BO_WR);

	    BEGIN_RING(chan, sifm, NV04_SCALED_IMAGE_FROM_MEMORY_COLOR_CONVERSION, 9);
	    OUT_RING  (chan, NV04_SCALED_IMAGE_FROM_MEMORY_COLOR_CONVERSION_TRUNCATE);
	    OUT_RING  (chan, nv04_scaled_image_format(src->format));
	    OUT_RING  (chan, NV04_SCALED_IMAGE_FROM_MEMORY_OPERATION_SRCCOPY);
	    OUT_RING  (chan, 0);
	    OUT_RING  (chan, sub_h << 16 | sub_w);
	    OUT_RING  (chan, 0);
	    OUT_RING  (chan, sub_h << 16 | sub_w);
	    OUT_RING  (chan, 1 << 20);
	    OUT_RING  (chan, 1 << 20);

	    BEGIN_RING(chan, sifm, NV04_SCALED_IMAGE_FROM_MEMORY_SIZE, 4);
	    OUT_RING  (chan, sub_h << 16 | sub_w);
	    OUT_RING  (chan, src_pitch |
			     NV04_SCALED_IMAGE_FROM_MEMORY_FORMAT_ORIGIN_CENTER |
			     NV04_SCALED_IMAGE_FROM_MEMORY_FORMAT_FILTER_POINT_SAMPLE);
	    OUT_RELOCl(chan, src_bo, src->offset + cy * src_pitch +
			     cx * src->texture->block.size, NOUVEAU_BO_GART |
			     NOUVEAU_BO_VRAM | NOUVEAU_BO_RD);
	    OUT_RING  (chan, 0);
	  }
	}

	return 0;
}

static int
nv04_surface_copy_m2mf(struct nv04_surface_2d *ctx,
		       struct pipe_surface *dst, int dx, int dy,
		       struct pipe_surface *src, int sx, int sy, int w, int h)
{
	struct nouveau_channel *chan = ctx->nvws->channel;
	struct nouveau_grobj *m2mf = ctx->m2mf;
	struct nouveau_bo *src_bo = ctx->nvws->get_bo(ctx->buf(src));
	struct nouveau_bo *dst_bo = ctx->nvws->get_bo(ctx->buf(dst));
	unsigned src_pitch = ((struct nv04_surface *)src)->pitch;
	unsigned dst_pitch = ((struct nv04_surface *)dst)->pitch;
	unsigned dst_offset = dst->offset + dy * dst_pitch +
	                      dx * dst->texture->block.size;
	unsigned src_offset = src->offset + sy * src_pitch +
	                      sx * src->texture->block.size;

	WAIT_RING (chan, 3 + ((h / 2047) + 1) * 9);
	BEGIN_RING(chan, m2mf, NV04_MEMORY_TO_MEMORY_FORMAT_DMA_BUFFER_IN, 2);
	OUT_RELOCo(chan, src_bo,
		   NOUVEAU_BO_GART | NOUVEAU_BO_VRAM | NOUVEAU_BO_RD);
	OUT_RELOCo(chan, dst_bo,
		   NOUVEAU_BO_GART | NOUVEAU_BO_VRAM | NOUVEAU_BO_WR);

	while (h) {
		int count = (h > 2047) ? 2047 : h;

		BEGIN_RING(chan, m2mf, NV04_MEMORY_TO_MEMORY_FORMAT_OFFSET_IN, 8);
		OUT_RELOCl(chan, src_bo, src_offset,
			   NOUVEAU_BO_VRAM | NOUVEAU_BO_GART | NOUVEAU_BO_RD);
		OUT_RELOCl(chan, dst_bo, dst_offset,
			   NOUVEAU_BO_VRAM | NOUVEAU_BO_GART | NOUVEAU_BO_WR);
		OUT_RING  (chan, src_pitch);
		OUT_RING  (chan, dst_pitch);
		OUT_RING  (chan, w * src->texture->block.size);
		OUT_RING  (chan, count);
		OUT_RING  (chan, 0x0101);
		OUT_RING  (chan, 0);

		h -= count;
		src_offset += src_pitch * count;
		dst_offset += dst_pitch * count;
	}

	return 0;
}

static int
nv04_surface_copy_blit(struct nv04_surface_2d *ctx, struct pipe_surface *dst,
		       int dx, int dy, struct pipe_surface *src, int sx, int sy,
		       int w, int h)
{
	struct nouveau_channel *chan = ctx->nvws->channel;
	struct nouveau_grobj *surf2d = ctx->surf2d;
	struct nouveau_grobj *blit = ctx->blit;
	struct nouveau_bo *src_bo = ctx->nvws->get_bo(ctx->buf(src));
	struct nouveau_bo *dst_bo = ctx->nvws->get_bo(ctx->buf(dst));
	unsigned src_pitch = ((struct nv04_surface *)src)->pitch;
	unsigned dst_pitch = ((struct nv04_surface *)dst)->pitch;
	int format;

	format = nv04_surface_format(dst->format);
	if (format < 0)
		return 1;

	WAIT_RING (chan, 12);
	BEGIN_RING(chan, surf2d, NV04_CONTEXT_SURFACES_2D_DMA_IMAGE_SOURCE, 2);
	OUT_RELOCo(chan, src_bo, NOUVEAU_BO_VRAM | NOUVEAU_BO_RD);
	OUT_RELOCo(chan, dst_bo, NOUVEAU_BO_VRAM | NOUVEAU_BO_WR);
	BEGIN_RING(chan, surf2d, NV04_CONTEXT_SURFACES_2D_FORMAT, 4);
	OUT_RING  (chan, format);
	OUT_RING  (chan, (dst_pitch << 16) | src_pitch);
	OUT_RELOCl(chan, src_bo, src->offset, NOUVEAU_BO_VRAM | NOUVEAU_BO_RD);
	OUT_RELOCl(chan, dst_bo, dst->offset, NOUVEAU_BO_VRAM | NOUVEAU_BO_WR);

	BEGIN_RING(chan, blit, 0x0300, 3);
	OUT_RING  (chan, (sy << 16) | sx);
	OUT_RING  (chan, (dy << 16) | dx);
	OUT_RING  (chan, ( h << 16) |  w);

	return 0;
}

static void
nv04_surface_copy(struct nv04_surface_2d *ctx, struct pipe_surface *dst,
		  int dx, int dy, struct pipe_surface *src, int sx, int sy,
		  int w, int h)
{
	unsigned src_pitch = ((struct nv04_surface *)src)->pitch;
	unsigned dst_pitch = ((struct nv04_surface *)dst)->pitch;
	int src_linear = src->texture->tex_usage & NOUVEAU_TEXTURE_USAGE_LINEAR;
	int dst_linear = dst->texture->tex_usage & NOUVEAU_TEXTURE_USAGE_LINEAR;

	assert(src->format == dst->format);

	/* Setup transfer to swizzle the texture to vram if needed */
	if (src_linear && !dst_linear && w > 1 && h > 1) {
		nv04_surface_copy_swizzle(ctx, dst, dx, dy, src, sx, sy, w, h);
		return;
	}

	/* NV_CONTEXT_SURFACES_2D has buffer alignment restrictions, fallback
	 * to NV_MEMORY_TO_MEMORY_FORMAT in this case.
	 */
	if ((src->offset & 63) || (dst->offset & 63) ||
	    (src_pitch & 63) || (dst_pitch & 63) ||
	    debug_get_bool_option("NOUVEAU_NO_COPYBLIT", FALSE)) {
		nv04_surface_copy_m2mf(ctx, dst, dx, dy, src, sx, sy, w, h);
		return;
	}

	nv04_surface_copy_blit(ctx, dst, dx, dy, src, sx, sy, w, h);
}

static void
nv04_surface_fill(struct nv04_surface_2d *ctx, struct pipe_surface *dst,
		  int dx, int dy, int w, int h, unsigned value)
{
	struct nouveau_channel *chan = ctx->nvws->channel;
	struct nouveau_grobj *surf2d = ctx->surf2d;
	struct nouveau_grobj *rect = ctx->rect;
	struct nouveau_bo *dst_bo = ctx->nvws->get_bo(ctx->buf(dst));
	unsigned dst_pitch = ((struct nv04_surface *)dst)->pitch;
	int cs2d_format, gdirect_format;

	cs2d_format = nv04_surface_format(dst->format);
	assert(cs2d_format >= 0);

	gdirect_format = nv04_rect_format(dst->format);
	assert(gdirect_format >= 0);

	WAIT_RING (chan, 16);
	BEGIN_RING(chan, surf2d, NV04_CONTEXT_SURFACES_2D_DMA_IMAGE_SOURCE, 2);
	OUT_RELOCo(chan, dst_bo, NOUVEAU_BO_VRAM | NOUVEAU_BO_WR);
	OUT_RELOCo(chan, dst_bo, NOUVEAU_BO_VRAM | NOUVEAU_BO_WR);
	BEGIN_RING(chan, surf2d, NV04_CONTEXT_SURFACES_2D_FORMAT, 4);
	OUT_RING  (chan, cs2d_format);
	OUT_RING  (chan, (dst_pitch << 16) | dst_pitch);
	OUT_RELOCl(chan, dst_bo, dst->offset, NOUVEAU_BO_VRAM | NOUVEAU_BO_WR);
	OUT_RELOCl(chan, dst_bo, dst->offset, NOUVEAU_BO_VRAM | NOUVEAU_BO_WR);

	BEGIN_RING(chan, rect, NV04_GDI_RECTANGLE_TEXT_COLOR_FORMAT, 1);
	OUT_RING  (chan, gdirect_format);
	BEGIN_RING(chan, rect, NV04_GDI_RECTANGLE_TEXT_COLOR1_A, 1);
	OUT_RING  (chan, value);
	BEGIN_RING(chan, rect,
		   NV04_GDI_RECTANGLE_TEXT_UNCLIPPED_RECTANGLE_POINT(0), 2);
	OUT_RING  (chan, (dx << 16) | dy);
	OUT_RING  (chan, ( w << 16) |  h);
}

void
nv04_surface_2d_takedown(struct nv04_surface_2d **pctx)
{
	struct nv04_surface_2d *ctx;

	if (!pctx || !*pctx)
		return;
	ctx = *pctx;
	*pctx = NULL;

	nouveau_notifier_free(&ctx->ntfy);
	nouveau_grobj_free(&ctx->m2mf);
	nouveau_grobj_free(&ctx->surf2d);
	nouveau_grobj_free(&ctx->swzsurf);
	nouveau_grobj_free(&ctx->rect);
	nouveau_grobj_free(&ctx->blit);
	nouveau_grobj_free(&ctx->sifm);

	FREE(ctx);
}

struct nv04_surface_2d *
nv04_surface_2d_init(struct nouveau_winsys *nvws)
{
	struct nv04_surface_2d *ctx = CALLOC_STRUCT(nv04_surface_2d);
	struct nouveau_channel *chan = nvws->channel;
	unsigned handle = 0x88000000, class;
	int ret;

	if (!ctx)
		return NULL;

	ret = nouveau_notifier_alloc(chan, handle++, 1, &ctx->ntfy);
	if (ret) {
		nv04_surface_2d_takedown(&ctx);
		return NULL;
	}

	ret = nouveau_grobj_alloc(chan, handle++, 0x0039, &ctx->m2mf);
	if (ret) {
		nv04_surface_2d_takedown(&ctx);
		return NULL;
	}

	BEGIN_RING(chan, ctx->m2mf, NV04_MEMORY_TO_MEMORY_FORMAT_DMA_NOTIFY, 1);
	OUT_RING  (chan, ctx->ntfy->handle);

	if (chan->device->chipset < 0x10)
		class = NV04_CONTEXT_SURFACES_2D;
	else
		class = NV10_CONTEXT_SURFACES_2D;

	ret = nouveau_grobj_alloc(chan, handle++, class, &ctx->surf2d);
	if (ret) {
		nv04_surface_2d_takedown(&ctx);
		return NULL;
	}

	BEGIN_RING(chan, ctx->surf2d,
			 NV04_CONTEXT_SURFACES_2D_DMA_IMAGE_SOURCE, 2);
	OUT_RING  (chan, chan->vram->handle);
	OUT_RING  (chan, chan->vram->handle);

	if (chan->device->chipset < 0x10)
		class = NV04_IMAGE_BLIT;
	else
		class = NV12_IMAGE_BLIT;

	ret = nouveau_grobj_alloc(chan, handle++, class, &ctx->blit);
	if (ret) {
		nv04_surface_2d_takedown(&ctx);
		return NULL;
	}

	BEGIN_RING(chan, ctx->blit, NV04_IMAGE_BLIT_DMA_NOTIFY, 1);
	OUT_RING  (chan, ctx->ntfy->handle);
	BEGIN_RING(chan, ctx->blit, NV04_IMAGE_BLIT_SURFACE, 1);
	OUT_RING  (chan, ctx->surf2d->handle);
	BEGIN_RING(chan, ctx->blit, NV04_IMAGE_BLIT_OPERATION, 1);
	OUT_RING  (chan, NV04_IMAGE_BLIT_OPERATION_SRCCOPY);

	ret = nouveau_grobj_alloc(chan, handle++, NV04_GDI_RECTANGLE_TEXT,
				  &ctx->rect);
	if (ret) {
		nv04_surface_2d_takedown(&ctx);
		return NULL;
	}

	BEGIN_RING(chan, ctx->rect, NV04_GDI_RECTANGLE_TEXT_DMA_NOTIFY, 1);
	OUT_RING  (chan, ctx->ntfy->handle);
	BEGIN_RING(chan, ctx->rect, NV04_GDI_RECTANGLE_TEXT_SURFACE, 1);
	OUT_RING  (chan, ctx->surf2d->handle);
	BEGIN_RING(chan, ctx->rect, NV04_GDI_RECTANGLE_TEXT_OPERATION, 1);
	OUT_RING  (chan, NV04_GDI_RECTANGLE_TEXT_OPERATION_SRCCOPY);
	BEGIN_RING(chan, ctx->rect,
			 NV04_GDI_RECTANGLE_TEXT_MONOCHROME_FORMAT, 1);
	OUT_RING  (chan, NV04_GDI_RECTANGLE_TEXT_MONOCHROME_FORMAT_LE);

	switch (chan->device->chipset & 0xf0) {
	case 0x00:
	case 0x10:
		class = NV04_SWIZZLED_SURFACE;
		break;
	case 0x20:
		class = NV20_SWIZZLED_SURFACE;
		break;
	case 0x30:
		class = NV30_SWIZZLED_SURFACE;
		break;
	case 0x40:
	case 0x60:
		class = NV40_SWIZZLED_SURFACE;
		break;
	default:
		/* Famous last words: this really can't happen.. */
		assert(0);
		break;
	}

	ret = nouveau_grobj_alloc(chan, handle++, class, &ctx->swzsurf);
	if (ret) {
		nv04_surface_2d_takedown(&ctx);
		return NULL;
	}

	switch (chan->device->chipset & 0xf0) {
	case 0x10:
	case 0x20:
		class = NV10_SCALED_IMAGE_FROM_MEMORY;
		break;
	case 0x30:
		class = NV30_SCALED_IMAGE_FROM_MEMORY;
		break;
	case 0x40:
	case 0x60:
		class = NV40_SCALED_IMAGE_FROM_MEMORY;
		break;
	default:
		class = NV04_SCALED_IMAGE_FROM_MEMORY;
		break;
	}

	ret = nouveau_grobj_alloc(chan, handle++, class, &ctx->sifm);
	if (ret) {
		nv04_surface_2d_takedown(&ctx);
		return NULL;
	}

	ctx->nvws = nvws;
	ctx->copy = nv04_surface_copy;
	ctx->fill = nv04_surface_fill;
	return ctx;
}
