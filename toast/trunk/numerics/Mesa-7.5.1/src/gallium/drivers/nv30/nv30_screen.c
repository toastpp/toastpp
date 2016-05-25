#include "pipe/p_screen.h"
#include "util/u_simple_screen.h"

#include "nv30_context.h"
#include "nv30_screen.h"

#define NV30TCL_CHIPSET_3X_MASK 0x00000003
#define NV34TCL_CHIPSET_3X_MASK 0x00000010
#define NV35TCL_CHIPSET_3X_MASK 0x000001e0

static const char *
nv30_screen_get_name(struct pipe_screen *pscreen)
{
	struct nv30_screen *screen = nv30_screen(pscreen);
	struct nouveau_device *dev = screen->nvws->channel->device;
	static char buffer[128];

	snprintf(buffer, sizeof(buffer), "NV%02X", dev->chipset);
	return buffer;
}

static const char *
nv30_screen_get_vendor(struct pipe_screen *pscreen)
{
	return "nouveau";
}

static int
nv30_screen_get_param(struct pipe_screen *pscreen, int param)
{
	switch (param) {
	case PIPE_CAP_MAX_TEXTURE_IMAGE_UNITS:
		return 16;
	case PIPE_CAP_NPOT_TEXTURES:
		return 0;
	case PIPE_CAP_TWO_SIDED_STENCIL:
		return 1;
	case PIPE_CAP_GLSL:
		return 0;
	case PIPE_CAP_S3TC:
		return 0;
	case PIPE_CAP_ANISOTROPIC_FILTER:
		return 1;
	case PIPE_CAP_POINT_SPRITE:
		return 1;
	case PIPE_CAP_MAX_RENDER_TARGETS:
		return 2;
	case PIPE_CAP_OCCLUSION_QUERY:
		return 1;
	case PIPE_CAP_TEXTURE_SHADOW_MAP:
		return 1;
	case PIPE_CAP_MAX_TEXTURE_2D_LEVELS:
		return 13;
	case PIPE_CAP_MAX_TEXTURE_3D_LEVELS:
		return 10;
	case PIPE_CAP_MAX_TEXTURE_CUBE_LEVELS:
		return 13;
	case PIPE_CAP_TEXTURE_MIRROR_CLAMP:
		return 0;
	case PIPE_CAP_TEXTURE_MIRROR_REPEAT:
		return 1;
	case PIPE_CAP_MAX_VERTEX_TEXTURE_UNITS:
		return 0;
	case NOUVEAU_CAP_HW_VTXBUF:
	case NOUVEAU_CAP_HW_IDXBUF:
		return 1;
	default:
		NOUVEAU_ERR("Unknown PIPE_CAP %d\n", param);
		return 0;
	}
}

static float
nv30_screen_get_paramf(struct pipe_screen *pscreen, int param)
{
	switch (param) {
	case PIPE_CAP_MAX_LINE_WIDTH:
	case PIPE_CAP_MAX_LINE_WIDTH_AA:
		return 10.0;
	case PIPE_CAP_MAX_POINT_WIDTH:
	case PIPE_CAP_MAX_POINT_WIDTH_AA:
		return 64.0;
	case PIPE_CAP_MAX_TEXTURE_ANISOTROPY:
		return 8.0;
	case PIPE_CAP_MAX_TEXTURE_LOD_BIAS:
		return 4.0;
	default:
		NOUVEAU_ERR("Unknown PIPE_CAP %d\n", param);
		return 0.0;
	}
}

static boolean
nv30_screen_surface_format_supported(struct pipe_screen *pscreen,
				     enum pipe_format format,
				     enum pipe_texture_target target,
				     unsigned tex_usage, unsigned geom_flags)
{
	if (tex_usage & PIPE_TEXTURE_USAGE_RENDER_TARGET) {
		switch (format) {
		case PIPE_FORMAT_A8R8G8B8_UNORM:
		case PIPE_FORMAT_R5G6B5_UNORM:
		case PIPE_FORMAT_Z24S8_UNORM:
		case PIPE_FORMAT_Z16_UNORM:
			return TRUE;
		default:
			break;
		}
	} else {
		switch (format) {
		case PIPE_FORMAT_A8R8G8B8_UNORM:
		case PIPE_FORMAT_A1R5G5B5_UNORM:
		case PIPE_FORMAT_A4R4G4B4_UNORM:
		case PIPE_FORMAT_R5G6B5_UNORM:
		case PIPE_FORMAT_L8_UNORM:
		case PIPE_FORMAT_A8_UNORM:
		case PIPE_FORMAT_I8_UNORM:
		case PIPE_FORMAT_A8L8_UNORM:
		case PIPE_FORMAT_Z16_UNORM:
		case PIPE_FORMAT_Z24S8_UNORM:
			return TRUE;
		default:
			break;
		}
	}

	return FALSE;
}

static struct pipe_buffer *
nv30_surface_buffer(struct pipe_surface *surf)
{
	struct nv30_miptree *mt = (struct nv30_miptree *)surf->texture;

	return mt->buffer;
}

static void
nv30_screen_destroy(struct pipe_screen *pscreen)
{
	struct nv30_screen *screen = nv30_screen(pscreen);
	struct nouveau_winsys *nvws = screen->nvws;

	nvws->res_free(&screen->vp_exec_heap);
	nvws->res_free(&screen->vp_data_heap);
	nvws->res_free(&screen->query_heap);
	nvws->notifier_free(&screen->query);
	nvws->notifier_free(&screen->sync);
	nvws->grobj_free(&screen->rankine);

	FREE(pscreen);
}

struct pipe_screen *
nv30_screen_create(struct pipe_winsys *ws, struct nouveau_winsys *nvws)
{
	struct nv30_screen *screen = CALLOC_STRUCT(nv30_screen);
	struct nouveau_stateobj *so;
	unsigned rankine_class = 0;
	unsigned chipset = nvws->channel->device->chipset;
	int ret, i;

	if (!screen)
		return NULL;
	screen->nvws = nvws;

	/* 2D engine setup */
	screen->eng2d = nv04_surface_2d_init(nvws);
	screen->eng2d->buf = nv30_surface_buffer;

	/* 3D object */
	switch (chipset & 0xf0) {
	case 0x30:
		if (NV30TCL_CHIPSET_3X_MASK & (1 << (chipset & 0x0f)))
			rankine_class = 0x0397;
		else
		if (NV34TCL_CHIPSET_3X_MASK & (1 << (chipset & 0x0f)))
			rankine_class = 0x0697;
		else
		if (NV35TCL_CHIPSET_3X_MASK & (1 << (chipset & 0x0f)))
			rankine_class = 0x0497;
		break;
	default:
		break;
	}

	if (!rankine_class) {
		NOUVEAU_ERR("Unknown nv3x chipset: nv%02x\n", chipset);
		return NULL;
	}

	ret = nvws->grobj_alloc(nvws, rankine_class, &screen->rankine);
	if (ret) {
		NOUVEAU_ERR("Error creating 3D object: %d\n", ret);
		return FALSE;
	}

	/* Notifier for sync purposes */
	ret = nvws->notifier_alloc(nvws, 1, &screen->sync);
	if (ret) {
		NOUVEAU_ERR("Error creating notifier object: %d\n", ret);
		nv30_screen_destroy(&screen->pipe);
		return NULL;
	}

	/* Query objects */
	ret = nvws->notifier_alloc(nvws, 32, &screen->query);
	if (ret) {
		NOUVEAU_ERR("Error initialising query objects: %d\n", ret);
		nv30_screen_destroy(&screen->pipe);
		return NULL;
	}

	ret = nvws->res_init(&screen->query_heap, 0, 32);
	if (ret) {
		NOUVEAU_ERR("Error initialising query object heap: %d\n", ret);
		nv30_screen_destroy(&screen->pipe);
		return NULL;
	}

	/* Vtxprog resources */
	if (nvws->res_init(&screen->vp_exec_heap, 0, 256) ||
	    nvws->res_init(&screen->vp_data_heap, 0, 256)) {
		nv30_screen_destroy(&screen->pipe);
		return NULL;
	}

	/* Static rankine initialisation */
	so = so_new(128, 0);
	so_method(so, screen->rankine, NV34TCL_DMA_NOTIFY, 1);
	so_data  (so, screen->sync->handle);
	so_method(so, screen->rankine, NV34TCL_DMA_TEXTURE0, 2);
	so_data  (so, nvws->channel->vram->handle);
	so_data  (so, nvws->channel->gart->handle);
	so_method(so, screen->rankine, NV34TCL_DMA_COLOR1, 1);
	so_data  (so, nvws->channel->vram->handle);
	so_method(so, screen->rankine, NV34TCL_DMA_COLOR0, 2);
	so_data  (so, nvws->channel->vram->handle);
	so_data  (so, nvws->channel->vram->handle);
	so_method(so, screen->rankine, NV34TCL_DMA_VTXBUF0, 2);
	so_data  (so, nvws->channel->vram->handle);
	so_data  (so, nvws->channel->gart->handle);
/*	so_method(so, screen->rankine, NV34TCL_DMA_FENCE, 2);
	so_data  (so, 0);
	so_data  (so, screen->query->handle);*/
	so_method(so, screen->rankine, NV34TCL_DMA_IN_MEMORY7, 1);
	so_data  (so, nvws->channel->vram->handle);
	so_method(so, screen->rankine, NV34TCL_DMA_IN_MEMORY8, 1);
	so_data  (so, nvws->channel->vram->handle);

	for (i=1; i<8; i++) {
		so_method(so, screen->rankine, NV34TCL_VIEWPORT_CLIP_HORIZ(i), 1);
		so_data  (so, 0);
		so_method(so, screen->rankine, NV34TCL_VIEWPORT_CLIP_VERT(i), 1);
		so_data  (so, 0);
	}

	so_method(so, screen->rankine, 0x220, 1);
	so_data  (so, 1);

	so_method(so, screen->rankine, 0x03b0, 1);
	so_data  (so, 0x00100000);
	so_method(so, screen->rankine, 0x1454, 1);
	so_data  (so, 0);
	so_method(so, screen->rankine, 0x1d80, 1);
	so_data  (so, 3);
	so_method(so, screen->rankine, 0x1450, 1);
	so_data  (so, 0x00030004);

	/* NEW */
	so_method(so, screen->rankine, 0x1e98, 1);
	so_data  (so, 0);
	so_method(so, screen->rankine, 0x17e0, 3);
	so_data  (so, fui(0.0));
	so_data  (so, fui(0.0));
	so_data  (so, fui(1.0));
	so_method(so, screen->rankine, 0x1f80, 16);
	for (i=0; i<16; i++) {
		so_data  (so, (i==8) ? 0x0000ffff : 0);
	}

	so_method(so, screen->rankine, 0x120, 3);
	so_data  (so, 0);
	so_data  (so, 1);
	so_data  (so, 2);

	so_method(so, screen->rankine, 0x1d88, 1);
	so_data  (so, 0x00001200);

	so_method(so, screen->rankine, NV34TCL_RC_ENABLE, 1);
	so_data  (so, 0);

	so_method(so, screen->rankine, NV34TCL_DEPTH_RANGE_NEAR, 2);
	so_data  (so, fui(0.0));
	so_data  (so, fui(1.0));

	so_method(so, screen->rankine, NV34TCL_MULTISAMPLE_CONTROL, 1);
	so_data  (so, 0xffff0000);

	/* enables use of vp rather than fixed-function somehow */
	so_method(so, screen->rankine, 0x1e94, 1);
	so_data  (so, 0x13);

	so_emit(nvws, so);
	so_ref(NULL, &so);
	nvws->push_flush(nvws, 0, NULL);

	screen->pipe.winsys = ws;
	screen->pipe.destroy = nv30_screen_destroy;

	screen->pipe.get_name = nv30_screen_get_name;
	screen->pipe.get_vendor = nv30_screen_get_vendor;
	screen->pipe.get_param = nv30_screen_get_param;
	screen->pipe.get_paramf = nv30_screen_get_paramf;

	screen->pipe.is_format_supported = nv30_screen_surface_format_supported;

	nv30_screen_init_miptree_functions(&screen->pipe);
	nv30_screen_init_transfer_functions(&screen->pipe);
	u_simple_screen_init(&screen->pipe);

	return &screen->pipe;
}
