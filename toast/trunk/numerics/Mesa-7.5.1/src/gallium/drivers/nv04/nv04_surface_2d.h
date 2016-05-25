#ifndef __NV04_SURFACE_2D_H__
#define __NV04_SURFACE_2D_H__

struct nv04_surface {
	struct pipe_surface base;
	unsigned pitch;
};

struct nv04_surface_2d {
	struct nouveau_winsys *nvws;
	struct nouveau_notifier *ntfy;
	struct nouveau_grobj *surf2d;
	struct nouveau_grobj *swzsurf;
	struct nouveau_grobj *m2mf;
	struct nouveau_grobj *rect;
	struct nouveau_grobj *blit;
	struct nouveau_grobj *sifm;

	struct pipe_buffer *(*buf)(struct pipe_surface *);

	void (*copy)(struct nv04_surface_2d *, struct pipe_surface *dst,
		     int dx, int dy, struct pipe_surface *src, int sx, int sy,
		     int w, int h);
	void (*fill)(struct nv04_surface_2d *, struct pipe_surface *dst,
		     int dx, int dy, int w, int h, unsigned value);
};

struct nv04_surface_2d *
nv04_surface_2d_init(struct nouveau_winsys *nvws);

void
nv04_surface_2d_takedown(struct nv04_surface_2d **);

#endif
