#ifndef __SH_H
#define __SH_H

typedef struct { double x, y; } vec_t, *vec;
typedef struct { int len, alloc; vec v; } poly_t, *poly;

poly poly_clip(poly sub, poly clip);
void poly_free(poly p);

#endif // !__SH_H
