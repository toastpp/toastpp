/*
 * Copyright 2008 Corbin Simpson <MostAwesomeDude@gmail.com>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * on the rights to use, copy, modify, merge, publish, distribute, sub
 * license, and/or sell copies of the Software, and to permit persons to whom
 * the Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice (including the next
 * paragraph) shall be included in all copies or substantial portions of the
 * Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NON-INFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHOR(S) AND/OR THEIR SUPPLIERS BE LIABLE FOR ANY CLAIM,
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
 * USE OR OTHER DEALINGS IN THE SOFTWARE. */

#ifndef R300_EMIT_H
#define R300_EMIT_H

#include "util/u_math.h"

#include "r300_context.h"
#include "r300_cs.h"
#include "r300_screen.h"
#include "r300_state_inlines.h"

void r300_emit_blend_state(struct r300_context* r300,
                           struct r300_blend_state* blend);

void r300_emit_blend_color_state(struct r300_context* r300,
                                 struct r300_blend_color_state* bc);

void r300_emit_dsa_state(struct r300_context* r300,
                         struct r300_dsa_state* dsa);

void r300_emit_fragment_shader(struct r300_context* r300,
                               struct r300_fragment_shader* fs);

void r500_emit_fragment_shader(struct r300_context* r300,
                               struct r500_fragment_shader* fs);

void r300_emit_fb_state(struct r300_context* r300,
                        struct pipe_framebuffer_state* fb);

void r300_emit_rs_state(struct r300_context* r300, struct r300_rs_state* rs);

void r300_emit_rs_block_state(struct r300_context* r300,
                              struct r300_rs_block* rs);

void r300_emit_sampler(struct r300_context* r300,
                       struct r300_sampler_state* sampler, unsigned offset);

void r300_emit_scissor_state(struct r300_context* r300,
                             struct r300_scissor_state* scissor);

void r300_emit_texture(struct r300_context* r300,
                       struct r300_texture* tex, unsigned offset);

void r300_emit_vertex_format_state(struct r300_context* r300);

void r300_emit_vertex_shader(struct r300_context* r300,
                             struct r300_vertex_shader* vs);

void r300_emit_viewport_state(struct r300_context* r300,
                              struct r300_viewport_state* viewport);

void r300_flush_textures(struct r300_context* r300);

/* Emit all dirty state. */
void r300_emit_dirty_state(struct r300_context* r300);

#endif /* R300_EMIT_H */
