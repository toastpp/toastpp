/* 
 * Copyright © 2009 Corbin Simpson
 * All Rights Reserved.
 * 
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sub license, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NON-INFRINGEMENT. IN NO EVENT SHALL THE COPYRIGHT HOLDERS, AUTHORS
 * AND/OR ITS SUPPLIERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE 
 * USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 * The above copyright notice and this permission notice (including the
 * next paragraph) shall be included in all copies or substantial portions
 * of the Software.
 */
/*
 * Authors:
 *      Corbin Simpson <MostAwesomeDude@gmail.com>
 */
#ifndef RADEON_DRM_H
#define RADEON_DRM_H

#include "pipe/p_screen.h"

#include "util/u_memory.h"

#include "state_tracker/drm_api.h"

#include "radeon_buffer.h"
#include "radeon_r300.h"
#include "radeon_winsys_softpipe.h"

struct pipe_screen* radeon_create_screen(int drmFB,
					 struct drm_create_screen_arg *arg);

struct pipe_context* radeon_create_context(struct pipe_screen* screen);

boolean radeon_buffer_from_texture(struct pipe_texture* texture,
                                   struct pipe_buffer** buffer,
                                   unsigned* stride);

struct pipe_buffer* radeon_buffer_from_handle(struct pipe_screen* screen,
                                              const char* name,
                                              unsigned handle);

boolean radeon_handle_from_buffer(struct pipe_screen* screen,
                                  struct pipe_buffer* buffer,
                                  unsigned* handle);

boolean radeon_global_handle_from_buffer(struct pipe_screen* screen,
                                         struct pipe_buffer* buffer,
                                         unsigned* handle);

#endif
