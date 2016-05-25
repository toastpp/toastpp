/**************************************************************************
 * 
 * Copyright 2007 Tungsten Graphics, Inc., Cedar Park, Texas.
 * All Rights Reserved.
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sub license, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 * 
 * The above copyright notice and this permission notice (including the
 * next paragraph) shall be included in all copies or substantial portions
 * of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NON-INFRINGEMENT.
 * IN NO EVENT SHALL TUNGSTEN GRAPHICS AND/OR ITS SUPPLIERS BE LIABLE FOR
 * ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 * 
 **************************************************************************/

#ifndef P_INLINES_H
#define P_INLINES_H

#include "p_context.h"
#include "p_defines.h"
#include "p_screen.h"


#ifdef __cplusplus
extern "C" {
#endif


/**
 * Convenience wrappers for screen buffer functions.
 */

static INLINE struct pipe_buffer *
pipe_buffer_create( struct pipe_screen *screen,
                    unsigned alignment, unsigned usage, unsigned size )
{
   return screen->buffer_create(screen, alignment, usage, size);
}

static INLINE struct pipe_buffer *
pipe_user_buffer_create( struct pipe_screen *screen, void *ptr, unsigned size )
{
   return screen->user_buffer_create(screen, ptr, size);
}

static INLINE void *
pipe_buffer_map(struct pipe_screen *screen,
                struct pipe_buffer *buf,
                unsigned usage)
{
   if(screen->buffer_map_range) {
      unsigned offset = 0;
      unsigned length = buf->size;

      /* XXX: Actually we should be using/detecting DISCARD
       * instead of assuming that WRITE implies discard */
      if((usage & PIPE_BUFFER_USAGE_CPU_WRITE) &&
         !(usage & PIPE_BUFFER_USAGE_DISCARD))
         usage |= PIPE_BUFFER_USAGE_CPU_READ;

      return screen->buffer_map_range(screen, buf, offset, length, usage);
   }
   else
      return screen->buffer_map(screen, buf, usage);
}

static INLINE void
pipe_buffer_unmap(struct pipe_screen *screen,
                  struct pipe_buffer *buf)
{
   screen->buffer_unmap(screen, buf);
}

static INLINE void *
pipe_buffer_map_range(struct pipe_screen *screen,
                struct pipe_buffer *buf,
                unsigned offset,
                unsigned length,
                unsigned usage)
{
   assert(offset < buf->size);
   assert(offset + length <= buf->size);
   assert(length);
   if(screen->buffer_map_range)
      return screen->buffer_map_range(screen, buf, offset, length, usage);
   else
      return screen->buffer_map(screen, buf, usage);
}

static INLINE void
pipe_buffer_flush_mapped_range(struct pipe_screen *screen,
                               struct pipe_buffer *buf,
                               unsigned offset,
                               unsigned length)
{
   assert(offset < buf->size);
   assert(offset + length <= buf->size);
   assert(length);
   if(screen->buffer_flush_mapped_range)
      screen->buffer_flush_mapped_range(screen, buf, offset, length);
}

static INLINE void
pipe_buffer_write(struct pipe_screen *screen,
                  struct pipe_buffer *buf,
                  unsigned offset, unsigned size,
                  const void *data)
{
   uint8_t *map;
   
   assert(offset < buf->size);
   assert(offset + size <= buf->size);
   assert(size);

   map = pipe_buffer_map_range(screen, buf, offset, size, 
                               PIPE_BUFFER_USAGE_CPU_WRITE | 
                               PIPE_BUFFER_USAGE_FLUSH_EXPLICIT);
   assert(map);
   if(map) {
      memcpy(map + offset, data, size);
      pipe_buffer_flush_mapped_range(screen, buf, offset, size);
      pipe_buffer_unmap(screen, buf);
   }
}

static INLINE void
pipe_buffer_read(struct pipe_screen *screen,
                 struct pipe_buffer *buf,
                 unsigned offset, unsigned size,
                 void *data)
{
   uint8_t *map;
   
   assert(offset < buf->size);
   assert(offset + size <= buf->size);
   assert(size);

   map = pipe_buffer_map_range(screen, buf, offset, size, PIPE_BUFFER_USAGE_CPU_READ);
   assert(map);
   if(map) {
      memcpy(data, map + offset, size);
      pipe_buffer_unmap(screen, buf);
   }
}


#ifdef __cplusplus
}
#endif

#endif /* P_INLINES_H */
