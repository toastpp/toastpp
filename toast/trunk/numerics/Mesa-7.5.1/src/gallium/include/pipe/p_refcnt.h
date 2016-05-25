/**************************************************************************
 * 
 * Copyright 2009 VMware, Inc.
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

#ifndef P_REFCNT_H
#define P_REFCNT_H


#include "p_defines.h"
#include "p_atomic.h"


#ifdef __cplusplus
extern "C" {
#endif


struct pipe_reference
{
   struct pipe_atomic count;
};


static INLINE void
pipe_reference_init(struct pipe_reference *reference, unsigned count)
{
   p_atomic_set(&reference->count, count);
}


static INLINE bool
pipe_is_referenced(struct pipe_reference *reference)
{
   return p_atomic_read(&reference->count) != 0;
}


/**
 * Set 'ptr' to point to 'reference' and update reference counting.
 * The old thing pointed to, if any, will be unreferenced first.
 * 'reference' may be NULL.
 */
static INLINE bool
pipe_reference(struct pipe_reference **ptr, struct pipe_reference *reference)
{
   bool destroy = FALSE;

   if(*ptr != reference) {
      /* bump the reference.count first */
      if (reference) {
         assert(pipe_is_referenced(reference));
         p_atomic_inc(&reference->count);
      }
   
      if (*ptr) {
         assert(pipe_is_referenced(*ptr));
         if (p_atomic_dec_zero(&(*ptr)->count)) {
            destroy = TRUE;
         }
      }
   
      *ptr = reference;
   }

   return destroy;
}


#ifdef __cplusplus
}
#endif

#endif /* P_REFCNT_H */
