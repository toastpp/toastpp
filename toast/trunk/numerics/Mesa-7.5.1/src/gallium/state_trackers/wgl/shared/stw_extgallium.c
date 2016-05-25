/**************************************************************************
 *
 * Copyright 2008 Tungsten Graphics, Inc., Cedar Park, Texas.
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


#include "pipe/p_screen.h"
#include "stw_public.h"
#include "stw_device.h"
#include "stw_winsys.h"

#ifdef DEBUG
#include "trace/tr_screen.h"
#include "trace/tr_context.h"
#endif


struct pipe_screen * APIENTRY
wglGetGalliumScreenMESA(void)
{
   return stw_dev ? stw_dev->screen : NULL;
}


/* XXX: Unify with stw_create_layer_context */
struct pipe_context * APIENTRY
wglCreateGalliumContextMESA(void)
{
   struct pipe_screen *screen = NULL;
   struct pipe_context *pipe = NULL;

   if(!stw_dev)
      return NULL;

   screen = stw_dev->screen;

#ifdef DEBUG
   /* Unwrap screen */
   if(stw_dev->trace_running)
      screen = trace_screen(screen)->screen;
#endif

   pipe = stw_dev->stw_winsys->create_context( screen );
   if (pipe == NULL)
      goto no_pipe;

#ifdef DEBUG
   /* Wrap context */
   if(stw_dev->trace_running)
      pipe = trace_context_create(stw_dev->screen, pipe);
#endif

   return pipe;

no_pipe:
   return NULL;
}
