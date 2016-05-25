/**************************************************************************
 * 
 * Copyright 2007 Tungsten Graphics, Inc., Bismarck, ND., USA
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NON-INFRINGEMENT. IN NO EVENT SHALL
 * THE COPYRIGHT HOLDERS, AUTHORS AND/OR ITS SUPPLIERS BE LIABLE FOR ANY CLAIM,
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE 
 * USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 * The above copyright notice and this permission notice (including the
 * next paragraph) shall be included in all copies or substantial portions
 * of the Software.
 * 
 * 
 **************************************************************************/

/*
 * Authors:
 *   Keith Whitwell
 *   Brian Paul
 */


#include "xm_api.h"

#undef ASSERT
#undef Elements

#include "pipe/internal/p_winsys_screen.h"
#include "pipe/p_format.h"
#include "pipe/p_context.h"
#include "pipe/p_inlines.h"
#include "util/u_math.h"
#include "util/u_memory.h"
#include "softpipe/sp_winsys.h"
#include "softpipe/sp_texture.h"

#include "xlib.h"

/**
 * Subclass of pipe_buffer for Xlib winsys.
 * Low-level OS/window system memory buffer
 */
struct xm_buffer
{
   struct pipe_buffer base;
   boolean userBuffer;  /** Is this a user-space buffer? */
   void *data;
   void *mapped;

   XImage *tempImage;
#ifdef USE_XSHM
   int shm;
   XShmSegmentInfo shminfo;
#endif
};


/**
 * Subclass of pipe_winsys for Xlib winsys
 */
struct xmesa_pipe_winsys
{
   struct pipe_winsys base;
/*   struct xmesa_visual *xm_visual; */
#ifdef USE_XSHM
   int shm;
#endif
};



/** Cast wrapper */
static INLINE struct xm_buffer *
xm_buffer( struct pipe_buffer *buf )
{
   return (struct xm_buffer *)buf;
}


/**
 * X Shared Memory Image extension code
 */
#ifdef USE_XSHM
#define XSHM_ENABLED(b) ((b)->shm)
#else
#define XSHM_ENABLED(b) 0
#endif

#ifdef USE_XSHM

static volatile int mesaXErrorFlag = 0;

/**
 * Catches potential Xlib errors.
 */
static int
mesaHandleXError(Display *dpy, XErrorEvent *event)
{
   (void) dpy;
   (void) event;
   mesaXErrorFlag = 1;
   return 0;
}


static GLboolean alloc_shm(struct xm_buffer *buf, unsigned size)
{
   XShmSegmentInfo *const shminfo = & buf->shminfo;

   shminfo->shmid = shmget(IPC_PRIVATE, size, IPC_CREAT|0777);
   if (shminfo->shmid < 0) {
      return GL_FALSE;
   }

   shminfo->shmaddr = (char *) shmat(shminfo->shmid, 0, 0);
   if (shminfo->shmaddr == (char *) -1) {
      shmctl(shminfo->shmid, IPC_RMID, 0);
      return GL_FALSE;
   }

   shminfo->readOnly = False;
   return GL_TRUE;
}


/**
 * Allocate a shared memory XImage back buffer for the given XMesaBuffer.
 */
static void
alloc_shm_ximage(struct xm_buffer *b, struct xmesa_buffer *xmb,
                 unsigned width, unsigned height)
{
   /*
    * We have to do a _lot_ of error checking here to be sure we can
    * really use the XSHM extension.  It seems different servers trigger
    * errors at different points if the extension won't work.  Therefore
    * we have to be very careful...
    */
   int (*old_handler)(Display *, XErrorEvent *);

   b->tempImage = XShmCreateImage(xmb->xm_visual->display,
                                  xmb->xm_visual->visinfo->visual,
                                  xmb->xm_visual->visinfo->depth,
                                  ZPixmap,
                                  NULL,
                                  &b->shminfo,
                                  width, height);
   if (b->tempImage == NULL) {
      b->shm = 0;
      return;
   }


   mesaXErrorFlag = 0;
   old_handler = XSetErrorHandler(mesaHandleXError);
   /* This may trigger the X protocol error we're ready to catch: */
   XShmAttach(xmb->xm_visual->display, &b->shminfo);
   XSync(xmb->xm_visual->display, False);

   if (mesaXErrorFlag) {
      /* we are on a remote display, this error is normal, don't print it */
      XFlush(xmb->xm_visual->display);
      mesaXErrorFlag = 0;
      XDestroyImage(b->tempImage);
      b->tempImage = NULL;
      b->shm = 0;
      (void) XSetErrorHandler(old_handler);
      return;
   }

   b->shm = 1;
}

#endif /* USE_XSHM */



/* Most callbacks map direcly onto dri_bufmgr operations:
 */
static void *
xm_buffer_map(struct pipe_winsys *pws, struct pipe_buffer *buf,
              unsigned flags)
{
   struct xm_buffer *xm_buf = xm_buffer(buf);
   xm_buf->mapped = xm_buf->data;
   return xm_buf->mapped;
}

static void
xm_buffer_unmap(struct pipe_winsys *pws, struct pipe_buffer *buf)
{
   struct xm_buffer *xm_buf = xm_buffer(buf);
   xm_buf->mapped = NULL;
}

static void
xm_buffer_destroy(struct pipe_buffer *buf)
{
   struct xm_buffer *oldBuf = xm_buffer(buf);

   if (oldBuf->data) {
#ifdef USE_XSHM
      if (oldBuf->shminfo.shmid >= 0) {
         shmdt(oldBuf->shminfo.shmaddr);
         shmctl(oldBuf->shminfo.shmid, IPC_RMID, 0);
         
         oldBuf->shminfo.shmid = -1;
         oldBuf->shminfo.shmaddr = (char *) -1;
      }
      else
#endif
      {
         if (!oldBuf->userBuffer) {
            align_free(oldBuf->data);
         }
      }

      oldBuf->data = NULL;
   }

   free(oldBuf);
}


/**
 * Display/copy the image in the surface into the X window specified
 * by the XMesaBuffer.
 */
static void
xlib_softpipe_display_surface(struct xmesa_buffer *b,
                              struct pipe_surface *surf)
{
   XImage *ximage;
   struct softpipe_texture *spt = softpipe_texture(surf->texture);
   struct xm_buffer *xm_buf = xm_buffer(spt->buffer);
   static boolean no_swap = 0;
   static boolean firsttime = 1;

   if (firsttime) {
      no_swap = getenv("SP_NO_RAST") != NULL;
      firsttime = 0;
   }

   if (no_swap)
      return;

#ifdef USE_XSHM
   if (XSHM_ENABLED(xm_buf) && (xm_buf->tempImage == NULL)) {
      assert(surf->texture->block.width == 1);
      assert(surf->texture->block.height == 1);
      alloc_shm_ximage(xm_buf, b, spt->stride[surf->level] /
                       surf->texture->block.size, surf->height);
   }
#endif

   ximage = (XSHM_ENABLED(xm_buf)) ? xm_buf->tempImage : b->tempImage;
   ximage->data = xm_buf->data;

   /* display image in Window */
#ifdef USE_XSHM
   if (XSHM_ENABLED(xm_buf)) {
      XShmPutImage(b->xm_visual->display, b->drawable, b->gc,
                   ximage, 0, 0, 0, 0, surf->width, surf->height, False);
   } else
#endif
   {
      /* check that the XImage has been previously initialized */
      assert(ximage->format);
      assert(ximage->bitmap_unit);

      /* update XImage's fields */
      ximage->width = surf->width;
      ximage->height = surf->height;
      ximage->bytes_per_line = spt->stride[surf->level];

      XPutImage(b->xm_visual->display, b->drawable, b->gc,
                ximage, 0, 0, 0, 0, surf->width, surf->height);
   }
}


static void
xm_flush_frontbuffer(struct pipe_winsys *pws,
                     struct pipe_surface *surf,
                     void *context_private)
{
   /*
    * The front color buffer is actually just another XImage buffer.
    * This function copies that XImage to the actual X Window.
    */
   XMesaContext xmctx = (XMesaContext) context_private;
   xlib_softpipe_display_surface(xmctx->xm_buffer, surf);
   xmesa_check_and_update_buffer_size(xmctx, xmctx->xm_buffer);
}



static const char *
xm_get_name(struct pipe_winsys *pws)
{
   return "Xlib";
}


static struct pipe_buffer *
xm_buffer_create(struct pipe_winsys *pws, 
                 unsigned alignment, 
                 unsigned usage,
                 unsigned size)
{
   struct xm_buffer *buffer = CALLOC_STRUCT(xm_buffer);
#ifdef USE_XSHM
   struct xmesa_pipe_winsys *xpws = (struct xmesa_pipe_winsys *) pws;

   buffer->shminfo.shmid = -1;
   buffer->shminfo.shmaddr = (char *) -1;

   if (xpws->shm && (usage & PIPE_BUFFER_USAGE_PIXEL) != 0) {
      buffer->shm = xpws->shm;

      if (alloc_shm(buffer, size)) {
         buffer->data = buffer->shminfo.shmaddr;
         buffer->shm = 1;
      }
   }
#endif

   pipe_reference_init(&buffer->base.reference, 1);
   buffer->base.alignment = alignment;
   buffer->base.usage = usage;
   buffer->base.size = size;

   if (buffer->data == NULL) {
      /* align to 16-byte multiple for Cell */
      buffer->data = align_malloc(size, max(alignment, 16));
   }

   return &buffer->base;
}


/**
 * Create buffer which wraps user-space data.
 */
static struct pipe_buffer *
xm_user_buffer_create(struct pipe_winsys *pws, void *ptr, unsigned bytes)
{
   struct xm_buffer *buffer = CALLOC_STRUCT(xm_buffer);
   pipe_reference_init(&buffer->base.reference, 1);
   buffer->base.size = bytes;
   buffer->userBuffer = TRUE;
   buffer->data = ptr;
#ifdef USE_XSHM
   buffer->shm = 0;
#endif

   return &buffer->base;
}


static struct pipe_buffer *
xm_surface_buffer_create(struct pipe_winsys *winsys,
                         unsigned width, unsigned height,
                         enum pipe_format format,
                         unsigned usage,
                         unsigned *stride)
{
   const unsigned alignment = 64;
   struct pipe_format_block block;
   unsigned nblocksx, nblocksy;

   pf_get_block(format, &block);
   nblocksx = pf_get_nblocksx(&block, width);
   nblocksy = pf_get_nblocksy(&block, height);
   *stride = align(nblocksx * block.size, alignment);

   return winsys->buffer_create(winsys, alignment,
                                usage,
                                *stride * nblocksy);
}


/*
 * Fence functions - basically nothing to do, as we don't create any actual
 * fence objects.
 */

static void
xm_fence_reference(struct pipe_winsys *sws, struct pipe_fence_handle **ptr,
                   struct pipe_fence_handle *fence)
{
}


static int
xm_fence_signalled(struct pipe_winsys *sws, struct pipe_fence_handle *fence,
                   unsigned flag)
{
   return 0;
}


static int
xm_fence_finish(struct pipe_winsys *sws, struct pipe_fence_handle *fence,
                unsigned flag)
{
   return 0;
}



static struct pipe_winsys *
xlib_create_softpipe_winsys( void )
{
   static struct xmesa_pipe_winsys *ws = NULL;

   if (!ws) {
      ws = CALLOC_STRUCT(xmesa_pipe_winsys);

      /* Fill in this struct with callbacks that pipe will need to
       * communicate with the window system, buffer manager, etc. 
       */
      ws->base.buffer_create = xm_buffer_create;
      ws->base.user_buffer_create = xm_user_buffer_create;
      ws->base.buffer_map = xm_buffer_map;
      ws->base.buffer_unmap = xm_buffer_unmap;
      ws->base.buffer_destroy = xm_buffer_destroy;

      ws->base.surface_buffer_create = xm_surface_buffer_create;

      ws->base.fence_reference = xm_fence_reference;
      ws->base.fence_signalled = xm_fence_signalled;
      ws->base.fence_finish = xm_fence_finish;

      ws->base.flush_frontbuffer = xm_flush_frontbuffer;
      ws->base.get_name = xm_get_name;
   }

   return &ws->base;
}


static struct pipe_screen *
xlib_create_softpipe_screen( void )
{
   struct pipe_winsys *winsys;
   struct pipe_screen *screen;

   winsys = xlib_create_softpipe_winsys();
   if (winsys == NULL)
      return NULL;

   screen = softpipe_create_screen(winsys);
   if (screen == NULL)
      goto fail;

   return screen;

fail:
   if (winsys)
      winsys->destroy( winsys );

   return NULL;
}


static struct pipe_context *
xlib_create_softpipe_context( struct pipe_screen *screen,
                              void *context_private )
{
   struct pipe_context *pipe;
   
   pipe = softpipe_create(screen);
   if (pipe == NULL)
      goto fail;

   pipe->priv = context_private;
   return pipe;

fail:
   /* Free stuff here */
   return NULL;
}

struct xm_driver xlib_softpipe_driver = 
{
   .create_pipe_screen = xlib_create_softpipe_screen,
   .create_pipe_context = xlib_create_softpipe_context,
   .display_surface = xlib_softpipe_display_surface
};



