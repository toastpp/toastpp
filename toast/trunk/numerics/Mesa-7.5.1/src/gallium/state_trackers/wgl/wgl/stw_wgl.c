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

#include <windows.h>

#include "util/u_debug.h"
#include "shared/stw_public.h"
#include "stw_wgl.h"


WINGDIAPI BOOL APIENTRY
wglCopyContext(
   HGLRC hglrcSrc,
   HGLRC hglrcDst,
   UINT mask )
{
   return stw_copy_context( (UINT_PTR)hglrcSrc, 
                            (UINT_PTR)hglrcDst, 
                            mask );
}

WINGDIAPI HGLRC APIENTRY
wglCreateContext(
   HDC hdc )
{
   return wglCreateLayerContext(hdc, 0);
}

WINGDIAPI HGLRC APIENTRY
wglCreateLayerContext(
   HDC hdc,
   int iLayerPlane )
{
   return (HGLRC) stw_create_layer_context( hdc, iLayerPlane );
}

WINGDIAPI BOOL APIENTRY
wglDeleteContext(
   HGLRC hglrc )
{
   return stw_delete_context( (UINT_PTR)hglrc );
}


WINGDIAPI HGLRC APIENTRY
wglGetCurrentContext( VOID )
{
   return (HGLRC)stw_get_current_context();
}

WINGDIAPI HDC APIENTRY
wglGetCurrentDC( VOID )
{
   return stw_get_current_dc();
}

WINGDIAPI BOOL APIENTRY
wglMakeCurrent(
   HDC hdc,
   HGLRC hglrc )
{
   return stw_make_current( hdc, (UINT_PTR)hglrc );
}


WINGDIAPI BOOL APIENTRY
wglSwapBuffers(
   HDC hdc )
{
   return stw_swap_buffers( hdc );
}


WINGDIAPI BOOL APIENTRY
wglSwapLayerBuffers(
   HDC hdc,
   UINT fuPlanes )
{
   return stw_swap_layer_buffers( hdc, fuPlanes );
}

WINGDIAPI PROC APIENTRY
wglGetProcAddress(
    LPCSTR lpszProc )
{
   return stw_get_proc_address( lpszProc );
}


WINGDIAPI int APIENTRY
wglChoosePixelFormat(
   HDC hdc,
   CONST PIXELFORMATDESCRIPTOR *ppfd )
{
   if (ppfd->nSize != sizeof( PIXELFORMATDESCRIPTOR ) || ppfd->nVersion != 1)
      return 0;
   if (ppfd->iPixelType != PFD_TYPE_RGBA)
      return 0;
   if (!(ppfd->dwFlags & PFD_DRAW_TO_WINDOW))
      return 0;
   if (!(ppfd->dwFlags & PFD_SUPPORT_OPENGL))
      return 0;
   if (ppfd->dwFlags & PFD_DRAW_TO_BITMAP)
      return 0;
   if (ppfd->dwFlags & PFD_SUPPORT_GDI)
      return 0;
   if (!(ppfd->dwFlags & PFD_STEREO_DONTCARE) && (ppfd->dwFlags & PFD_STEREO))
      return 0;

   return stw_pixelformat_choose( hdc, ppfd );
}

WINGDIAPI int APIENTRY
wglDescribePixelFormat(
   HDC hdc,
   int iPixelFormat,
   UINT nBytes,
   LPPIXELFORMATDESCRIPTOR ppfd )
{
   return stw_pixelformat_describe( hdc, iPixelFormat, nBytes, ppfd );
}

WINGDIAPI int APIENTRY
wglGetPixelFormat(
   HDC hdc )
{
   return stw_pixelformat_get( hdc );
}

WINGDIAPI BOOL APIENTRY
wglSetPixelFormat(
   HDC hdc,
   int iPixelFormat,
   const PIXELFORMATDESCRIPTOR *ppfd )
{
   if (ppfd->nSize != sizeof( PIXELFORMATDESCRIPTOR ))
      return FALSE;

   return stw_pixelformat_set( hdc, iPixelFormat );
}


WINGDIAPI BOOL APIENTRY
wglUseFontBitmapsA(
   HDC hdc,
   DWORD first,
   DWORD count,
   DWORD listBase )
{
   (void) hdc;
   (void) first;
   (void) count;
   (void) listBase;

   assert( 0 );

   return FALSE;
}

WINGDIAPI BOOL APIENTRY
wglShareLists(
   HGLRC hglrc1,
   HGLRC hglrc2 )
{
   return stw_share_lists( (UINT_PTR)hglrc1, (UINT_PTR)hglrc2);;
}

WINGDIAPI BOOL APIENTRY
wglUseFontBitmapsW(
   HDC hdc,
   DWORD first,
   DWORD count,
   DWORD listBase )
{
   (void) hdc;
   (void) first;
   (void) count;
   (void) listBase;

   assert( 0 );

   return FALSE;
}

WINGDIAPI BOOL APIENTRY
wglUseFontOutlinesA(
   HDC hdc,
   DWORD first,
   DWORD count,
   DWORD listBase,
   FLOAT deviation,
   FLOAT extrusion,
   int format,
   LPGLYPHMETRICSFLOAT lpgmf )
{
   (void) hdc;
   (void) first;
   (void) count;
   (void) listBase;
   (void) deviation;
   (void) extrusion;
   (void) format;
   (void) lpgmf;

   assert( 0 );

   return FALSE;
}

WINGDIAPI BOOL APIENTRY
wglUseFontOutlinesW(
   HDC hdc,
   DWORD first,
   DWORD count,
   DWORD listBase,
   FLOAT deviation,
   FLOAT extrusion,
   int format,
   LPGLYPHMETRICSFLOAT lpgmf )
{
   (void) hdc;
   (void) first;
   (void) count;
   (void) listBase;
   (void) deviation;
   (void) extrusion;
   (void) format;
   (void) lpgmf;

   assert( 0 );

   return FALSE;
}

WINGDIAPI BOOL APIENTRY
wglDescribeLayerPlane(
   HDC hdc,
   int iPixelFormat,
   int iLayerPlane,
   UINT nBytes,
   LPLAYERPLANEDESCRIPTOR plpd )
{
   (void) hdc;
   (void) iPixelFormat;
   (void) iLayerPlane;
   (void) nBytes;
   (void) plpd;

   assert( 0 );

   return FALSE;
}

WINGDIAPI int APIENTRY
wglSetLayerPaletteEntries(
   HDC hdc,
   int iLayerPlane,
   int iStart,
   int cEntries,
   CONST COLORREF *pcr )
{
   (void) hdc;
   (void) iLayerPlane;
   (void) iStart;
   (void) cEntries;
   (void) pcr;

   assert( 0 );

   return 0;
}

WINGDIAPI int APIENTRY
wglGetLayerPaletteEntries(
   HDC hdc,
   int iLayerPlane,
   int iStart,
   int cEntries,
   COLORREF *pcr )
{
   (void) hdc;
   (void) iLayerPlane;
   (void) iStart;
   (void) cEntries;
   (void) pcr;

   assert( 0 );

   return 0;
}

WINGDIAPI BOOL APIENTRY
wglRealizeLayerPalette(
   HDC hdc,
   int iLayerPlane,
   BOOL bRealize )
{
   (void) hdc;
   (void) iLayerPlane;
   (void) bRealize;

   assert( 0 );

   return FALSE;
}
