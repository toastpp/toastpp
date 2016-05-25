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


#include "main/imports.h"
#include "main/mtypes.h"
#include "main/bufferobj.h"

#include "st_inlines.h"
#include "st_context.h"
#include "st_cb_bufferobjects.h"

#include "pipe/p_context.h"
#include "pipe/p_defines.h"
#include "pipe/p_inlines.h"



/* Pixel buffers and Vertex/index buffers are handled through these
 * mesa callbacks.  Framebuffer/Renderbuffer objects are
 * created/managed elsewhere.
 */



/**
 * There is some duplication between mesa's bufferobjects and our
 * bufmgr buffers.  Both have an integer handle and a hashtable to
 * lookup an opaque structure.  It would be nice if the handles and
 * internal structure where somehow shared.
 */
static struct gl_buffer_object *
st_bufferobj_alloc(GLcontext *ctx, GLuint name, GLenum target)
{
   struct st_buffer_object *st_obj = ST_CALLOC_STRUCT(st_buffer_object);

   if (!st_obj)
      return NULL;

   _mesa_initialize_buffer_object(&st_obj->Base, name, target);

   return &st_obj->Base;
}



/**
 * Deallocate/free a vertex/pixel buffer object.
 * Called via glDeleteBuffersARB().
 */
static void
st_bufferobj_free(GLcontext *ctx, struct gl_buffer_object *obj)
{
   struct st_buffer_object *st_obj = st_buffer_object(obj);

   if (st_obj->buffer) 
      pipe_buffer_reference(&st_obj->buffer, NULL);

   _mesa_free(st_obj);
}



/**
 * Replace data in a subrange of buffer object.  If the data range
 * specified by size + offset extends beyond the end of the buffer or
 * if data is NULL, no copy is performed.
 * Called via glBufferSubDataARB().
 */
static void
st_bufferobj_subdata(GLcontext *ctx,
		     GLenum target,
		     GLintptrARB offset,
		     GLsizeiptrARB size,
		     const GLvoid * data, struct gl_buffer_object *obj)
{
   struct st_buffer_object *st_obj = st_buffer_object(obj);

   if (offset >= st_obj->size || size > (st_obj->size - offset))
      return;

   st_cond_flush_pipe_buffer_write(st_context(ctx), st_obj->buffer,
				   offset, size, data);
}


/**
 * Called via glGetBufferSubDataARB().
 */
static void
st_bufferobj_get_subdata(GLcontext *ctx,
                         GLenum target,
                         GLintptrARB offset,
                         GLsizeiptrARB size,
                         GLvoid * data, struct gl_buffer_object *obj)
{
   struct st_buffer_object *st_obj = st_buffer_object(obj);

   if (offset >= st_obj->size || size > (st_obj->size - offset))
      return;

   st_cond_flush_pipe_buffer_read(st_context(ctx), st_obj->buffer,
				  offset, size, data);
}


/**
 * Allocate space for and store data in a buffer object.  Any data that was
 * previously stored in the buffer object is lost.  If data is NULL,
 * memory will be allocated, but no copy will occur.
 * Called via glBufferDataARB().
 */
static void
st_bufferobj_data(GLcontext *ctx,
		  GLenum target,
		  GLsizeiptrARB size,
		  const GLvoid * data,
		  GLenum usage, 
		  struct gl_buffer_object *obj)
{
   struct st_context *st = st_context(ctx);
   struct pipe_context *pipe = st->pipe;
   struct st_buffer_object *st_obj = st_buffer_object(obj);
   unsigned buffer_usage;

   st_obj->Base.Size = size;
   st_obj->Base.Usage = usage;
   
   switch(target) {
   case GL_PIXEL_PACK_BUFFER_ARB:
   case GL_PIXEL_UNPACK_BUFFER_ARB:
      buffer_usage = PIPE_BUFFER_USAGE_PIXEL;
      break;
   case GL_ARRAY_BUFFER_ARB:
      buffer_usage = PIPE_BUFFER_USAGE_VERTEX;
      break;
   case GL_ELEMENT_ARRAY_BUFFER_ARB:
      buffer_usage = PIPE_BUFFER_USAGE_INDEX;
      break;
   default:
      buffer_usage = 0;
   }

   pipe_buffer_reference( &st_obj->buffer, NULL );

   st_obj->buffer = pipe_buffer_create( pipe->screen, 32, buffer_usage, size );

   if (!st_obj->buffer) {
      _mesa_error(ctx, GL_OUT_OF_MEMORY, "glBufferDataARB");
      return;
   }

   st_obj->size = size;

   if (data)
      st_no_flush_pipe_buffer_write(st_context(ctx), st_obj->buffer, 0,
				    size, data);
}


/**
 * Called via glMapBufferARB().
 */
static void *
st_bufferobj_map(GLcontext *ctx, GLenum target, GLenum access,
                 struct gl_buffer_object *obj)
{
   struct st_buffer_object *st_obj = st_buffer_object(obj);
   GLuint flags;

   switch (access) {
   case GL_WRITE_ONLY:
      flags = PIPE_BUFFER_USAGE_CPU_WRITE;
      break;
   case GL_READ_ONLY:
      flags = PIPE_BUFFER_USAGE_CPU_READ;
      break;
   case GL_READ_WRITE:
      /* fall-through */
   default:
      flags = PIPE_BUFFER_USAGE_CPU_READ | PIPE_BUFFER_USAGE_CPU_WRITE;
      break;      
   }

   obj->Pointer = st_cond_flush_pipe_buffer_map(st_context(ctx),
						st_obj->buffer,
						flags);
   if(obj->Pointer) {
      obj->Offset = 0;
      obj->Length = obj->Size;
   }
   return obj->Pointer;
}



/**
 * Called via glMapBufferRange().
 */
static void *
st_bufferobj_map_range(GLcontext *ctx, GLenum target, 
                       GLintptr offset, GLsizeiptr length, GLbitfield access,
                       struct gl_buffer_object *obj)
{
   struct pipe_context *pipe = st_context(ctx)->pipe;
   struct st_buffer_object *st_obj = st_buffer_object(obj);
   GLuint flags = 0;
   char *map;

   if (access & GL_MAP_WRITE_BIT)
      flags |= PIPE_BUFFER_USAGE_CPU_WRITE;

   if (access & GL_MAP_READ_BIT)
      flags |= PIPE_BUFFER_USAGE_CPU_READ;

   if (access & GL_MAP_FLUSH_EXPLICIT_BIT)
      flags |= PIPE_BUFFER_USAGE_FLUSH_EXPLICIT;
   
   /* ... other flags ...
    */

   if (access & MESA_MAP_NOWAIT_BIT)
      flags |= PIPE_BUFFER_USAGE_DONTBLOCK;

   assert(offset >= 0);
   assert(length >= 0);
   assert(offset < obj->Size);
   assert(offset + length <= obj->Size);

   map = obj->Pointer = pipe_buffer_map_range(pipe->screen, st_obj->buffer, offset, length, flags);
   if(obj->Pointer) {
      obj->Offset = offset;
      obj->Length = length;
      map += offset;
   }
   
   return map;
}


static void
st_bufferobj_flush_mapped_range(GLcontext *ctx, GLenum target, 
                                GLintptr offset, GLsizeiptr length,
                                struct gl_buffer_object *obj)
{
   struct pipe_context *pipe = st_context(ctx)->pipe;
   struct st_buffer_object *st_obj = st_buffer_object(obj);

   /* Subrange is relative to mapped range */
   assert(offset >= 0);
   assert(length >= 0);
   assert(offset + length <= obj->Length);
   
   pipe_buffer_flush_mapped_range(pipe->screen, st_obj->buffer, 
                                  obj->Offset + offset, length);
}


/**
 * Called via glUnmapBufferARB().
 */
static GLboolean
st_bufferobj_unmap(GLcontext *ctx, GLenum target, struct gl_buffer_object *obj)
{
   struct pipe_context *pipe = st_context(ctx)->pipe;
   struct st_buffer_object *st_obj = st_buffer_object(obj);

   pipe_buffer_unmap(pipe->screen, st_obj->buffer);
   obj->Pointer = NULL;
   obj->Offset = 0;
   obj->Length = 0;
   return GL_TRUE;
}


void
st_init_bufferobject_functions(struct dd_function_table *functions)
{
   functions->NewBufferObject = st_bufferobj_alloc;
   functions->DeleteBuffer = st_bufferobj_free;
   functions->BufferData = st_bufferobj_data;
   functions->BufferSubData = st_bufferobj_subdata;
   functions->GetBufferSubData = st_bufferobj_get_subdata;
   functions->MapBuffer = st_bufferobj_map;
   functions->MapBufferRange = st_bufferobj_map_range;
   functions->FlushMappedBufferRange = st_bufferobj_flush_mapped_range;
   functions->UnmapBuffer = st_bufferobj_unmap;
}
