/**
 * Many similar implementations exist. See for example libwsbm
 * or the linux kernel include/atomic.h
 *
 * No copyright claimed on this file.
 *
 */

#ifndef P_ATOMIC_H
#define P_ATOMIC_H

#include "p_compiler.h"
#include "p_defines.h"

#ifdef __cplusplus
extern "C" {
#endif


/* Favor OS-provided implementations.
 *
 * Where no OS-provided implementation is available, fall back to
 * locally coded assembly, compiler intrinsic or ultimately a
 * mutex-based implementation.
 */
#if (defined(PIPE_SUBSYSTEM_WINDOWS_DISPLAY) || \
     defined(PIPE_SUBSYSTEM_WINDOWS_MINIPORT))
#define PIPE_ATOMIC_OS_UNLOCKED 
#elif (defined(PIPE_CC_MSVC) && defined(PIPE_SUBSYSTEM_WINDOWS_USER))
#define PIPE_ATOMIC_OS_MS_INTERLOCK
#elif (defined(PIPE_CC_MSVC) && defined(PIPE_ARCH_X86))
#define PIPE_ATOMIC_ASM_MSVC_X86                
#elif (defined(PIPE_CC_GCC) && defined(PIPE_ARCH_X86))
#define PIPE_ATOMIC_ASM_GCC_X86
#elif defined(PIPE_CC_GCC)
#define PIPE_ATOMIC_GCC_INTRINSIC
#else
#define PIPE_ATOMIC_MUTEX                       
#endif



#if defined(PIPE_ATOMIC_ASM_GCC_X86)

#define PIPE_ATOMIC "GCC x86 assembly"

struct pipe_atomic {
   int32_t count;
};

#define p_atomic_set(_v, _i) ((_v)->count = (_i))
#define p_atomic_read(_v) ((_v)->count)


static INLINE boolean
p_atomic_dec_zero(struct pipe_atomic *v)
{
   unsigned char c;

   __asm__ __volatile__("lock; decl %0; sete %1":"+m"(v->count), "=qm"(c)
			::"memory");

   return c != 0;
}

static INLINE void
p_atomic_inc(struct pipe_atomic *v)
{
   __asm__ __volatile__("lock; incl %0":"+m"(v->count));
}

static INLINE void
p_atomic_dec(struct pipe_atomic *v)
{
   __asm__ __volatile__("lock; decl %0":"+m"(v->count));
}

static INLINE int32_t
p_atomic_cmpxchg(struct pipe_atomic *v, int32_t old, int32_t _new)
{
   return __sync_val_compare_and_swap(&v->count, old, _new);
}
#endif



/* Implementation using GCC-provided synchronization intrinsics
 */
#if defined(PIPE_ATOMIC_GCC_INTRINSIC)

#define PIPE_ATOMIC "GCC Sync Intrinsics"

struct pipe_atomic {
   int32_t count;
};

#define p_atomic_set(_v, _i) ((_v)->count = (_i))
#define p_atomic_read(_v) ((_v)->count)


static INLINE boolean
p_atomic_dec_zero(struct pipe_atomic *v)
{
   return (__sync_sub_and_fetch(&v->count, 1) == 0);
}

static INLINE void
p_atomic_inc(struct pipe_atomic *v)
{
   (void) __sync_add_and_fetch(&v->count, 1);
}

static INLINE void
p_atomic_dec(struct pipe_atomic *v)
{
   (void) __sync_sub_and_fetch(&v->count, 1);
}

static INLINE int32_t
p_atomic_cmpxchg(struct pipe_atomic *v, int32_t old, int32_t _new)
{
   return __sync_val_compare_and_swap(&v->count, old, _new);
}
#endif



/* Unlocked version for single threaded environments, such as some
 * windows kernel modules.
 */
#if defined(PIPE_ATOMIC_OS_UNLOCKED) 

#define PIPE_ATOMIC "Unlocked"

struct pipe_atomic
{
   int32_t count;
};

#define p_atomic_set(_v, _i) ((_v)->count = (_i))
#define p_atomic_read(_v) ((_v)->count)
#define p_atomic_dec_zero(_v) ((boolean) --(_v)->count)
#define p_atomic_inc(_v) ((void) (_v)->count++)
#define p_atomic_dec(_v) ((void) (_v)->count--)
#define p_atomic_cmpxchg(_v, old, _new) ((_v)->count == old ? (_v)->count = (_new) : (_v)->count)

#endif


/* Locally coded assembly for MSVC on x86:
 */
#if defined(PIPE_ATOMIC_ASM_MSVC_X86)

#define PIPE_ATOMIC "MSVC x86 assembly"

struct pipe_atomic
{
   int32_t count;
};

#define p_atomic_set(_v, _i) ((_v)->count = (_i))
#define p_atomic_read(_v) ((_v)->count)

static INLINE boolean
p_atomic_dec_zero(struct pipe_atomic *v)
{
   int32_t *pcount = &v->count;
   unsigned char c;

   __asm {
      mov       eax, [pcount]
      lock dec  dword ptr [eax]
      sete      byte ptr [c]
   }

   return c != 0;
}

static INLINE void
p_atomic_inc(struct pipe_atomic *v)
{
   int32_t *pcount = &v->count;

   __asm {
      mov       eax, [pcount]
      lock inc  dword ptr [eax]
   }
}

static INLINE void
p_atomic_dec(struct pipe_atomic *v)
{
   int32_t *pcount = &v->count;

   __asm {
      mov       eax, [pcount]
      lock dec  dword ptr [eax]
   }
}

static INLINE int32_t
p_atomic_cmpxchg(struct pipe_atomic *v, int32_t old, int32_t _new)
{
   int32_t *pcount = &v->count;
   int32_t orig;

   __asm {
      mov ecx, [pcount]
      mov eax, [old]
      mov edx, [_new]
      lock cmpxchg [ecx], edx
      mov [orig], eax
   }

   return orig;
}
#endif


#if defined(PIPE_ATOMIC_OS_MS_INTERLOCK)

#define PIPE_ATOMIC "MS userspace interlocks"

#include <windows.h>

struct pipe_atomic
{
   volatile long count;
};

#define p_atomic_set(_v, _i) ((_v)->count = (_i))
#define p_atomic_read(_v) ((_v)->count)

static INLINE boolean
p_atomic_dec_zero(struct pipe_atomic *v)
{
   return InterlockedDecrement(&v->count) == 0;
}

static INLINE void
p_atomic_inc(struct pipe_atomic *v)
{
   InterlockedIncrement(&v->count);
}

static INLINE void
p_atomic_dec(struct pipe_atomic *v)
{
   InterlockedDecrement(&v->count);
}

static INLINE int32_t
p_atomic_cmpxchg(struct pipe_atomic *v, int32_t old, int32_t _new)
{
   return InterlockedCompareExchange(&v->count, _new, old);
}

#endif



#if defined(PIPE_ATOMIC_MUTEX)

#define PIPE_ATOMIC "mutex-based fallback"

#include "pipe/p_thread.h"

/**
 * This implementation should really not be used.
 * Add an assembly port instead. It may abort and
 * doesn't destroy used mutexes.
 */

struct pipe_atomic {
   pipe_mutex mutex;
   int32_t count;
};

static INLINE void
p_atomic_set(struct pipe_atomic *v, int32_t i)
{
   pipe_mutex_init(v->mutex);
   pipe_mutex_lock(v->mutex);
   v->count = i;
   pipe_mutex_unlock(v->mutex);
}

static INLINE int32_t
p_atomic_read(struct pipe_atomic *v)
{
   int32_t ret;

   pipe_mutex_lock(v->mutex);
   ret = v->count;
   pipe_mutex_unlock(v->mutex);
   return ret;
}

static INLINE void
p_atomic_inc(struct pipe_atomic *v)
{
   pipe_mutex_lock(v->mutex);
   ++v->count;
   pipe_mutex_unlock(v->mutex);
}

static INLINE void
p_atomic_dec(struct pipe_atomic *v)
{
   pipe_mutex_lock(v->mutex);
   --v->count;
   pipe_mutex_unlock(v->mutex);
}

static INLINE boolean
p_atomic_dec_zero(struct pipe_atomic *v)
{
   boolean ret;

   pipe_mutex_lock(v->mutex);
   ret = (--v->count == 0);
   pipe_mutex_unlock(v->mutex);
   return ret;
}

static INLINE int32_t
p_atomic_cmpxchg(struct pipe_atomic *v, int32_t old, int32_t _new)
{
   int32_t ret;

   pipe_mutex_lock(v->mutex);
   ret = v->count;
   if (ret == old)
      v->count = _new;
   pipe_mutex_unlock(v->mutex);

   return ret;
}

#endif 


#ifndef PIPE_ATOMIC
#error "No pipe_atomic implementation selected"
#endif



#ifdef __cplusplus
}
#endif

#endif /* P_ATOMIC_H */
