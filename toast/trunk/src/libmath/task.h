// -*-C++-*-

#ifndef __TASK_H
#define __TASK_H

#ifdef TOAST_THREAD

#include <pthread.h>

#define NUMTHREAD 2


#ifdef TOAST_PARALLEL
#define INITMUTEX(m) pthread_mutex_t m = PTHREAD_MUTEX_INITIALIZER
#define STATICINITMUTEX(m) static pthread_mutex_t m = PTHREAD_MUTEX_INITIALIZER
#define MUTEXLOCK(m) pthread_mutex_lock (m)
#define MUTEXUNLOCK(m) pthread_mutex_unlock (m)
#else
#define INITMUTEX(m)
#define STATICINITMUTEX(m)
#define MUTEXLOCK(m)
#define MUTEXUNLOCK(m)
#endif

typedef struct {
    int proc;
    int np;
    void *data;
} task_data;

MATHLIB void Task_Init (int nth);

class MATHLIB Task {
public:
    Task() {}

    static int nProcessor();
    // number of processors available

    static void SetThreadCount (int _nthread) { nthread = _nthread; }
    // set the number of threads to use for multiprocessing tasks.
    // Default is nProcessor()

    static int GetThreadCount () { return nthread; }
    // return current thread count setting

    static double GetThreadCPUTiming () { return ttime; }
    // returns user time [sec] spent by threads inside Multiprocess().
    // Note that the interpretation of this time is platform-dependent:
    // For Linux, this is the time spent by the master thread. For Sun and
    // SGI, it is the sum of all threads.

    static double GetThreadWallTiming () { return wtime; }
    // Wall clock (real) time spent inside Multiprocess()

    static void Multiprocess (void (*func)(task_data*), void *data, int np = 0);
    // run function 'func' in parallel. User data 'data' are passed to
    // each instance of 'func'. 'np' is the number of threads to create. If
    // np==0 then nthread is used

    static bool IsMultiprocessing() { return is_multiprocessing; }

    inline static void UserMutex_lock() { 
		int res = pthread_mutex_lock (&user_mutex);
		dASSERT(!res, "Mutex could not be locked: %d", res);
	}
    inline static void UserMutex_unlock() {
		int res = pthread_mutex_unlock (&user_mutex);
		dASSERT(!res, "Mutex could not be unlocked: %d", res);
	}

private:
    static int nthread;
    static double ttime; // accumulated cpu time spent multiprocessing
    static double wtime; // accumulated wall clock time spent multiprocessing
    static pthread_mutex_t user_mutex;
    static bool is_multiprocessing;  // we are currently in Multiprocess
};

// ===========================================================================
// class ThreadPool

typedef struct tpool_work {
    void (*routine)(void*,int,int); // task
    void *arg;                      // task arguments
    int idx0, idx1;                 // task index range
    int id;                         // task no. in sequence
    int *counter;                   // pointer to sequence task counter
    struct tpool_work *next;
} tpool_work_t;

typedef struct tpool {
    int num_threads;                // number of threads
    int queue_size;                 // current queue size
    pthread_t *threads;             // array of worker threads
    tpool_work_t *queue_head, *queue_tail;

    pthread_mutex_t queue_lock;
    pthread_cond_t queue_not_empty;
    pthread_cond_t thread_done;
} tpool_t;

class ThreadPool {
public:
    ThreadPool (int num_threads);
    // Construct a pool with `num_threads' threads

    void ProcessSequence (void (*routine)(void*,int,int), void *arg,
        int from, int to, int grain=1);
    // Calls `routine(arg,i)' for a sequence from <= i < to of indices
    // `grain' defines the granularity, i.e. the number of indices assigned
    // to each task
    // Function returns when complete sequence is processed

    inline void LockUserMutex() { pthread_mutex_lock (&user_lock); }
    inline void UnlockUserMutex() { pthread_mutex_unlock (&user_lock); }

private:
    tpool_t *tpool;              // pool properties
    pthread_mutex_t user_lock;   // user-space mutex
};


// ===========================================================================
// class ThreadPool2

typedef struct {
    void(*task)(int,void*);
    void *context;
    pthread_mutex_t mutex;  // general-purpose mutex to be used by workers
} THREAD_GLOBAL;

typedef struct {
    int nth;                      // thread index
    pthread_t thread;             // thread handle
    pthread_mutex_t wakeup_mutex; // locked by worker during task processing
    pthread_mutex_t done_mutex;
    pthread_cond_t wakeup_cond;
    pthread_cond_t done_cond;
    bool wakeup;
    bool done;
    THREAD_GLOBAL *tg;           // pointer to global pool data
} THREAD_DATA;

class ThreadPool2 {
public:
    ThreadPool2 (int num_threads);
    ~ThreadPool2 ();
    static void Initialise (int num_threads);
    static ThreadPool2 *Pool() { return g_tpool2; }

    void MutexLock() { pthread_mutex_lock (&tg.mutex); }
    void MutexUnlock() { pthread_mutex_unlock (&tg.mutex); }
    inline int NumThread() const { return nthread; }

    void Invoke (void(*func)(int,void*), void *context);

private:
    int nthread;                   // number of threads
    THREAD_GLOBAL tg;              // global pool data
    THREAD_DATA *td;               // array of worker threads
    static ThreadPool2 *g_tpool2;

    static void *tpool_thread (void *context);
};

#ifndef __TASK_CC
#ifndef OLD_PARALLEL
void TPool_Init (int nt);
extern ThreadPool *g_tpool;
#endif
#endif

#endif // TOAST_THREAD
#endif // !__TASK_H
