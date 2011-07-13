// -*-C++-*-

#ifndef __TASK_H
#define __TASK_H

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

class Task {
public:
    Task() {}

    static int nProcessor();
    // number of processors available

    static void SetThreadCount (int _nthread) { nthread = _nthread; }
    // set the number of threads to use for multiprocessing tasks.
    // Default is nProcessor()

    static int GetThreadCount () { return nthread; }
    // return current thread count setting

    static double GetThreadTiming () { return ttime; }
    // returns user time [sec] spent by threads inside Multiprocess().
    // Note that the interpretation of this time is platform-dependent:
    // For Linux, this is the time spent by the master thread. For Sun and
    // SGI, it is the sum of all threads.

    static void Multiprocess (void *func(task_data*), void *data, int np = 0);
    // run function 'func' in parallel. User data 'data' are passed to
    // each instance of 'func'. 'np' is the number of threads to create. If
    // np==0 then nthread is used


private:
    static int nthread;
    static double ttime; // accumulated time spent multiprocessing
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

#ifndef __TASK_CC
#ifndef OLD_PARALLEL
void TPool_Init (int nt);
extern ThreadPool *g_tpool;
#endif
#endif

#endif // !__TASK_H
