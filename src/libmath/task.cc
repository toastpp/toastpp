#define __TASK_CC
//#define DEBUG_THREAD

#include <sys/types.h>
#include <sys/sysinfo.h>
#if defined (sun)
#include <thread.h>
#include <unistd.h>
#elif defined (__sgi)
#include <sys/sysmp.h>
#endif
#include <iostream.h>
#include <stdio.h>
#include <sys/times.h>
#include "task.h"

// comment this out to autodetect number of processors
// #define FIXEDNUMPROC 2

// initialisation of static members

int Task::nthread = Task::nProcessor();
double Task::ttime = 0.0;

int Task::nProcessor ()
{
#if defined (FIXEDNUMPROC)
    cerr << "fixednumproc" << endl;
    return FIXEDNUMPROC;
#elif defined (__linux)
    cerr << "linux" << endl;
    return get_nprocs();
#elif defined (__sgi)
    cerr << "sgi" << endl;
    return sysmp (MP_NAPROCS);
#elif defined (sun)
    cerr << "sunos" << endl;
    return sysconf(_SC_NPROCESSORS_ONLN);
#else
    cerr << "unknown" << endl;
    return 1;
#endif
}

void Task::Multiprocess (void *func(task_data*), void *data, int np)
{
    struct tms tm1, tm2;

    if (!np) np = nthread;
    int p;
    task_data *td = new task_data[np];
    pthread_t *thread = new pthread_t[np];

    cerr << "Multiprocess: Branch into " << np << " threads\n";
    times (&tm1);

    // create worker threads
    for (p = 1; p < np; p++) {
        td[p].proc = p;
	td[p].np   = np;
	td[p].data = data;
	pthread_create (thread+p, NULL, (void*(*)(void*))func, (void*)(td+p));
    }

    // let the master thread do some work
    td[0].proc = 0;
    td[0].np   = np;
    td[0].data = data;
    func (td);

    // wait for workers to finish
    for (p = 1; p < np; p++)
        pthread_join (thread[p], NULL);

    times (&tm2);
    ttime += (double)(tm2.tms_utime - tm1.tms_utime)/(double)CLOCKS_PER_SEC;
    cerr << "Multiprocess: Threads joined\n";

    delete []td;
    delete []thread;
}

// ===========================================================================
// class ThreadPool

// this is the loop for worker threads waiting for tasks to arrive
// in the queue

void tpool_thread (tpool_t *tpool)
{
    tpool_work_t *my_workp;

    for (;;) {
        pthread_mutex_lock (&tpool->queue_lock);
	while (!tpool->queue_size)
	    pthread_cond_wait (&tpool->queue_not_empty, &tpool->queue_lock);
	my_workp = tpool->queue_head;
	tpool->queue_size--;
	if (tpool->queue_size == 0)
	    tpool->queue_head = tpool->queue_tail = NULL;
	else
	   tpool->queue_head = my_workp->next;

#ifdef DEBUG_THREAD
	printf ("Waking thread for task %3d of sequence %x\n",
		my_workp->id, my_workp->counter);
#endif

	pthread_mutex_unlock (&tpool->queue_lock);
	(*(my_workp->routine))(my_workp->arg, my_workp->idx0, my_workp->idx1);
	pthread_mutex_lock (&tpool->queue_lock);
	(*my_workp->counter)--;

#ifdef DEBUG_THREAD
	printf ("Thread finished task   %3d of sequence %x (%d remaining)\n",
		my_workp->id, my_workp->counter, *my_workp->counter);
#endif
	delete my_workp;
	pthread_cond_broadcast (&tpool->thread_done);
	pthread_mutex_unlock (&tpool->queue_lock);
    }
}

ThreadPool::ThreadPool (int num_threads)
{
    tpool = new tpool_t;
    tpool->num_threads = num_threads;
    tpool->queue_size = 0;
    tpool->queue_head = NULL;
    tpool->queue_tail = NULL;
    tpool->threads = new pthread_t[num_threads];
    if (pthread_mutex_init (&(tpool->queue_lock), NULL))
        cerr << "ThreadPool: pthread_mutex_init failed\n";
    if (pthread_mutex_init (&user_lock, NULL))
        cerr << "ThreadPool: pthread_mutex_init failed\n";
    if (pthread_cond_init (&(tpool->queue_not_empty), NULL))
        cerr << "ThreadPool: pthread_cond_init failed\n";
    if (pthread_cond_init (&(tpool->thread_done), NULL))
        cerr << "ThreadPool: pthread_cond_init failed\n";

    // create worker threads
    for (int i = 0; i < num_threads; i++) {
        if (pthread_create (tpool->threads+i, NULL,
			    (void*(*)(void*))tpool_thread, (void*)tpool))
	    cerr << "TreadPool: pthread_create failed\n";
    }

#ifdef DEBUG_THREAD
    cerr << "Created thread pool with " << num_threads << " threads\n";
#endif
}

void ThreadPool::ProcessSequence (void (*routine)(void*,int,int), void *arg,
    int from, int to, int grain)
{
    tpool_work_t *workp;
    int n, *task_count;

    pthread_mutex_lock (&tpool->queue_lock);
    if (grain <= 0) grain = 1;
    n = (to-from+grain-1)/grain;
    if (n <= 0) { // nothing to do
        pthread_mutex_unlock (&tpool->queue_lock);
	return;
    } else {
        task_count = new int; // counts tasks remaining for this sequence
	*task_count = n;
    }

    // add loops to queue
    for (int i = 0; i < n; i++) {
        workp = new tpool_work_t;
	workp->routine = routine;
	workp->arg     = arg;
	if ((workp->idx1 = (workp->idx0 = i*grain+from) + grain) > to)
	    workp->idx1 = to;
	workp->id      = i;
	workp->counter = task_count;
	workp->next    = NULL;
	if (tpool->queue_tail) tpool->queue_tail->next = workp;
	else                   tpool->queue_head = workp;
	tpool->queue_tail = workp;
    }

    tpool->queue_size += n;
    if (tpool->queue_size == n) // wake up dormant threads
        pthread_cond_broadcast (&tpool->queue_not_empty);

#ifdef DEBUG_THREAD
    cerr << "Added sequence " << task_count << " (" << n
	 << " tasks) to queue, queue size now " << tpool->queue_size << endl;
#endif // DEBUG_THREAD

    // wait until sequence counter is down to zero
    while (*task_count)
        pthread_cond_wait (&tpool->thread_done, &tpool->queue_lock);

#ifdef DEBUG_THREAD
    cerr << "Sequence finished, queue size now " << tpool->queue_size << endl;
    cerr << "Returning from ProcessSequence\n";
#endif

    delete task_count;
    pthread_mutex_unlock (&tpool->queue_lock);
}

#ifndef OLD_PARALLEL
ThreadPool *g_tpool = 0;

void TPool_Init (int nt)
{
    g_tpool = new ThreadPool (nt);
}

#endif
