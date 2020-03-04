

/*************************************************************************/
/*                                                                       */
/*  Copyright (c) 1994 Stanford University                               */
/*                                                                       */
/*  All rights reserved.                                                 */
/*                                                                       */
/*  Permission is given to use, copy, and modify this software for any   */
/*  non-commercial purpose as long as this copyright notice is not       */
/*  removed.  All other uses, including redistribution in whole or in    */
/*  part, are forbidden without prior written permission.                */
/*                                                                       */
/*  This software is provided with absolutely no warranty and no         */
/*  support.                                                             */
/*                                                                       */
/*************************************************************************/

#ifndef _Memory_H
#define _Memory_H 1

#include "defs.h"
#include "particle.h"
#include "box.h"

#define MAX_LOCKS 2048



#include <pthread.h>
#include <stdlib.h>
#include <semaphore.h>
#include <assert.h>
#if !(defined PAGE_SIZE)
#define PAGE_SIZE 4096
#endif
#define __MAX_THREADS__ 256

extern pthread_t __tid__[__MAX_THREADS__];
extern unsigned __threads__;
extern pthread_mutex_t __intern__;
void *our_malloc(size_t size, char * file, unsigned line);


typedef struct _G_Mem g_mem;

struct _G_Mem
{
   pthread_mutex_t io_lock;
   pthread_mutex_t mal_lock;
   pthread_mutex_t single_lock;
   pthread_mutex_t count_lock;
   int count;
   pthread_mutex_t (lock_array)[MAX_LOCKS];
   
struct { pthread_mutex_t bar_mutex; pthread_cond_t bar_cond; unsigned bar_teller; } synch;

   volatile int *i_array;
   volatile double *d_array;
   real f_array[MAX_PROCS][NUM_DIM_POW_2];
   real max_x;
   real min_x;
   real max_y;
   real min_y;
   int id;
};
extern g_mem *G_Memory;

typedef struct _Local_Memory local_memory;
struct _Local_Memory {
#if !(defined NO_PADDING)
   int pad_begin[PAD_SIZE];
#endif
   box *B_Heap;
   int Index_B_Heap;
   int Max_B_Heap;

   particle **Particles;
   int Num_Particles;
   int Max_Particles;

   box *Childless_Partition;
   box *Parent_Partition[MAX_LEVEL];
   int Max_Parent_Level;

   box *Local_Grid;
   real Local_X_Max;
   real Local_X_Min;
   real Local_Y_Max;
   real Local_Y_Min;

   int Total_Work;
   int Min_Work;
   int Max_Work;

   int Time_Step;
   double Time;
   unsigned int init_done_times;
   time_info Timing[MAX_TIME_STEPS];
   
#if !(defined NO_PADDING)
   int pad_end[PAD_SIZE];
#endif
};
extern local_memory Local[MAX_PROCS];

extern void InitGlobalMemory();

#endif /* _Memory_H */

/* Generated from ../Changed/memory.H */
