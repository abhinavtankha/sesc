

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

#include "defs.h"
#include "memory.h"



#include <pthread.h>
#include <stdlib.h>
#include <semaphore.h>
#include <assert.h>
#if !(defined PAGE_SIZE)
#define PAGE_SIZE 4096
#endif
#define __MAX_THREADS__ 256

pthread_t __tid__[__MAX_THREADS__];
unsigned __threads__=0;
pthread_mutex_t __intern__;
void *our_malloc(size_t size, char * file, unsigned line) { return malloc(size); }


g_mem *G_Memory;
local_memory Local[MAX_PROCS];

/*
 *  InitGlobalMemory ()
 *
 *  Args : none.
 *
 *  Returns : nothing.
 *
 *  Side Effects : Allocates all the global storage for G_Memory.
 *
 */
void
InitGlobalMemory ()
{
   int i;

   G_Memory = (g_mem *) our_malloc(sizeof(g_mem),__FILE__,__LINE__);;
   G_Memory->i_array = (int *) our_malloc(Number_Of_Processors * sizeof(int),__FILE__,__LINE__);;
   G_Memory->d_array = (double *) our_malloc(Number_Of_Processors
					 * sizeof(double),__FILE__,__LINE__);;
   if (G_Memory == NULL) {
      printf("Ran out of global memory in InitGlobalMemory\n");
      exit(-1);
   }
   G_Memory->count = 0;
   G_Memory->id = 0;
   {pthread_mutex_init(&(G_Memory->io_lock),NULL);};
   {pthread_mutex_init(&(G_Memory->mal_lock),NULL);};
   {pthread_mutex_init(&(G_Memory->single_lock),NULL);};
   {pthread_mutex_init(&(G_Memory->count_lock),NULL);};
   { int i; for(i = 0; i < (MAX_LOCKS); i++) pthread_mutex_init(&((G_Memory->lock_array)[i]), NULL); };
   {
pthread_mutex_init(&((G_Memory->synch).bar_mutex), NULL);
pthread_cond_init(&((G_Memory->synch).bar_cond), NULL);
(G_Memory->synch).bar_teller=0;
};
   G_Memory->max_x = -MAX_REAL;
   G_Memory->min_x = MAX_REAL;
   G_Memory->max_y = -MAX_REAL;
   G_Memory->min_y = MAX_REAL;
}



/* Generated from ../Source/memory.C */
