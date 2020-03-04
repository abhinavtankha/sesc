

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

/*************************************************************************/
/*                                                                       */
/*  SPLASH Ocean Code                                                    */
/*                                                                       */
/*  This application studies the role of eddy and boundary currents in   */
/*  influencing large-scale ocean movements.  This implementation uses   */
/*  dynamically allocated four-dimensional arrays for grid data storage. */
/*                                                                       */
/*  Command line options:                                                */
/*                                                                       */
/*     -nN : Simulate NxN ocean.  N must be (power of 2)+2.              */
/*     -pP : P = number of processors.  P must be power of 2.            */
/*     -eE : E = error tolerance for iterative relaxation.               */
/*     -rR : R = distance between grid points in meters.                 */
/*     -tT : T = timestep in seconds.                                    */
/*     -s  : Print timing statistics.                                    */
/*     -o  : Print out relaxation residual values.                       */
/*     -h  : Print out command line options.                             */
/*                                                                       */
/*  Default: OCEAN -n130 -p1 -e1e-7 -r20000.0 -t28800.0                  */
/*                                                                       */
/*  NOTE: This code works under both the FORK and SPROC models.          */
/*                                                                       */
/*************************************************************************/



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


#define DEFAULT_N      258
#define DEFAULT_P        1
#define DEFAULT_E        1e-7
#define DEFAULT_T    28800.0
#define DEFAULT_R    20000.0
#define UP               0
#define DOWN             1
#define LEFT             2
#define RIGHT            3
#define UPLEFT           4
#define UPRIGHT          5
#define DOWNLEFT         6
#define DOWNRIGHT        7
#if !(defined NO_PADDING)
#define PAGE_SIZE     4096
#else
#define PAGE_SIZE        0
#endif
#include <stdio.h>
#include <math.h>
#include <time.h>

struct multi_struct {
   double err_multi;
} *multi;

struct global_struct {
   int id;
   int starttime;
   int trackstart;
   double psiai;
   double psibi;
} *global;

double ****psi;
double ****psim;
double ***psium;
double ***psilm;
double ***psib;
double ***ga;
double ***gb;
double ****work1;
double ***work2;
double ***work3;
double ****work4;
double ****work5;
double ***work6;
double ****work7;
double ****temparray;
double ***tauz;
double ***oldga;
double ***oldgb;
double *f;
double ****q_multi;
double ****rhs_multi;

struct locks_struct {
   pthread_mutex_t idlock;
   pthread_mutex_t psiailock;
   pthread_mutex_t psibilock;
   pthread_mutex_t donelock;
   pthread_mutex_t error_lock;
   pthread_mutex_t bar_lock;
} *locks;

struct bars_struct {
   
struct { pthread_mutex_t bar_mutex; pthread_cond_t bar_cond; unsigned bar_teller; } iteration;

   
struct { pthread_mutex_t bar_mutex; pthread_cond_t bar_cond; unsigned bar_teller; } gsudn;

   
struct { pthread_mutex_t bar_mutex; pthread_cond_t bar_cond; unsigned bar_teller; } p_setup;
 
   
struct { pthread_mutex_t bar_mutex; pthread_cond_t bar_cond; unsigned bar_teller; } p_redph;
 
   
struct { pthread_mutex_t bar_mutex; pthread_cond_t bar_cond; unsigned bar_teller; } p_soln;
 
   
struct { pthread_mutex_t bar_mutex; pthread_cond_t bar_cond; unsigned bar_teller; } p_subph;
 
   
struct { pthread_mutex_t bar_mutex; pthread_cond_t bar_cond; unsigned bar_teller; } sl_prini;

   
struct { pthread_mutex_t bar_mutex; pthread_cond_t bar_cond; unsigned bar_teller; } sl_psini;

   
struct { pthread_mutex_t bar_mutex; pthread_cond_t bar_cond; unsigned bar_teller; } sl_onetime;

   
struct { pthread_mutex_t bar_mutex; pthread_cond_t bar_cond; unsigned bar_teller; } sl_phase_1;

   
struct { pthread_mutex_t bar_mutex; pthread_cond_t bar_cond; unsigned bar_teller; } sl_phase_2;

   
struct { pthread_mutex_t bar_mutex; pthread_cond_t bar_cond; unsigned bar_teller; } sl_phase_3;

   
struct { pthread_mutex_t bar_mutex; pthread_cond_t bar_cond; unsigned bar_teller; } sl_phase_4;

   
struct { pthread_mutex_t bar_mutex; pthread_cond_t bar_cond; unsigned bar_teller; } sl_phase_5;

   
struct { pthread_mutex_t bar_mutex; pthread_cond_t bar_cond; unsigned bar_teller; } sl_phase_6;

   
struct { pthread_mutex_t bar_mutex; pthread_cond_t bar_cond; unsigned bar_teller; } sl_phase_7;

   
struct { pthread_mutex_t bar_mutex; pthread_cond_t bar_cond; unsigned bar_teller; } sl_phase_8;

   
struct { pthread_mutex_t bar_mutex; pthread_cond_t bar_cond; unsigned bar_teller; } sl_phase_9;

   
struct { pthread_mutex_t bar_mutex; pthread_cond_t bar_cond; unsigned bar_teller; } sl_phase_10;

   
struct { pthread_mutex_t bar_mutex; pthread_cond_t bar_cond; unsigned bar_teller; } error_barrier;

} *bars;

void subblock();
void slave();
int log_2(int);
void printerr(char *);

int nprocs = DEFAULT_P;
double h1 = 1000.0;
double h3 = 4000.0;
double h = 5000.0;
double lf = -5.12e11;
double res = DEFAULT_R;
double dtau = DEFAULT_T;
double f0 = 8.3e-5;
double beta = 2.0e-11;
double gpr = 0.02;
int im = DEFAULT_N;
int jm;
double tolerance = DEFAULT_E;
double eig2;
double ysca;
int jmm1;
double pi;
double t0 = 0.5e-4 ;
double outday0 = 1.0;
double outday1 = 2.0;
double outday2 = 2.0;
double outday3 = 2.0;
double factjacob;
double factlap;
int numlev;
int *imx;
int *jmx;
double *lev_res;
double *lev_tol;
double maxwork = 10000.0;

struct Global_Private {
#if !(defined NO_PADDING)
  char pad[PAGE_SIZE];
#endif
  int *rel_num_x;
  int *rel_num_y;
  int *eist;     
  int *ejst;     
  int *oist;     
  int *ojst;     
  int *rlist;    
  int *rljst;    
  int *rlien;    
  int *rljen;    
  int rownum;
  int colnum;
  int neighbors[8];
  double multi_time;
  double total_time;
} *gp;

double *i_int_coeff;
double *j_int_coeff;
int xprocs;
int yprocs;
int *xpts_per_proc;
int *ypts_per_proc;
int minlevel;
int do_stats = 0;
int do_output = 0;

int main(argc, argv)

int argc;
char *argv[];

{
   int i;
   int j;
   int k;
   double work_multi;
   int my_num;
   int x_part;
   int y_part;
   int d_size;
   int itemp;
   int jtemp;
   double procsqrt;
   FILE *fileptr;
   int iindex;
   int temp = 0;
   char c;
   double min_total;
   double max_total;
   double avg_total;
   double min_multi;
   double max_multi;
   double avg_multi;
   double min_frac;
   double max_frac;
   double avg_frac;
   int ch;
   extern char *optarg;
   unsigned int computeend;
   unsigned int start;

   {long time(); (start) = time(0);}

   while ((ch = getopt(argc, argv, "n:p:e:r:t:soh")) != -1) {
     switch(ch) {
     case 'n': im = atoi(optarg);
               if (log_2(im-2) == -1) {
                 printerr("Grid must be ((power of 2)+2) in each dimension\n");
                 exit(-1);
               }
               break;
     case 'p': nprocs = atoi(optarg);
               if (nprocs < 1) {
                 printerr("P must be >= 1\n");
                 exit(-1);
               }
               if (log_2(nprocs) == -1) {
                 printerr("P must be a power of 2\n");
                 exit(-1);
               }
               break;
     case 'e': tolerance = atof(optarg); break;
     case 'r': res = atof(optarg); break;
     case 't': dtau = atof(optarg); break;
     case 's': do_stats = !do_stats; break;
     case 'o': do_output = !do_output; break;
     case 'h': printf("Usage: OCEAN <options>\n\n");
               printf("options:\n");
               printf("  -nN : Simulate NxN ocean.  N must be (power of 2)+2.\n");
               printf("  -pP : P = number of processors.  P must be power of 2.\n");
               printf("  -eE : E = error tolerance for iterative relaxation.\n");
               printf("  -rR : R = distance between grid points in meters.\n");
               printf("  -tT : T = timestep in seconds.\n");
               printf("  -s  : Print timing statistics.\n");
               printf("  -o  : Print out relaxation residual values.\n");
               printf("  -h  : Print out command line options.\n\n");
               printf("Default: OCEAN -n%1d -p%1d -e%1g -r%1g -t%1g\n",
                       DEFAULT_N,DEFAULT_P,DEFAULT_E,DEFAULT_R,DEFAULT_T);
               exit(0);
               break;
     }
   }

   {__tid__[__threads__++]=pthread_self();} 

   jm = im;
   printf("\n");
   printf("Ocean simulation with W-cycle multigrid solver\n");
   printf("    Processors                         : %1d\n",nprocs);
   printf("    Grid size                          : %1d x %1d\n",im,jm);
   printf("    Grid resolution (meters)           : %0.2f\n",res);
   printf("    Time between relaxations (seconds) : %0.0f\n",dtau);
   printf("    Error tolerance                    : %0.7g\n",tolerance);
   printf("\n");

   xprocs = 0;
   yprocs = 0;
   procsqrt = sqrt((double) nprocs);
   j = (int) procsqrt;
   while ((xprocs == 0) && (j > 0)) {
     k = nprocs / j;
     if (k * j == nprocs) {
       if (k > j) {
         xprocs = j;
         yprocs = k;
       } else {
         xprocs = k;
         yprocs = j;
       }
     }
     j--;
   }
   if (xprocs == 0) {
     printerr("Could not find factors for subblocking\n");
     exit(-1);
   }  

   minlevel = 0;
   itemp = 1;
   jtemp = 1;
   numlev = 0;
   minlevel = 0;
   while (itemp < (im-2)) {
     itemp = itemp*2;
     jtemp = jtemp*2;
     if ((itemp/yprocs > 1) && (jtemp/xprocs > 1)) {
       numlev++;
     }
   }  
   
   if (numlev == 0) {
     printerr("Must have at least 2 grid points per processor in each dimension\n");
     exit(-1);
   }

   imx = (int *) our_malloc(numlev*sizeof(int),__FILE__,__LINE__);;
   jmx = (int *) our_malloc(numlev*sizeof(int),__FILE__,__LINE__);;
   lev_res = (double *) our_malloc(numlev*sizeof(double),__FILE__,__LINE__);;
   lev_tol = (double *) our_malloc(numlev*sizeof(double),__FILE__,__LINE__);;
   i_int_coeff = (double *) our_malloc(numlev*sizeof(double),__FILE__,__LINE__);;
   j_int_coeff = (double *) our_malloc(numlev*sizeof(double),__FILE__,__LINE__);;
   xpts_per_proc = (int *) our_malloc(numlev*sizeof(int),__FILE__,__LINE__);;
   ypts_per_proc = (int *) our_malloc(numlev*sizeof(int),__FILE__,__LINE__);;

   imx[numlev-1] = im;
   jmx[numlev-1] = jm;
   lev_res[numlev-1] = res;
   lev_tol[numlev-1] = tolerance;

   for (i=numlev-2;i>=0;i--) {
     imx[i] = ((imx[i+1] - 2) / 2) + 2;
     jmx[i] = ((jmx[i+1] - 2) / 2) + 2;
     lev_res[i] = lev_res[i+1] * 2;
   }

   for (i=0;i<numlev;i++) {
     xpts_per_proc[i] = (jmx[i]-2) / xprocs;
     ypts_per_proc[i] = (imx[i]-2) / yprocs;
   }  
   for (i=numlev-1;i>=0;i--) {
     if ((xpts_per_proc[i] < 2) || (ypts_per_proc[i] < 2)) {
       minlevel = i+1;
       break;
     }
   }    
 
   for (i=0;i<numlev;i++) {
     temp += imx[i];
   }
   temp = 0;
   j = 0;
   for (k=0;k<numlev;k++) {
     for (i=0;i<imx[k];i++) {
       j++;
       temp += jmx[k];
     }
   }

   d_size = nprocs*sizeof(double ***);
   psi = (double ****) our_malloc(d_size,__FILE__,__LINE__);;
   psim = (double ****) our_malloc(d_size,__FILE__,__LINE__);;
   work1 = (double ****) our_malloc(d_size,__FILE__,__LINE__);;
   work4 = (double ****) our_malloc(d_size,__FILE__,__LINE__);;
   work5 = (double ****) our_malloc(d_size,__FILE__,__LINE__);;
   work7 = (double ****) our_malloc(d_size,__FILE__,__LINE__);;
   temparray = (double ****) our_malloc(d_size,__FILE__,__LINE__);;

   d_size = 2*sizeof(double **);
   for (i=0;i<nprocs;i++) {
     psi[i] = (double ***) our_malloc(d_size,__FILE__,__LINE__);;
     psim[i] = (double ***) our_malloc(d_size,__FILE__,__LINE__);;
     work1[i] = (double ***) our_malloc(d_size,__FILE__,__LINE__);;
     work4[i] = (double ***) our_malloc(d_size,__FILE__,__LINE__);;
     work5[i] = (double ***) our_malloc(d_size,__FILE__,__LINE__);;
     work7[i] = (double ***) our_malloc(d_size,__FILE__,__LINE__);;
     temparray[i] = (double ***) our_malloc(d_size,__FILE__,__LINE__);;
   }

   d_size = nprocs*sizeof(double **);
   psium = (double ***) our_malloc(d_size,__FILE__,__LINE__);;
   psilm = (double ***) our_malloc(d_size,__FILE__,__LINE__);;
   psib = (double ***) our_malloc(d_size,__FILE__,__LINE__);;
   ga = (double ***) our_malloc(d_size,__FILE__,__LINE__);;
   gb = (double ***) our_malloc(d_size,__FILE__,__LINE__);;
   work2 = (double ***) our_malloc(d_size,__FILE__,__LINE__);;
   work3 = (double ***) our_malloc(d_size,__FILE__,__LINE__);;
   work6 = (double ***) our_malloc(d_size,__FILE__,__LINE__);;
   tauz = (double ***) our_malloc(d_size,__FILE__,__LINE__);;
   oldga = (double ***) our_malloc(d_size,__FILE__,__LINE__);;
   oldgb = (double ***) our_malloc(d_size,__FILE__,__LINE__);;

   gp = (struct Global_Private *) our_malloc((nprocs+1)*sizeof(struct Global_Private),__FILE__,__LINE__);;
   for (i=0;i<nprocs;i++) {
     gp[i].rel_num_x = (int *) our_malloc(numlev*sizeof(int),__FILE__,__LINE__);;
     gp[i].rel_num_y = (int *) our_malloc(numlev*sizeof(int),__FILE__,__LINE__);;
     gp[i].eist = (int *) our_malloc(numlev*sizeof(int),__FILE__,__LINE__);;
     gp[i].ejst = (int *) our_malloc(numlev*sizeof(int),__FILE__,__LINE__);;
     gp[i].oist = (int *) our_malloc(numlev*sizeof(int),__FILE__,__LINE__);;
     gp[i].ojst = (int *) our_malloc(numlev*sizeof(int),__FILE__,__LINE__);;
     gp[i].rlist = (int *) our_malloc(numlev*sizeof(int),__FILE__,__LINE__);;
     gp[i].rljst = (int *) our_malloc(numlev*sizeof(int),__FILE__,__LINE__);;
     gp[i].rlien = (int *) our_malloc(numlev*sizeof(int),__FILE__,__LINE__);;
     gp[i].rljen = (int *) our_malloc(numlev*sizeof(int),__FILE__,__LINE__);;
     gp[i].multi_time = 0;
     gp[i].total_time = 0;
   }

   subblock();

   x_part = (jm - 2)/xprocs + 2;
   y_part = (im - 2)/yprocs + 2;

   d_size = x_part*y_part*sizeof(double) + y_part*sizeof(double *);

   global = (struct global_struct *) our_malloc(sizeof(struct global_struct),__FILE__,__LINE__);;  
   for (i=0;i<nprocs;i++) {
     psi[i][0] = (double **) our_malloc(d_size,__FILE__,__LINE__);;
     psi[i][1] = (double **) our_malloc(d_size,__FILE__,__LINE__);;
     psim[i][0] = (double **) our_malloc(d_size,__FILE__,__LINE__);;
     psim[i][1] = (double **) our_malloc(d_size,__FILE__,__LINE__);;
     psium[i] = (double **) our_malloc(d_size,__FILE__,__LINE__);;
     psilm[i] = (double **) our_malloc(d_size,__FILE__,__LINE__);;
     psib[i] = (double **) our_malloc(d_size,__FILE__,__LINE__);;
     ga[i] = (double **) our_malloc(d_size,__FILE__,__LINE__);;
     gb[i] = (double **) our_malloc(d_size,__FILE__,__LINE__);;
     work1[i][0] = (double **) our_malloc(d_size,__FILE__,__LINE__);;
     work1[i][1] = (double **) our_malloc(d_size,__FILE__,__LINE__);;
     work2[i] = (double **) our_malloc(d_size,__FILE__,__LINE__);;
     work3[i] = (double **) our_malloc(d_size,__FILE__,__LINE__);;
     work4[i][0] = (double **) our_malloc(d_size,__FILE__,__LINE__);;
     work4[i][1] = (double **) our_malloc(d_size,__FILE__,__LINE__);;
     work5[i][0] = (double **) our_malloc(d_size,__FILE__,__LINE__);;
     work5[i][1] = (double **) our_malloc(d_size,__FILE__,__LINE__);;
     work6[i] = (double **) our_malloc(d_size,__FILE__,__LINE__);;
     work7[i][0] = (double **) our_malloc(d_size,__FILE__,__LINE__);;
     work7[i][1] = (double **) our_malloc(d_size,__FILE__,__LINE__);;
     temparray[i][0] = (double **) our_malloc(d_size,__FILE__,__LINE__);;
     temparray[i][1] = (double **) our_malloc(d_size,__FILE__,__LINE__);;
     tauz[i] = (double **) our_malloc(d_size,__FILE__,__LINE__);;
     oldga[i] = (double **) our_malloc(d_size,__FILE__,__LINE__);;
     oldgb[i] = (double **) our_malloc(d_size,__FILE__,__LINE__);;
   }
   f = (double *) our_malloc(im*sizeof(double),__FILE__,__LINE__);;

   multi = (struct multi_struct *) our_malloc(sizeof(struct multi_struct),__FILE__,__LINE__);;

   d_size = numlev*sizeof(double **);
   if (numlev%2 == 1) {         /* To make sure that the actual data 
                                   starts double word aligned, add an extra
                                   pointer */
     d_size += sizeof(double **);
   }
   for (i=0;i<numlev;i++) {
     d_size += ((imx[i]-2)/yprocs+2)*((jmx[i]-2)/xprocs+2)*sizeof(double)+
              ((imx[i]-2)/yprocs+2)*sizeof(double *);
   }

   d_size *= nprocs;

   if (nprocs%2 == 1) {         /* To make sure that the actual data 
                                   starts double word aligned, add an extra
                                   pointer */
     d_size += sizeof(double ***);
   }

   d_size += nprocs*sizeof(double ***);
   q_multi = (double ****) our_malloc(d_size,__FILE__,__LINE__);;
   rhs_multi = (double ****) our_malloc(d_size,__FILE__,__LINE__);;

   locks = (struct locks_struct *) our_malloc(sizeof(struct locks_struct),__FILE__,__LINE__);;
   bars = (struct bars_struct *) our_malloc(sizeof(struct bars_struct),__FILE__,__LINE__);;

   {pthread_mutex_init(&(locks->idlock),NULL);}
   {pthread_mutex_init(&(locks->psiailock),NULL);}
   {pthread_mutex_init(&(locks->psibilock),NULL);}
   {pthread_mutex_init(&(locks->donelock),NULL);}
   {pthread_mutex_init(&(locks->error_lock),NULL);}
   {pthread_mutex_init(&(locks->bar_lock),NULL);}

   {
pthread_mutex_init(&((bars->iteration).bar_mutex), NULL);
pthread_cond_init(&((bars->iteration).bar_cond), NULL);
(bars->iteration).bar_teller=0;
}
   {
pthread_mutex_init(&((bars->gsudn).bar_mutex), NULL);
pthread_cond_init(&((bars->gsudn).bar_cond), NULL);
(bars->gsudn).bar_teller=0;
}
   {
pthread_mutex_init(&((bars->p_setup).bar_mutex), NULL);
pthread_cond_init(&((bars->p_setup).bar_cond), NULL);
(bars->p_setup).bar_teller=0;
} 
   {
pthread_mutex_init(&((bars->p_redph).bar_mutex), NULL);
pthread_cond_init(&((bars->p_redph).bar_cond), NULL);
(bars->p_redph).bar_teller=0;
} 
   {
pthread_mutex_init(&((bars->p_soln).bar_mutex), NULL);
pthread_cond_init(&((bars->p_soln).bar_cond), NULL);
(bars->p_soln).bar_teller=0;
} 
   {
pthread_mutex_init(&((bars->p_subph).bar_mutex), NULL);
pthread_cond_init(&((bars->p_subph).bar_cond), NULL);
(bars->p_subph).bar_teller=0;
} 
   {
pthread_mutex_init(&((bars->sl_prini).bar_mutex), NULL);
pthread_cond_init(&((bars->sl_prini).bar_cond), NULL);
(bars->sl_prini).bar_teller=0;
}
   {
pthread_mutex_init(&((bars->sl_psini).bar_mutex), NULL);
pthread_cond_init(&((bars->sl_psini).bar_cond), NULL);
(bars->sl_psini).bar_teller=0;
}
   {
pthread_mutex_init(&((bars->sl_onetime).bar_mutex), NULL);
pthread_cond_init(&((bars->sl_onetime).bar_cond), NULL);
(bars->sl_onetime).bar_teller=0;
}
   {
pthread_mutex_init(&((bars->sl_phase_1).bar_mutex), NULL);
pthread_cond_init(&((bars->sl_phase_1).bar_cond), NULL);
(bars->sl_phase_1).bar_teller=0;
}
   {
pthread_mutex_init(&((bars->sl_phase_2).bar_mutex), NULL);
pthread_cond_init(&((bars->sl_phase_2).bar_cond), NULL);
(bars->sl_phase_2).bar_teller=0;
}
   {
pthread_mutex_init(&((bars->sl_phase_3).bar_mutex), NULL);
pthread_cond_init(&((bars->sl_phase_3).bar_cond), NULL);
(bars->sl_phase_3).bar_teller=0;
}
   {
pthread_mutex_init(&((bars->sl_phase_4).bar_mutex), NULL);
pthread_cond_init(&((bars->sl_phase_4).bar_cond), NULL);
(bars->sl_phase_4).bar_teller=0;
}
   {
pthread_mutex_init(&((bars->sl_phase_5).bar_mutex), NULL);
pthread_cond_init(&((bars->sl_phase_5).bar_cond), NULL);
(bars->sl_phase_5).bar_teller=0;
}
   {
pthread_mutex_init(&((bars->sl_phase_6).bar_mutex), NULL);
pthread_cond_init(&((bars->sl_phase_6).bar_cond), NULL);
(bars->sl_phase_6).bar_teller=0;
}
   {
pthread_mutex_init(&((bars->sl_phase_7).bar_mutex), NULL);
pthread_cond_init(&((bars->sl_phase_7).bar_cond), NULL);
(bars->sl_phase_7).bar_teller=0;
}
   {
pthread_mutex_init(&((bars->sl_phase_8).bar_mutex), NULL);
pthread_cond_init(&((bars->sl_phase_8).bar_cond), NULL);
(bars->sl_phase_8).bar_teller=0;
}
   {
pthread_mutex_init(&((bars->sl_phase_9).bar_mutex), NULL);
pthread_cond_init(&((bars->sl_phase_9).bar_cond), NULL);
(bars->sl_phase_9).bar_teller=0;
}
   {
pthread_mutex_init(&((bars->sl_phase_10).bar_mutex), NULL);
pthread_cond_init(&((bars->sl_phase_10).bar_cond), NULL);
(bars->sl_phase_10).bar_teller=0;
}
   {
pthread_mutex_init(&((bars->error_barrier).bar_mutex), NULL);
pthread_cond_init(&((bars->error_barrier).bar_cond), NULL);
(bars->error_barrier).bar_teller=0;
}

   link_all();

   multi->err_multi = 0.0;
   i_int_coeff[0] = 0.0;
   j_int_coeff[0] = 0.0;
   for (i=0;i<numlev;i++) {
     i_int_coeff[i] = 1.0/(imx[i]-1);
     j_int_coeff[i] = 1.0/(jmx[i]-1);
   }

/* initialize constants and variables

   id is a global shared variable that has fetch-and-add operations
   performed on it by processes to obtain their pids.   */

   global->id = 0;
   global->psibi = 0.0;
   pi = atan(1.0);
   pi = 4.*pi;

   factjacob = -1./(12.*res*res);
   factlap = 1./(res*res);
   eig2 = -h*f0*f0/(h1*h3*gpr);

   jmm1 = jm-1 ;
   ysca = ((double) jmm1)*res ;

   im = (imx[numlev-1]-2)/yprocs + 2;
   jm = (jmx[numlev-1]-2)/xprocs + 2;

   for (i=1;i<nprocs;i++) {
     {
pthread_mutex_lock(&__intern__);
assert(__threads__<__MAX_THREADS__);
pthread_create(&(__tid__[__threads__++]), NULL, (void*(*)(void *))(slave), NULL);
pthread_mutex_unlock(&__intern__);
}  
   }

   if (do_output) {
     printf("                       MULTIGRID OUTPUTS\n");
   }

   slave();
   {int aantal; for(aantal=nprocs-1;aantal>0;aantal--) pthread_join(__tid__[--__threads__], NULL);}
   {long time(); (computeend) = time(0);}

   printf("\n");
   printf("                       PROCESS STATISTICS\n");
   printf("                  Total          Multigrid         Multigrid\n");
   printf(" Proc             Time             Time            Fraction\n");
   printf("    0   %15.0f    %15.0f        %10.3f\n",
          gp[0].total_time,gp[0].multi_time,
          gp[0].multi_time/gp[0].total_time);

   if (do_stats) {
     min_total = max_total = avg_total = gp[0].total_time;
     min_multi = max_multi = avg_multi = gp[0].multi_time;
     min_frac = max_frac = avg_frac = gp[0].multi_time/gp[0].total_time;
     for (i=1;i<nprocs;i++) {
       if (gp[i].total_time > max_total) {
         max_total = gp[i].total_time;
       }
       if (gp[i].total_time < min_total) {
         min_total = gp[i].total_time;
       }
       if (gp[i].multi_time > max_multi) {
         max_multi = gp[i].multi_time;
       }
       if (gp[i].multi_time < min_multi) {
         min_multi = gp[i].multi_time;
       }
       if (gp[i].multi_time/gp[i].total_time > max_frac) {
         max_frac = gp[i].multi_time/gp[i].total_time;
       }
       if (gp[i].multi_time/gp[i].total_time < min_frac) {
         min_frac = gp[i].multi_time/gp[i].total_time;
       }
       avg_total += gp[i].total_time;
       avg_multi += gp[i].multi_time;
       avg_frac += gp[i].multi_time/gp[i].total_time;
     }
     avg_total = avg_total / nprocs;
     avg_multi = avg_multi / nprocs;
     avg_frac = avg_frac / nprocs;
     for (i=1;i<nprocs;i++) {
       printf("  %3d   %15.0f    %15.0f        %10.3f\n",
	      i,gp[i].total_time,gp[i].multi_time,
	      gp[i].multi_time/gp[i].total_time);
     }
     printf("  Avg   %15.0f    %15.0f        %10.3f\n",
            avg_total,avg_multi,avg_frac);
     printf("  Min   %15.0f    %15.0f        %10.3f\n",
            min_total,min_multi,min_frac);
     printf("  Max   %15.0f    %15.0f        %10.3f\n",
            max_total,max_multi,max_frac);
   }
   printf("\n");

   global->starttime = start;
   printf("                       TIMING INFORMATION\n");
   printf("Start time                        : %16d\n",
           global->starttime);
   printf("Initialization finish time        : %16d\n",
           global->trackstart);
   printf("Overall finish time               : %16d\n",
           computeend);
   printf("Total time with initialization    : %16d\n",
           computeend-global->starttime);
   printf("Total time without initialization : %16d\n",
           computeend-global->trackstart);
   printf("    (excludes first timestep)\n");
   printf("\n");

   {exit(0);}

   return 0;
}

int log_2(number)

int number;

{
  int cumulative = 1;
  int out = 0;
  int done = 0;

  while ((cumulative < number) && (!done) && (out < 50)) {
    if (cumulative == number) {
      done = 1;
    } else {
      cumulative = cumulative * 2;
      out ++;
    }
  }

  if (cumulative == number) {
    return(out);
  } else {
    return(-1);
  }
}

void printerr(s)

char *s;

{
  fprintf(stderr,"ERROR: %s\n",s);
}


/* Generated from ../Changed/main.C */
