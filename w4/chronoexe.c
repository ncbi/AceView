/*  File: chronoexe.c
 *  Author: Jean Thierry-Mieg (mieg@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1992
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: General Timing statistics
 *              This package is self contained since any function         
 *              called from here cannot be timed.
 *              If you compile without flag CHRONO all calls to the     
 *              chrono package just disappear out of the code.          
 *              See chrono.h for how this works.
 *              The basis of this code is the UNIX getrusage() routine,
 *              this is _not_ a POSIX call and may not be available on
 *              all platforms.

 *                                                             *
 *  2 routines are public :                                    *
 *     chrono, chronoReturn                                    *
 *                                                             *
 *     chrono("procedure") switches the timing to "procedure"  *
 *     chronoReturn ends it.                                   *
 *     chrono.h defines chrono = chronoSwitch or void          *
 *              and     chronoReturn = chronoDoReturn or void  *
 *     If you compile without flag CHRONO all calls to the     *
 *     chrono package just disappear out of the code.          *

 * Exported functions: see chrono.h
 * HISTORY:
 * Last edited: Dec  4 11:36 1998 (fw)
 * * Oct 13 14:37 1998 (edgrif): Removed ACEDB defines, replaced with
 *              function calls.
 * * Nov 27 12:36 1995 (mieg)
 * Created: Mon Jun 15 14:42:32 1990 (mieg)
 *-------------------------------------------------------------------
 */
/* $Id: chronoexe.c,v 1.1 2007/03/26 22:47:27 mieg Exp $ */

#include "regular.h"
#include "array.h"
#include "chrono.h"

/* Specify platforms to run chrono on...       */
/* Basically, you can't run this on much, I assume this is 
   because it uses getrusage() to get the timings and this is
   not a POSIX call.              */

/* If we are on the correct platform then define the chrono stuff.  */
/*                                                                  */

#define CHRONO

#ifdef CHRONO
#include "chrono_.h"
#include "freeout.h"

#include <sys/time.h>
#include <sys/resource.h>
#include <sys/types.h>

/* subroutines to return info about system and real time */

typedef struct rusage *Timer_ ;
static struct rusage previous[1], current[1] ;

static BOOL chronoIsRunning = FALSE ;

static int last = 0, level = 0 ;
static Stack chronoStack = 0 ;
static Array chronoArray = 0 ;

/**********************************************/
/* This is not portable at all, we should be using clock ticks to do this,   */
/* see a Unix manual....                                                     */
static double MICRO = .000001 ;

static double convert (struct timeval *tval) /* conversion to seconds */
{
  return (double) (tval->tv_sec + MICRO*tval->tv_usec) ;
}

/******************************************/

static void chronoGo(int i)
{
  Chrono *ch ;

  getrusage (RUSAGE_SELF,current) ;

  ch = arrayp (chronoArray, i, Chrono) ;

  ch->tSys +=  convert (&current->ru_stime) 
    - convert (&previous->ru_stime) ;
  ch->tUser += convert (&current->ru_utime) 
    - convert (&previous->ru_utime) ;
  
  ch = arrayp (chronoArray, 0, Chrono) ;
  *previous = *current ;
 }

/******************************************/

static int chronoOrder (const void *a, const void *b)
{
  const Chrono *ca = (const Chrono *) a, *cb = (const Chrono *) b ;
  /*
    double x = ca->tUser - cb->tUser ;

    if (x < 0) return -1 ;
    else if (x > 0) return 1 ;
  */
  return strcmp(ca->prok, cb->prok) ;
}

/******************************************/
/*********** Public Routines **************/
/******************************************/
/******************************************/

BOOL chronoStart(void)
{
  Chrono *ch ;
  
  chronoArray = arrayReCreate (chronoArray, 64, Chrono) ;
  
  ch = arrayp (chronoArray, 0, Chrono) ;
  ch->prok = "Chrono" ;
  getrusage (RUSAGE_SELF,previous) ;
  
  chronoStack = stackReCreate(chronoStack, 64) ;
  chronoIsRunning = TRUE ;
  last = 0 ;
  return TRUE ;
}

/******************************************/

void chronoStop(void)
{
  chronoIsRunning = FALSE ;
  last = 0 ;
  chronoStack = stackReCreate(chronoStack, 64) ;
}

/******************************************/

static Array chronoDoReport (void)
{
  Array aa = 0 ;
  
  if (chronoArray)
    {
      aa = arrayCopy (chronoArray) ;
      arraySort (aa, chronoOrder) ;
    }
  return aa ;
} /* chronoDoReport */

/******************************************/

void chronoReport (void)
{ 
  Chrono *ch ; 
  double tSysTotal = 0, tUserTotal = 0;
  Array aa = chronoDoReport() ;
  int i = aa ? arrayMax(aa) : 0 ; 

  ch = arrayp (aa, 0, Chrono) - 1 ;
  while(ch++, i--)
    {
      tSysTotal += ch->tSys ;
      tUserTotal += ch->tUser ;
    }
  freeOutf("// Total time : %.2f s  system,  %.2f user  level = %d\n", 
	   tSysTotal, tUserTotal, level) ;
  
  freeOutf ("//                     # of calls      System        %%           User       %%\n") ;

  for (i = 0 ; i < arrayMax(aa) ; i++)
    {
      ch = arrayp (aa, i, Chrono) ;
      
      freeOutf ("// %22s  %6d  %8.2f s  %5d %%     %8.2f s %5d %%\n",
		ch->prok, ch->call,
		ch->tSys,  (int)(100 * (ch->tSys)/(.001+tSysTotal)) ,
		ch->tUser, (int)( 100 * (ch->tUser)/(.001+tUserTotal))) ;
    }
  return ;
} /* chronoReport  */

/******************************************/

/* aliased to chrono in chrono.h */
void chronoSwitch(char *cp)
{
  int i = 0, n ;
  Chrono *ch = 0 ;
  
  if (!chronoIsRunning) 
    return ;
  
  n = arrayMax(chronoArray) ;
  if (n > 0)
    for (i = 0, ch = arrayp (chronoArray, 0, Chrono) ; 
	 ch->prok != cp && i < n ; ch++, i++) ;
  if (i == n)
    {
      ch = arrayp (chronoArray, n, Chrono) ;
      ch->prok = cp ;
    }
  ch->call ++ ;
  chronoGo(last) ;
  level++ ;
  push(chronoStack, last, int) ;
  last = i ;
}

/*****************************************/

/* aliased to chronoReturn in chrono.h */
void chronoDoReturn(void)
{
  if (!chronoIsRunning) 
    return ;
  
  chronoGo(last) ;
  level-- ;
  last = pop(chronoStack, int) ;
}

/*****************************************/
/*****************************************/
#else

BOOL chronoStart (void) { return FALSE ; }
void chronoStop (void) { return ; }
void chronoReport (void) { return ; }

void chronoSwitch(char *cp) { return ; }
void chronoDoReturn(void) { return ; }

#endif  /* CHRONO */
