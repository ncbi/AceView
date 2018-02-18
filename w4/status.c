/*  File: status.c
 *  Author: Jean Thierry-Mieg (mieg@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1991
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
 **  Status  of the ACeDB program.                            **
 * Exported functions:
 **  acedbstatus and dosomething                               
     dosomething is for a single user DOS system                
 * HISTORY:
 * Last edited: Dec  3 14:43 1998 (edgrif)
 * * Dec  3 14:43 1998 (edgrif): Change calls to new interface to aceversion.
 * Created: Fri Nov 29 12:20:30 1991 (mieg)
 *-------------------------------------------------------------------
 */

/* $Id: status.c,v 1.4 2015/09/18 22:05:16 mieg Exp $ */
 
#include "acedb.h"
#include "aceversion.h"					    /* For version reporting. */
#include "version.h"
#include "disk.h"    /*for diskavail() and saveAll*/
#include "block.h"   /*for blockavail()*/
#include "lex.h"
#include "key.h"
#include "bs.h"
#include "cache.h"
#include "cachedisp.h"		/* for cacheShow */
#include "graph.h"
#include "display.h"
#include "pick.h"
#include "session.h"
#include "a.h"			/* for aStatus */


/******************************************************************/

void acedbstatus(void);


extern int lexHashMax(int t) ;

/******************************************************************/

static int NbMaxArray = 0 ;

static void arrayDoShow (void) { arrayReport(NbMaxArray) ; }

static void arrayDoMark(void) { NbMaxArray = arrayReportMark() ; }

/******************************************************************/

extern mysize_t _stklen, stackused(void) ;
#ifndef __MSDOS__
#define coreleft()     0
#define _stklen        0
#endif

static Graph statusGraph = 0;
static void lexDetails(void);

static MENUOPT statusMenu[]={
        {graphDestroy,"Quit"},
        {help,"Help"},
	{graphPrint, "Print"},
	{acedbstatus,"Update"},
	{lexDetails,"Number of Objects per class"},
#if !defined(MACINTOSH)
	{arrayDoShow, "Array report on stderr"},
	{arrayDoMark, "Array Mark"},
#endif
#ifdef JUNK /* need to recompile blocksubs with define TESTBLOCK */
        {blockshow,"About Cache 1"},
#endif
	{cacheShow, "About Cache 2"},
        {0,0} };

/**********************************************************/

static void lexDetails(void)
{ int i , line = 3 ;
  graphClear() ;
  graphPop() ;
  graphTextBounds (100,10) ;
  for(i=0; i<256;i++)
    {
      switch (pickType(KEYMAKE(i,0)))
	{
	case 'B':
	  if (lexMax(i) < 3)
	    continue ;
	  break ;
	default: 
	  if (lexMax(i) < 2)
	    continue ;
	  break ;
	}
      line++ ;
      graphText("Class :",4, line) ;
      graphText(pickClass2Word(i),12, line) ;
      graphText(messprintf("%6d entries    %8d char  %8d hashed keys",
			   lexMax(i), vocMax(i), lexHashMax (i)), 25, line) ;
    }
  graphTextBounds (100,(line+5)) ;
  graphRedraw() ;
}

/******************************************/

static void statusDestroy(void)
{
  statusGraph = 0 ;
}

/**********************************************************/
extern int NII, NIC, NWC, NWF, NWG  ;
extern int blockMax (void) ;

void acedbstatus (void)
{
 unsigned long keynum,keynum2, vocnum,vocspace,hspace,tspace ;
 mysize_t  aUsed, aMade, aAlloc, aReal ;

 unsigned long int
   dused, dfree, plus, minus, ndread, ndwrite ;
 int
   n,
   bsused, bsalloc, btused, btalloc,bsmemory,
   nmessalloc, messalloctot,
   bused, bpinned, bfree, bmodif,
   cacheused, cachel, cachek, cachemod,   
   nc1r, nc1w, ncc1r, ncc1w,naread,nawrite,
   nc2read,nc2write,nccr,nccw,
   level ;
 int line ;
 char buffer [2] ;
#ifdef __MSDOS__
 unsigned long corelft=coreleft() ;
 unsigned int stacklft = _stklen - stackused();
#endif

 diskavail (&dfree,&dused, &plus, &minus, &ndread, &ndwrite) ;
 lexavail (&vocnum,&keynum,&keynum2,&vocspace, &hspace,&tspace);
 BSstatus (&bsused, &bsalloc, &n) ; bsmemory = n ;
 BTstatus (&btused, &btalloc, &n) ; bsmemory += n ;
 arrayStatus (&aMade, &aUsed, &aAlloc, &aReal) ;
 aStatus (&naread, &nawrite) ;
 nmessalloc = messAllocStatus (&messalloctot) ;
 blockavail(&bused,&bpinned,&bfree,&bmodif) ; 
 blockStatus (&nc1r, &nc1w, &ncc1r, &ncc1w) ;
 cacheStatus(&cacheused,&cachel,&cachek,&cachemod,&nc2read,&nc2write,&nccr,&nccw);
 if(graphActivate(statusGraph))
   { graphPop() ;
     graphClear();
     graphTextBounds (100,15) ;
   }
 else
   { statusGraph =  displayCreate(DtStatus) ;
     graphMenu (statusMenu);
     graphRegister (DESTROY, statusDestroy) ;
     graphTextBounds (100,15) ;
     graphColor (BLACK) ;
   }
 line = 1 ;


 graphText (messprintf ("Program %s, %s", messGetErrorProgram(), utAppGetCompileDate()), 2, line++) ;
 graphText (messprintf ("Using %s, %s", aceGetVersionString(), aceGetLinkDateString()), 2, line++) ;
 line++ ;
 graphText (messprintf ("Data directory %s, release %d-%d",
			sessionFilName ("", 0, 0),
			thisSession.mainDataRelease,
			thisSession.subDataRelease), 2, line++) ;
 line++ ;
 graphText(messprintf("Session : %d, User %s,  Write Access %s ",
		      thisSession.session, thisSession.name,
		      (isWriteAccess() ? "Yes" : "No" )  ) ,
                      12,line++);
 graphText(messprintf
	   ("Global Address %d,  Index: %d", 
	    thisSession.gAddress, bIndexVersion(-1)), 12, line++) ;
 
 level = freesettext(buffer,"") ;
 graphText(messprintf("Stream level %d", level),12,line++) ;
 freeclose(level) ;
 line++ ;
 graphText
   (messprintf
    ("Disk: %ld blocks used,  %ld blocks available.",
                     dused,dfree),2,line++) ;
 graphText
   (messprintf
    ("        : %ld blocks allocated and %ld freed in this session.",
                     plus, minus),2,line++);
 line++ ;
 graphText(messprintf(
   "Lexiques: %lu classes, %lu kb allocated", vocnum,tspace),2,line++);
 graphText(messprintf(
   " %lu keys, %lu in lex2, %lu kb of names %lu hashing keys",
          keynum,keynum2,vocspace/1024, hspace),3,line++);
 graphText(messprintf("Hash info: nii = %d, nic = %d, nwc = %d, nwf = %d, nwg = %d",
		       NII  , NIC  , NWC  , NWF  , NWG) , 2, line++) ;

 line++ ;
 blockavail(&bused,&bpinned,&bfree,&bmodif);

 graphText (messprintf("Cache 1: %d kb", blockMax()), 2, line++) ;
  graphText(messprintf(
  "usage : %d obj served, %d saved, %d read from disk, %d saved to disk",
               nc1r, nc1w, ncc1r, ncc1w),
            4,line++);
 graphText(messprintf(
  "Current content: %d blocks used, %d modified, %d locked, %d free",
               bused,bmodif,bpinned,bfree),
            4, line++) ;
 graphText(messprintf(
		   "Usage : %d obj served, %d saved, %d read from disk, %d saved to disk",
		   nc1r, nc1w, ncc1r, ncc1w),
	4,line++);

 graphText(messprintf("Cache 2: %d slots.%d kbytes allocated", cacheused,bsmemory/1024),2,line++);
 graphText(messprintf("Current content: %d obj locked, %d obj known, %d mod.",
                     cachel,cachek,cachemod), 4, line++);
 graphText(messprintf(
		   "Usage : %d obj served, %d saved, %d read from cache1, %d saved to cache1",
		   nc2read, nc2write, nccr, nccw),
	4,line++);
 graphText(messprintf(
    "%d/%d BS, %d/%d BT nodes used/allocated, total allocation: %d kbytes",
		      bsused,bsalloc,btused,btalloc,bsmemory/1024), 4, line++) ;

 line++ ;
 graphText("Cache statistics: each layer serves the known objects or asks the lower layer", 2, line++) ;
 graphText(messprintf ("Array: %d read from cache1, %d saved to cache1", naread, nawrite), 4, line++) ;
 graphText(messprintf ("Cache2: %d served, %d saved; %d read from cache1, %d saved to cache1", nc2read, nc2write, nccr, nccw), 4, line++) ;
 graphText(messprintf ("Cache1: %d served, %d saved; %d read from disk, %d saved to disk", nc1r, nc1w, ncc1r, ncc1w), 4, line++) ;
 graphText(messprintf ("Disk blocks: %d served, %d saved", ndread, ndwrite), 4, line++) ;

 line++ ;
 graphText(messprintf
   ("Arrays (Including the lexiques):"),
	   2, line++) ;
 graphText(messprintf
   ("%d made, %d current, %d Kb allocated, %d Kb used (0: unknown)",
		      aMade, aUsed, aAlloc/1024, aReal/1024),
	   4, line++) ;

 line ++ ;
 graphText ("Overall Memory Usage (caches 1 & 2 + arrays + various)", 2, line++) ;
 if (messalloctot)
   graphText(messprintf("Messalloc: %d blocks, %d kb allocated", 
			nmessalloc, messalloctot/1024),
	     4, line++);
 else
   graphText(messprintf("Messalloc: %d blocks", nmessalloc), 5, line++);
 line++ ;

#ifdef __MSDOS__
 line ++ ;
 graphText(messprintf("Available memory  %lu ", corelft/1024),
	   2,line++);
 graphText(messprintf("Available Stack   %u", stacklft), 2, line++);
#endif
  
 graphTextBounds (100,++line) ;
 graphRedraw() ;
}
 
/*********************************************/
/*********************************************/

