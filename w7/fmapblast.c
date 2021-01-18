 /*  File: fmapblast.c
 *  Author: Richard Durbin (rd@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1992
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: Interface to BLAST, Oligo selection program
 * HISTORY:
 * Last edited: Dec 11 16:42 1998 (fw)
 * * Jul 16 10:09 1998 (edgrif): Introduce private header fmap_.h
 * Created: Tue Oct 7 1996 (mieg)
 *-------------------------------------------------------------------
 */

/* $Id: fmapblast.c,v 1.2 2020/05/30 16:50:33 mieg Exp $ */

#include "fmap_.h"

#include "display.h"
#include "call.h"
#include "dna.h"
#include "systags.h"
#include "query.h"
#include "session.h"
#include "parse.h"

static Graph blastControlGraph = 0 ;
static void fMapBlastControl (void) ;
static void fMapBlastControlDraw (void) ;
static LOOK BLASTGET (char *title) ;
static void popBlastControl(void) ;
static void fMapBlastX_pir (void) ;
static void fMapBlastX_nr (void) ;
static void fMapBlastN (void) ;

static BOOL running = FALSE ;
static int MSPupper = 50 ; /* limit for conserved hits, see /usr/local/bin/MSPcrunch -C help */

/***************************************************************************************/

void fMapBlastInit (void)
{
  /*   LOOK look = 0 ; 

    if (!graphAssFind ((void*)MAP_LOOK, &look) ||
      look->magic != MAGIC) 
    return ;
    */
  fMapBlastControl () ;
}

/***********************************************************************/
/***************** actual calls to blast *******************************/
/***********************************************************************/

/* i just need to parse the returned file 
   and the look buffer which contains the alignment info 
*/
static void fMapBlastEnd (FILE *f, void *lk)
{
  KEY key = 0 ;
  Stack s = (Stack) lk ;
  char *cp ;

  sessionGainWriteAccess () ;
  if (f && stackExists (s))
    { 
      stackCursor (s, 0) ;
      cp = stackNextText (s) ;
      lexaddkey (cp, &key, _VSequence) ;
      cp = stackNextText (s) ;
      parseBuffer (cp, 0) ;
      
      parseKeepGoing = TRUE ;	 
      parseFile (f, 0, 0) ;  
      parseFile (f, 0, 0) ; /* do it twice because alias system is bugged */
      parseKeepGoing = FALSE ;
    }
  stackDestroy (s) ;
  if ( fMapActivateGraph())
    { graphUnMessage() ;
      displayPreserve () ;
    }
  if (key) display (key, 0, 0) ;
  running = FALSE ;  
  popBlastControl() ;
  fMapBlastControlDraw () ;
}

static void fMapBlastAny (char *pgm, char *db)
{
  Array dna, colors ; 
  Stack lk = 0 ;
  KEY seqKey = 0, father = 0 ;
  int from, to, origin, i ;
  char *tmpName ;
  FILE *fil ;
  LOOK look = BLASTGET("fMapBlastAny") ;
   
  if (!look) return ;
  if (!fMapActive (&dna, &colors, &seqKey, 0))
    { messout ("First select a Sequence window") ;
    return ;
    }
  if (!graphCheckEditors(blastControlGraph, TRUE))
    return ;
  fMapFindZone (look, &from, &to, &origin) ;
  
  if (!(fil =  filtmpopen (&tmpName, "w")))
    return ;
  fMapFindZoneFather (look, from, to, &father, &origin) ;
  dnaDumpFastA (dna, from, to-1, 
		messprintf ("%s_%d-%d", name(father), 
			    origin+1, origin + to - from), 
		fil, 0) ;
  filclose (fil) ;
  lk = stackCreate (50) ;
  pushText (lk,	messprintf ("%s_%d-%d", name(father), 
			    origin+1, origin + to - from)) ;
  pushText(lk,  messprintf(
   "Sequence %s_%d-%d\nSource %s\n\nSequence %s\nSubsequence %s_%d-%d %d %d \n\n",
   name(father), origin+1, origin + to - from, name(father), 
   name(father), name(father),  origin+1, origin + to - from, 
   origin+1, origin + to - from)) ;

  i = stackMark (lk) ;
  pushText (lk, messprintf("%s %s %s %d",tmpName, pgm, db, MSPupper)) ;
  if ( externalAsynchroneCommand 
       ("blast_search", stackText (lk, i), lk, fMapBlastEnd))
    { displayPreserve() ;
      running = TRUE ;
      popBlastControl() ;
      fMapBlastControlDraw () ;
    }
}

static void fMapBlastX_pir (void)
{ fMapBlastAny ("blastx", "pir") ;
}

static void fMapBlastX_nr (void)
{ fMapBlastAny ("blastx", "nr") ;
}

static void fMapBlastN (void)
{ fMapBlastAny ("blastn", "nr") ;
}

/**************************************************************************/
#ifdef JUNK
static void fMapBlastMail (void)
{
  Array dna, colors ; 
  void  *look ;
  KEY	seqKey ;
  Stack s = stackCreate(50) ;
  OBJ obj ;
  KEY titleKey = 0 ;
  int from, to, origin, i ;
  char buf[2] ;
  DNACPTGET("fMapBlastMail") ;
  
  if(!fMapActive(&dna,&colors,&seqKey, &look))
    { messout("First select a dna window or the Search Active KeySet button") ;
      return ;
    }
  fMapFindZone (look, &from, &to, &origin) ;
  
  if ((obj = bsCreate(seqKey)))
    { 
      bsGetKey (obj, _Title, &titleKey) ;
      bsDestroy (obj) ;
    }
 
  pushText(s, titleKey ? name(titleKey) : "acedb_submission") ;
  catText(s,"  ") ;
  
  buf[1] = 0 ;
  for (i = from ; i < to ;)
    { buf[0] = dnaDecodeChar[((int)arr(dna, i++, char)) & 0xff] ;
      catText(s,buf) ;
    }
  catText(s," &\n") ;
 
  if (!callScript ("blast_mailer", stackText(s,0)))
    messout ("I sent your request to the blast mailer") ;
  else
    messout ("Failed to find script to mail to blast") ;
  stackDestroy(s) ;
}

/**************************************************************************/

static void dnacptFastamailMail (void)
{
  Array dna, colors ; 
  void  *look ;
  KEY	seqKey ;
  Stack s = stackCreate(50) ;
  OBJ obj ;
  KEY titleKey = 0 ;
  int from, to, origin, i ;
  char buf[2] ;
  DNACPTGET("dnacptFastaMail") ;
  

  if(!fMapActive(&dna,&colors,&seqKey, &look))
    { messout("First select a dna window or the Search Active KeySet button") ;
      return ;
    }
  fMapFindZone (look, &from, &to, &origin) ;
  
  if ((obj = bsCreate(seqKey)))
    { 
      bsGetKey (obj, _Title, &titleKey) ;
      bsDestroy (obj) ;
    }
 
  pushText(s, titleKey ? name(titleKey) : "acedb_submission") ;
  catText(s,"  ") ;
  buf[1] = 0 ;
  for (i = from ; i < to ;)
    { buf[0] = dnaDecodeChar[((int)arr(dna, i++, char)) & 0xff] ;
      catText(s,buf) ;
    }
  catText(s,"  &\n") ;
 
  if (!callScript ("fastamail_mailer", stackText(s,0)))
    messout ("I sent your request to the fasta mailer") ;
  else
    messout ("failed to find script to mail to fasta") ;
  stackDestroy(s) ;
}

#endif /* JUNK */

/***********************************************************************/
/****************** Blast control graph ********************************/
/***********************************************************************/

static int selecting = 0 ;

static void popBlastControl(void)
{
  if (graphActivate (blastControlGraph))
    graphPop() ;
}

static LOOK BLASTGET (char *title)
{
  LOOK look ; 
  void *vp ;
  
  if (!graphActivate (blastControlGraph))
    { blastControlGraph = 0 ; return 0 ; }
  graphCheckEditors(blastControlGraph, TRUE) ;

  if (selecting == 1) selecting = 0 ;
  fMapBlastControlDraw () ;
  if (!fMapActive (0,0,0,&vp)) 
    { messout("Sorry, no active sequence, I cannot proceed") ;
      return 0 ;
    }
  graphUnMessage () ;
  look = (LOOK) vp ;
  graphActivate (look->graph) ;

  return look ;
}

/********************
commented out, seems useless
static void blastRedraw (void)
{
  LOOK look = BLASTGET ("blastRedraw") ;
  if (!look) return ;

  fMapDraw(look,0) ;   
  popBlastControl() ;
  fMapBlastControlDraw () ;
}


******************************************************/

static MENUOPT blastControlMenu[] =
{ { graphDestroy, "Quit" },
  { help, "Help" },
  { 0, 0 }
} ;

/*****************************************************************/

static BOOL MSPlimitcheck (int n)
{ return (n <= 80 && n >= 30) ;
}

static void fMapBlastControlDraw (void)
{ 
  int line = 3;
  Graph old = graphActive () ;


  if (!graphActivate (blastControlGraph))
    return ;
  graphClear () ;
  graphText ("A direct interface to BLAST",
	     2, line++) ;
  graphText (" the homology search program",
	     2, line++) ;
  graphText (" from the NCBI",
	     2, line++) ;
  line++ ;

  if (running)
    {
      graphText ("I sent your request to blast", 3, line++) ; 
      graphText ("on the server defined in the file $ACEDB/wscripts/blast-search", 3, line++) ; 
      line += 2 ;
      graphText ("It takes around 5 minutes to search a 5 kb", 3, line++) ;
      graphText ("fragment against the nr database", 3, line++) ;
      graphText ("During the search you can keep working", 3, line++) ;
      graphText ("but you should not edit the searched sequence", 3, line++) ;
      line += 2 ;
    }
  else
    {
      graphText ("It takes around 5 minutes to search a 5 kb", 3, line++) ;
      graphText ("fragment against the nr database", 3, line++) ;
      graphText ("During the search you can keep working", 3, line++) ;
      graphText ("but you should not edit the searched sequence", 3, line++) ;

      line++ ;
      graphText ("If all works", 3, line++) ;
      graphText ("   a new window will pop out with the results", 3, line++) ;
      graphText ("If not", 3, line++) ;
      graphText ("   check the script wscripts/blast_search", 3, line++) ;
      graphText ("This is a new tool, please send comments", 3, line++) ;
      line += 3 ;
      graphButton ("BlastX/PIR:", fMapBlastX_pir, 2,line) ;  /* line += 2.2 ; */
      graphButton ("BlastX/nr:", fMapBlastX_nr, 15,line) ;  /* line += 2.2 ; */
      graphButton ("BlastN/nr:", fMapBlastN, 30,line) ;  line += 2.2 ;
      graphIntEditor ("Minimal score [30-80]", &MSPupper, 18, line++, MSPlimitcheck) ; 
      graphText ("Compare the 6-frame translations against", 4,line++) ; 
      graphText ("the protein database PIR", 4,line) ; 
      line += 2.2 ;
      /*
      graphButton("tBlastX: translate 6 frames and blast against the 6 translations of the dna database nr", fMapBlastN, 2,line) ; line += 2.2 ;  
      */
    }

  graphButton ("Help", help, 3,line) ; 
  graphButton ("Comments", acedbMailComments, 13,line) ; line += 2 ;
  if (!running) graphButton ("Cancel", graphDestroy, 3,line) ;
  graphRedraw () ;
  graphActivate (old) ;

}

static void fMapBlastControlDestroy (void)
{
  blastControlGraph = 0 ;
}

static void fMapBlastControl (void)
{
  Graph old = graphActive () ;
  
  if (graphActivate (blastControlGraph))
    { graphPop () ;
    return ;
    }
  
  blastControlGraph = graphCreate (TEXT_SCROLL, "BLAST parameters", 
				 0, 0, 0.6, 0.5) ;
  graphHelp("BLAST") ;
  graphRegister (DESTROY, fMapBlastControlDestroy) ;
  graphMenu (blastControlMenu) ;
  
  fMapBlastControlDraw () ;
  graphActivate (old) ;
}

/***************************************************************************************/
/***************************************************************************************/
