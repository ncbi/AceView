/*  File: banner.c
 *  Author: Extracted from xacemain.c by R. Bruskiewich
 *          xacemain.c by Jean Thierry-Mieg (mieg@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1991
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@crbm.cnrs-mop.fr
 *
 * Description: extern String array for introductory banner text
 *              which can be modified and used in a flexible manner
 * SCCS: @(#)banner.c	1.2 10/7/96
 * HISTORY:
 * Last edited: Dec 21 10:16 1998 (edgrif)
 * * Dec  1 15:39 1998 (edgrif): Increase intial array size for bannerMainStrings,
 *              change title of banner to show both acedb version/build date &
 *              calling applications build date.
 * * Nov 26 10:19 1998 (fw): 
 *  -   thisSession defined locally and exported from DLL
 * * Jun  3 19:28 1996 (rd)
 * Created: Jun 28 00:20 1995 (rbrusk)
 *-------------------------------------------------------------------
 */

#include "acedb.h"
#include "aceversion.h"
#include "session.h"		/* for thisSession.subDataRelease */
#include "version.h"				  /* For getting program compile date. */



#ifdef RD_PROPOSAL_9812

/*------- These comment lines stop the compiler moaning... ---------------
From rd Thu Dec 17 16:00:32 1998
To: Danielle et jean Thierry-Mieg <mieg@crbm.cnrs-mop.fr>
Subject: Re:  banner

I think in trying to be specific the banner tends to become unbalanced and out
of date.  Also it is too long to put on stdout - it takes up more than one
screen.  I would prefer something 10 lines long.  However, the last thing I want
to do is get into a debate with you about attributing everyone correctly.  What
about:

****** xace - ACEDB Version 4_6, Release INTERNAL, Link date Nov 2 1998 ******
Developers include: Jean Thierry-Mieg (mieg@ncbi.nlm.nih.gov), Richard Durbin
(rd@sanger.ac.uk), Friedemann Wobus, Ed Griffiths, Michel Potdevin, Danielle
Thierry-Mieg, Lincoln Stein, Simon Kelley, Richard Bruskiewich, Sam Cartinhour,
Ian Longden, Erik Sonnhammer, Phil Green, Otto Ritter, Frank Eeckman, Detlef
Wolf, Cyrus Harmon, John McCarthy, Martin Senger, Petr Kocab, John Morris, Gary
Aochi, Mike Holman.
Available under Gnu Public License via http://www.sanger.ac.uk/Software/Acedb/ 
*******************************************************************************

Add anyone you want.  I guess we should continue to have two separate lists for
the graphical and non-graphical (kernel only) versions.  If you want to put up a
mirror of our site, or something comparable, or add anything to it, please
suggest.

In case this upsets anyone there, dont take it as seriously thought out - I am
sleep deprived and have lots of other things going on.
------------------------------------------------------------------------*/

#endif


/* Somehow this seemed to have disappeared from this file ? I've put it back.*/
#define ADD(z) array(result, arrayMax(result), char*) = strnew (z,0)



ACEDB_FUNC_DEF Array bannerMainStrings (char *program, BOOL isGraphical, BOOL isAcembly)
{ 
  Array result = arrayCreate (64, char*) ;

  /* Header to show application compile date and acedb version.              */
  ADD(messprintf ("\n**** Program %s,  %s ****",
		  program, utAppGetCompileDate())) ;

  ADD(messprintf ("**** Using  %s,  %s ****\n",
		  aceGetVersionString(), 
		  aceGetLinkDateString())) ;


  if (isAcembly)
    {
      ADD("** Acembly: Sequence Assembly and Edition system") ;
      ADD("   by  Danielle et Jean Thierry-Mieg") ;
      ADD("       (CNRS, France) mieg@ncbi.nlm.nih.gov") ;
      ADD("   Built over the acedb data manager, developed by ") ;
      ADD("     Thierry-Mieg and Durbin (1989-...)") ;
    }
  else
    {
      ADD("Code by: Jean Thierry-Mieg (CNRS, France) mieg@ncbi.nlm.nih.gov") ;
      ADD("         Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk") ;
    }

#if defined(WIN32)
  ADD("") ;
  ADD("Windows-32 adaptation by:") ;
  ADD("   Richard Bruskiewich (UBC, Canada) rbrusk@octogene.medgen.ubc.ca") ;
  ADD("   funded by the Canadian Genome Analysis and Technology") ;
  ADD("   Program (CGAT) of the Medical Research Council of Canada.") ;
  ADD("") ;
#endif

  if (isGraphical)
    {
      ADD("with ideas and code contributed from many laboratories including:") ;
      ADD("   Montpellier:") ;
      ADD("\t Danielle Thierry-Mieg on the user interface") ;
      ADD("\t Michel Potdevin on more-info, biblio and debugging") ;
      ADD("   NCBI:") ;
      ADD("\t Mark Sienkiewicz on AceC") ;
      ADD("   Sanger:") ;
      ADD("\t Simon Kelley (Sanger Centre, UK) srk@sanger.ac.uk") ;
      ADD("\t Erik Sonnhammer for the blixem and dotter packages,") ;
      ADD("\t Ed Griffiths on code design") ;
      ADD("\t Friedemann Wobus on AQL and the graphic package") ;
      ADD("\t Ian Longden on graphics applications") ;
      ADD("   Cold Spring Harbour:") ;
      ADD("\t Lincoln Stein on AcePerl,") ;
      ADD("   DKFZ:") ;
      ADD("\t Otto Ritter on large data sets and the attach mechanism,") ;
      ADD("\t Detlef Wolf for quality control") ;
      ADD("\t Petr Kocab and Martin Senger for rpc layer;") ;
      ADD("   LBL:") ;
      ADD("\t Franck Eeckman for the Mac port ") ;
      ADD("\t John McCarthy on metadata and documentation,") ;
      ADD("\t Gary Aochi for the query interface") ;
      ADD("   UC Berkeley Drosophila Genome group: ") ;
      ADD("\t Cyrus Harmon also for the Mac port ") ;
      ADD("\t Suzana Lewis on the graphic toolbox") ;
#if !defined(WIN32)
      ADD("University of British Columbia, Canada:") ;
      ADD("\t Richard Bruskiewich on the Microsoft Windows port.") ;
      ADD("") ;
#endif
      ADD("   and also:") ;
      ADD("\t Phil Green (genefinder), Jaime Prilusky (Mosaic help),") ;
      ADD("\t Sam Cartinhour and John Morris on documentation.") ;
      ADD(" This list is only tentative, we apologise for omissions, and wish to ") ;
      ADD(" stress that these people and others have also helped fix all sorts") ;
      ADD(" of details that cannot be listed individually. ") ;
      ADD(" Thank you et merci.") ;
    }

  ADD("") ;
  ADD("You may redistribute this program and database subject to the") ;
  ADD("conditions in the accompanying copyright file.  Anyone interested in") ;
  ADD("maintaining an up-to-date version should contact one of the authors") ;
  ADD("at the above email addresses.") ;
  ADD("") ;

  return result ;
} /* bannerMainStrings */

/************************************************************/

ACEDB_FUNC_DEF Array bannerDataInfoStrings (void)
     /* must be called after sessionInit, so thisSession is set */
{ 
  Array result = arrayCreate (32, char*) ;
  FILE *fil ;
  
  if (filName ("wspec/datainfo","wrm","r") &&
      (fil = filopen("wspec/datainfo","wrm","r")))
    {
    int level = freesetfile (fil, messprintf ("%d %d %d", 
					      aceGetVersion(), aceGetRelease(),
					      thisSession.subDataRelease)) ;
    while (freecard(level))
      ADD(freepos()) ;
    }

  return result ;
} /* bannerDataInfoStrings */

/************************************************************/

ACEDB_FUNC_DEF void bannerWrite (Array strings)
{ 
  int i ;
  char *banner_line;

    for (i = 0 ; i < arrayMax(strings) ; ++i)
      {
	banner_line = arr(strings, i, char*);

	if (!getenv("ACEDB_NO_BANNER"))
	  printf ("%s\n", banner_line) ;
	messfree (banner_line);
      }

  arrayDestroy (strings) ;
} /* bannerWrite */

/**************** end of file *****************/
