/*  File: colours.h
 *  Author: Fred Wobus (fw@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1998
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@crbm.cnrs-mop.fr
 *
 * Description: enum for colur definitions
 *              they should be avilable for even for non-graphical
 *              programs that don't include graph.h
 * Exported functions:
 *              enum Colour {..};
 * HISTORY:
 * Last edited: Nov 19 15:52 1998 (fw)
 * Created: Thu Nov 19 15:49:30 1998 (fw)
 *-------------------------------------------------------------------
 */

#ifndef _COLOURS_H
#define _COLOURS_H


/* ACEDB "TRANSPARENT" conflicts with a WIN32 definition... */
#if defined(WIN32) && defined(TRANSPARENT)
#define WIN32_TRANSPARENT TRANSPARENT  // for Windows CDC::SetBkMode() ? 
#undef TRANSPARENT  /* as defined in Windows: wingdi.h include or elsewhere */
#endif

  /* These colors match those declared in systags, they must appear in the same order */	
  /* These colors match those declared in systags, they must appear in the same order */
  /* open mainWindow->admin->palette (in w2) for a reasonable order */
  /* WHITE must be first */

enum Colour    {WHITE, BLACK, LIGHTGRAY, DARKGRAY,
                RED, GREEN, BLUE,
                YELLOW, CYAN, MAGENTA,
		LIGHTRED, LIGHTGREEN, LIGHTBLUE,
		DARKRED, DARKGREEN, DARKBLUE,
		PALERED, PALEGREEN, PALEBLUE,
		PALEYELLOW, PALECYAN, PALEMAGENTA,
		BROWN, ORANGE, PALEORANGE,
		PURPLE, VIOLET, PALEVIOLET,
		GRAY, PALEGRAY,
		CERISE, MIDBLUE,
                LIGHTMAGENTA, LIGHTCYAN, DARKVIOLET, LAVANDER, 
		
		BLUE1, BLUE2, BLUE3, BLUE4,
		BLUE5, BLUE6, BLUE7, BLUE8,

		GREEN1, GREEN2, GREEN3, GREEN4,
		GREEN5, GREEN6, GREEN7, GREEN8,

		RED1, RED2, RED3, RED4,
		RED5, RED6, RED7, RED8,

                LIGHTVIOLET,
		DARKCYAN,
		LIGHTORANGE,
		NUM_TRUECOLORS,
                TRANSPARENT,	/* pseudocolour only for background */
		FORECOLOR,	/* pseudocolor to force box->fcol after graphColor */
		BACKCOLOR	/* pseudocolor to force box->bcol after graphColor */
               } ;


#endif /* _COLOURS_H */
