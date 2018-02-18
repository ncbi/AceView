/*  File: sysclass.h
 *  Author: Ed Griffiths (edgrif@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1998
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.crbm.cnrs-mop.fr
 *
 * Description:
 *   Before editing this file, please contact us.
 *
 *	This file is read only at compile time
 *
 *    It holds the enumeration of the classes needed by the acedb kernel.
 *    The additional application classes are listed in wspec/classes.wrm
 *    but do not need a hard choice of their number
 *
 *    Each #define line correspond to a class
 *    These numbers must NOT be repeated nor modified
 *    They are explicitelly repeated in pickDefineSystemClasses.
 *    
 *    Class properties, formerly defined in this file, are
 *    now in wspec/options.wrm 
 *
 *    The total number of MainClasses cannot exceed 255
 *
 * Exported functions:
 * HISTORY:
 * Last edited: Oct 22 13:35 1998 (edgrif)
 * * Oct 22 13:33 1998 (edgrif): Added this comment header. Added function
 *              decs for the sysclass.c functions.
 * Created: Thu Oct 22 13:27:58 1998 (edgrif)
 *-------------------------------------------------------------------
 */

/* $Id: sysclass.h,v 1.5 2007/11/13 00:29:11 mieg Exp $ */



/* The kernel explicitelly relies on classes 0, 1, 2, 3 order during bootstrap */
/* 24: The main classes, used in KEYMAKE, only 256 allowed */
/* number 24 is chosen for backward compatibility to earlier code releases */

#ifndef SYSCLASS_DEF
#define SYSCLASS_DEF

#define _VSystem          0 
#define _VGlobal          1 
#define _VSession         2 
#define _VVoc             3 

#define _VDisplay        23
#define _VMainClasses    24
    

/* from here on, the class numbers can be chosen dynamically,
   we list here simply those class names which appear in the kernel code
   they are allocated and initialised dynamically in lexsubs.c

   IMPORTANT: 
   the files wspec/sysclass.wrm and wspec/classes.wrm are no longer needed
   However, in release 3-* we still search them for consistency with datasets
   created by earlier subreleases. 
*/

extern int
  _VBat, _VKeySet, _VCalcul, 
  _VClass,  /* the class hierarchy and class aliases */
  _VModel,  /* the model of all the main classes and the subtypes */
  _VText, _VLongText, _VImage, 
  _VComment, _VUserSession,   	/* must be 1 + _VComment, seems useless needs checking */
  _VQuery, _VConstraint, _VTable, _VTableResult , _VJade, _VFicheView, _VView,
  _VBinary ;


Stack sysModels (void) ;
void sysClassInit (void) ;
void sysClassOptions(void) ;	/* class options for system classes */
void sysClassDisplayTypes(void) ; /* display types for system classes */
void lexDefineSystemTags (void) ;
void classHardInit (char *titre) ;


/************ end of file **********/
#endif /* SYSCLASS_DEF */

