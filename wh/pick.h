/*  File: pick.h
 *  Author: Ed Griffiths (edgrif@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1998
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.crbm.cnrs-mop.fr
 *
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Dec  4 14:48 1998 (fw)
 * * Oct 22 14:11 1998 (edgrif): Added this header, added dec for pickRegisterConstraints
 * Created: Thu Oct 22 14:11:20 1998 (edgrif)
 *-------------------------------------------------------------------
 */

/* $Id: pick.h,v 1.6 2007/03/26 22:48:18 mieg Exp $ */

#ifndef DEFINE_PICK_H
#define DEFINE_PICK_H

#include "regular.h"
#include "keyset.h"

typedef  struct { 
                 int 	name, alias ;
		 char	visible ;
		 char	type ;
                 KEY	displayKey ;
                 KEY	mapKey ; /* Arun's SYNMAP maptype */
		 BOOL	isCaseSensitive ;
			/*
			* isCaseSensitive means that objects in this class
			* have case sensitive names.  For most classes, the
			* object names are not case sensitive.
			*/
 		 BOOL 	updateNames;
			/*
			* updateNames means that whenever the name of the
			* object is written, we want to modify the name
			* to match the name specified.  This means that
			* an object named "arf" that is updated as "Arf"
			* would have it's name changed to use the new
			* capitalization.
			*/
		 BOOL	protected ;  /* do not parse or addkey */
		 unsigned char mask ;
		 KEY    model ;
		 Array	conditions ;
		 void*	constraints ;
		 Stack	filters ;
                 BOOL	Xref ;
                 BOOL	known ;
		 BOOL	sybase ; /* experimental sybase storage */
		 KEY    tag, symbol ;
		 float	width, height ;
                 BOOL   private ; /* private to kernel, do not dump, invisibleKey */
                 KEY    classe ;  /* corresponding entry in class Class */
	       } PICKLIST ;

extern PICKLIST pickList[] ;
extern FREEOPT pickVocList[] ;	/* in picksubs.c */

void pickPreInit (void);
void pickInit (void);


                           /* returns # of class Word */
                           /* used by parse and readmodels */
int pickWord2Class (const char *word);
char *pickClass2Word(int classe) ;
                            /* returns name of class */
                   /* used by lexaddkey to initialise the vocs */

int pickMakeClass (const char* cp);
void pickCheckConsistency(void);


#define pickXref(t)               (t<255 && pickList[(t) & 255].Xref)
#define pickVisible(t)               (t<255 && ace_lower(pickList[(t) & 255].visible) == 'v')

#define pickCheck(t)               (pickList[(t) & 255].check)
#define pickKnown(t)               (pickList[(t) & 255].known)

#define pickDisplayKey(key)      pickList[(class(key)) & 255 ].displayKey
#define pickMapKey(key)      pickList[(class(key)) & 255 ].mapKey
#define pickType(key)             pickList[(class(key)) & 255 ].type
#define pickModel(key)             pickList[(class(key)) & 255 ].model
#define pickCaseSensitive(key)    pickList[(class(key)) & 255 ].isCaseSensitive
#define pickExpandTag(key)        (KEYKEY(key) ?  \
				   pickList[(class(key)) & 255 ].tag : 0)

 /* kernel function */
void pickSetClassOption (const char *nom, char type,  char v, 
			 const char *disp, BOOL protected, BOOL private,
		  	 BOOL case_sensitive, BOOL update_names );


/* register size for next display create (now in whooks/class.c) */
void pickRememberDisplaySize (const char *display);

void pickRegisterConstraints (void) ;

BOOL pickIsA (KEY *classp, unsigned char *maskp) ;

int superClass (KEY classe) ;  /* class < 256 and class in ?Class both recognised */
Array pickIsComposite(KEY classe) ;

#endif

