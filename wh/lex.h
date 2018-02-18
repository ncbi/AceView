/*  File: lex.h
 *  Author: Fred Wobus (fw@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1998
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@crbm.cnrs-mop.fr
 *
 * Description: public function header for lexsubs4.c
 *              for everything dealing with the lexique
 * HISTORY:
 * Last edited: Nov 23 18:00 1998 (fw)
 * Created: Mon Nov 23 12:06:50 1998 (fw)
 *-------------------------------------------------------------------
 */

/* $Id: lex.h,v 1.8 2007/03/21 21:03:14 mieg Exp $ */

#ifndef _LEX_H
#define _LEX_H

/************************************************************/

#include "disk.h"		/* for DISK type */

/************************************************************/
int   lexMax(int t) ;  /* Numbers of entries in the lexique */
int   lexi2Max(int t) ;  /* Numbers of cached entries of the t lexique */
int   vocMax(int t) ;
void  lex2clear (KEY key) ; /* call this function to clear the lex2 array */

int   iskey(KEY kk);
                      /* Returns 0 if kk is unknown,
                       * 1 if kk is just a vocabulary entry
                       * 2 if kk corresponds to a full object
                       */
BOOL  iskeyold(KEY kk);
                      /* TRUE if key is on disk or in first cache */
BOOL  lexReClass(KEY key,KEY *kp, int t);

BOOL  lexIsInClass (KEY key, KEY classe) ;  /* Valid also for sub classes */
BOOL  lexClassKey (char *cp,KEY *kp) ; /* Recognises class:name */
BOOL  lexword2key (const char *cp,KEY *kp, KEY classe);
                           /*given a word *cp, sets its key *kp */
                           /* search done in t lexique */
                           /*returns TRUE if found, FALSE if not */

BOOL  nextName(KEY key, const char **cpp)  ;
                           /* Runs along tha aliases of key,
			    * on first call *cpp should be null
			    */  
BOOL  lexAlias(KEY *keyp, const char *newName, BOOL doAsk, BOOL keepBothNames) ;
                 /* lexCleanUp newName and modify the lexique
		    return TRUE if newName is accepted
		    FALSE if error or newName corresponds to some
		    other object
		    */
int lexIsAlias (KEY key) ;  /* returns 0:no, 1:may_be, 2:yes, key is alias of another key */
KEY   lexAliasOf(KEY key) ;  /* Returns the end of the alias list */
BOOL lexDestroyAlias (KEY key); /* removes an aliased key and 
				   relinks the alias list */

BOOL  lexIsKeyVisible(KEY key);
/* returns false if key does not exist, has aliasstatus or emptystatus, 
   or is a KEYKEY 0 model */
/* cannot use lexGetStatus for this as it follows the alias list */

BOOL  lexaddkey(const char *cp,KEY *kp, KEY classe);
                  /*add to the t lexique the word *cp */
               /*returns TRUE if added, FALSE if known, crashes if Full */

char  *lexcleanup(const char *cp, AC_HANDLE h) ;   /* Remove leadind and trailing spaces */

int   lexstrcmp (const char *a, const char *b);
   /*returns 1 if *a>*b  or 0 or -1 as strcmp but is case unsensitive*/
   /*This is used to better the presentation, yet avoid doubles*/
BOOL lexNext(KEY classe, KEY *p); /*gives the next KEY in table n*/
KEY lexLastKey(int classe) ;

void lexkill(KEY key);      /*deallocates a key*/
void lexshow (void) ;   /* RMD 6/19/90, for debugging */
void lexInit(void);            /*Called at beginning of program*/
void lexRead(void);            /*Called at beginning of program*/
void lexSave(void) ;
void lexavail(unsigned long *vocnum,
	      unsigned long *lex1num,
	      unsigned long *lex2num,
	      unsigned long *vocspace,
              unsigned long *hspace,
	      unsigned long *totalspace
	      ) ;

#define LOCKSTATUS	(unsigned char)1  /* bits 1 and 2 are reset to 0 by lexiread */
#define CALCULSTATUS	(unsigned char)2  /* all bits are reset by lexAlias */
#define EMPTYSTATUS	(unsigned char)4
#define ALIASSTATUS	(unsigned char)8
#define TOUCHSTATUS	(unsigned char)16   /* key touched this session */
#define ALPHASTATUS	(unsigned char)32   /* reserved for the alphabetic sorter */
#define SERVERSTATUS    (unsigned char)64   /* reserved for xclient */
#define ISALIASSTATUS   (unsigned char)128  /* this key is the alias of some other key */
#define LEXPRIVATESTATUS (unsigned char) 9  /* bits 1 8 are protected in lexsetstatus */

void lexSetStatus(KEY key, unsigned char c) ;
void lexUnsetStatus(KEY key, unsigned char c) ;
void lexClearClassStatus(KEY classe, unsigned char c) ;
unsigned char lexGetStatus(KEY key) ;
void lexunlock(KEY key) ; /* these 2 check for double locking */
BOOL lexlock(KEY key) ;
#define lexiskeylocked(key) (lexGetStatus(key) & LOCKSTATUS)

void lexSetIsA(KEY key, unsigned char c);

void lexShutDown (void) ;

/************ in lexalpha.c ****************/

void lexAlphaMark (int c) ;
void lexAlphaMakeAll (void) ;
void lexAlphaSaveAll (void) ;
void lexAlphaClear (void) ;
int  lexMaxVisible (int c) ;

void lexAlphaShutDown (void) ;
DISK lexDisk(KEY key) ;
void lexmark(int tt) ; /* kernel use only */
KEYSET lexAlphaClassList (BOOL redo) ;

BOOL lexIsGoodName(char *name) ;

#endif /*  _LEX_H */

/***************************** eof ****************************/
 

