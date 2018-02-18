/*  File: query.h
 *  Author: Jean Thierry-Mieg (mieg@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1999
 * -------------------------------------------------------------------
 * Acedb is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 * or see the on-line version at http://www.gnu.org/copyleft/gpl.txt
 * -------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@crbm.cnrs-mop.fr
 *
 * Description: include for the query-submodule of the ace-kernel
 * Exported functions:
 *              all public functions from queryexe.c
 * HISTORY:
 * Last edited: Aug  4 09:10 2000 (edgrif)
 * Created: Thu Feb  4 10:47:23 1999 (fw)
 *-------------------------------------------------------------------
 */

/* $Id: query.h,v 1.9 2014/05/28 22:30:24 mieg Exp $ */

#ifndef ACEDB_QUERY_H
#define ACEDB_QUERY_H
 
#include "acedb.h"
#include "bs.h"

#ifndef CONDITION_DEFINED
typedef void* COND ;
#endif

BOOL  condCheckSyntax (const char *cp) ;   /* checks the syntax of a query */
BOOL condConstruct (const char *text, COND* cdp) ;
void  condDoDestroy(COND cond) ;
#define condDestroy(_cond) ((_cond) ? (condDoDestroy(_cond),(_cond)=0,TRUE) : FALSE)

KEYSET queryParametrized (KEYSET oldSet, const char *text, const char *parms) ;
#define query(ks,t) queryParametrized(ks,t,"")

KEYSET queryLocalParametrized (KEYSET ks1, const char *text, const char *parms);

KEYSET queryKeyParametrized  (KEY key, const char *text, const char *parms) ;
#define queryKey(k,t) queryKeyParametrized(k,t,"")

KEYSET queryKeyTransitiveClosure (KEY key, KEY tag)  ;

KEYSET queryGrep (KEYSET active, const char *text) ;
void queryCreate (void) ;

    /* Perfoms query on single obj, if cp == 0, reuses preceeding condition */
/* force localization in obj */
BOOL queryFindLocalise (OBJ obj, KEY key, const char *cp) ; 
BOOL queryFindLocalise3 (COND cond, OBJ *objp, KEY key) ;  

/* do not force localization in obj */
BOOL queryFind3 (COND cond, OBJ *objp, KEY key) ;
BOOL queryFind21 (OBJ *objp, KEY key, const char *cp, const char *tname) ;  /* set name of Text key */
BOOL queryFind31 (COND cond, OBJ *objp, KEY key, const char *tname) ;

BOOL queryIsAInit(Stack s, Array a) ;

BOOL queryCheckConstraints (OBJ obj, KEY key, void *v) ;
void queryClearConstraints (void) ;
BOOL queryConstraintsInit (const char *text, int classe) ;
BOOL queryCalcul (const char *buf, int *ip, float *fp) ; /* buf = [arith expression], ip or fp must be non zero */

/* specific function for acembly */
typedef KEYSET (*WEBQUERYFUNC) (KEYSET ks0, char *text) ;

/* REGEXP package */
/* 0: bad query, !NULL: use in queryRegExpRun, then call queryRegExpDoFree */
typedef void QueryRegExp ;
QueryRegExp *queryRegExpCreate (const char *regExpPattern, BOOL getPos, AC_HANDLE h) ;
#define queryRegExpDestroy(_compiledQuery) messfree(_compiledQuery)
/* 0: not found, >0 position */
int queryRegExpFind (const char *data, QueryRegExp *vRegP, BOOL getPos) ;
/* all in one function */
int queryRegExpMatch (const char *data, const char *regExpPattern, BOOL getPos) ; 
   /* 0: not found, > 0 found [at position if getPos==TRUE], -1 pattern is wrong */

BOOL dnaSearch_textQuery_textDna (const char *query, const char *data) ;
KEYSET OWQProduct2Gene (KEYSET ks, BOOL force) ;

#endif
