/*  File: bql_.h
 *  Author: Richard Durbin (rd@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1996
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
 * 	Richard Durbin (MRC LMB, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * SCCS: $Id: bql.h,v 1.10 2020/06/21 12:55:24 mieg Exp $
 * Description: bql is a new lighter implementation of aql
 *   suppressing many complexities of the original implementation
 *   removing the dependency on lex/yacc
 *   unifyilng aql and table maker
 *     Not all the refinments of the AQL language have been kept
 *   but it is doubtful that they ever worked and certainly did not
 *   follow the available documentation

 * Exported functions:
 * HISTORY:
 * Last edited: 
 * Created: Fri Mar 18 13:46:46 EDT 2016 (mieg)
 *-------------------------------------------------------------------
 */

#ifndef BQL_H
#define BQL_H

#include "ac.h"

/****************************************************/
/* public data structure for an BQL query object :- */
/****************************************************/

typedef struct bqlStruct BQL;	/* opaque type for public interface */
typedef struct bqlIterStruct BQL_ITER;	/* opaque type for public interface */

/********************************/
/* public exported functions :- */
/********************************/

/******************************************************************/
/*                    bqlCreate                                   */
/*                                                                */
/* (1) create an BQL object on the objectHandle (if given)        */
/* (2) parse the query string to set up query object (set error)  */
/* (3) pre-process the query to check semantics	(set error)       */
/*                                                                */
/*  dump_out - if non-Null, is the stream on which to output the  */
/*             resulting table if no sorting is required. If this */
/*             is done (as the query progresses) no output table  */
/*             is produced. This is so that we can output very    */
/*             large tables without having to accumulate them in  */
/*             memory.                                            */
/*                                                                */
/*  beauty - an output control character, as used in the acedb    */
/*           command loop. This is used only if dump_out is       */
/*           being written to directly.                           */
/*                                                                */
/* debugLevel - activates debug output during parsing/execution - */
/*    0 - no debug info (use it for all except test programs)     */
/*    1 - only the parse-structure output                         */
/*    2 - all intermediate parsetrees                             */
/*    3 - all info during query execution                         */
/*                                                                */
/******************************************************************/

BQL *bqlCreate (BOOL debug, AC_HANDLE h)  ;
BQL *bqlDbCreate (AC_DB db, BOOL debug, AC_HANDLE h)  ;
BOOL bqlParse (BQL *bql, const char *query, BOOL acedbQuery) ;
int bqlRun (BQL *bql, KEYSET activeKeyset, KEYSET resultKeyset) ; /* returns number of lines */
void bqlMaxLine (BQL *bql, int maxLine) ; /* will block the calculation of the top iteration */

BOOL bqlExport (BQL *bql, char beauty, BOOL showTitle, int beginLine, int maxLine, ACEOUT ao) ;
AC_TABLE bqlResults (BQL *bql) ;

BQL_ITER *bqlIterCreate (const char *bqlQuery, KEYSET ksIn, const char **errors, AC_HANDLE h) ;

/* To manipulate large tables you may optionally
 * process table one KEY at a time
 * 1: iter = bqlIterCreate ("bql query", ksIn, &errors, h) ;
 * 2: while ((mini_table = bqlIterTable (iter)) 
 *        { do_something_with_the_mini_table ; }
 *     The mini table corresponds to the next key producing some data
 *     returns NULL when no more key produces a non empty mini_table
 *     DO NOT free the mini_table, it is freed internally on each iteration
 * 3: Call ac_free (iter) OR ac_free (h) ;
 */
BQL_ITER *bqlIterCreate (const char *bqlQuery, KEYSET ksIn, const char **errors, AC_HANDLE h) ;
AC_TABLE bqlIterTable (BQL_ITER *bqlIter) ;

const char *bqlError (BQL *bql) ;
/******************************************************************/
/*                BQL Errorhandling                               */
/*                                                                */
/* After calling bqlCreate and bqlExecute, one can check if an    */
/* error ocurred during processing, and take appropriate action.  */
/* Note: both those strings might consist of multiple lines       */
/*                                                                */
/******************************************************************/

/************************* end of file ****************************/

#endif /* BQL_H */

 
 
