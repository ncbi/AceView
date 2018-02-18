/*  File: aql_.h
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
 * SCCS: $Id: aql.h,v 1.7 2007/05/19 02:49:05 mieg Exp $
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Nov  9 10:54 1999 (fw)
 * * Aug  5 17:19 1998 (fw): made public AQL type opaque
 * * Aug  4 14:14 1998 (fw): completed public interface specs
 * Created: Mon Oct 21 23:09:29 1996 (rd)
 *-------------------------------------------------------------------
 */

#ifndef AQL_H
#define AQL_H

#include "table.h"
#include <stdarg.h>

/****************************************************/
/* public data structure for an AQL query object :- */
/****************************************************/

typedef struct AqlStruct *AQL;	/* opaque type for public interface */

/********************************/
/* public exported functions :- */
/********************************/

/******************************************************************/
/*                    aqlCreate                                   */
/*                                                                */
/* (1) create an AQL object on the objectHandle (if given)        */
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
/* The variable-length parameter list is a NULL-pointer           */
/* terminated list of (char*) that defines the names of the       */
/* context variables (table and scalar) that are going to be used */
/*     aql = aqlCreate(cp, debug, "@active", 0);                  */
/* Registers the table-variable name @active, the subsequent call */
/* to aqlExecute has to provide this name again, but then also    */
/* pass the actual value.                                         */
/*                                                                */
/******************************************************************/
AQL   aqlCreate (const char *queryText,
                 ACEOUT dump_out,
                 char beauty,
                 int parseDebugLevel,
                 int evalDebugLevel,
                 ...);

/* mieg: 2007_05_14 
 * passed with argc = -1, will loop as it used to
 * with argc>=0, looping stops after argc calls
 * i add this controled called becuase on 64 bits machine 
 * we crahs while looping on  va_arg(args, )
 */

AQL   aqlCreate2 (const char *queryText,
                 ACEOUT dump_out,
                 char beauty,
                 int parseDebugLevel,
                 int evalDebugLevel,
		  int argc,  
                 ...);



/******************************************************************/
/*                    aqlExecute                                  */
/*                                                                */
/* (1) runs the internal pre-processed query structure against    */
/*     the database (set errors if necessary)                     */
/* (2) create a table object on the tableHandle (if given)        */
/*     and returns the results in that table.                     */
/*     If no tableHandle is given the table has to be explicitly  */
/*     destroyed (using tableDestroy()) after its content is no   */
/*     longer needed.                                             */
/*                                                                */
/* If a non-NULL pointer to a variable of type KEYSET is passed   */
/* in, a keyset for the result is created (if appropriate) using  */
/* the handle supplied in the following argument.                 */
/*                                                                */
/* The variable-length parameter list is a NULL-pointer           */
/* terminated list of pairs of names and values that bind actual  */
/* values to the context variables, e.g.                          */
/*	results = aqlExecute(aql, th, "Table @active", lastT, 0); */
/*                                                                */
/* The different types of variables are bound in pairs like these */
/*         "Table @t", TABLE *t					  */
/*         "Int $i", int i					  */
/*         "Float $f", float i					  */
/*         "Text $s", char *s					  */
/*         "DateType $t", mytime_t t				  */
/*                                                                */
/* NOTE: tace uses tableCreateFromKeySet() to pass in a table     */
/*       called @active as the active keyset                      */
/******************************************************************/
TABLE* aqlExecute(AQL aqlObject,
                  AC_HANDLE tableHandle,
                  KEYSET *result_keyset,
                  AC_HANDLE result_keyset_handle,
                  ...);



/******************************************************************/
/*                    aqlDestroy                                  */
/*                                                                */
/* Cleans up the storage allocated inside the AQL object.         */
/* Then destroys the actual object itself.                        */
/*                                                                */
/******************************************************************/
void   aqlDestroy(AQL aqlObject);


/******************************************************************/
/*                AQL Errorhandling                               */
/*                                                                */
/* After calling aqlCreate and aqlExecute, one can check if an    */
/* error ocurred during processing, and take appropriate action.  */
/*                                                                */
/* aqlIsError         - true/false whether there was an error     */
/* aqlGetErrorNumber  - get the number of the error               */
/* aqlGetErrorMessage - short string with an error message        */
/* aqlGetErrorReport  - longer string with a pointer to the       */
/*                      erroneous part of the query string        */
/*                                                                */
/* Note: both those strings might consist of multiple lines       */
/*                                                                */
/******************************************************************/
BOOL  aqlIsError         (AQL aql);
int   aqlGetErrorNumber  (AQL aql);
char *aqlGetErrorMessage (AQL aql);
char *aqlGetErrorReport  (AQL aql);


/******************************************************************/
/*                    aqlTable                                    */
/*                                                                */
/* Single entry call to the AQL system, no AQL object is returned,*/
/* no context variables can be registered, the errors are written */
/* to error_out.                                                  */
/* Every memory taken up during processing is cleaned up.         */
/* The returned table object is allocated upon the passed handle, */
/* or if NULL upon handle0, it can be destroyed by tableDestroy() */
/* or by messfree on the original handle after the results are no */
/* longer needed.                                                 */
/*                                                                */
/******************************************************************/
TABLE* aqlTable (char *queryText, ACEOUT error_out, AC_HANDLE handle);


/******************************************************************/
/*                    aqlQuick                                    */
/*                                                                */
/* Even quicker access to the AQL functionality, the query is     */
/* processed and the results or any errors are written to freeOut */
/* using the standard acedb format with a TAB delimiter.          */
/* This function is useful for debugging in gdb -                 */
/*   (gdb) print aqlQuick("select ...")                           */
/* to get an answer during a gdb session.                         */
/*                                                                */
/******************************************************************/
void aqlQuick (char *queryText);


/************************* end of file ****************************/

#endif /* AQL_H */

 
 
