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
 * SCCS: $Id: aql_.h,v 1.3 2007/03/21 21:03:14 mieg Exp $
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Feb 28 14:11 1999 (fw)
 * * Aug  5 17:20 1998 (fw): completed with new interface specs
 * Created: Mon Oct 21 23:09:29 1996 (rd)
 *-------------------------------------------------------------------
 */

#ifndef AQL__H
#define AQL__H

/************************************************************/

#include "aql.h"		/* to import the opaque public type */

/************************************************************/

#include "dict.h"
#include "keyset.h"
#include "table.h"
#include <setjmp.h>

/************************************************************/
/******************** DATA STRUCTURES ***********************/
/************************************************************/

/* forward declarations of struct types */
typedef struct AqlNodeStruct	 AqlNode ;
typedef struct AqlQueryStruct	 AqlQuery ;
typedef struct AqlScopeStruct	 AqlScope ;
typedef struct AqlLocStruct	 AqlLoc ;
typedef struct AqlCheckEnvStruct AqlCheckEnv ;

/****************** variant record for node->value **********/

typedef union { 
  KEY k ;			/* Key */
  KEY g ;			/* Tag */
  int i ;			/* Integer */
  float f ;			/* Float */
  char *s ;			/* String */
  mytime_t t ;			/* DateType */
  BOOL b ;			/* Boolean */
  KEYSET K ;			/* KeySet */
  TABLE *T ;			/* Table */
} AqlType ;

/********************* scopeflag for node->flags ******************/
  
#define AQL_SWITCH_SCOPE    0x00001000

/********************* enum types for node->type ********************/

typedef enum {			/* left      right     nxt   name      number    value   */

  nTABLE_ASSIGN=1,		/* table	       tab   name  ->  tvar       .T     */
  nVAR_ASSIGN,			/* table     tablFunc  tab   name                 .i/f/t */
  nTABLE_VAR,			/*		       tab   name  ->  tvar       .T     */
  nTABLE_OP,			/* table     table     tab 	           	  .T     */
  nTABLE_SFW,			/* fromwhere select    tab                        .T     */
  nTABLE_SFW_ALL,		/* fromwhere	       tab	           	       4 */
  nTABLE_ORDER,			/* table     sort		           	  .T     */

  nSORT_FIELD,			/* qual			sf	       fieldnum		 */
  nSORT_FIELD_NAME,		/* qual			sf   fieldname		         */
  nSORT_FIELD_EXPR,		/* expr				   		         */
  nSORT_QUALIFIER,		/*	                               op-val	         */

  nSELECT_FIELD,		/* expr                 +    fieldname 		      10 */

  nFROM_LOC,			/* loc	     bool	+    name  ->  var	         */
  nFROM_TABLE,			/* table     bool	+    name  ->  var	         */

  nCLASS,			/*			     classname        -> .k    8 */
  nOBJECT,			/* class_expr name_expr      obj-name  class-key .k  2,9 */

  nLOC_VAR,			/*                      +    name  ->  var	       1 */
  nLOC_TABLE_FIELD,		/*			+              field	     2,7 */
  nLOC_TABLE_FIELD_NAME,	/* 			+    fieldname		     2,7 */
  nLOC_LOCAL_TAG,		/*			+      	                 .k    3 */
  nLOC_LOCAL_TAG_NAME,		/*			+    name	           3,5,6 */
  nLOC_LOCAL_POS,		/*			+     	       offset	     3,6 */
  nLOC_LOCAL_FIELD,		/*                      +              fieldnum      3,6 */
  nLOC_LOCAL_FIELD_NAME,	/*                      +    name                  3,5,6 */
  nLOC_FOLLOW_TAG,		/*			+	                 .k    3 */
  nLOC_FOLLOW_TAG_NAME,		/*			+    name	             3,5 */
  nLOC_FOLLOW_POS,		/*			+	       offset	       3 */

  nLOC_METHOD,			/* loc                       methodname                  */
  nEXPR_TABLE_FUNC,		/* table                                          .i/f/s/t/b */
  nEXPR_TABLE_FUNC_FIELD,	/* table                               fieldnum   .i/f/s/t/b */
  nEXPR_TABLE_FUNC_FIELD_NAME,	/* table                                          .i/f/s/t/b */
  nEXPR_OP,			/* expr       expr           fieldname            .i/f/s/t/b */

  nBOOL_EXISTS,			/* loc						         */
  nBOOL_EXISTS_TAG,		/* loc						         */
  nBOOL_COMPARISON,		/* expr       expr 				         */
  nBOOL_NOT,			/* bool						         */
  nBOOL_OP,			/* bool       bool				         */

  nVAR,				/*                                               .i/f/t  */
  nTEXT,			/*						 .s      */
  nINT,				/*                                               .i      */
  nFLOAT,			/*                                               .f      */
  nDATE,			/*                                               .t      */
  nKEY,				/*                                               .k      */
  nBOOL 		        /*						 .b      */
} AqlNodeType ;

/* key for above:
y  1: can be followed in initial parse by nLOC_LOCAL* or nLOC_FOLLOW*
      after check1 should only be followed in FROM position, by one item only
y  2: can be followed in initial parse by nLOC_LOCAL* or nLOC_FOLLOW*
      after check1 should be only in FROM position, not followed by anything
y  3: see 1
y  4: converted to sfw_all in check1
y  5: converted to un-named varieties in check2
y  6: LOCAL_TAG_NAME, LOCAL_POS converted to FIELD, FIELD_NAME in where nec. check1
y  7: converted to FROM_TABLE and LOC_LOCAL_FIELD(_NAME) in check1
y  8: number set to key in check2
y  9: in check2, class-key, obj-name and .k calculated if possible (if constant
      epxressions), else calculation is deferred.
y 10: absorbed away during check1.  Field names left because no further conflict in
      expressions.
  12: convert EXISTS into a flag on the Loc
  13: merge declarations when building them
  14: problem of intentional copies of variables - currently forbid!
  15: lazy opening of objects
y  16: remove [0], insert [1], and exand out [k] (k > 1)
*/

/* more comments by fw

   '+' in the nxt column means that the nxt-pointer is used for a linked list
   of nodes of the same type. Thereby nodes are chained together at the 
   same level in the tree.

   TABLE_SFW->name is initially 000 for a table with three columns.
   The 0's are later replaced with characters dentoting the type of
   the values in each column
*/


/******************* operators for node->op ***************/
		  
typedef enum {
  oUNION = 1, 
  oINTERSECT, 
  oDIFF,
  oASC, 
  oDESC,
  oEQ, 
  oNE,
  oLIKE,
  oLT, 
  oGT, 
  oLE, 
  oGE,       
  oNOT, 
  oOR, 
  oAND, 
  oXOR,
  oYEARDIFF, 
  oMONTHDIFF, 
  oWEEKDIFF, 
  oDAYDIFF, 
  oHOURDIFF, 
  oMINDIFF, 
  oSECDIFF,       
  oUMINUS, 
  oPLUS, 
  oMINUS, 
  oTIMES, 
  oDIVIDE,       
  oCOUNT, 
  oFIRST, 
  oLAST, 
  oSUM, 
  oMIN, 
  oMAX, 
  oAVG,
  oABS,
  oMOD
} AqlOpType ;

/************** enum source types for node->loc->source ********/

typedef enum { 
  IS_CLASS = 1, 
  IS_ROW, 
  IS_COLUMN, 
  IS_OBJECT, 
  IS_LOCAL_TAG, 
  IS_FOLLOW_TAG, 
  IS_LOCAL_POS, 
  IS_FOLLOW_POS,
  IS_METHOD
} AqlLocSourceType ;


/****************************************************************/
/************************ The Structures ************************/
/****************************************************************/
/* ATTENTION, mieg april 2005
* Because they are handled via handleAlloc/Free, it is mandatory
* to start each struct on a basic type, not on another allocated
* memory like an Array or and AC_*
*/

/********************* AQL *******************************/

typedef struct _AqlStruct {
  BOOL         IsError;		        /* flag to signal an error (need to clean up error string's storage) */
  AqlQuery    *query;			/* the query object, allocated on ->handle */
  jmp_buf     *errorJmpEnv;	        /* exception handling for errors */
  int          errorNumber;		/* set if error occurs during processing */
  char        *errorMessage;	        /* string passed into aqlError() - 
					   allocated (on ->handle) and filled, if error occurs */
  char        *errorReport;	        /* fancy error report to point out location of error in query string - 
					   allocated (on ->handle) and filled, if error occurs */
  int          debugLevel;	        /* debugging level
					   0 - no debug info (default)
					   1 - only the pretty output at the end
					   2 - print the intermediate parsetrees
					   3 - all info during query execution (slows it down) */
  AqlCheckEnv *checkEnv;	        /* all stuff needed during check1, alloc'ed on query->handle */

  AC_HANDLE handle;		        /* contains the actual object for query/errorJmpEnv/errorMessage and Report/checkEnv */
} AqlStruct;

/******************* AqlCheckEnv ***************************/
/************* environment during aqlCheck1() **************/

struct AqlCheckEnvStruct {
  int       declCount ;	/* to generate unique names 
			   for new var-declarations, like "_1", "_2" etc.. */
  AqlScope *currentLocScope ;
  AqlScope *currentTableScope ;
  AqlScope *currentScalarScope ;
  
  AqlNode **pCurrDecl;
  AqlNode **pCurrWhere ;
  AqlNode  *currFromLoc ;
  
  /* to keep currs when changing scope */
  Stack     declStack;	/* of pCurrDecl  : **AqlNode */
  Stack     whereStack;	/* of pCurrWhere : **AqlNode */
  Stack     fromStack ;	/* of curFromLoc : *AqlNode  */
  
} ;

/********************* AqlQuery *******************************/

struct AqlQueryStruct {
  int          dummy ;          /* to be sure we do not have an handlealloc error */
  char	       *text ;		/* text of query (all whitespace cleared), allocated on ->handle */
  AqlNode      *root ;		/* root-node of the parsetree */
  AqlScope     *locScope ;	/* for locator identifiers */
  AqlScope     *tableScope ;	/* for table variable identifiers */
  AqlScope     *scalarScope ;	/* for scalar variable identifiers */
  AC_HANDLE	handle ;	/* handle in which all allocation for a query is done, 
				   its destruction clears all memory allocated during processing,
				   including the ->text, the parsetree and the ->scopes within this object */
} ;

/********************* AqlNode ***************************/

struct AqlNodeStruct {
  int          dummy ;          /* to be sure we do not have an handlealloc error */
  AqlNodeType	type ;

  AqlNode	*nxt, *left, *right ;
  AqlOpType	op ;
  char		*name ;
  int		number ;
  AqlScope	*scope ;
  AqlLoc	*loc ;
  int		flags ;		/* used to switch scope */

/* next block for values */

  BOOL		isEmpty ;	/* set by evalExpr() */
  AqlType	value ;
  DICT          *resultTableDict;
  char		vtype ;	 /* type of the value associated with this node
			    'k' = KEY of object
			    'g' = KEY of tag
			    'i' = Integer
			    'f' = Float
			    't' = DateTime
			    'b' = Boolean
			    's' = String
			  */

  int		pos ;		/* for error reporting */
} ;

/********************** AqlLoc ***********************/

#ifndef DEFINE_OBJ
#define OBJ void*
#define BSMARK void*
#endif
struct AqlLocStruct {
  char		type ;		/* only relevant if not isEmpty */
  AqlType	value ;		/* only relevant if not isEmpty */
  BOOL		isEmpty ;	/* TRUE if the value is NULL */
  BOOL		mustExist ;	/* cannot be ->isEmpty (no NULL values possible) */
  OBJ		trueObj;	/* the object that is bs-created (from loc->value.k), on locs of type KEY */
  OBJ		useObj ;	/* more lightweight locs that look inside objects
				 * (to dereference or find tags) use this variable which 
				 * gets copied from trueObj every time the KEYtype loc iterates to a new value
				 */
  Stack         markStack;
  BSMARK	mark ;		/* remembers the position in the object of the loc, 
				 * initially positioned on trueObj, then as we start iterating through
				 * the object (using useObj) this marks the position on useObj */
  AqlNode	*definition ;	/* points back to AqlNode x for this variable, i.e. node->loc->definition == node */
  BOOL		isContext;	/* locs of context vars need to get past check1 without having a definition node yet */
  AqlLocSourceType	source ;
  int		column ;	/* derived from row-variable LOC_VARs in aqlrun.c
				 *     column is node->left->nxt->number - 1
				 * also use for FOLLOW_POS and LOCAL_POS, where it is 0 ("obj->tag[0]") or 1 ("obj->tag[1][1]")*/
  KEY		key ;		/* class, tag etc. */
  TABLE		*table ;	/* source for row-variables */
} ;


/***************** AqlScope ************************/

/* var scope changes at SFW boundaries */
/* no predetermined scope boundary for tables - shadowing allowed at any point */

struct AqlScopeStruct {
  int          dummy ;          /* to be sure we do not have an handlealloc error */
  DICT		*dict ;		/* of node->name */
  AqlScope	*parent ;
  Array		loc ;		/* of AqlLoc*    */
} ;

/************************************************************/
/*********************** FUNCTIONS  *************************/
/************************************************************/

/***** aqlparse.[ly]:the lex/yacc parser *******/

void aqlParse (AqlStruct *aql);

/***** aqlcheck.c: routines to check structure and precompile *******/

void      aqlCheck1 (AqlStruct   *aql);
void      aqlTraverse (AqlNode   *p, 
		       AqlStruct *aql,
		       BOOL (*pre)(AqlNode*,AqlStruct*), 
		       void (*post)(AqlNode*,AqlStruct*)) ;
AqlQuery* aqlQueryCreate (const char *text, AC_HANDLE objectHandle);
AqlScope* aqlScopeCreate (AC_HANDLE handle);

/***** aqlrun.c: database dependent checking and evaluation *****/

void   aqlCheck2 (AqlStruct *aql) ;
TABLE* aqlEval   (AqlStruct *aql) ;

/**** aqlerror.c: error handling routines using exceptions ******/

void aqlError (AqlStruct *aql, int errNo, int errPos, char *format, ...);
/*
void aqlWarning (int err, int pos, char *format, ...) ;
void aqlReportError (AqlQuery *query, int err) ;
*/

/***** aqldebug.c: debugging routines using freeout() *******/

void aqlQueryOut (AqlQuery *x) ;
void aqlNodeOut (AqlNode *x, int level, char *rel) ;
void aqlNodeOutPretty (AqlNode *x, int ival) ;
void aqlLocOut (char   *text, AqlLoc *theLocator);
void aqlNodeValueOut (AqlNode *x) ;
char* aqlNodeTypeName (AqlNodeType inType);
char* aqlOpTypeName (AqlOpType inType);
char* aqlLocSourceTypeName (AqlLocSourceType inType);

#endif /* AQL__H */

/********************** end of file *************************/
/************************************************************/
 
 
 
