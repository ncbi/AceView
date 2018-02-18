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
 * SCCS: $Id: aql_.h,v 1.4 2007/09/14 00:49:17 mieg Exp $
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Sep  2 15:19 2004 (edgrif)
 * * Aug  5 17:20 1998 (fw): completed with new interface specs
 * Created: Mon Oct 21 23:09:29 1996 (rd)
 *-------------------------------------------------------------------
 */

#ifndef ACEDB_AQL__H
#define ACEDB_AQL__H

/************************************************************/

#include <wh/acedb.h>
#include <wh/aceiotypes.h>
#include <wh/dict.h>
#include <wh/table.h>
#include <wh/bs.h>
#include <wh/aql.h>					    /* to import the opaque public type */

/************************************************************/
/******************** DATA STRUCTURES ***********************/
/************************************************************/

/* forward declarations of struct types */
typedef struct AqlNodeStruct	 AqlNode ;
typedef struct AqlQueryStruct	 AqlQuery ;
typedef struct AqlScopeStruct	 AqlScope ;
typedef struct AqlLocStruct	 AqlLoc ;
typedef struct AqlCheckEnvStruct AqlCheckEnv ;
typedef struct AqlResultConsumerStruct AqlResultConsumer;

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



typedef enum {NULLVAL=-1, KNOWNVAL=-10, FUZZYVAL=2} AqlFuzzyType ;

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
  IS_METHOD,
  IS_EXPR
} AqlLocSourceType ;

/***************** dispositions for result consumption *****************/
typedef enum {
    dOutputAlmostImmediately,
    dTableForSorting,
    dJustCountThem
} AqlResultDisposition;

/****************************************************************/
/************************ The Structures ************************/
/****************************************************************/


/********************* AQL *******************************/

struct AqlStruct {
  AqlQuery    *query;			/* the query object, allocated on ->handle */
  jmp_buf     *errorJmpEnv;	        /* exception handling for errors */
  BOOL         IsError;		        /* flag to signal an error (need to clean up error string's storage) */
  int          errorNumber;		/* set if error occurs during processing */
  char        *errorMessage;	        /* string passed into aqlError() - 
					   allocated (on ->handle) and filled, if error occurs */
  char        *errorReport;	        /* fancy error report to point out location of error in query string - 
					   allocated (on ->handle) and filled, if error occurs */
  int          parseDebugLevel;	        /* debugging level
					   0 - no debug info (default)
					   1 - only the pretty output at the end
					   2 - print the intermediate parsetrees
					   3 - all info during query execution (slows it down) */
  int          evalDebugLevel;	        /* debugging level
					   0 - no debug info (default)
					   1 - only the pretty output at the end
					   2 - print the intermediate parsetrees
					   3 - all info during query execution (slows it down) */
  ACEOUT       debug_out;
  AqlCheckEnv *checkEnv;	        /* all stuff needed during check1, alloc'ed on query->handle */

  AC_HANDLE handle;		        /* contains the actual object for query/errorJmpEnv/errorMessage and Report/checkEnv */




    ACEOUT dump_out;
    char beauty;
    KEYSET *result_keyset;
    AC_HANDLE result_keyset_handle;



} ;

/******************* AqlCheckEnv ***************************/
/************* environment during aqlCheck1() **************/

struct AqlCheckEnvStruct {
  AqlNode **pCurrDecl;
  AqlNode **pCurrWhere;
  AqlNode  *currFromLoc;
  
  AqlNode *currParent;

  /* to keep currs when changing scope */
  Stack     declStack;		/* of pCurrDecl  : **AqlNode */
  Stack     whereStack;		/* of pCurrWhere : **AqlNode */
  Stack     fromStack;		/* of curFromLoc : *AqlNode  */
  
  int       declCount;		/* to generate unique names 
				 * for new var-declarations, like "_1", "_2" etc.. */

  /**************************************/
  /* The next 3 scopes must not be allocated upon the checkEnv->handle
   * as scopes in aql->query struct will continue to refer to them.
   * They are allocated on the aql->query->handle to outlast the
   * lifetime of this object */
  AqlScope *currentLocScope;
  AqlScope *currentTableScope;
  AqlScope *currentScalarScope;
} ;

/********************* AqlQuery *******************************/

struct AqlQueryStruct {
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
  AqlNodeType	type ;

  AqlNode	*nxt, *left, *right, *up;
  AqlOpType	op ;
  char		*name ;
  int		number ;
  AqlScope	*scope ;
  AqlLoc	*loc ;
  int		flags ;		/* used to switch scope */

/* next block for values */

  AqlFuzzyType	isValue ;	/* FALSE - value is NULL
				 * TRUE  - use the value
				 * FUZZY - can't say yet - need to re-run */
  AqlType	value ;
  char		vtype ;	 /* type of the value associated with this node
			    'k' = KEY of object
			    'g' = KEY of tag
			    'i' = Integer
			    'f' = Float
			    't' = DateTime
			    'b' = Boolean
			    's' = String
                'T' = Table
                '0' = uninitialized
			  */

  AqlResultConsumer *resultConsumer;
  int		pos ;		/* for error reporting */
} ;

/********************** AqlLoc ***********************/

struct AqlLocStruct {
  AqlFuzzyType	isValue ;	/* NULL/FUZZY/TRUE */
  char		type ;		/* only relevant if isValue is TRUE */
  AqlType	value ;		/* only relevant if isValue is TRUE */
  BOOL		mustExist ;	/* cannot be ->isValue=NULL (no NULL values possible) */
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
				 *     column is node->left->nxt->number - 1 */
  KEY		key ;		/* class, tag etc. */

  /* for row-variables (locs of FROM_TABLE nodes) */
  TABLE		*table ;	/* the table over which to iterate */
  int           row;		/* current row number in the table */
} ;


/***************** AqlScope ************************/

/* var scope changes at SFW boundaries */
/* no predetermined scope boundary for tables - shadowing allowed at any point */

struct AqlScopeStruct {
  DICT		*dict ;		/* of node->name */
  AqlScope	*parent ;
  Array		loc ;		/* of AqlLoc*    */
} ;

/***************** AqlResultConsumer *****************/

struct AqlResultConsumerStruct {
    DICT    *resultTableDict;
    Stack    dictStringStack;
    int      nRows;
    BOOL     sorting;
    AqlNode *sortNode;
    AqlType	 value;
    char     vtype ;
};

/************************************************************/
/*********************** FUNCTIONS  *************************/
/************************************************************/

/***** aqlparse.[ly]:the lex/yacc parser *******/

void aqlParse (AQL aql);

/***** aqlcheck.c: routines to check structure and precompile *******/

void      aqlCheck1 (AQL aql);
void      aqlTraverse (AqlNode *p, 
                       AQL      aql,
                       BOOL   (*pre)(AqlNode*,AQL,int), 
                       void   (*post)(AqlNode*,AQL,int),
                       int level) ;
AqlQuery* aqlQueryCreate (const char *text, AC_HANDLE objectHandle);
AqlScope* aqlScopeCreate (AC_HANDLE handle);

/***** aqlrun.c: database dependent checking and evaluation *****/

void   aqlCheck2 (AQL aql) ;
TABLE* aqlEval   (AQL aql) ;

/**** aqlerror.c: error handling routines using exceptions ******/

void aqlError (AQL aql, int errNo, int errPos, char *format, ...);

/***** aqldebug.c: debugging routines using ACEOUT *******/

void aqlQueryOut (AQL aql, ACEOUT debug_out, AqlQuery *query);
void aqlNodeOut (ACEOUT debug_out, AqlNode *inNode, int level, char *rel, BOOL suppress_next);
void aqlLocOut (ACEOUT debug_out, char *text, AqlLoc *theLocator);
void aqlNodeValueOut (ACEOUT debug_out, AqlNode *inNode) ;
void aqlNodeOutPretty (AQL aql, ACEOUT debug_out, AqlNode *inNode, BOOL doBrackets, int indentLevel) ;
void aqlCheckEnvOut(ACEOUT debug_out,
                    char *text,
                    AqlCheckEnv *env,
                    int level);

char* aqlNodeTypeName (AqlNodeType inType);
char* aqlOpTypeName (AqlOpType inType);
char* aqlLocSourceTypeName (AqlLocSourceType inType);

#if defined(IBM) || defined (CYGWIN)
/* predeclare lex.yy.c fns */
void yyerror (char *s);
#endif

#endif /* !ACEDB_AQL__H */

/********************** end of file *************************/
/************************************************************/
 
 
 
