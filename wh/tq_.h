/*  File: tq_.h
 *  Author: Richard Durbin (rd@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1996
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * SCCS: $Id: tq_.h,v 1.2 2007/03/21 21:03:14 mieg Exp $
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Nov  5 22:02 1996 (rd)
 * Created: Mon Oct 21 23:09:29 1996 (rd)
 *-------------------------------------------------------------------
 */

#include "table.h"
#include "dict.h"

/******************** DATA STRUCTURES ***********************/
/************************************************************/

typedef union { KEY k ; KEY g ; int i ; float f ;
		char *s ; mytime_t t ; BOOL b ; } TqType ;

typedef struct TqQueryStruct	TqQuery ;
typedef struct TqTableStruct	TqTable ;
typedef struct TqOrderStruct	TqOrder ;
typedef struct TqSFWStruct	TqSFW ;
typedef struct TqFieldStruct	TqField ;
typedef struct TqDeclStruct	TqDecl ;
typedef struct TqLocStruct	TqLoc ;
typedef struct TqExprStruct	TqExpr ;
typedef struct TqBoolStruct	TqBool ;

/********************* query *******************************/

struct TqQueryStruct {
  TABLE *initialActive ;
  char *text ;			/* text of query */
  AC_HANDLE handle ;
  TqTable *table ;
} ;

/********* structure for error reporting **********/

typedef struct TqTextPositionStruct {
  int begin, end, point ;
} TqTextPosition ;

/* Use this when there are errors to report the offending region
   of the query text. point can be used for an operator position
   if there is one. e.g. report:
Query type mismatch in: "7 + `1996-10"
                around:    ^
*/

/********************* table_expr ***************************/

typedef enum { tOR = 1, tAND, tDIFF } TableOp ;

struct TqTableStruct {
  TABLE  *active ;
  DICT   *idDict ;		/* of identifiers */
  BOOL    isActive ;
  TqSFW  *sfw ;			/* select from where statement */
  TableOp  op ;
    TqTable *left, *right ;
  KEY    key ;
    char* keyName ;
  BOOL   isAssign ;
  TqOrder *order ;
/* next three used during evaluation */
  TABLE  *result ;		/* result table */
  char   *idType ;		/* type of each identifier */
  TqType *idValue ;		/* current values of each identifier */
  BOOL   *idExists ;		/* is current value non-null */
  TqTextPosition pos ;
  TqTable *nxt ;
} ;

/* One and only one of isActive, sfw, op, key must non-zero.
   If key is non-zero: if left is non-zero then assign.  
   An id dict must be maintained for the scope of the current table, 
   for nesting operations in tableFuncs (see in expr below).
   The field names are contained in the result table.
*/

/********************** sort_criterion **********************/

typedef enum { ASC = 1, DESC } OrderType ;

struct TqOrderStruct {
  char *fieldName ;
  OrderType type ;
  TqOrder *nxt ;
  TqTextPosition pos ;
} ;

/********************** SELECT FROM WHERE *******************/

struct TqSFWStruct {
  TqField *select ;
  TqDecl *from ;
  TqBool *where ;
} ;

/********* field **********/

struct TqFieldStruct {
  char *name ;
  TqExpr *expr ;
  TqField *nxt ;
  TqTextPosition pos ;
} ;

/********* identifier declaration *********/

/* deal with all complex locator stuff here */

struct TqLocStruct {
  char *idString ;
  char *obClassName ;
    char *obName ;
  int relPos ;
  char* localTagName ;
  char* followTagName ;
  TqLoc *nxt ;
} ;

/* Only one of localTag, followTag can be non-zero.  If neither then
   relpos should be non-zero (else declaration is vacuous).
*/

struct TqDeclStruct {
  char *idName ;
  BOOL isActive ;
  char *className ;
  TqLoc *loc ;
  TqBool *cond ;
  TqDecl *nxt ;
  TqTextPosition pos ;
} ;

/* One and only one of isActive, className, loc must be non-zero.
*/

/* On initial parsing, the WHERE condition is placed in sfw.where.
   Before evaluation, it should be broken up into maximal pieces 
   that can be evaluated after each decl is made, based on what 
   identifiers are required.  These are then transfered to the 
   decl->cond positions.  We could allow cond's to be interspersed
   with decls, giving syntax more like Jean's.
*/

/*********************** bool_expr **************************/

typedef enum { bNOT = 1, bOR, bAND, bXOR } BoolOp ;

typedef enum { uEQ = 1, kkEQ, ksEQ, skEQ, ssEQ, iiEQ, ifEQ, fiEQ, ffEQ,
	       uNE, kkNE, ksNE, skNE, ssNE, iiNE, ifNE, fiNE, ffNE,
	       uLT, iiLT, ifLT, fiLT, ffLT,
	       uLE, iiLE, ifLE, fiLE, ffLE,
	       uGE, iiGE, ifGE, fiGE, ffGE,
	       uGT, iiGT, ifGT, fiGT, ffGT
	     } Comparator ;

struct TqBoolStruct {
  BoolOp op ;
  Comparator comp ;
  int idExists ;
  TqLoc *locExists ;
  union { TqBool *bool ; TqExpr *expr ; } left, right ;
/* used during evaluation */
  BOOL exists ;
  BOOL value ;
  TqTextPosition pos ;
} ;

/* one and only one of op and comp is non_zero */

/*********************** expr *******************************/

typedef enum { kID = 1, kNAME, kCLASS,           /* operators on objects */
	       uUMINUS, iUMINUS, fUMINUS,
	       uPLUS,  uMINUS,  uTIMES,  uDIVIDE,
	       iiPLUS, iiMINUS, iiTIMES, iiDIVIDE,
	       ifPLUS, ifMINUS, ifTIMES, ifDIVIDE,
	       fiPLUS, fiMINUS, fiTIMES, fiDIVIDE,
	       ffPLUS, ffMINUS, ffTIMES, ffDIVIDE,
	       dYEARDIFF, dMONTHDIFF, dWEEKDIFF, dDAYDIFF, /* on dates */
		 dHOURDIFF, dMINDIFF, dSECDIFF
	     } ExprOp ;

typedef enum { tCOUNT = 1, tFIRST, tLAST, tSUM, tMIN, tMAX, tAVG } TableFunc ;

struct TqExprStruct {
  int id ; 
  TqLoc *loc ;
  TableFunc tableFunc ;
    TqTable *table ;
  ExprOp op ;
    TqExpr *left, *right ;
/* used during evaluation, but may be filled in earlier */
  char type ;			/* one of k,t,i,f,s,d */
  TqType value ;
  BOOL exists ;			/* is value non-null */
  TqTextPosition pos ;
} ;

/* at most one of id, tableFunc, op is non-zero
   if all are zero then we have a literal
   if a more complex locator is used than a single identifier, 
     the parser must generate an implicit identifier declaration
*/

/*********************** FUNCTIONS  *************************/
/************************************************************/

TqQuery *tqParse (char *text) ;

/***** tqdebug.c: debugging routines using freeout() *******/

void tqOutLoc (TqLoc *x) ;
void tqOutExpr (TqExpr *e, TqTable *currTable, int ival) ;
void tqOutBool (TqBool *x, TqTable *currTable, int ival) ;
void tqOutDecl (TqDecl *x, TqTable *currTable, int ival) ;
void tqOutField (TqField *x, TqTable *currTable, int ival) ;
void tqOutSFW (TqSFW *x, TqTable *currTable, int ival) ;
void tqOutTable (TqTable *x, int ival) ;
void tqOutQuery (TqQuery *x) ;

/******************** end of file *************************/


 
