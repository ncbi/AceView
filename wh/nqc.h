/*  Last edited: Dec  4 14:48 1998 (fw) */
/* $Id: nqc.h,v 1.1.1.1 2002/07/19 20:23:19 sienkiew Exp $ */

/* 
 This grammar is taken from the document bnf.4,received from rd  oct 22

 I made a few modif:
 $= set definition, i think all set operations must be $prefixed
   it is better for the parser and also for the human user
 {} i use this notation to create a cleant_list
   of course, these were the conventions of the old system
   but i d like to preserve them.
 * had a double meaning, i denote the default identifier as .>
 : had a double meaning, i denote the field name tokenizer as ::
 DATE i add this function to parse date_literrals as DATE("date_litteral")
   in old system, to parse anything as date creates problem because nearly
   anything can possibly be a date. At the same time, the "" prevents
   parsing of embedded . and -
 ID  i am not sure i understand it, but i do not mind
 -> [] # .   I think they do differ and have this precedence order

 This definition of the pass1 is very terse, as usual in acedb, but already
 quite complete, fully parenthesize the expression, and checks bindings.
 The pass2 still has to check that the resulting tree is meaningful
  for example that you only multiply numbers by numbers
 The pass3 should be model dependant and check for existence of class names
  and the possibility to cast the acedb_data into the requested types
 */

#ifndef DEF_NQ_H  
#define DEF_NQ_H

#include "regular.h"


/* PASS1 of the query compiler runs as follows:
   - Fundamental types (int, float, "text") are recognize,
   - Existing () {} [] are tokenised recursivelly
   - The right most strongest operator, according to the following
    precedence list, is found,
   - It is bound to its neighbours according to its BINDING type
   - An error is reported if the binding is not optional but we hit 
    a closed parenthese or a ; or the end
   - bound to its neighbours and () the compound
   - PASS1 ends when the strongest operator is weaker than LAST_ATOM
*/

/* BINDINGS are used to parenthesize the expression */
typedef enum 
{ 
  NONE = 0 ,  /* and atom */
  LEFT,       /* must gobble one expr on its left,  unusual,  example: ASC */
  RIGHT,      /* must gobble one expr on its right, usual,    example: NOT */
  BOTH,       /* must gobble on both side,          frequent, example: MULTIPLY */
  LEFT_OPT,   /* may on left, must on right,        ackward,  example: MINUS */
  RIGHT_OPT,  /* may on right,                      unusual,  example: FROM */
  BOTH_OPT,   /* may on left and rigth,             peculiar, example: SEMICOL */
  FUNCTION,   /* must gobble a () on its right,     usual,    example: COUNT */
  FUNCTION_2  /* must gobble a (,) on its right,    usual,    example: YEARDIFF */
} NQ_BINDING ;

/* PRECEDENCE:
    ATTENTION:
   To allow C-compiler checks and the usage of fast switch() 
   contructions, all the operators are declared in an enum

   THIS enum DEFINES PRECEDENCE
   example PLUS is left of MULTIPLY because
         3 + 2 * 7 == 3 + (2 * 7)
*/

typedef enum
{
   /*  0 */ NULLQOP = 0,  /* Null operator */
  /*-------------------------------------------------------------------------*/ 
   /*  1 */ SUBEXP,      /* a round  () paranthese pair */
            SUBEXP_C,    /* a curly  {} paranthese pair */
            SUBEXP_S,    /* a square [] paranthese pair */

   /*  4 */ TEXT,        /* a double quoted expression or a non recognised string */
            FLOAT, INTEGER, DATETYPE,   /* a real or int number */

            VARIABLE, TAG, CLASSTYPE, CALCUL,
            ID, THIS, NAME, CLASSE, NOW,
            ACTIVE, FIELD,

	    LAST_ATOM,   /* pseudo operator, separates singlets from true operators */
  /*------------------- stop pass 1 recursion here -------------------------*/ 
  /*-------------------------------------------------------------------------*/ 
            SEMICOL,     /* Separates full queries. Prevents  bindings */ 
            SETDEFINE,   /* affects a table to a table name */
            ORDER_BY,
  /*-------------------------------------------------------------------------*/ 
            SETMINUS, SETOR, SETXOR, SETAND, /* set operators */

  /* 9 */   SELECT, FROM, WHERE, 
	    COMA,
            ASC, DSC,
            AS, CLASS,

            COLUMN,                     /* don t understand col */


  /* 23 */  OR, XOR, AND, NO, NOT, EXISTS,    /* logical operators */

  /* 32 */  LT, LE, GT, GE, NEQ, EQ,       /* comparators */

  /* 38 */  PLUS, MINUS, MULTIPLY, DIVIDE, LIKE, /* arithmetics */

	    COUNT, FIRST, LAST, SUM, MIN, MAX, AVG,

            DATE, FIRSTD, LASTD, 
            YEARDIFF , MONTHDIFF , WEEKDIFF, DAYDIFF , HOURDIFF , MINDIFF , SECDIFF,
          
            OBJECT,
            ARROW, ARRAY, HASH, DOT,
            DOLLAR,
  /*-------------------------------------------------------------------------*/ 
  /* 62 */  LAST_QOP         /* pseudo operator, must come last */
} QOP ;

/* To define the user interface, 
 * the same operators are now listed, with their eventual
 * abbreviation and full name, and the groups they belong to.
 *
 * Although it may, for legibility,
 * if this list is not in the same order as the enum
 * it DOES NOT DEFINE the PRECEDENCE order, the enum does 
 *
 * Note that some operators, for internal use, have no interface.
 * For example a NUMBER or a SUBEXP
 * their names starts with and underscore,
 *   _names are used in report messages, but cannot be parsed.
 *
 * The only extra tokenisation rule is the definition of
 *   string litteral: "" protected
 *   int litteral:  freeint_acceptable
 *   float litteral: freefloat_acceptable
 *
 * I beleive that date litteral should be protected in a special
 * way, so i added a DATE function acting on a "" portected string
 */ 


typedef struct fullop 
{ int qop; char *symbol; char *name ; NQ_BINDING binding ;} NQL ;

#ifndef _define_op_list_  /* a trick to declare the op list only once */
extern Array nqArray ;
#else

struct fullop opListDeclaration [] = 
{
  /* query separator */
  SEMICOL,    ";",  "_SEMI_COLUMN",       BOTH_OPT, 
  SETDEFINE, "$=",  "_SET_DEFINE",        BOTH,

  /* list expressions */

  SELECT,    0,  "SELECT",             RIGHT ,
  FROM,      0,  "FROM",               BOTH,
  WHERE,     0,  "WHERE",              BOTH ,
  ORDER_BY,  0,  "ORDER_BY",           RIGHT ,

  COMA,    ",",  "_COMA",              BOTH,

  ASC,       0,  "ASC",                LEFT,
  DSC,       0,  "DSC",                LEFT,

  AS,        0,  "AS",                 BOTH,
  CLASS,     0,  "CLASS",              RIGHT,
  
  /* i do not understand this operator */
  COLUMN,  ":",  "_COLUMN",            LEFT_OPT,

  /* set operators, possibly redundant, but allows type checking  */
  SETMINUS,"$-",  "SETMINUS",          BOTH,
  SETOR,   "$|",  "SETOR",             BOTH,
  SETXOR,  "$^",  "SETXOR",            BOTH,
  SETAND,  "$&",  "SETAND",            BOTH,

  /* boolean operators */
  OR,      "|",  "OR",                 BOTH,
  XOR,     "^",  "XOR",                BOTH,
  AND,     "&",  "AND",                BOTH,
  NO,        0,  "NO",                 RIGHT,
  NOT,     "!",  "NOT",                RIGHT,
  EXISTS, "!!",  "EXISTS",             RIGHT,

  /* comparators */
  LT,      "<",  "LT",                 BOTH,
  LE,      "<=", "LE",                 BOTH,
  GT,      ">",  "GT",                 BOTH,
  GE,      ">=", "GE",                 BOTH,
  EQ,      "==", "EQ",                 BOTH,
  NEQ,     "!=", "NEQ",                BOTH,

  /* numeric operators */
  PLUS,    "+",  "_PLUS",              LEFT_OPT,
  MINUS,   "-",  "_MINUS",             LEFT_OPT,  /* i do not understand - non-numeric-expre */
  MULTIPLY,"*",  "_MULTIPLY",          BOTH,
  DIVIDE,  "/",  "_DIVIDE",            BOTH,
  LIKE,    "~",  "_LIKE",              BOTH,

  /* function without parentheses */

  /* list functions, take 1 arg */
  COUNT,   0,    "COUNT",              FUNCTION,
  FIRST,   0,    "LAST",               FUNCTION,
  LAST,    0,    "LAST",               FUNCTION,
  SUM,     0,    "SUM",                FUNCTION,
  MIN,     0,    "MIN",                FUNCTION,
  MAX,     0,    "MAX",                FUNCTION,
  AVG,     0,    "AVG",                FUNCTION,

  /* time difference functions, take 2 args */
  DATE,     0,   "DATE",               FUNCTION,  /** i do not like date litterals **/
  FIRSTD,   0,   "LAST",               FUNCTION,
  LASTD,    0,   "LAST",               FUNCTION,
  YEARDIFF, 0,  "YEARDIFF",            FUNCTION_2,
  MONTHDIFF,0,  "MONTHDIFF",           FUNCTION_2,
  WEEKDIFF, 0,  "WEEKDIFF",            FUNCTION_2,
  DAYDIFF,  0,  "DAYDIFF",             FUNCTION_2,
  HOURDIFF, 0,  "HOURDIFF",            FUNCTION_2,
  MINDIFF,  0,  "MINDIFF",             FUNCTION_2,
  SECDIFF,  0,  "SECDIFF",             FUNCTION_2,

  /* locators */
  ARROW,  "->",  "_ARROW",             BOTH,
  ARRAY,     0,  "_ARRAY",             BOTH, /* automatically added left of [] */
  HASH,    "#",  "_HASH",              BOTH,
  DOT,     ".",  "_DOT",               BOTH,

  /* pseudo sets */
  ACTIVE,    0,  "ACTIVE",             NONE,
  DOLLAR,   "$", "_DOLLAR",            RIGHT,

  /* pseudo variables */
  OBJECT,    0, "OBJECT",              FUNCTION_2,
  THIS,    "_", "_THIS",               NONE,   /* was * in bnf.4, conflict with mutiply */
  FIELD,  "::", "_FIELD",              RIGHT,  /* was : in bnf.4, conflict otheyr : */

  /* pseudo atoms */
  ID,        0,  "ID",                 NONE,
  NAME,      0,  "NAME",               NONE,
  CLASSE,    0,  "CLASSE",             NONE,    /******* CLASSE/CLASS conflict ***/
  NOW,       0,  "NOW",                NONE,

  /* From here on, these operators are recognised by the tokenizer
   * outside of the grammar
   */

  SUBEXP,   0,    "__()__",            NONE, 
  SUBEXP_C, 0,    "__{}__",            NONE, 
  SUBEXP_S, 0,    "__[]__",            NONE, 

  /* basic types */
  TEXT,     0,    "_TEXT",             NONE,  /* double quote delimited */
  FLOAT,    0,    "_FLOAT",            NONE,  /* freeint recognized */
  INTEGER,  0,    "_INTEGER",          NONE,  /* freefloat recognized */
  DATETYPE, 0,    "_DATETYPE",         NONE,  /* freedate recognized */

  /* these operators are for the conveneience of the parser */
  VARIABLE, 0,    "_VARIABLE",         NONE,    /* a recognized delared variable */
  TAG,      0,    "_ACEDB_TAG",        NONE,    /* a recognised acedb tag */
  CLASSTYPE,0,    "_ACEDB_CLASS",      NONE,    /* a recognised acedb class-name */
  CALCUL,   0,    "_CALCUL",           NONE,    /* a sub exp casted to numerics: redundant ? */
  LAST_ATOM,0,    "__last_atom_",      NONE,
  LAST_QOP, 0,    "__last_qop_",       NONE /* must come last */ 
} ;
#endif

typedef struct nqNode *NQ;
struct nqNode
{ QOP qop ;
  NQ up, left, right ; /* tree structure */
  KEY tag ; KEY classe ; int xi; float xf ; Array xa ;
  int mark ;   /* mark in the associated stack */
} ;

typedef struct nq 
{ void*   magic;		/* == &MAGIC */
  Stack stack ;
  NQ nq ;
} *NQCOND ;


/* nqTools */

BOOL nqInit(void) ;

void nqDump(Stack s, NQ nq, int level, BOOL isLeft) ;
BOOL nqCheckParentheses (char *text) ; /* verify that "" () {} [] are balanced */

void nqDestroy(NQ c) ;  /* Recursively destroys a tree of nqs */

NQCOND NQCONDalloc (void) ;      /* self managed calloc */
void NQCONDfree (NQCOND nqCond) ;
void NQCONDstatus (int *used, int *alloc) ; /* report */
void  nqCondDestroy (NQCOND nqCond) ;

NQ NQalloc (void)  ;     /* self managed calloc */
void NQfree (NQ nq) ;
void NQstatus (int *used, int *alloc) ;
char *nqName(QOP qop) ;
NQ_BINDING nqBinds (QOP qop) ;


QOP nqTokenize (char *text, char **item) ;
int  nqTokenizeDate (char *text, QOP *newQop) ;
int  nqTokenizeNumber (char *text, QOP *newQop) ;
int  nqTokenizeSymbol (char *cp, BOOL searchNames, QOP *newQop) ;
BOOL nqPass1 (char *text, NQCOND* cdp) ;
BOOL nqPass2 (NQCOND* cdp) ;
BOOL nqPass3 (NQCOND* cdp) ;

#endif

