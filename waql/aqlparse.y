/*  File: aqlparse.y
 *  Author: Stefan Wiesmann and Richard Durbin (wiesmann,rd@sanger.ac.uk)
 *          and Fred Wobus (fw@sanger.ac.uk)
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
 * SCCS: $Id: aqlparse.y,v 1.4 2015/01/25 04:52:50 mieg Exp $
 * Description: yacc grammar for AQL
 * Exported functions:
 *              aqlParse()
 * HISTORY:
 * Last edited: Nov  9 10:47 1999 (fw)
 * * Sep  7 17:12 1999 (fw): substantial addition to the grammar -
 *              'basic_decl' can now be 'expr' for a much more
 *              general kind of declaration
 * * May  3 11:25 1999 (fw): added 'class <cName>."<templ>"' notation
 * * Mar  2 10:07 1999 (fw): all identifiers referring to names in the database
 *              can now be quoted. select a->"min" from a in class "class"
 *              is now possible. All user chosen identifiers have to be non-keywords still.
 * * Feb 26 10:17 1999 (fw): merged bool_expr: into expr: so all boolean
 *              evaluation is now just an expression value of type BOOL
 * * Aug 11 13:48 1998 (fw): introduction of a default row-variable locator
                             for 'select all @t where :1 = "xx"' type queries.
 * * Aug 10 16:53 1998 (fw): made class() locFunc possible
 * * Jul 10 11:55 1998 (fw): SORT_QUALIFIER now sets ->op to ASC/DESC, rather than number
 * * Jul 9 10:23 1998 (fw): changed alias syntax from Id:(expr) to Id::expr
 * Created: Tue Oct 29 00:01:02 1996 (rd)
 *-------------------------------------------------------------------
 */
%{

#include "acedb.h"
#include "aceio.h"
#include "waql/aql_.h"

/**************************************************************/

static mytime_t timeCreateFromString (char *dateString);

/********* globals for this module - NOT threadsafe **********/

static AQL aql_L;		/* current aql object to work with */

static int tokPos[1024];

/**************************************************************/

  /* macros for creating new nodes, they should really used some
     standard interface, like aqlcheck.c:makeNewNode()  */
  /* create a node of type t */
#define zMake(t)			AqlNode* z = halloc(sizeof(AqlNode),aql_L->query->handle); yyval.Ptr = z; z->type = t; z->vtype = '0'

  /* used to chain multiples of select_fields, sortfields or table_exprs together */
#define zzFinal(x)			AqlNode* zz = (AqlNode*)x; yyval.Ptr = zz; while (zz->nxt) zz = zz->nxt

  /* used to chain local field dereferencer together, e.g. a->paper->author->name */
#define zzFinalMake(x,t)		zzFinal(x); zz->nxt = (AqlNode*)halloc(sizeof(AqlNode),aql_L->query->handle); zz = zz->nxt; zz->type = t; zz->vtype = '0'
#define posSet(x,t)			(x)->pos = tokPos[t]

/**************************************************************/

%}
 
%union  {
  char*          String;          
  int		 Int;
  float		 Float;
  int		 Command;
  void*		 Ptr;		
}

%token <String>  Identifier
%token <Int>     Number
%token <Float>   FloatLiteral
%token <String>  StringLiteral
%token <String>  DateLiteral
%token <Int>     Comparator
%token <Int>	 TableFunc 
%token <Int>	 LocatorFunc 
%token <Int>	 ExprFunc 
%token <Int>	 ExprExprFunc
%token <Int>	 Ordering 
%token <Command> SELECT FROM WHERE ALL ORDER BY
%token <Command> TRUEtok FALSEtok
%token <Command> AS OBJECT EXISTS EXISTS_TAG CLASS
%token <Command> ARROW ASSIGN DOUBLE_COLON IN
%token <Command> NOW TODAY
%token <Command> AND OR XOR
%token <Command> UNION INTERSECT DIFF 

%left ';'
%left ASSIGN 
%nonassoc ORDER BY
%left UNION
%left INTERSECT DIFF
%left OR XOR
%left AND
%left NOT
%left Comparator
%left '+' '-'
%left '*' '/'
%left '%'
%nonassoc UMINUS
%nonassoc EXISTS
%nonassoc EXISTS_TAG

%%

start:        query                             { aql_L->query->root = $<Ptr>1; $<Ptr>$ = $<Ptr>1; }
             ;

query:        table_expr			{ $<Ptr>$ = $<Ptr>1; }
	    | query ';' query			{ zzFinal($<Ptr>1); zz->nxt = $<Ptr>3; aql_L->query->root = $<Ptr>1; }
	    | '@' Identifier ASSIGN table_expr  { zMake(nTABLE_ASSIGN); z->left = $<Ptr>4; z->name = $<String>2; $<Ptr>$ = z; posSet(z,ASSIGN); }
            | '$' Identifier ASSIGN expr	{ zMake(nVAR_ASSIGN); z->left = $<Ptr>4;  z->name = $<String>2; $<Ptr>$ = z; posSet(z,ASSIGN); }
	   ;


table_expr:   safe_table			{ $<Ptr>$ = $<Ptr>1; }
	    | SELECT fieldlist          	{ zMake(nTABLE_SFW); z->right = $<Ptr>2; posSet(z,SELECT); }
	    | SELECT fieldlist FROM fwlist	{ zMake(nTABLE_SFW); z->left = $<Ptr>4; z->right = $<Ptr>2; posSet(z,SELECT); }
	    | SELECT ALL fwlist		        { zMake(nTABLE_SFW_ALL); z->left = $<Ptr>3; posSet(z,SELECT); } /* select all class XXX  -or-
														 * select all object(..) -or-
														 * select all @table */
	    | SELECT ALL FROM fwlist	        { zMake(nTABLE_SFW_ALL); z->left = $<Ptr>4; posSet(z,SELECT); }
	    | table_expr ORDER			{ zMake(nTABLE_ORDER); z->left = $<Ptr>1; z->op = oASC; posSet(z,ORDER); }
            | table_expr ORDER Ordering		{ zMake(nTABLE_ORDER); z->left = $<Ptr>1; z->op = $<Int>3; posSet(z,ORDER); }
	    | table_expr ORDER BY sortlist	{ zMake(nTABLE_ORDER); z->left = $<Ptr>1; z->right = $<Ptr>4; $<Ptr>$ = z; posSet(z,ORDER); }
            | table_expr UNION table_expr	{ zMake(nTABLE_OP); z->op = oUNION; z->left = $<Ptr>1; z->right = $<Ptr>3; posSet(z,UNION); }
            | table_expr INTERSECT table_expr	{ zMake(nTABLE_OP); z->op = oINTERSECT; z->left = $<Ptr>1; z->right = $<Ptr>3; posSet(z,INTERSECT); }
            | table_expr DIFF table_expr 	{ zMake(nTABLE_OP); z->op = oDIFF; z->left = $<Ptr>1; z->right = $<Ptr>3; posSet(z,DIFF); }
	   ;


safe_table:   '(' table_expr ')'		{ $<Ptr>$ = $<Ptr>2; }
	    | '@' Identifier			{ zMake(nTABLE_VAR); z->name = $<String>2; posSet(z,'@'); }
	     ;		

sortlist:     sort_criterion			{ $<Ptr>$ = $<Ptr>1; }
            | sortlist ',' sort_criterion	{ zzFinal($<Ptr>1); zz->nxt = $<Ptr>3; }
             ;

sort_criterion: ':' Identifier		{ zMake(nSORT_FIELD_NAME); z->name = $<String>2; posSet(z,':'); }
              | ':' TableFunc		{ zMake(nSORT_FIELD_NAME); z->name = strnew(aqlOpTypeName($<Int>2), aql_L->query->handle); posSet(z,':'); }
              | ':' Number		{ zMake(nSORT_FIELD); z->number = $<Int>2; posSet(z,':'); }
	      | ':' Identifier Ordering	{ zMake(nSORT_FIELD_NAME); z->name = $<String>2;
					  z->left = (AqlNode*)halloc(sizeof(AqlNode),aql_L->query->handle);
	      				  z->left->type = nSORT_QUALIFIER; z->left->op = $<Int>3; posSet(z,':'); }
	      | ':' Number  Ordering	{ zMake(nSORT_FIELD); z->number = $<Int>2; 
	                                  z->left = (AqlNode*)halloc(sizeof(AqlNode),aql_L->query->handle);
					  z->left->type = nSORT_QUALIFIER; z->left->op = $<Int>3; posSet(z,':'); }
/*	      | ':' expr		{ zMake(nSORT_FIELD_EXPR); z->left = $<Ptr>2; posSet(z,':'); }*/
               ; 

fieldlist:    field			{ $<Ptr>$ = $<Ptr>1; }
	    | fieldlist ',' field	{ zzFinal($<Ptr>1); zz->nxt = $<Ptr>3; }
             ;

field:	      expr			{ zMake(nSELECT_FIELD); z->left = $<Ptr>1; }
	    | Identifier DOUBLE_COLON  expr { zMake(nSELECT_FIELD); z->name = $<String>1; z->left = $<Ptr>3; posSet(z,DOUBLE_COLON); }
             ;

fwlist:       fw			{ $<Ptr>$ = $<Ptr>1; }
            | fwlist ',' fw             { zzFinal($<Ptr>1); zz->nxt = $<Ptr>3; }
             ;
                
fw:	      basic_decl		{ AqlNode *z = $<Ptr>1; z->name = "_DEF"; $<Ptr>$ = z; }
            | basic_decl WHERE expr	{ AqlNode *z = $<Ptr>1; z->name = "_DEF"; z->right = $<Ptr>3; $<Ptr>$ = z; }
            | Identifier IN basic_decl	{ AqlNode *z = $<Ptr>3; z->name = $<String>1; $<Ptr>$ = z; }
            | Identifier IN basic_decl WHERE expr { AqlNode *z = $<Ptr>3; z->name = $<String>1; z->right = $<Ptr>5; $<Ptr>$ = z; }
            | basic_decl AS Identifier	{ AqlNode *z = $<Ptr>1; z->name = $<String>3; $<Ptr>$ = z; }
            | basic_decl AS Identifier WHERE expr { AqlNode *z = $<Ptr>1; z->name = $<String>3; z->right = $<Ptr>5; $<Ptr>$ = z; }
            | basic_decl Identifier	{ AqlNode *z = $<Ptr>1; z->name = $<String>2; $<Ptr>$ = z; }	
            | basic_decl Identifier WHERE expr { AqlNode *z = $<Ptr>1; z->name = $<String>2; z->right = $<Ptr>4; $<Ptr>$ = z; }	
             ;

basic_decl:   expr			{ zMake(nFROM_LOC); z->left = $<Ptr>1; z->pos = z->left->pos; }
	    | safe_table		{ zMake(nFROM_TABLE); z->left = $<Ptr>1; z->pos = z->left->pos; }
            | CLASS Identifier		{ zMake(nFROM_LOC); z->left = (AqlNode*)halloc(sizeof(AqlNode),aql_L->query->handle);
	    				  z = z->left; z->type = nCLASS; z->name = $<String>2; posSet(z,Identifier); }
            | CLASS Identifier '.' StringLiteral { zMake(nFROM_LOC); z->left = (AqlNode*)halloc(sizeof(AqlNode),aql_L->query->handle);
	    				  z = z->left; z->type = nCLASS; z->name = $<String>2; posSet(z,Identifier); 
	                                  z->nxt = (AqlNode*)halloc(sizeof(AqlNode),aql_L->query->handle);
					  z->nxt->type = nTEXT; z->nxt->isValue = KNOWNVAL; z->nxt->value.s = $<String>4; z->nxt->vtype = 's';
	                                  }
            | CLASS StringLiteral	{ zMake(nFROM_LOC); z->left = (AqlNode*)halloc(sizeof(AqlNode),aql_L->query->handle);
	    				  z = z->left; z->type = nCLASS; z->name = $<String>2; posSet(z,StringLiteral); }
            | CLASS StringLiteral '.' StringLiteral { zMake(nFROM_LOC); z->left = (AqlNode*)halloc(sizeof(AqlNode),aql_L->query->handle);
	    				  z = z->left; z->type = nCLASS; z->name = $<String>2; posSet(z,StringLiteral); 
	                                  z->nxt = (AqlNode*)halloc(sizeof(AqlNode),aql_L->query->handle);
					  z->nxt->type = nTEXT; z->nxt->isValue = KNOWNVAL; z->nxt->value.s = $<String>4; z->nxt->vtype = 's';
	                                  }
	     ;

locator:      Identifier		{ zMake(nLOC_VAR); z->name = $<String>1; posSet(z,Identifier); }
     	    | OBJECT '(' expr ',' expr ')' { zMake(nOBJECT); z->left = $<Ptr>3; z->right = $<Ptr>5; posSet(z,'('); }
            | locator ARROW Identifier	{ zzFinalMake($<Ptr>1,nLOC_FOLLOW_TAG_NAME); zz->name = $<String>3; posSet(zz,ARROW); }
            | locator ARROW StringLiteral { zzFinalMake($<Ptr>1,nLOC_FOLLOW_TAG_NAME); zz->name = $<String>3; posSet(zz,ARROW); }
            | locator ARROW Number      { zzFinalMake($<Ptr>1,nLOC_FOLLOW_POS); zz->number = $<Int>3; posSet(zz,ARROW); }
            | locator '[' Identifier ']'{ zzFinalMake($<Ptr>1,nLOC_LOCAL_TAG_NAME); zz->name = $<String>3; posSet(zz,'['); }
            | locator '[' StringLiteral ']'{ zzFinalMake($<Ptr>1,nLOC_LOCAL_TAG_NAME); zz->name = $<String>3; posSet(zz,'['); }
            | locator '[' Number ']'	{ zzFinalMake($<Ptr>1,nLOC_LOCAL_POS); zz->number = $<Int>3; posSet(zz,'['); }
            | locator ':' Identifier	{ zzFinalMake($<Ptr>1,nLOC_LOCAL_FIELD_NAME); zz->name = $<String>3; posSet(zz,':'); }
            | locator ':' StringLiteral	{ zzFinalMake($<Ptr>1,nLOC_LOCAL_FIELD_NAME); zz->name = $<String>3; posSet(zz,':'); }
            | locator ':' Number	{ zzFinalMake($<Ptr>1,nLOC_LOCAL_FIELD); zz->number = $<Int>3; posSet(zz,':'); }
            | locator '.' CLASS		{ zzFinalMake($<Ptr>1,nLOC_METHOD); zz->name = $<String>3; posSet(zz,'.');
                                          /* by matching the keyword 'class' we can write loc.class as a valid method. All other
					   * method names are matched by the Identifier rule below */
	                                }
	    | locator '.' Identifier	{ zzFinalMake($<Ptr>1,nLOC_METHOD); zz->name = $<String>3; posSet(zz,'.'); }

	    | ARROW Identifier		{ zMake(nLOC_VAR); z->name = "_DEF"; 
					  { zzFinalMake(z,nLOC_FOLLOW_TAG_NAME); zz->name = $<String>2; posSet(zz,Identifier); }}
	    | ARROW StringLiteral	{ zMake(nLOC_VAR); z->name = "_DEF"; 
					  { zzFinalMake(z,nLOC_FOLLOW_TAG_NAME); zz->name = $<String>2; posSet(zz,StringLiteral); }}
	    | ARROW Number              { zMake(nLOC_VAR); z->name = "_DEF"; 
					  { zzFinalMake(z,nLOC_FOLLOW_POS); zz->number = $<Int>2; posSet(zz,Number); }}
            | '[' Identifier ']'	{ zMake(nLOC_VAR); z->name = "_DEF"; 
					  { zzFinalMake(z,nLOC_LOCAL_TAG_NAME); zz->name = $<String>2; posSet(zz,Identifier); }}
            | '[' StringLiteral ']'	{ zMake(nLOC_VAR); z->name = "_DEF"; 
					  { zzFinalMake(z,nLOC_LOCAL_TAG_NAME); zz->name = $<String>2; posSet(zz,StringLiteral); }}
	    | '[' Number ']'		{ zMake(nLOC_VAR); z->name = "_DEF"; 
					  { zzFinalMake(z,nLOC_LOCAL_POS); zz->number = $<Int>2; posSet(zz,Number); }}
            | ':' Identifier		{ zMake(nLOC_VAR); z->name = "_DEF";
					  { zzFinalMake(z,nLOC_LOCAL_FIELD_NAME); zz->name = $<String>2; posSet(zz,Identifier); }}
            | ':' Number		{ zMake(nLOC_VAR); z->name = "_DEF";
					  { zzFinalMake(z,nLOC_LOCAL_FIELD); zz->number = $<Int>2; posSet(zz,Number); }}
            | '.' CLASS			{ zMake(nLOC_VAR); z->name = "_DEF";
					  { zzFinalMake(z,nLOC_METHOD); zz->name = $<String>2; posSet(zz,'.'); }}
	    | '.' Identifier		{ zMake(nLOC_VAR); z->name = "_DEF";
					  { zzFinalMake(z,nLOC_METHOD); zz->name = $<String>2; posSet(zz,'.'); }}
 	    | safe_table ':' Identifier	{ zMake(nLOC_TABLE_FIELD_NAME); z->left = $<Ptr>1; z->name = $<String>3; posSet(z,':'); }
	    | safe_table ':' Number	{ zMake(nLOC_TABLE_FIELD); z->left = $<Ptr>1; z->number = $<Int>3; posSet(z,':'); }
             ;     

expr:         locator          		{ $<Ptr>$ = $<Ptr>1; }
            | '$' Identifier		{ zMake(nVAR); z->name = $<String>2; posSet(z,'$'); }
            | Number       		{ zMake(nINT); z->value.i = $<Int>1; posSet(z,Number); }
            | FloatLiteral		{ zMake(nFLOAT); z->value.f = $<Float>1; posSet(z,FloatLiteral); }
            | StringLiteral		{ zMake(nTEXT); z->value.s = $<String>1; posSet(z,StringLiteral); }
	    | date                      { $<Ptr>$ = $<Ptr>1; }

	    | TableFunc safe_table	{ zMake(nEXPR_TABLE_FUNC); z->op = $<Int>1; z->left = $<Ptr>2; posSet(z,'(');}
	    | TableFunc safe_table ':' Number     { zMake(nEXPR_TABLE_FUNC_FIELD); 
                                                    z->op = $<Int>1; z->left = $<Ptr>2; z->number = $<Int>4; posSet(z,':');}
	    | TableFunc safe_table ':' Identifier { zMake(nEXPR_TABLE_FUNC_FIELD_NAME); 
                                                    z->op = $<Int>1; z->left = $<Ptr>2; z->name = $<String>4; posSet(z,':');}
            | ExprExprFunc '(' expr ',' expr ')' { zMake(nEXPR_OP); z->op = $<Int>1; z->left = $<Ptr>3; z->right = $<Ptr>5; posSet(z,'('); }
            | ExprFunc '(' expr ')'	{ zMake(nEXPR_OP); z->op = $<Int>1; z->left = $<Ptr>3; posSet(z,'('); }
            | expr '+' expr		{ zMake(nEXPR_OP); z->op = oPLUS; z->left = $<Ptr>1; z->right = $<Ptr>3; posSet(z,'+'); }
            | expr '-' expr		{ zMake(nEXPR_OP); z->op = oMINUS; z->left = $<Ptr>1; z->right = $<Ptr>3; posSet(z,'-'); }
            | expr '*' expr		{ zMake(nEXPR_OP); z->op = oTIMES; z->left = $<Ptr>1; z->right = $<Ptr>3; posSet(z,'*'); }
            | expr '/' expr		{ zMake(nEXPR_OP); z->op = oDIVIDE; z->left = $<Ptr>1; z->right = $<Ptr>3; posSet(z,'/'); }
            | '-' expr %prec UMINUS	{ zMake(nEXPR_OP); z->op = oUMINUS; z->left = $<Ptr>2; posSet(z,'-'); }
            | expr '%' expr             { zMake(nEXPR_OP); z->op = oMOD; z->left = $<Ptr>1; z->right = $<Ptr>3; posSet(z,'%'); }
            | TRUEtok			{ zMake(nBOOL); z->value.b = TRUE; posSet(z,TRUEtok); }
            | FALSEtok			{ zMake(nBOOL); z->value.b = FALSE; posSet(z,FALSEtok); }
            | EXISTS  locator		{ zMake(nBOOL_EXISTS); z->left = $<Ptr>2; posSet(z,EXISTS); }
            | EXISTS_TAG  locator	{ zMake(nBOOL_EXISTS_TAG); z->left = $<Ptr>2; posSet(z,EXISTS_TAG); }
            | EXISTS  '(' locator ')'	{ zMake(nBOOL_EXISTS); z->left = $<Ptr>3; posSet(z,EXISTS); }
            | EXISTS_TAG  '(' locator ')' { zMake(nBOOL_EXISTS_TAG); z->left = $<Ptr>3; posSet(z,EXISTS_TAG); }
            | expr Comparator expr	{ zMake(nBOOL_COMPARISON); z->op = $<Int>2; z->left = $<Ptr>1; z->right = $<Ptr>3; posSet(z,Comparator); }
            | NOT expr		{ zMake(nBOOL_NOT); z->left = $<Ptr>2; posSet(z,NOT); }
            | expr AND expr	{ zMake(nBOOL_OP); z->op = oAND; z->left = $<Ptr>1; z->right = $<Ptr>3; posSet(z,AND); }
            | expr OR expr	{ zMake(nBOOL_OP); z->op = oOR; z->left = $<Ptr>1; z->right = $<Ptr>3; posSet(z,OR); }
            | expr XOR expr	{ zMake(nBOOL_OP); z->op = oXOR; z->left = $<Ptr>1; z->right = $<Ptr>3; posSet(z,XOR); }
            | '(' expr ')'		{ $<Ptr>$ = $<Ptr>2; }
             ;

date:	      NOW			{ zMake(nDATE); z->value.t = timeParse ("now"); posSet(z,NOW); }
            | TODAY			{ zMake(nDATE); z->value.t = timeParse ("today"); posSet(z,TODAY); }
            | DateLiteral		{ zMake(nDATE); z->value.t = timeCreateFromString ($<String>1); posSet(z,DateLiteral); }
             ;  
%%

/************ entry point: aqlParse() ***********************/

static char *lexString;
static int   lexPos;
int yyparse ();

void aqlParse (AQL aql)
{

#if defined(LINUX_JUNK) 

  extern FILE *yyin;
  extern FILE *yyout;
  
  /* this will fool the flex generated lexer into thinking that we 
   * have an input file stream. It avoids the default setting to 'stdin'. */
  yyin = (FILE*)1; yyout = (FILE*)1;
#endif

  aql_L = aql;		/* make current object global to this file (not threadsafe)*/

  lexString = strnew (aql_L->query->text, aql_L->query->handle); /* let the lexer work on a copy */
			       /* global char * pointer lexString is modified during parsing */
  lexPos = 0;


  yyparse ();			/* return code ignored - error handle by aqlError longjmp */

  return;
} /* aqlParse */

static mytime_t timeCreateFromString (char *dateString)
{
  mytime_t returnTimeElement;
  
  returnTimeElement = timeParse (dateString);
  if (!returnTimeElement)
    aqlError (aql_L, 702, lexPos-1, "Bad date/time literal");

  return returnTimeElement;
}

/***************** lex code ************************/

#include "lex.yy.c"

/***************** end of file **********************/
