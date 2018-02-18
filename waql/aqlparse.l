%{
/*  File: aqlparse.l
 *  Author: Stefan Wiesmann and Richard Durbin (wiesmann,rd@sanger.ac.uk)
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
 * Description: lex for AQL
 * Exported functions:
 * HISTORY:
 * Last edited: Sep 13 15:57 1999 (fw)
 * * Mar  1 13:19 1999 (fw): fixed lex/flex incompatability for Linux
 * * Feb 24 10:26 1999 (fw): reintroduced exist_tag (interpreted as tag[0])
 * * Jan 11 11:26 1999 (fw): replaced order_by with order and by
 *                           removed keyword 'exists_tag' (now done with tag[0])
 * * Aug 17 16:13 1998 (fw): ~ changed to 'like'
 * * Aug 14 18:11 1998 (fw): added ~ comparator for string wildcard matching
 * * Aug 10 11:55 1998 (fw): removed keyword ACTIVE (now replaced by context vars in tace)
 * * Jul 27 15:33 1998 (fw): changed error code for invalid character from 700 to 701
 * Created: Tue Oct 29 00:01:02 1996 (rd)
 *-------------------------------------------------------------------
 */

/* $Id: aqlparse.l,v 1.12 2015/08/11 22:15:01 mieg Exp $ */

#if !defined(LINUX) && !defined(OPTERON) && !FreeBSD && !defined(__CYGWIN__) && !defined(DARWIN)

  /* for system that use the AT&T lex and its derivatives
   * we can redefine the input method for the lexer quite
   * easily - by replacing the macros input/unput
   */

#undef input
#undef unput
  /* redefine them using a function that accesses the global string */
  static char input (void) { ++lexPos ; return *lexString++ ; }
  static void unput (char c) { --lexPos ; *--lexString = c ; }
  
#else  /* LINUX || FreeBSD  or CYGWIN */
  
  /* LINUX uses flex, which doesn't provide the above functionality
   * the fact that input/unput are macros was never defined in the
   * POSIX standard. The standard define setting the global 'yyin' to
   * a user supplied FILE* stream as the only way to change the input
   * from stdin to another stream. So flex rigerously implements the standard.
   * We have to do a lot of trickery to convince flex not to work on a FILE* stream.
   */

#define YY_SKIP_YYWRAP		/* yywrap() is used to switch to the next input file 
				 *  if it's finished with one. We're not dealing with
				 *  multiple inputs and ignore this feature.
				 * We need to link with -lfl (flex library) to provide
				 *  the yywrap stub.
				 */
  

#define YY_NEVER_INTERACTIVE 1	/* flex has a notion of interactive buffers
				 *  (streams from tty devices such as stdin) and
				 *  non-interactive buffers (streams from fopen on a file)
				 * To find out which it is, the function 'isatty' is
				 *  is called. When this compiler flag is set, we
				 *  make sure it doesn't need to call this function
				 *  at all. We need to prevent it from calling file-ops
				 *  on yyin, because we set it to (FILE*)1 in aqlParse()
				 *  to fool flex into thinking we have a stream (which prevents
				 *  it from defaulting to stdin)
				 */

  /* redefine the macro which is used by flex to get at another character
   * from the ficticous input stream. We have set yyin to (FILE*)1 and
   * now we let flex get the input from out global input string instead.
   * the redefinition is based on an example given on the Linux flex man-page.
   */
#define YY_INPUT(buf,result,max_size) \
               { \
               int c = *lexString++; ++lexPos ; \
               result = (c == '\0') ? YY_NULL : (buf[0] = c, 1); \
               }

#endif /* LINUX */

  /****************************************************/

#define token(x)	tokPos[x] = lexPos - yyleng ; return x

%}

letter			[A-Za-z]
digit			[0-9]
id			{letter}({letter}|{digit}|"_")*
number			{digit}+

%%

\"([^"]|\\["\n])*\"     { yytext[strlen((const char*)&yytext[0])-1] = 0 ;
                          yylval.String = strnew ((char*)&yytext[1], aql_L->query->handle) ;
                          token(StringLiteral) ;
			  /* match string literals before all keyword, 
			   * so we can quote keywords and match them as tagnames etc. */
                        }

select 			{ token(SELECT) ; }
from  			{ token(FROM) ; }
where   		{ token(WHERE) ; }
all                     { token(ALL); }
as   			{ token(AS) ; }
in   			{ token(IN) ; }
class                   { yylval.String = strnew ("class", aql_L->query->handle); /* for EXPR_LOC_FUNC */
                          token(CLASS) ; }
object  		{ token(OBJECT) ; }
exists   		{ token(EXISTS) ; }
exists_tag   		{ token(EXISTS_TAG) ; }
now			{ token(NOW) ; }
today			{ token(TODAY) ; }
not  			{ token(NOT) ; }
and  			{ token(AND) ; }
or   			{ token(OR) ; }
xor   			{ token(XOR) ; }
union  			{ token(UNION) ; }
intersect		{ token(INTERSECT) ; }
except 			{ token(DIFF) ; }
diff 			{ token(DIFF) ; }
order   		{ token(ORDER) ; }
sort	   		{ token(ORDER) ; }
by                      { token(BY) ; }
true   			{ token(TRUEtok) ; }
false  			{ token(FALSEtok) ; }

yeardiff		{ yylval.Int = oYEARDIFF ; token(ExprExprFunc) ; }
monthdiff		{ yylval.Int = oMONTHDIFF ; token(ExprExprFunc) ; }
weekdiff		{ yylval.Int = oWEEKDIFF ; token(ExprExprFunc) ; }
daydiff			{ yylval.Int = oDAYDIFF ; token(ExprExprFunc) ; }
hourdiff		{ yylval.Int = oHOURDIFF ; token(ExprExprFunc) ; }
mindiff			{ yylval.Int = oMINDIFF ; token(ExprExprFunc) ; }
secdiff			{ yylval.Int = oSECDIFF ; token(ExprExprFunc) ; }

count			{ yylval.Int = oCOUNT ; token(TableFunc) ; }
first			{ yylval.Int = oFIRST ; token(TableFunc) ; }
last			{ yylval.Int = oLAST ; token(TableFunc) ; }
min			{ yylval.Int = oMIN ; token(TableFunc) ; }
max			{ yylval.Int = oMAX ; token(TableFunc) ; }
sum			{ yylval.Int = oSUM ; token(TableFunc) ; }
avg			{ yylval.Int = oAVG ; token(TableFunc) ; }

abs			{ yylval.Int = oABS ; token(ExprFunc) ; }
modulo			{ yylval.Int = oMOD ; token(ExprExprFunc) ; }
mod			{ yylval.Int = oMOD ; token(ExprExprFunc) ; }

asc   			{ yylval.Int = oASC ; token(Ordering) ; }
desc   	 		{ yylval.Int = oDESC ; token(Ordering) ; }

like                    { yylval.Int = oLIKE ; token(Comparator) ; }

{id}			{ yylval.String = strnew ((char*)&yytext[0], aql_L->query->handle) ;
			  token(Identifier) ;
			}

[\@\$\;\:\.\(\)\[\]\+\-\*\/\%\,\_\|]	token((int)yytext[(int)0]) ;

":="			{ token(ASSIGN) ; }
"->"			{ token(ARROW) ; }
"::"			{ token(DOUBLE_COLON) ; }

"="                     { yylval.Int = oEQ ; token(Comparator) ; }
"<"                     { yylval.Int = oLT ; token(Comparator) ; }
">"                     { yylval.Int = oGT ; token(Comparator) ; }
"!="                    { yylval.Int = oNE ; token(Comparator) ; }
"<="                    { yylval.Int = oLE ; token(Comparator) ; }
">="                    { yylval.Int = oGE ; token(Comparator) ; }

{number}                { sscanf ((const char*)&yytext[0], "%d", &yylval.Int) ;
                          token(Number) ;
                        }

[0-9]+[eE][+-]?[0-9]+		|
[0-9]+"."[0-9]*[eE][+-]?[0-9]+	|
"."[0-9]*[eE][+-]?[0-9]+	|
[0-9]+"."[0-9]*			|
"."[0-9]+		{ sscanf ((const char*)&yytext[0], "%f", &yylval.Float) ;
			  token(FloatLiteral) ;
                        }


"`"[12][90][0-9][0-9]"-"[01][0-9]"-"[0-3][0-9]"_"[0-2][0-9]":"[0-5][0-9]":"[0-5][0-9]	|	/* YYYY-MM-DD_HH:MM:SS */	
"`"[12][90][0-9][0-9]"-"[01][0-9]"-"[0-3][0-9]"_"[0-2][0-9]":"[0-5][0-9]		|	/* YYYY-MM-DD_HH:MM */	
"`"[12][90][0-9][0-9]"-"[01][0-9]"-"[0-3][0-9]"_"[0-2][0-9]				|	/* YYYY-MM-DD_HH */	
"`"[12][90][0-9][0-9]"-"[01][0-9]"-"[0-3][0-9]						|	/* YYYY-MM-DD */	
"`"[12][90][0-9][0-9]"-"[01][0-9]							|	/* YYYY-MM */	
"`"[12][90][0-9][0-9]	{ yylval.String = strnew ((char*)&yytext[1], aql_L->query->handle) ;    /* YYYY */	
			  token(DateLiteral) ;
			}



[ \t\r\f\n]+              ;  /* whitespace */

.                       aqlError (aql_L, 701, lexPos == 0 ? 1 : lexPos-1, "syntax error: invalid character") ;

%%

int yyerror (char *s)
{
  /* Syntax error */
  aqlError (aql_L, 700, lexPos ? lexPos-1 : 1, s);  /* aql is a static in aqlparse.y - the current AQL object */
  return 1 ;
}


/********************** end of file *******************/
