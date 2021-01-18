/*
* cbeau - C beautifier.
*
* usage: cbeau [ -s <number> ] [ -d ] [ files ]
*
*	The optional -s <number> is the number of columns to indent for each
*	indent level.  default is -2.
*
*       -d debug mode, stop if there is a { or a } inside a ""
*
*	If file names are given, all files are contcatenated and beautified
*	as one.  If no file name is given, standard input is beautified.
*
*	The beautified output is on standard output.
*
*/
#include <stdlib.h>
#include <ctype.h>
#include <stdio.h>
#include <string.h>

static int space_per_level=2;

/*
static int string_wrap_desire = 80 ;
static int string_wrap_force = 200;
static int no_tabs = 0;
*/
static int debug = 0;
static int verbose = 1 ;
static char word [2048] ; /* word buffer */

/*
* These are global to preserve them from one call to beautify
* to the next.
*
* level is the current indent level.
*
* needspace is true if we saw a space in the input, and no space has yet
* been put in the output at this point.  The next non-space character to
* come along will cause a space to be output.
*
* lastspace is true if the last character we output was a space.  This means
* we can (really must) suppress further consecutive spaces.
*


static int level=0;
static int lastspace=0;
*/

static int neededSpace=0;
static int c0 ;
static FILE *input = 0 ;

#define GETC(_c) _c=getc(input);if (_c == EOF) return 0 

static void showWord (void)
{
  if (word [0]) /* consume word */
    {
      char *cp = word - 1 ;
      while (*++cp) putchar (*cp) ; /* avoid ++ in putchar which may be a macro */
      c0 = *(cp - 1) ;
      word [0] = 0 ;
    }
}

static void noSpace (void)
{
  showWord () ;
  neededSpace=0;
}

static void space (void)
{
  {  /* lay out */
    int i = neededSpace ;  
    while (i-- > 0)
      putchar (' ') ;
  }
  showWord () ;
  neededSpace = 1 ;
}

static int lineBreak (int level)
{ 
  if (word[0]) space () ;
  neededSpace = 2 * level ;
  if (c0 != '\n')
    {
      putchar ('\n') ;
      c0 = '\n' ;
    }
  return c0 ;
} /* lineBreak */

static int cPut (int c)
{
  space () ;
  putchar (c) ;

  c0 = c ;
  return c0 ;
}

static int cGet (void)
{
  int c = getc(input) ;
  if (c == EOF) /* do not transform EOF to int ! */
    return 0 ;
  switch (c)
    {
    case 0:
      return 0 ;
    case '\\':  /* universal protection of next char */
      putchar(c) ; 
      c = getc(input) ;
      if (c == 0 || c == EOF) /* do not transform EOF to int ! */
	return 0 ;
      putchar(c) ; 
      return '\\' ; /* take no action in calling function */
    }
  return c ;
} /* cGet */

static int beautify (FILE *myInput)
{
  int c, cTmp, nw = 0, inWord = 0, nBlock = 0 ;
  int blocks [200] ;
  int    level = 0 ,
    /* nQuote = 0 ,    depth inide  ' ' signs */
    nDoubleQuote = 0 ,   /* depth inide " " signs */
    nPar = 0 ;  /* depth inide () signs */
    
  word [0] = 0 ;
  c0 = ' ' ;
  input = myInput ;

  for (;;)
    {
      c = cGet () ;
      if (inWord) inWord = 1 ; /* convenient way to reset nw = 0 */
      
      switch (c)
	{
	case 0:   /* end of file */
	  lineBreak (0) ;
	  return 0 ; 
	case '\\':  /* already done, do nothing */
	  c0 = 0 ;
	  continue ; 
	case '"':  /* export as is the whole string */
	case '\'':
	  space () ;
	  cTmp = c ;
	  nDoubleQuote = 1 ;
	  putchar(c) ; 
	  while (nDoubleQuote)
	    {
	      c = cGet () ;
	      switch (c)
		{
		case 0:   /* end of file */
		  return 0 ; 
		case '"':
		case '\'':
		  if (c == cTmp)
		    {
		      nDoubleQuote = 0 ;
		      break ;
		    }
		default:
		  putchar(c) ; 
		  break ;
		}
	    }
	  noSpace () ;
	  c0 = cPut(c) ;  /* ends by exporting the final double quote */
	  break ;
	case '{':
	  level++ ;
	  if (c0 != '\n')
	    lineBreak (level) ;
	  cPut(c) ; 
	  if (verbose)
	    printf (" /* C_B_block_start %d */", blocks [level] = nBlock++) ;
	  c0 = lineBreak (level) ;
	  break ;
	case '}':
	  lineBreak (--level) ;
	  cPut(c) ; 
	  if (verbose && level >= 0) 
	    printf (" /* C_B_block_end %d level=%d */",  blocks [level+1], level) ;
	  c0 = lineBreak (level) ;
	  if (level == 0) putchar ('\n') ;
	  break ;
	case '(':
	  nPar++ ;
	  c0 = cPut(c) ;
	  noSpace () ;
	  break ;
	case ')':
          nPar-- ;
	  noSpace () ;
	  c0 = cPut(c) ; 
	  break ;
	case ':': 
	  noSpace () ;
	  c0 = cPut(c) ; 
	  lineBreak (level) ;
	  break ;
	case '#':
	  lineBreak (0) ;
	  c0 = cPut(c) ; 
	  break ;
	case ';':
	  cPut(c) ;
	  if (!nPar) c0 = lineBreak (level) ; 
	  break ;
	case '\n':
	  c0 = lineBreak (level) ; 
	  break ;
	case ' ':
	  if (word[0]) space () ;
	  break ;
	default:
	  c0 = c ;
	  inWord = 2 ;
	  word [nw++] = c ;
	  word [nw] = 0 ;
	  break ;
	}
      if (inWord < 2)
	{ inWord = 0 ; nw = 0 ; }
    }
  return 0;
}

static void usage (void)
{
  fprintf (stderr,
	   "Usage: cbeau [ -s <number> ] [ -d ] [ files ]\n"
	   "  -s <number>:number of spce per level of identation, default = 2\n"
	   "  -d debug mode\n"
	   "  -f If file names are given, all files are contcatenated and beautified as one."
	   "     else beautify stdin\n\n"
	   ) ;
  exit (1) ;
}

/********************************************************************/
/********************************************************************/
/* utility to find and consume an argument on the unix command line */
/* get the existence of the arg */

int getArg (int *argcp, char **argv, char *arg)
{
  int i = *argcp ;

  for (i = 1 ; i < *argcp ; i++)
    if (argv[i] && !strcmp (arg, argv[i]))
      {
	for (;i < *argcp - 1 ; i++)
	  argv[i] = argv[i+1] ;
	*argcp -= 1 ;
	return  1 ;
      }
  return 0 ;
}

/*************************************************************************************/
/* get the value of the arg */
char* getArgV (int *argcp, char **argv, char *arg)
{
  int i = *argcp ;
  char *cp = 0 ;

  for (i = 1 ; i < *argcp - 1 ; i++)
    if (argv[i] && !strcmp (arg, argv[i]))
      {
	cp = argv[i+1] ;
	for (;i < *argcp - 2 ; i++)
	  argv[i] = argv[i+2] ;
	*argcp -= 2 ;
	break ;	
      }
  return cp ;
}

/*******************************************************************/
/*******************************************************************/

int main(int argc,char **argv)
{
FILE *f;
char *cp=0;

if((cp=getArgV(&argc,argv,"-p")))
space_per_level=atoi(cp);

debug=getArg(&argc,argv,"-d");

 if (argc == 1) usage() ;
if(argv[0])
{
while(argv[0])
{
f=fopen(argv[0],"r");
if(!f)
{
perror(argv[0]);
return 1;
}
else
{
beautify(f);
fclose(f);
}
argv++;
}
}
else
beautify(stdin);

putchar('\n');

return 0;
}

