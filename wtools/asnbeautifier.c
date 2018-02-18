/*
* asnb - ASN.1 beautifier.
*
* usage: asnb [ -<number> ] [ -s ] [ -d ] [ files ]
*
*	The optional -<number> is the number of columns to indent for each
*	indent level.  default is -2.
*
*	-s means only use spaces for indenting.  By default, tabs are used
*	where possible.
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

int space_per_level=2 ;

int string_wrap_desire = 120 ;
int string_wrap_force = 200 ;

int no_tabs = 0 ;
int debug = 0 ;

/****************************************************************/

static void indent(int n)
{
  n = n * space_per_level ;
  if (! no_tabs)
    while (n >= 8)
      {
	putchar('\t') ;
	n -= 8 ;
      }
  while (n > 0)
    {
      putchar(' ') ;
      n-- ;
    }
}

/****************************************************************/
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
*/

static int level=0 ;
static int needspace=0 ;
static int lastspace=0 ;

/****************************************************************/

static int beautify(FILE *input)
{
  int c ;
  int qcount ;
  
  for ( ; ;)
    {
      c=getc(input) ; 
      if (c == EOF)
	goto end ;
      switch (c)
	{
	case '"':
	  /*
	   * quoted strings - Acedb ASN.1 export does not generate 
	   * strings with embedded quotes.  I don't know if ASN1
	   * has an escape mechanism, but I don't need it even
	   * if it does.
	   */
	  if (needspace)
	    putchar(' ') ;
	  putchar(c) ;
	  qcount=0 ;
	  for ( ; ;)
	    {
	      c=getc(input) ; 
	      if (c == EOF)
		goto end ;
	      if (c == '\n')
		{
		  putchar('~') ;
		  putchar('\n') ;
		  qcount=0 ;
		  continue ;
		}
	      if (c == '"')
		{
		  putchar('"') ;
		  break ;
		}
	      if (c == '{' || c == '}')
		{
		  putchar(c) ;
		  if (debug) { printf (" ERROR %c inside \"\"\n", c) ; exit (1) ;}
		  continue ;
		}
	      
	      if ( isspace(c) && (qcount > string_wrap_desire) )
		{
		  putchar('\n') ;
		  qcount=0 ;
		}
	      if ( qcount > string_wrap_force )
		{
		  putchar(c) ;
		  putchar('\n') ;
		  qcount=0 ;
		  continue ;
		}
	      putchar(c) ;
	      qcount++ ;
	    }
	  needspace=1 ;
	  lastspace=0 ;
	  break ;
	case '{':
	  if (!lastspace)
	    putchar(' ') ;
	  putchar('{') ;
	  putchar('\n') ;
	  indent(++level) ;
	  needspace=0 ;
	  lastspace=1 ;
	  break ;
	case '}':
	  putchar('\n') ;
	  indent(--level) ;
	  putchar('}') ;
	  needspace=1 ;
	  lastspace=0 ;
	  break ;
	case ',':
	  putchar(c) ;
	  putchar('\n') ;
	  indent(level) ;
	  needspace=0 ;
	  lastspace=1 ;
	  break ;
	case ' ':
	case '\t':
	case '\n':
	case '\f':
	  if (! lastspace)
	    needspace=1 ;
	  break ;
	default:
	  if (needspace && ! lastspace)
	    putchar(' ') ; 
	  putchar(c) ;
	  needspace=0 ;
	  lastspace=0 ;
	  break ;
	}
      
    }
 end:
  return 0 ;
}

/****************************************************************/

static int exitcode=0 ;

int main(int argc,char **argv)
{
  FILE *f ;
  
  if (argv[0]) argv++ ;
  
  while (argv[0] && argv[0][0] == '-')
    {
      switch (argv[0][1])
	{
	case '0': case '1': case '2': case '3': case '4':
	case '5': case '6': case '7': case '8': case '9':
	  space_per_level = atoi(argv[0]+1) ;
	  break ;
	case 's':
	  no_tabs = 1 ;
	  break ;
	case 'd':
	  debug = 1 ;
	  break ;
	default:
	  fprintf(stderr,"unknown option %s\n",argv[0]) ;
	  exitcode=1 ;
	}
      argv++ ;
    }
  
  if (argv[0])
    {
      while (argv[0])
	{
	  f=fopen(argv[0],"r") ;
	  if (!f)
	    {
	      perror(argv[0]) ;
	      exitcode=1 ;
	    }
	  else
	    {
	      beautify(f) ;
	      fclose(f) ;
	    }
	  argv++ ;
	}
    }
  else
    beautify(stdin) ;
  
  putchar('\n') ;
  
  return exitcode ;
}

/****************************************************************/
/****************************************************************/
/****************************************************************/
