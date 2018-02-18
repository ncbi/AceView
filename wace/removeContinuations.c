/*  Last edited: Aug  7 19:35 1994 (rd) */

/* $Id: removeContinuations.c,v 1.1.1.1 2002/07/19 20:23:11 sienkiew Exp $ */

/* simple filter to remove .ace file continuation lines */

#include <stdio.h>

void main (int argc, char **argv)
{
  char c, oldc ;

  if (feof (stdin))
    return ;
  oldc = getchar() ;
  while (!feof (stdin))
    { c = getchar() ;
      if (oldc == '\\' && c == '\n')
	{ while (!feof (stdin))
	    { oldc = getchar () ;
	      if (oldc != ' ' && oldc != '\t')
		break ;
	    }
	}
      else
	{ putchar (oldc) ;
	  oldc = c ;
	}
    }
  putchar (oldc) ;
}
