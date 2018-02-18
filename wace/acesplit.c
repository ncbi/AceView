/*  Last edited: Jan 23 18:21 1996 (mieg) */
/* stand alone, will split big ace files in a naive way */

/*   $Id: acesplit.c,v 1.2 2015/04/04 00:26:59 mieg Exp $   */
#include "ac.h"

int main (int argc, char **argv)
{ int i = 0, n = 0, nn = 0, nf = 0 ;
  char buf[10000], name[120], *cp ;
  FILE *out = 0 ;

  if (argc < 2)
    { fprintf (stderr, "Usage split prefix, splits stdin into n files of max 10^5 lines called prefix.n \n") ;
      exit(1) ;
    }

  strncpy (name, argv[1], 100) ;
  cp = name + strlen(name) ;
  sprintf(cp, ".%d", ++nf) ;
  out = fopen (name, "w") ;
  if (!out)
    { fprintf (stderr, "cannot ope output file\n") ;
      exit(1) ;
    }

  while (fgets(buf, 9999, stdin))
    { n++ ; nn++ ;
      if ((i = strlen(buf)) > 9990)
	{ fprintf(stderr, "Input too long: %d bytes in line %d\n", i, n) ;
	  exit(1) ;
	}
      
      if (nn > 500000 && i == 1)
	{ fclose (out) ;
	  sprintf(cp, ".%d", ++nf) ;
	  out = fopen (name, "w") ;
	  if (!out)
	    { fprintf (stderr, "cannot ope output file %s\n", name) ;
	      exit (1) ;
	    }
	  nn = 0 ;
	}
      else
	fputs(buf, out) ;
    }
  fprintf (stderr, "done, read %d lines, created %d files upto %s\n", n, nf, name) ;
  return 0 ;
}
