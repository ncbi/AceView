/*  Last edited: Feb 22 21:38 1996 (rd) */
/* SCCS $Id: gd2gif.c,v 1.1.1.1 2002/07/19 20:23:17 sienkiew Exp $ */
#include "stdio.h"
#include "gd.h"

void main (int argc, char **argv)
{ 
  FILE *in, *out ;
  gdImagePtr im ;
  char buf[256] ;

  if (argc != 2)
    { fprintf (stderr, "Usage: gd2gif <name>\n"
	       "  converts <name>.gd to <name>.gif\n") ;
      exit (-1) ;
    }
  strcpy (buf, argv[1]) ;
  strcat (buf, ".gd") ;
  if (!(in = fopen (buf, "r")))
    { fprintf (stderr, "failed to open %s\n", buf) ;
      exit (-1) ;
    }
  if (!(im = gdImageCreateFromGd (in)))
    { fprintf (stderr, "failed to get image from %s\n", buf) ;
      exit (-1) ;
    }
  fclose (in) ;
  strcpy (buf, argv[1]) ;
  strcat (buf, ".gif") ;
  if (!(out = fopen (buf, "w")))
    { fprintf (stderr, "failed to open %s\n", buf) ;
      exit (-1) ;
    }
  gdImageGif (im, out) ;
  fclose (in) ;
}
