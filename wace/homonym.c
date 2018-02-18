/*  Last edited: Sep  8 17:00 1994 (srk) */

/* $Id: homonym.c,v 1.3 2016/01/26 04:02:50 mieg Exp $ */

  /* The pupose of this stand alone is to resolve the orthograph conflicts
     in author names. It takes as  arguments, a homonym file of the form

correct_name AA =
alias1 posibly with several words =
alias2 no double quotes expected 
dummyentry
newname =
newalias

   and reads from an ace file.
   It will act only on the author entries of the ace file and, on these
   will replace the occurences of alias by the correct name 
   */
#define ARRAY_NO_CHECK
#include "regular.h"
#include "array.h"

static Stack ss = 0 ;
static Array alias = 0 ;
static int line = 0, transfo = 0 , nnew = 0 ;

typedef struct {int bad; int good; } ALIAS ;

/****************************************************************/

static void scanInput(void)
{  char *cp , cutter ; register int i , n ;
   FILE *new = fopen("newnames","w") ;
   
   if (!new)
     messcrash("Can't open the newname file") ;
  freespecial ("\n\"\\") ;
  fprintf(stderr,"I begin the scan line %d\n", line) ;
  while(line++ , freeread(stdin))
    { cp = freeword() ;
      if(!cp)
	{ printf("\n") ; continue ;}
      if (   strcmp("Author",cp)
          && strcmp("Mapper",cp)
          && strcmp("From_Author",cp)
          && strcmp("Representative",cp)
	  )  /* Just copy the whole line to the output */
	{ freeback() ;
	  printf("%s", freewordcut("\n",&cutter)) ;
	  printf("\n") ;
	  continue ;
	}
          /* Found author, do check */
      printf("%s : \"",cp) ;
      freenext() ; freestep(':') ; freenext() ;
      cp = freeword() ;
      if(cp)
	{ i = arrayMax(alias) ;
	  while(i--)
	    if(!strcmp(cp,stackText(ss,arr(alias,i,ALIAS).good)))
	      goto found;
	    else if( (n = arr(alias,i,ALIAS).bad) &&
		      !strcmp(cp,stackText(ss,n)))
	      { cp = stackText(ss,arr(alias,i,ALIAS).good) ;
		transfo++ ;
		goto found ;
	      }
	  nnew++ ;
	  fprintf(new,"%s\n",cp) ;
	found: 
	  printf("%s", cp) ;
	}
      printf("\"\n") ;
    }
   filclose(new) ;
}


/****************************************************************/

static BOOL readRules(char *name)
{
  FILE *f = fopen(name,"r") ;
  int nA = 0 , n ;
  char *cp , *cq , cutter ;
  BOOL found = FALSE ;
  int lastnA = 0 ;

  if(!f)
    { fprintf(stdout, "Cannot open the homonym file %s", name) ;
      return FALSE ;
    }

  ss = stackCreate(1000) ;
  pushText(ss,"dummy") ;
  alias = arrayCreate(100, ALIAS) ;

  freespecial ("\n\t\"\\") ;

  while (freeread(f))
    if((cp = freewordcut("\n",&cutter)))
     { cq = cp + strlen(cp) - 1 ;
       while(cq > cp && *cq == ' ') cq--;
       if(*cq == '=')
	 { *cq = ' ' ;
           found = TRUE ;
	   while(cq > cp && *cq == ' ') cq--;
	   *++cq = 0 ;
	   n = stackMark(ss) ;
	   pushText(ss,cp) ;
           array(alias, nA++,ALIAS).bad = n ;
	 }
       else
	 if (found)
	   { found = FALSE ;
	     *++cq = 0 ;
	     n = stackMark(ss) ;
	     pushText(ss,cp) ;
	     for ( ;lastnA <nA; lastnA++)
	       array(alias, lastnA,ALIAS).good = n ;
	   }
	 else
	   { array(alias, nA++,ALIAS).bad = 0 ;
	     *++cq = 0 ;
	     n = stackMark(ss) ;
	     pushText(ss,cp) ;
	     for ( ;lastnA <nA; lastnA++)
	       array(alias, lastnA,ALIAS).good = n ;
	   }
     }
  fprintf(stderr,"Found %d alias to eliminate\n",nA) ;
  filclose(f) ;
  return nA > 0 ? TRUE : FALSE ;
}

/****************************************************************/

int main (int argc, char **arg)
{
  if (argc != 2)
    { fprintf (stderr,
	       "Usage : homonym homoFile < input > output \n") ;
      return 1  ;
    }
 
  freeinit() ;
  if(readRules(arg[1]))
    scanInput() ; 

  fprintf(stderr,"\n done, scanned %d lines, did %d transformations\n found %d new names ",
	  line, transfo, nnew ) ;
  return 0 ;
}

/****************************************************************/
/****************************************************************/



