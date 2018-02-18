#include <stdio.h>
#include <string.h>
#include <stdlib.h>


void
puttext(char *t)
{
for (;;)
	{
	switch (*t)
		{
	case '\0': 	return;
	case '"': 	break;
	case '\t': 	putchar(' '); 
			break;
	default: 	putchar(*t);
		}
	t++;
	}
}

void
putdna(char *dna)
{
int n;
n = 0;
while (*dna)
	{
	if (*dna == '"')
		{
		dna++;
		continue;
		}
	putchar(*dna);
	dna++;
	n++;
	if (n >= 80)
		{
		printf("\n");
		n=0;
		}
	}
/* do not need a \n - there is one in the input string */
}


int main (int argc, char **argv)
{
  char b[1000000];
  char *cp, *cq, *dna;
  char dbchar, fmt;
  char gene[256], newname[256],product[256], title[1024];
  if (! argv[1])
    { fprintf(stderr,"must give dbchar\n"); exit(1); }
  dbchar=argv[1][0];
  
  if (! argv[2])
    { fprintf(stderr,"must give format\n"); exit(1); }
  fmt=argv[2][0];
  /* gene newname product dna */
  /* gene locusid product dna */

  if (fmt == 'd')
    while (fgets(b,sizeof(b),stdin))
      {
	cp = cq = b ;
	/* jump comments */
	if (*cp == '\n' || *cp == '/')
	  continue ;
	/* get gene name */
	cq = strstr(cp, "\t") ;
	if (! cq) continue ;
	*cq = 0 ;
	strncpy (gene, cp, 255) ;
	cp = cq = cq+1 ;
	/* get linkname (newname or locusid) */
	cq = strstr(cp, "\t") ;
	if (! cq) continue ;
	*cq = 0 ;
	strncpy (newname, cp, 255) ;
	cp = cq = cq+1 ;
	/* get product */
	cq = strstr(cp, "\t") ;
	if (! cq) continue ;
	*cq = 0 ;
	strncpy (product, cp, 255) ;
	cp = cq = cq+1 ;
	/* get title */
	cq = strstr(cp, "\t") ;
	if (! cq) continue ;
	*cq = 0 ;
	strncpy (title, cp, 1023) ;
	cp = cq = cq+1 ;
	/* get dna */
	dna = cp ;
	if (!*dna)
	  continue ;
	printf(">%c",dbchar) ;
	puttext(gene) ; putchar('|') ;
	puttext(newname) ; putchar('|') ;
	puttext(product) ; putchar('|') ;
	puttext(title) ; putchar('\n') ;
	putdna(dna);
      }
  else
    {
      fprintf(stderr,"must give 'p'(sorry:mieg undefined it 2004-10-14) or 'd' for format\n");
      return 1;
    }
  return 0;
}
