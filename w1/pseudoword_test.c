#include <ctype.h>

#include "regular.h"

int main()
{
char b[100], buffer [80] ;
int c,d;
int x;
enum pseudoword_phoneme_set ph;
char *s;

ph = pseudoword_english;

 while (fgets(b,99,stdin))
	{
	s = b;
	if (*s == '=')
		{
		s++;
		switch (*s)
			{
		case 'e':  ph = pseudoword_english; s++; break;
		case 'j':  ph = pseudoword_japanese; s++; break;
			}
		s++;
		}
	
	if (*s == '/')
		{
		c=1;
		d=500;
		if (sscanf(s+1,"%d%d",&c,&d) == 1)
			{ d = c + 100 ; }
		for (x=c; x < d; x++)
			{
			printf("%4d %-10s ",x,generate_pseudoword(buffer,x,ph));
			if (x % 4 == 0) printf("\n");
			}
		printf("\n");
		}
	else if (*s == '!')
		{
		c = 0; d = 10000000;
		if (sscanf(s+1,"%d%d",&c,&d) == 1)
			{ d = c; c = 0; }
		for (x=c; x<d; x++)
			{
			char *s;
			s = generate_pseudoword(buffer,x,ph);
			if (decode_pseudoword(s,ph) != x) { printf("%d %s\n",x,s); abort(); }
			if (x % 100000 == 0) printf("at %d\n",x);
			}
		printf("done\n");
		}
	else if (isdigit(s[0]))
		{
		printf("english  %s\n",generate_pseudoword(buffer,atoi(s),pseudoword_english));
		printf("japanese %s\n",generate_pseudoword(buffer,atoi(s),pseudoword_japanese));
		}
	else
		{
		printf("english  %d\n",decode_pseudoword(s,pseudoword_english));
		printf("japanese %d\n",decode_pseudoword(s,pseudoword_japanese));
		}
	}
return 0;
}
