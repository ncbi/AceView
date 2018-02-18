#include    "../wac/ac.h"
#include    "vtxt.h"

/*******************************************************************/
/*******************************************************************/

/* match to template with wildcards.   Authorized wildchars are * ? #
     ? represents any single char
     * represents any set of chars
   case sensitive.   Example: *Nc*DE# fits abcaNchjDE23

   returns 0 if not found
           1 + pos of first sigificant match (i.e. not a *) if found
*/

static int doPickMatch (char *cp, char *tp, BOOL caseSensitive)
{
  char *c=cp, *t=tp;
  char *ts = 0, *cs = 0, *s = 0 ;
  int star=0;

  while (TRUE)
    switch(*t)
      {
      case '\0':
 /*
        return (!*c ? ( s ? 1 + (s - cp) : 1) : 0) ;
*/
	if(!*c)
	  return  ( s ? 1 + (s - cp) : 1) ;
	if (!star)
	  return 0 ;
        /* else not success yet go back in template */
	t=ts; c=cs+1;
	if(ts == tp) s = 0 ;
	break ;
      case '?':
	if (!*c)
	  return 0 ;
	if(!s) s = c ;
        t++ ;  c++ ;
        break;
      case '*':
        ts=t;
        while( *t == '?' || *t == '*')
          t++;
        if (!*t)
          return s ? 1 + (s-cp) : 1 ;
	if (!caseSensitive)
	  {
	    while (ace_upper(*c) != ace_upper(*t))
	      if(*c)
		c++;
	      else
		return 0 ;
	  }
	else
	  {
	    while (*c != *t)
	      if(*c)
		c++;
	      else
		return 0 ;
	  }
        star=1;
        cs=c;
	if(!s) s = c ;
        break;
      default:
	if (!caseSensitive)
	  {      
	    if (ace_upper(*t++) != ace_upper(*c++))
	      { if(!star)
		return 0 ;
	      t=ts; c=cs+1;
	      if(ts == tp) s = 0 ;
	      }
	    else
	      if(!s) s = c - 1 ;
	  }
	else
	  {
	    if (*t++ != *c++)
	      { if(!star)
		return 0 ;
	      t=ts; c=cs+1;
	      if(ts == tp) s = 0 ;
	      }
	    else
	      if(!s) s = c - 1 ;
	  }
        break;
      }
} /* pickMatch */

int pickMatch (char *cp,char *tp)
{ return doPickMatch (cp,tp,FALSE) ; }
int pickMatchCaseSensitive (char *cp,char *tp, BOOL caseSensitive)
{ return doPickMatch (cp,tp,caseSensitive) ; }

