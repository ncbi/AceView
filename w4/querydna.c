#include "regular.h"
#include "query.h"
#include <regex.h>


static void queryRegExpFinalise (void *queryRegExp)
{
  if (queryRegExp)
    regfree ((regex_t *) queryRegExp) ; 
}

/* 0: bad query, !NULL: use in queryRegExpFind, then call queryRegExpDoFree */
QueryRegExp *queryRegExpCreate (const char *regExpPattern, BOOL getPos, AC_HANDLE h)
{
  int nn ;
  regex_t *queryRegExp = halloc (sizeof (regex_t), h) ;

  blockSetFinalise (queryRegExp, queryRegExpFinalise) ;
  
  /* man regcomp for details */
  if (getPos)
    nn = regcomp (queryRegExp, regExpPattern, REG_EXTENDED | REG_ICASE) ;
  else
    nn = regcomp (queryRegExp, regExpPattern, REG_EXTENDED | REG_ICASE | REG_NOSUB ) ;
  if (nn) /* error */
    {
      messout ("queryRegExpMatch: sorry, I cannot understand regexp: %s", regExpPattern) ;
      if (queryRegExp)
	messfree (queryRegExp) ;
    }
  return queryRegExp ;
}

/* 0: not found, >0 position */
int queryRegExpFind (const char *data, QueryRegExp *vRegP, BOOL getPos)
{
  BOOL ok = FALSE ;
  int nn ;
  regex_t *regp = 0 ;
  regmatch_t pmatch ;

  /* man regcomp for details */
  if (vRegP)
    {
      regp = (regex_t *) vRegP ;
      if (data)
	{
	  if (getPos)
	    nn = regexec (regp, data, 1, &pmatch, 0) ;
	  else
	    nn = regexec (regp, data, 1, 0, 0) ;
	  if (nn) { /* report not found */ ; ok = 0 ; }
	  else  /* match found */ 
	    { 	     
	      if (getPos)
		ok = 1 + pmatch.rm_so ;
	      else
		ok = 1 ;
	    }
	}
      else
	ok = 0 ;
    }
  else
    ok = -1 ;
  return ok ;
}

int queryRegExpMatch (const char *data, const char *regExpPattern, BOOL getPos)
{
  int ok = 0 ;
  QueryRegExp *regp = 0 ;

  if ((regp = queryRegExpCreate (regExpPattern, getPos, 0)))
    {
      ok = queryRegExpFind (data, regp, getPos) ; /* 0: not found, >0 position */
      queryRegExpDestroy (regp) ;
    }
  else
    ok = -1 ; /* query is wrong */
  return ok ;
}


int queryRegExpMatch2 (const char *data, const char *regExpPattern, BOOL getPos)
{
  BOOL ok = FALSE ;
  int nn ;
  regex_t preg ;
  regmatch_t pmatch ;

  /* man regcomp for details */
  if (getPos)
    nn = regcomp (&preg, regExpPattern, REG_EXTENDED | REG_ICASE) ;
  else
    nn = regcomp (&preg, regExpPattern, REG_EXTENDED | REG_ICASE | REG_NOSUB ) ;
  if (!nn) /* no error */
    {
      if (data)
	{
	  if (getPos)
	    nn = regexec (&preg, data, 1, &pmatch, 0) ;
	  else
	    nn = regexec (&preg, data, 1, 0, 0) ;
	  if (nn) { /* report not found */ ; ok = 0 ; }
	  else  /* match found */ 
	    { 	     
	      if (getPos)
		ok = 1 + pmatch.rm_so ;
	      else
		ok = 1 ;
	    }
	}
      else
	ok = 0 ;
    }
  else
    {
      ok = -1 ;
      if (! data)
	messout ("queryRegExpMatch: sorry, I cannot understand regexp: %s", regExpPattern) ;
    }
  regfree (&preg) ;
  return ok ;
}
