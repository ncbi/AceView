#include <stdio.h>

#include <waceview/avlib.h>
#include <wh/regular.h>
#include "mytime.h"

struct av_db_desc *
get_db_desc(char *dbname, char **msg)
{
  FILE           *f;
  char            b[200], b1[100], b2[100], *s;
  static char return_message[500];
  struct av_db_desc desc, *ret;

  /*
  * Instead of writing all the possible database names and
  * parameters inline here, they are in the file db_config.
  * By using the current directory, the db_config file could
  * be fetched through the web server, but nobody knows it is
  * here and it doesn't contain anything important.
  */

  f = fopen("db_config", "r");
  if (!f) {
    sprintf(return_message, 
      "Database server for %s is currently unavailable - try again later", dbname);
    *msg = return_message;
    fprintf(stderr,"%s open failed\n", logPrefix("db_config"));
    return NULL;
  }

  /*
  * db_config lines are: 
  * ^# or blank - comment
  *
  * an entry is a line containing:
  *
  * dbname as passed to the cgi 
  * database type ("h" for human, "w" for worm)
  * hostname that runs gifacemblyserver 
  * rpc number to reach server 
  * species name
  * title 
  *
  * '_' characters in species name and title are converted to space. 
  * These fields exist to be displayed to the user.
  *
  */
  while (fgets(b, sizeof(b), f)) {
    if (b[0] == '#')
      continue;
    if (sscanf(b, "%s", b1) != 1)
      continue;
    if (strcmp(b1, dbname) == 0) {
      /*
      * ok, this line is for the database that we want - see
      * what it says
      */
      int             pgm;
      if (sscanf(b, "%s %c %s %s %s %s %s", desc.db_name, desc.db_species_char, desc.db_version, b2, b1, desc.db_species_title, desc.db_title) != 7)
        continue;

      if (strcmp(b2,"-") == 0) {
	sprintf(return_message,"database %s discontinued %s\n",b2,desc.db_title);
        *msg = return_message;
        fprintf(stderr,"%s\n", logPrefix("discontinued"));
	return NULL;
      }

      pgm = atoi(b1);
      sprintf(desc.db_access, "%s:%d", b2, pgm);

      for (s = desc.db_species_title; *s; s++)
      if (*s == '_')
        *s = ' ';

      for (s = desc.db_title; *s; s++)
      if (*s == '_')
        *s = ' ';

      ret = malloc(sizeof(desc));
      if (!ret) {
        sprintf(return_message, 
        "Database server for %s is currently unavailable - try again later", dbname);
	*msg = return_message;
        fprintf(stderr,"%s\n", logPrefix("get_db_desc_malloc"));
	return NULL;
      }
      memcpy(ret, &desc, sizeof(desc));
      return ret;
    }
  }

  /*
  * db not found in config file
  */
  sprintf(return_message,"Unknown database %s\n",dbname);
  *msg = return_message;
  fprintf(stderr,"%s %s\n",logPrefix("db_unknown"), dbname);
  return NULL;
}


char *pgetenv (char *a)
{
  char *s = getenv(a) ;
  if (!s)
    s = "(null)" ;
  return s ;
}

char *logPrefix (char *desc)
{
  char *s ;
  char *q, *rem, *p, *ref ;

  q = pgetenv("QUERY_STRING") ;
  rem = pgetenv("REMOTE_ADDR") ;
  p = pgetenv("PROXIED_IP") ;
  ref = pgetenv("HTTP_REFERER") ;
  s = messprintf("aceview: %d %s %s %s %s %s ", timeShowNow(),desc,q,rem,p,ref) ;
  return s ;
}
