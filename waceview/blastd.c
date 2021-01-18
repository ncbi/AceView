#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <sys/wait.h>
#include <time.h>
#include <sys/resource.h>
#include "../wtcp/tcp_connect.h"

/*
 * Author Mark Sienkiewicz, jan 2003
 *
 * blastall will happily produce incorrect output if _any_ _one_ of
 * the databases it is asked to search is missing...
 */

/*
 * The db list is going to be like "n.worm n.human n.ara ".  BIG is
 * more than enough characters to hold the whole list.
 */

#define BIG 200

static int check_a_database (char *s)
{
  char b[BIG+20] ;
  if (! *s)
    return 0 ;
  strcpy (b,s) ;
  strcat (b,".nhr") ;
  if (access (b,0) < 0)
    return 1 ;
  return 0 ;
}

static void check_dbnames (char *s)
{
  char b1[BIG] ;
  char b2[BIG] ;
  char b3[BIG] ;
  char b4[BIG] ;
  if (strlen (s) > BIG)
    goto error ;
  b1[0]=0 ;
  b2[0]=0 ;
  b3[0]=0 ;
  b4[0]=0 ;
  if (sscanf (s,"%s%s%s%s",b1,b2,b3,b4) == 0)
    goto error ;
  if (check_a_database (b1)) 
    goto error ;
  if (check_a_database (b2)) 
    goto error ;
  if (check_a_database (b3)) 
    goto error ;
  if (check_a_database (b4)) 
    goto error ;
  return ;
 error:
  printf ("error: database %s %s %s %s files missing\n", b1, b2, b3, b4) ;
  exit (1) ;
}



static int service (int s)
{
  char b[100] ;
  char *blasttype, *evalue, *dbnames, *t, *search_type ;
  int n ;
  time_t tyme ;
  char *time_str ;

  wait3 (NULL,WNOHANG,NULL) ;
  memset (b,0,sizeof (b)) ;
  n = read (s, b, sizeof (b)) ;
  if (n < 0)
    perror ("read") ;
  b[sizeof (b)-1] = 0 ;
  evalue=b ;
  search_type=strchr (b,' ') ;
  if (!search_type)
    return 0 ;
  *search_type++ = 0 ;
  dbnames = strchr (search_type,' ') ;
  *dbnames++ = 0  ;
  while ((t = strchr (dbnames, '_')))
    *t = ' ' ;
  if (*search_type == 'p')
    blasttype = "tblastn" ;
  else
    blasttype = "blastn" ;
  tyme=time (0) ;
  time_str=ctime (&tyme) ;
  time_str[24]=0 ;
  printf ("%s type %s e %s db %s\n",time_str,blasttype, evalue, dbnames) ;
  fflush (stdout)  ;
  switch (fork ())
    {
    case 0:	break ;
    case -1:
      printf ("error: unexpected error in server\n") ;
    default:
      return 0 ;
    }
  close (0) ;
  close (1) ;
  dup (s) ;
  dup (s) ;
  check_dbnames (dbnames) ;
  /*  /netopt/ncbi_tools64/ver0.0/ncbi/bin/blastall */
  execl ("echo" , "echo" , "hello vahan world" , NULL) ;
  execl ("/home/mieg/SERVER/BLAST/blastall",
	 "blastall",
	 "-F","\"m L\"",
	 "-e",evalue,
	 "-p",blasttype,
	 "-d",dbnames,
	 "-T",
	 NULL) ;
  /* arrive here  ONLY if execl failed */
  /* execlp will try the command on on your path components */
  execlp ("/home/mieg/SERVER/BLAST/blastall"
	 ,"blastall",
	 "-F","\"m L\"",
	 "-e",evalue,
	 "-p",blasttype,
	 "-d",dbnames,
	 "-T",
	 NULL) ;
  /* arrive here  ONLY if execlp failed */  
  write (1,"<h1>Fnternal error</h1>execl and execlp failed\n",47) ;
  exit (1) ;
  return 0 ;
}

int main (int argc, char **argv)
{
  int port = 3001 ;
  if (argv[1])
    port = atoi (argv[1]) ;
  daemon_main (port, 1, service) ;
  return 0 ;
}

