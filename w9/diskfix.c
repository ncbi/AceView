/*  File: diskfix.c
 *  Author: Richard Durbin (rd@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1993
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: to reset the session number and master block in
 	the super block!
 * Exported functions:
 * HISTORY:
 * Last edited: Oct 16 15:46 1998 (fw)
 * Created: Thu Aug 12 00:51:59 1993 (rd)
 *-------------------------------------------------------------------
 */

/* $Id: diskfix.c,v 1.2 2016/03/21 22:26:01 mieg Exp $ */

#include <fcntl.h>
#include <errno.h>
#include "acedb.h"
#include "disk_.h"
#include "regular.h"

#define BLKMX (BLOC_SIZE - sizeof(BLOCKHEADER) - sizeof(DISK)\
	       -sizeof(int) -sizeof(int) -sizeof(int))

typedef struct block /* sizeof(block) must be < BLOC_SIZE */
  {
    BLOCKHEADER  h ;
    DISK gAddress ;
    int mainRelease ,  subDataRelease, subCodeRelease ;
    char c[BLKMX] ;
  }
    BLOCK, *BLOCKP;   /*the transfer unit between disk and cache*/

/*******************************************/

int main (int argc, char **argv)
{
  BLOCK b ;
  char a, filename[1024] ;
  int fd ;
  myoff_t pos, pp ;

  if (argc != 3)
    messcrash ("Usage: diskfix session blockno\n") ;

  if (!getenv ("ACEDB"))
    messcrash ("Set ACEDB\n") ;
  sprintf (filename, "%s/database/blocks.wrm", getenv("ACEDB")) ;
  if ((fd = open (filename, O_WRONLY | O_BINARY)) == -1)
    messcrash ("Could not open %s\n", filename) ;
  else
    printf ("Opened %s\n", filename) ;

  b.h.disk = 1 ;
  b.h.session = atoi(argv[1]) ;
  b.gAddress = atoi(argv[2]) ;
  b.mainRelease = 2 ;

  printf ("Do you want session %d, lexi1 address %d? (y/n) ",
	  b.h.session, b.gAddress) ;
  scanf ("%c", &a) ;
  if (a == 'y' || a == 'Y')
    { pos = BLOC_SIZE ;
      if ((pp = lseek(fd, pos, SEEK_SET)) != pos)
	messcrash ("lseek failed: pos=%ld, returned=%ld\n",
		   (long) pos, (long) pp) ;
      if ((pp = write(fd,&b,BLOC_SIZE)) != BLOC_SIZE)
	messcrash ("Write failed: size=%d, returned=%d\n", 
		   BLOC_SIZE, pp) ;
      printf ("Write change succeeded\n") ;
    }
  return 0 ;
}
