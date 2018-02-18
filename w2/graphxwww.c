/*  File: graphxwww.c
 *  Author: Fred Wobus (fw@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1998
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * SCCS: %W% %G%
 * Description: function callMosaic to display a given link
      in an external browser. 

      At present this is implemented by callNetscape, which is
      using the functions in xnetscape_remote.c to drive
      netscape remotely in an UNIX-X11 environment.
 * Exported functions:
      BOOL callMosaic (char *link);
      BOOL callNetscape (char *link);
 * HISTORY:
 * Last edited: Jun 25 15:44 1999 (fw)
 * Created: Wed Oct 14 14:18:22 1998 (fw)
 *-------------------------------------------------------------------
 */

#if !(defined WIN32 || defined MACINTOSH)

#include "regular.h"
#include <X11/Xlib.h>

const char *expected_mozilla_version = "3.0"; /* extern'd by
						 xnetscape_remote.c */
const char *progname = 0;     /* extern'd by xnetscape_remote.c */

extern int mozilla_remote_commands (Display *dpy, 
				    Window window, 
				    char **commands);

BOOL graphWebBrowser(char *link)
/* start netscape navigator and remote control its display
   via x-atoms placed on the display and acted upon by an 
   active netscape process */
{
  char *command_buf[2];
  Display *dpy;
  char *dpy_string =0;
  int status;
  
  if (!link) return FALSE;
  if (strlen(link) == 0) return FALSE;

  dpy = XOpenDisplay(dpy_string);

  if(!dpy)
    return FALSE;

  progname = messGetErrorProgram();

  command_buf[0] = messalloc (10+strlen(link) + 1);
  command_buf[1] = "";

  if (*link == SUBDIR_DELIMITER)
    /* is it a filename path ? */
    sprintf(command_buf[0], "openFILE(%s)", link);
  else
    sprintf(command_buf[0], "openURL(%s)", link);

  /* try to submit the commands to an already opened netscape window */
  status = mozilla_remote_commands(dpy,	/* this display */
				   0, /* first netscape window */
				   command_buf); /* remote-commands */

  messfree (command_buf[0]);

  if (status != 0) 
    {
      /* we don't seem to be able to connect to a running 
       * netscape process. We could try and start a process
       * using a system call. However, the different ways
       * in which each individual user starts netscape
       * is impossible to guess. To it'd be best just
       * to ask him/her to start it how they usually do */
      messout("Netscape doesn't appear to be running.\n"
	      "Please start up Netscape in the usual way "
	      "and try again !");
    }


  return (status != 0);
} /* graphWebBrowser */


#elif !defined(WIN32) /* MACINTOSH ONLY, WIN32 callMosaic() implemented elsewhere */

BOOL graphWebBrowser (char *link)
{
  return FALSE ;
}

#endif /* MACINTOSH ONLY */

/********************** End of file **********************/
 
