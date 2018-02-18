/*  File: saceclient_.h
 *  Author: Ed Griffiths (edgrif@sanger.ac.uk) & Jean Thierry-Mieg (mieg@crbm.cnrs-mop.fr)
 *  Copyright (c) J Thierry-Mieg and R Durbin, 1999
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.crbm.cnrs-mop.fr
 *
 * Description: Private header for saceclient.
 * HISTORY:
 * Last edited: Apr 27 16:49 2000 (edgrif)
 * Created: Tue Sep 14 09:42:10 1999 (edgrif)
 * CVS info:   $Id: saceclient_.h,v 1.1.1.1 2002/07/19 20:23:35 sienkiew Exp $
 *-------------------------------------------------------------------
 */
#ifndef DEF_CLIENT_PRIV_H
#define DEF_CLIENT_PRIV_H


#include <wsocket/acesocket_.h>				    /* WHY DO WE NEED THIS ??????? */

enum {SOCKET_RETRIES = 100} ;


/* A set of routines to get userid/passwd/domain name stuff from the user    */
/* on the usual terminal stdin/stdout.                                       */

/* Prompt user for userid/passwd.                                            */
BOOL getUseridPasswd(char **userid, char **passwd) ;

/* Prompt user for                                                           */
char *getNewPasswd(char *userid, char *curr_passwd) ;

BOOL getUserUpdate(char **cp) ;

BOOL getDomainUpdate(char **cp) ;


/* These should be in the socket lib....                                     */
S_MSGState writeToSocket(int sock, void *handle, char *msg_type,
			 char *request) ;
S_MSGState readFromSocket(int sock, void *handle, char **msg_type_out,
			  char **answerp, int *lengthp) ;


extern BOOL debug_G ;					    /* Global debug flag for client code. */


#endif /* DEF_CLIENT_PRIV_H */
