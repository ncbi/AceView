/*  File: acesocket_.h
 *  Author: Ed Griffiths (edgrif@sanger.ac.uk) & Jean Thierry-Mieg (mieg@crbm.cnrs-mop.fr)
 *  Copyright (c) J Thierry-Mieg and R Durbin, 1999
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.crbm.cnrs-mop.fr
 *
 * Description: Private types/functions for acedb socket layer.
 * HISTORY:
 * Last edited: Jul 12 09:44 2000 (edgrif)
 * Created: Sat Jul 10 10:37:21 1999 (edgrif)
 * CVS info:   $Id: acesocket_.h,v 1.2 2017/02/15 20:39:41 mieg Exp $
 *-------------------------------------------------------------------
 */
#ifndef DEF_SOCKET_PRIV_H
#define DEF_SOCKET_PRIV_H

#include <sys/types.h>
#include <sys/socket.h>
#include <sys/ioctl.h>
#include <sys/time.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <netdb.h>
#include <errno.h>
#include <signal.h>

#include <wh/regular.h>
#include <wh/mytime.h>
#include <wsocket/servertransport.h>			    /* Defines interface required by acedb */
							    /* server. */


/* Messages:                                                                 */
/*                                                                           */

typedef struct S_MESSAGESTRUCT
{
  BOOL isNew ;                   /* controls sMessageSetRead/Write */
  BOOL doneHeader  ;
  struct _ACE_HEADER ah ;
  /* header handling */
  char hBuffer [ACE_HEADER_BYTES] ;
  int hBytesRequested ;
  int hBytesPending ;
  /* message handling */
  int mBytesRequested ;
  int mBytesPending ;
  char *readBuffer ;              /* filled by the client, known length */
  Stack writeStack ;              /* filled by the server, dynamic length */ 
} *S_MESSAGE ;




/* Messages are in one of these states while being sent, this is because even*/
/* without non-blocking I/O you can't guarantee sending the message in one   */
/* go if it is very big.                                                     */
typedef enum _S_MSGState {SMSG_WAIT, SMSG_DONE, SMSG_ERROR} S_MSGState ;


/* Clients:                                                                  */
/*                                                                           */
/* All clients must be in one of these states.                               */
typedef enum _ClientConnectState {CLCON_EMPTY, CLCON_CONN, CLCON_READWAIT,
				  CLCON_REQUEST, CLCON_WRITEWAIT, CLCON_CLOSE
                                 } ClientConnectState ;

/* Describes why a client is to be disconnected:                             */
typedef enum _ClientCloseReason {CL_NONE, CL_TIMEDOUT, CL_ERROR, CL_NORMAL} ClientCloseReason ;

/* client internal state                                                     */
typedef struct _ClientConnectStruct *S_CLIENT ;
typedef struct _ClientConnectStruct
{
  /* private descriptors */
  int magic ;
  ClientConnectState state ;
  ClientCloseReason close ;
  int fd ;						    /* socket number */
  mytime_t time ;					    /* Used to time out a client. */
  struct sockaddr_in socket ;
  

  S_MESSAGE message ;

  BOOL mustCallAce ;					    /* TRUE if ace must be called before
							       killing this client. */
  S_CLIENT previous, next ;				    /* client doubly linked list */
} ClientConnectStruct ;


/* message operations:   */

S_MESSAGE sMessageCreate(void) ;
void sMessageInit(S_MESSAGE msg) ;
void sMessageSetMagic(S_MESSAGE msg) ;
int  sMessageGetMagic(S_MESSAGE msg) ;
void sMessageSetRead(S_MESSAGE msg) ;
void sMessageSetWrite(S_MESSAGE msg) ;
S_MSGState sMessageSocketRead(int sockfd, S_MESSAGE message_p) ;
S_MSGState sMessageSocketWrite(int sockfd, S_MESSAGE message_p) ;
int sMessageGetMessage(S_MESSAGE msg, char **msg_text) ;
void sMessageDestroy(S_MESSAGE msg) ;
void directReply (S_MESSAGE message, char *text) ;



/* Socket definitions:                                                       */
/*                                                                           */
/* Some UNIX's have problems over various types for sockets, some are only   */
/* just being introduced. Here we do our own defines to help out with this.  */
/* (which I have copied from the sample code from Stevens "UNIX Network      */
/* Programming").                                                            */
/*                                                                           */
/* Dreadful hack for alphas which suddenly defined socklen_t somewhere       */
/* between levels  OSF1 V4.0 1229 & OSF1 V4.0 878                            */
/* And also for Suns which defined it somewhere between 5.5 & 5.7            */
#if     (defined(ALPHA) && !defined(SO_UMC))             \
     ||  defined (SGI)                                   \
     || (defined (SOLARIS) && !defined(_SOCKLEN_T))      \
     ||  defined(__CYGWIN__) 
typedef int  socklen_t ;
#endif

/* Following could be derived from SOMAXCONN in <sys/socket.h>, but many
   kernels still #define it as 5, while actually supporting many more */
#define	LISTENQ		1024				    /* 2nd argument to listen() */

/* absent on linux, definition copied from dec alpha */
#ifndef SHUT_RDWR
#define SHUT_WR       1					    /* Disables further send operations */
#define SHUT_RDWR     2					    /* Disables further send and receive */
							    /* operations */
#endif


/* Socket operations:                                                        */
/*                                                                           */
void Connect(int fd, const struct sockaddr *sa, socklen_t salen) ;
int aceSocketSelect(int nfds, fd_set *readfds, fd_set *writefds, fd_set *exceptfds, 
	   /*struct timeval*/ void *timeout) ;
int Fcntl(int fd, int cmd, int arg) ;

void Shutdown(int fd, int mode) ;

#endif /* DEF_SOCKET_PRIV_H */
