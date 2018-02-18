/*  File: acesocketlib.c
 *  Author: Ed Griffith (edgrif@sanger.ac.uk) & Jean Thierry-Mieg (mieg@crbm.cnrs-mop.fr)
 *  Copyright (c) J Thierry-Mieg and R Durbin, 1999
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.crbm.cnrs-mop.fr
 *
 * Description: Utility routines for the ace socket layer implementation.
 * Exported functions: See acesocket_.h
 * HISTORY:
 * Last edited: May 25 12:00 2000 (edgrif)
 * Created: Sat Jul 10 10:37:21 1999 (edgrif)
 * CVS info:   $Id: acesocketlib.c,v 1.2 2017/02/15 20:39:42 mieg Exp $
 *-------------------------------------------------------------------
 */

#include <wh/regular.h>
#include <wsocket/acesocket_.h>


/* Internal routines.                                                        */
static S_MSGState readSocket(int fd, char *readbuf, int *bytesPending) ;
static S_MSGState writeSocket(int fd, char *writebuf, int *bytesPending) ;


/**************************************************************/

/* This is for IPv4, some work needed to make it independent need to get     */
/* quite a lot of stephens code to make it protocol independent...           */
/* bad, bad, bad...need IPv6...as well...                                    */
/*                                                                           */
BOOL getClientInetAddressInfo(ClientConnect connect_struct,
			      char **dotted_name, char **host_name)
{
  BOOL name_found = FALSE ;
  struct sockaddr_in *cliaddr = (struct sockaddr_in *)&(connect_struct->socket) ;
  struct hostent *host_entry = NULL ;
  char *client_dotted_decimal = NULL ;

  /* Get dotted decimal name, this call cannot fail.                         */
  client_dotted_decimal = inet_ntoa(cliaddr->sin_addr) ;

  /* Now get the string form of the host name, which can fail.               */
  host_entry = gethostbyaddr((const char *)&(cliaddr->sin_addr), sizeof(struct in_addr),
			     cliaddr->sin_family) ;
  if (host_entry)
    {
      *host_name = strnew(host_entry->h_name, 0) ;
      *dotted_name = strnew(client_dotted_decimal, 0) ;
      name_found = TRUE ;
    }
  else
    {
      messdump("Lookup of client hostname translation from binary connection number failed, "
	       "dotted decimal host address of client was %s, "
	       "and reason was: %s",
	       client_dotted_decimal, messErrnoText(h_errno)) ;
    }

  return name_found ;
}



/************************************************************/

/* Return TRUE if we were started by inetd, FALSE otherwise. The function    */
/* only checks this the first time through, subsequent calls then just get   */
/* the stored boolean.                                                       */
/*                                                                           */
/* If we were started by inetd then file descriptors 0, 1 and 2 must all be  */
/* the same internet socket.                                                 */
/* This code determines this for file descriptors 0 to 2 by:                 */
/*                                                                           */
/*   1) checking that they are all sockets                                   */
/*   2) checking that the sockets are internet (AF_INET) sockets             */
/*   3) checking that all three refer to the same socket address/port number */
/*                                                                           */
/* Don't know if you can do much more than that.                             */
/*                                                                           */
/* Optionally we return the port number so that this knowledge is all        */
/* encapsulated in this routine.                                             */
/*                                                                           */
/* Note: on windows we don't support inetd so we always return FALSE here.   */
BOOL isInetdDaemon(int *port_out)
{
#ifdef __CYGWIN__
  return FALSE;
#else
  static BOOL first_time = TRUE ;
  static BOOL inetdStarted = FALSE ;
  static unsigned short port ;
  int fd ;
  struct sockaddr_in sockaddr ;
  int sockaddrlen ;
  unsigned long  addr ;


  if (first_time)
    {
      first_time = FALSE ;

      for (fd = 0, inetdStarted = TRUE ; fd <= 2 && inetdStarted == TRUE ; fd++)
	{
	  sockaddrlen = sizeof sockaddr;

	  /* A better call to use here would be isfdtype() which Posix.1g    */
	  /* defines should be in sys/stat.h, but most of our systems don't  */
	  /* have this yet...sigh... (see p.81 in Stephens Network Vol 1)    */
	  if (getsockname(fd, (struct sockaddr *) &sockaddr, &sockaddrlen) == -1)
	    {
#ifdef SOLARIS
	      /* for some reason Solaris 5.5 returns EINVAL, this disagrees  */
	      /* with the manpage...sigh. Maybe there is some problem on     */
	      /* solaris with doing a getsockname on a listening socket...   */
	      if (errno == EINVAL)
		inetdStarted = FALSE ;
	      else
#endif
	      if ((errno == EBADF) || (errno == ENOTSOCK)
		  || (errno == EOPNOTSUPP) || (errno == ENOENT))
		inetdStarted = FALSE ;
	      else
		  messcrash("Cannot determine if server started by inetd "
			    "because call to getsockname() failed: %s", messSysErrorText()) ;
	    }

	  if (inetdStarted)
	    {
	      if (sockaddr.sin_family != AF_INET)
		inetdStarted = FALSE ;
	    }

	  if (inetdStarted)
	    {
	      if (fd == 0)
		{
		  addr = ntohl(sockaddr.sin_addr.s_addr);
		  port = ntohs(sockaddr.sin_port);
		}
	      else
		{
		  if ((addr != ntohl(sockaddr.sin_addr.s_addr)) ||
		      (port != ntohs(sockaddr.sin_port))) 
		    {
		      inetdStarted = FALSE ;
		    }
		}
	    }
	}
    }

  /* return the listening port number for the daemon.                        */
  if (inetdStarted && port_out != NULL)
    *port_out = port ;

  return inetdStarted ;

#endif /* !__CYGWIN__ */
}
/************************************************************/
/* little utility to mimick a genuine server reply */

void directReply (S_MESSAGE message, char *text)
{
  Stack s = stackCreate (strlen(text)+1) ;
  pushText (s, text) ;  
  stackDestroy (message->writeStack) ;
  message->writeStack = s ;
  message->isNew = TRUE ;
}


/**************************************************************/

static void buf2hdr(char *buf, ACE_HEADER hh)
{
  memcpy(&(hh->swapMagic), buf, 4) ; buf += 4 ;
  memcpy(&(hh->length), buf, 4) ;  buf += 4 ;
  memcpy(&(hh->serverVersion), buf, 4) ; buf += 4 ;
  memcpy(&(hh->clientId), buf, 4) ; buf += 4 ;
  memcpy(&(hh->maxBytes), buf, 4) ; buf += 4 ;
  memcpy(&(hh->msgType), buf, ACESERV_MSGTYPE_BUFLEN) ;
}

/**************************************************************/

static void hdr2buf(ACE_HEADER hh, char *buf)
{
  memcpy(buf, &(hh->swapMagic), 4) ; buf += 4 ;
  memcpy(buf, &(hh->length), 4) ;    buf += 4 ;
  memcpy(buf, &(hh->serverVersion), 4) ; buf += 4 ;
  memcpy(buf, &(hh->clientId), 4) ; buf += 4 ;
  memcpy(buf, &(hh->maxBytes), 4) ; buf += 4 ;
  memcpy(buf, &(hh->msgType), ACESERV_MSGTYPE_BUFLEN) ;

}

/**************************************************************/

static void bufferDoSwap (char *buf)
{
  int ii , nn = ACE_HEADER_SWAPFIELDS ;
  int *ip = (int*) buf ;
  char *cp, *cq ;
  
  /* First int in buf is the magic tag which must not be swapped so skip it. */  
  /* Swop all the rest.                                                      */
  nn-- ;
  while (ip++, nn--)
    {
      ii = *ip ;
      cp = (char *)(&ii) ;
      cq = (char*)ip ;
      *cq = *(cp+3) ;
      *(cq+1) = *(cp+2) ;
      *(cq+2) = *(cp+1) ;
      *(cq+3) = *(cp) ;
    }

  return ;
}

/**************************************************************/

static BOOL bufferToHeader (char *buf, ACE_HEADER hh)
{
  int magic = *(int*)buf ;

  switch (magic)
    {
    case OK_MAGIC:
      break ;
    case SWAP_MAGIC:
      bufferDoSwap(buf) ;
      break ;
    default:
      return FALSE ;
    }
  buf2hdr(buf, hh) ;

  return TRUE ;
}

/**************************************************************/
static BOOL headerToBuffer (ACE_HEADER hh, char *buf)
{
  int magic = hh->swapMagic ;

  hdr2buf(hh, buf) ;

  switch (magic)
    {
    case SWAP_MAGIC:
      bufferDoSwap(buf) ;
      break ;
    case OK_MAGIC:
      break ;
    default:
      return FALSE ;
    }
  return TRUE ;
}


/**************************************************************/

int Fcntl(int fd, int cmd, int arg)
{
  int n = fcntl (fd, cmd, arg);
	
  if (n == -1)
    messcrash("fcntl error setting blocking state of socket, %s", messSysErrorText());

  return n ;
}

/**************************************************************/

void Connect(int fd, const struct sockaddr *sa, socklen_t salen)
{
  if (connect(fd, sa, salen) < 0)
    messcrash ("connect error %s", messSysErrorText());
}

/**************************************************************/

int aceSocketSelect(int nfds, fd_set *readfds, fd_set *writefds, fd_set *exceptfds,
	   /* struct timeval */ void  *timeout)
{
  int n;

  restart:
  n = select(nfds, readfds, writefds, exceptfds, timeout) ;
  
  if (n < 0)
    { 
      if (errno == EINTR)
	goto restart;
      else
	messcrash("Cannot select on socket %s",messSysErrorText()) ; 
    }

  return n ;		/* can return 0 on timeout */
}

/**************************************************************/

/* Use of this function is deprecated, should use Shutdown instead.       
void Close(int fd)
{
  if (close(fd) == -1)
    messcrash("Could not close socket connection, %s",messSysErrorText());
}
*/

/* Should use this for sockets...                                            */
/* We ignore ENOTCONN which means that the socket has probably been closed   */
/* by the client.                                                            */
void Shutdown(int fd, int mode)
{
  if (shutdown(fd, mode) < 0)
    {
      if (errno == ENOTCONN)
	  messout("Socket could not be shutdown because it is already disconnected.") ;
      else
	messcrash("Internal logic error, invalid option to shutdown call: %s", messSysErrorText()) ;
    }
  return ;
}

/**************************************************************/
/* Create a blank message header.                                            */
S_MESSAGE sMessageCreate(void)
{
  S_MESSAGE msg= (S_MESSAGE) messalloc(sizeof(struct S_MESSAGESTRUCT)) ;

  /* everything is initialized to zero by messalloc */  
  return msg ;
}

/* Set up a message to receive data from a socket.                           */
void sMessageSetRead(S_MESSAGE msg)
{
  if (msg->isNew)
    {
      msg->doneHeader = FALSE ;
      msg->hBytesRequested =  msg->hBytesPending = ACE_HEADER_BYTES ;
      memset (msg->hBuffer,0, ACE_HEADER_BYTES) ;
      msg->mBytesRequested =  msg->mBytesPending = 0 ;
      messfree (msg->readBuffer) ;
      stackDestroy (msg->writeStack) ; 
      msg->isNew = FALSE ;
    }
  return ;
}


/* Set up a message to be written out to a socket.                           */
/*                                                                           */
void sMessageSetWrite(S_MESSAGE msg)
{
 if (msg->isNew)
    {
      msg->doneHeader = FALSE ;
      msg->hBytesRequested =  msg->hBytesPending = ACE_HEADER_BYTES ;
      memset (msg->hBuffer,0, ACE_HEADER_BYTES) ;

      /* You need to be aware that stackMark gives the length of memory      */
      /* allocated to the current item on the stack, but this can be larger  */
      /* than the item...there is padding allowed to get alignment for       */
      /* stuffing ints etc. on to the stack. This can mess up attempts to    */
      /* find the length of text etc., it means that you have to know that   */
      /* something on the stack is text and do your own strlen, stack only   */
      /* returns its own internal length, not the length you originally      */
      /* supplied...                                                         */
      msg->ah.length =
	msg->mBytesPending = 
	msg->mBytesRequested = stackMark(msg->writeStack) ;


#ifdef ED_G_NEVER_INCLUDE_THIS_CODE
      /* Ok, I did this because I couldn't get the above to work but I think */
      /* my problems were caused by not appreciating the below stuff about   */
      /* pushTest and catText which are used non-obviously in the server     */
      /* code.                                                               */

      /* The above code now seems to work.                                   */

      /*                   YOU HAVE TO MAKE SURE THAT AFTER THE FIRST "PUSH" */
      /* ON TO THE STACK, ONLY CATTEXT IS USED AFTER THIS, NOT SO EASY       */
      /* BECAUSE THE ACECOMMANDEXECUTE DOES A PUSHTEXT FOR SOME COMMANDS     */
      /* BUT (YOU GUESSED IT) NOT OTHERS...AGGGHHHH SO IF YOU DO A PUSHTEXT  */
      /* AS WELL IT ALL BREAKS... see the code in serverace.c where it       */
      /* returns text...basically it does only "catext" once we have done    */
      /* an acecommand execute...aggghhh....actually its worse than this     */
      /* because the aceOut code which is used by acecommandexecute does     */
      /* NOT do a cattext....it does a memcpy basically....                  */
      msg->ah.length =
	msg->mBytesPending = 
	msg->mBytesRequested = strlen(stackText(msg->writeStack, 0)) + 1 ;
#endif /* ED_G_NEVER_INCLUDE_THIS_CODE */


      messfree (msg->readBuffer) ;  
      msg->isNew = FALSE ;
    }
 return ;
}



/* Set the magic, the client should set this on its first write, the server  */
/* should then copy the clients magic and use it for subsequent writes to    */
/* the client.                                                               */
void sMessageSetMagic(S_MESSAGE msg)
{
  msg->ah.swapMagic = OK_MAGIC ;
}

int  sMessageGetMagic(S_MESSAGE msg)
{
  return msg->ah.swapMagic ;
}

/* We will need this, maybe we already do...not sure of the interface to     */
/* ace...it should extract the length and message from ace...                */
/* N.B. no safe guards at the moment, could be used before message is all    */
/* read...                                                                   */
int sMessageGetMessage(S_MESSAGE msg, char **msgtext)
{
  *(msgtext) = msg->readBuffer ;
  return (msg->mBytesRequested) ;
}


/* Clean up a message of allocated resources                                 */
void sMessageDestroy(S_MESSAGE msg)
{ 
  if (msg)
    {
      messfree (msg->readBuffer) ;  
      stackDestroy (msg->writeStack) ;  
      messfree(msg) ;
    }
  messfree(msg) ;
}

/**************************************************************/

/* n.b. It is possible to peek socket buffers to see what is in there, but   */
/* this method gives extra complications, what if the client crashes between */
/* peek and read ?? what if the whole message is to big to be transmitted in */
/* one go, it's just not worth the hassle, hence we do straight reads.       */
S_MSGState sMessageSocketRead(int fd, S_MESSAGE req)
{
  S_MSGState state = SMSG_ERROR ;
  int sockFlags ;

  /* We can't just do a recv with MSG_DONTWAIT and MSG_PEEK because          */
  /* the former is not supported on many systems yet. So instead we make the */
  /* whole socket non-blocking, do the recv and then reset the socket.       */
  sockFlags = Fcntl(fd, F_GETFL, 0) ;
  Fcntl(fd, F_SETFL, sockFlags | O_NONBLOCK) ;

  /* Try to read the header, probably do this in one go.                     */
  if (req->doneHeader == FALSE)
    {
      /* First time for this message so allocate a buffer.                   */
      if (!req->hBytesPending)
	{
	  req->hBytesPending = req->hBytesRequested = ACE_HEADER_BYTES ;
	}

      state = readSocket(fd, req->hBuffer + (req->hBytesRequested - req->hBytesPending),
			 &(req->hBytesPending)) ;

      if (state == SMSG_DONE)
	{
	  if (!bufferToHeader(req->hBuffer, &(req->ah)))
	    {
	      state = SMSG_ERROR ;
	      messout("Error in socket header format from client.") ;
	    }
	  else
	    req->doneHeader = TRUE ;
	}
    }

  /* We have a header so try to read buffer (this may take several calls to  */
  /* this function).                                                         */
  if (req->doneHeader == TRUE)
    {
      /* First time for this message so allocate a buffer.                   */
      if (!req->readBuffer)
	{
	  req->readBuffer = messalloc (req->ah.length) ;
	  req->mBytesRequested = req->mBytesPending = req->ah.length ;
	}

      state = readSocket(fd, req->readBuffer + req->mBytesRequested - req->mBytesPending,
			 &(req->mBytesPending)) ;
    }

  /* Reset the socket to blocking.                                           */
  Fcntl(fd, F_SETFL, sockFlags) ;

  return state ;
}


/* Common read routine for handling reading from a socket and any resulting  */
/* errors. Note that it updates the the number of bytes left to read given   */
/* by bytesPending.                                                          */
static S_MSGState readSocket(int fd, char *readbuf, int *bytesPending)
{
  S_MSGState state = SMSG_ERROR ;
  int bytes ;

  bytes = read(fd, readbuf, *bytesPending) ;
  if (bytes < 0)
    {
      if (errno == EWOULDBLOCK || errno == EINTR)
	state = SMSG_WAIT ;
      else
	{
	  messdump("read error trying to read client data on socket,"
		   "connection terminated:  %s", messSysErrorText());
	  state = SMSG_ERROR ;
	}
    }
  else if (bytes == 0)
    {
      messdump("zero bytes returned on read, client has prematurely closed connection.") ;
      state = SMSG_ERROR ;
    }
  else
    {
      *bytesPending -= bytes ;
      
      if (*bytesPending < 0)
	{
	  messdump("read more bytes than number requested ") ;
	  state = SMSG_ERROR ;
	}
      else  if (*bytesPending == 0)
	state = SMSG_DONE ;
      else
	state = SMSG_WAIT ;
    }

  return state ;
}



/**************************************************************/

/* Jean I added code here to store any current handler for SIGPIPE, we don't */
/* know if the rest of the application will have installed one. We then      */
/* turn off signal handling for SIGPIPE so that we get EPIPE for a write to  */
/* a broken connection which we handle here. Then we turn back on any        */
/* existing signal handler.                                                  */
S_MSGState sMessageSocketWrite(int fd, S_MESSAGE req)
{
  S_MSGState state = SMSG_WAIT ;
  int sockFlags ;
  struct sigaction oursigpipe, oldsigpipe ;

  /* writes can deliver a SIGPIPE if the socket has been disconnected, by    */
  /* ignoring we will just receive errno = EPIPE.                            */
  oursigpipe.sa_handler = SIG_IGN ;
  sigemptyset(&oursigpipe.sa_mask) ;
  oursigpipe.sa_flags = 0 ;
  if (sigaction(SIGPIPE, &oursigpipe, &oldsigpipe) < 0)
    messcrash("Cannot set SIG_IGN for SIGPIPE for socket write operations") ;


  /* We can't just do a recv with MSG_DONTWAIT and MSG_PEEK because          */
  /* the former is not supported on many systems yet. So instead we make the */
  /* whole socket non-blocking, do the recv and then reset the socket.       */
  sockFlags = Fcntl(fd, F_GETFL, 0) ;
  Fcntl(fd, F_SETFL, sockFlags | O_NONBLOCK) ;


  /* OK, Let's try to write the header                                       */
  if (req->doneHeader == FALSE)
    {
      /* Do any byte swapping.                                               */
      headerToBuffer(&(req->ah), req->hBuffer) ;

      /* Now do the writes...                                                */
      state = writeSocket(fd, req->hBuffer + req->hBytesRequested - req->hBytesPending,
			  &(req->hBytesPending)) ;
      if (state == SMSG_DONE)
	req->doneHeader = TRUE ;
    }

  /* OK, header written, so try to write buffer.                             */
  if (req->doneHeader == TRUE)
    {
      if (req->mBytesPending <= 0)
	{
	  messdump("Logic error, caller is trying to write %d bytes to socket.", req->mBytesPending) ;
	  state = SMSG_ERROR ;
	}
      else
	{
	  state = writeSocket(fd,
			      stackText(req->writeStack, req->mBytesRequested - req->mBytesPending),
			      &(req->mBytesPending)) ;
	}
    }

  /* Reset the socket to blocking.                                           */
  Fcntl(fd, F_SETFL, sockFlags) ;

  /* Reset the old signal handler.                                           */
  if (sigaction(SIGPIPE, &oldsigpipe, NULL) < 0)
    messcrash("Cannot reset previous handler for signal SIGPIPE for socket write operations") ;

  return state ;
}


/* Common write routine for handling writing to a socket and any resulting   */
/* errors. Note that it updates the the number of bytes left to write given  */
/* by bytesPending.                                                          */
static S_MSGState writeSocket(int fd, char *writebuf, int *bytesPending)
{
  S_MSGState state = SMSG_ERROR ;
  int bytes ;

  bytes = write(fd, writebuf, *bytesPending) ;
  if (bytes < 0)
    {
      if (errno == EINTR || errno == EAGAIN)
	state = SMSG_WAIT ;
      else
	{
	  messdump("Cannot write message buffer to socket, %s", messSysErrorText()) ;
	  state = SMSG_ERROR ;
	}
    }
  else
    {
      *bytesPending -= bytes ;

      if (*bytesPending < 0)
	{
	  messdump("wrote more bytes than number requested ") ;
	  state = SMSG_ERROR ;
	}
      else  if (*bytesPending == 0)
	state = SMSG_DONE ;
      else 
	state = SMSG_WAIT ;
    }

  return state ;
}


/**************************************************************/
/**************************************************************/
