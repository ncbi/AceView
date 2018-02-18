/*  File: acesocket.c
 *  Author: Ed Griffith (edgrif@sanger.ac.uk) & Jean Thierry-Mieg (mieg@crbm.cnrs-mop.fr)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1992
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: Implements the TCP Socket communication layer that
 *              sits between acedb and the network. Handles non-blocking
 *              connections between client(s) and server, handles
 *              timeouts, socket errors. Keeps clients as a circularly
 *              linked list with clients being maintained as finite
 *              state machines (reading, writing, etc.).
 *
 * Exported functions:
 * HISTORY:
 * Last edited: Jan 22 16:54 2001 (edgrif)
 * * Mar 16 11:30 2000 (edgrif): Fix bug in timeout code so timeout
 *              is _since_ last socket communication.
 * Created: Wed July 8 1999 (mieg)
 * CVS info:   $Id: acesocket.c,v 1.1.1.1 2002/07/19 20:23:35 sienkiew Exp $
 *-------------------------------------------------------------------
 */

#include <wh/regular.h>
#include <wsocket/acesocket_.h>


/************************************************************/

static int  aceSocketCreate (int port) ;
static void setCleanEnv (void) ;
static int clientAccept (int fd, struct sockaddr *sa, socklen_t *salenptr) ;

static S_CLIENT sClientCreate (S_CLIENT* currentClientp) ;
static void sClientDestroy (S_CLIENT* currentClientp, S_CLIENT client) ;
static void closeAllClients (S_CLIENT* currentClientp, S_CLIENT shut_client,
			     fd_set allreads, fd_set allwrites, Stack shutDownMessage,
			     AceSocketServerRequestRoutine serverRequestFunc, void *serverData) ;
static BOOL timeOut(S_CLIENT client, int clientTimeOut, int *curr_timeout) ;
static void resetTime(S_CLIENT client) ;
static void checkServerReturn(AceSocketRequestTarget currTarget, AceSocketRequestType currRequest) ;



/************************************************************/

void aceSocketListen (int port, int serverTimeOut_in, int clientTimeOut_in,
		      AceSocketServerRequestRoutine serverRequestFunc, void *serverData) 
{
  int		maxfd, connfd ;
  int		nready ;
  fd_set	rset, wset, allreads, allwrites ;
  socklen_t	clilen;
  struct sockaddr_in cliaddr ;
  S_CLIENT currentClient = NULL , client = NULL ;
  int listenfd ;
  int time_out, serverTimeOut, clientTimeOut, shortest_timeout ;
  struct timeval tv, *tv_ptr ;
  int sockFlags ;
  AceSocketRequestTarget currTarget = ACESOCK_NONE ;
  AceSocketRequestType currRequest = ACESOCK_NULL ;
  BOOL shuttingDown = FALSE ;
  AceSocketRequestType serverRequest = ACESOCK_NULL ;
  Stack shutMessage = NULL ;


  /* Set timeouts for clients and the server itself.                         */
  /* If either is set to the special value zero, this implies an infinite    */
  /* timeout which should only be used for testing.                          */
  serverTimeOut = serverTimeOut_in >= ACESERV_MIN_TIMEOUT ? serverTimeOut_in : 600 ;
  clientTimeOut = clientTimeOut_in >= ACESERV_MIN_TIMEOUT ? clientTimeOut_in : 600 ;
  shortest_timeout = clientTimeOut ;

  /* Create a listening socket (non-blocking) on which to receive new clients*/
  listenfd = aceSocketCreate(port) ;
  sockFlags = Fcntl(listenfd, F_GETFL, 0) ;
  Fcntl(listenfd, F_SETFL, sockFlags | O_NONBLOCK) ;

  /* Set up the reading/writing socket sets for monitoring client connections*/
  maxfd = listenfd ;
  FD_ZERO(&allreads);
  FD_ZERO(&allwrites);
  FD_SET(listenfd, &allreads);

  /* Loop forever waiting for clients to connect and send requests, we exit  */
  /* when 'serverTimeOut' time has passed with no requests sent.             */
  while (TRUE)
    {

      /* Shutdown server/sockets in requested way.                           */
    shutdownCode:
      if (shuttingDown)
	{
	  /* Forced shutdown so kill all clients.                            */
	  if (serverRequest == ACESOCK_SHUTFORCE
	      || serverRequest == ACESOCK_TIMEDOUT)
	    {
	      closeAllClients(&currentClient, client, allreads, allwrites, shutMessage,
			      serverRequestFunc, serverData) ;
	    }

	  /* Once all clients are gone then do the shutdown.                 */
	  if (serverRequest == ACESOCK_SHUTFORCE
	      || serverRequest == ACESOCK_TIMEDOUT
	      || (!currentClient && serverRequest == ACESOCK_SHUTQUIESCE))
	    {	
	      /* Close listening socket and call server to do final exit.    */
	      /* n.b. only close the listening socket if we are _not_ inetd  */
	      /* started. Probably we could just always do a close because   */
	      /* it only closes the socket if the use count is 0.            */
	      if (!isInetdDaemon(NULL))
		Close(listenfd) ;
	      currTarget = ACESOCK_SERVER, currRequest = ACESOCK_EXIT ;
	      serverRequestFunc(serverData, client, 0, &currTarget, &currRequest, NULL, NULL) ;
	    }
	}


      /* Currently we only listen for reading or writing sockets to become   */
      /* ready....                                                           */
      /* ONE DAY WE MAY WISH TO LISTEN FOR EXCEPTIONS AS WELL...             */
      rset = allreads ;
      wset = allwrites ;

      /* Set the time out for the select call:                               */
      /*      no clients => user server time out                             */
      /*          client => use shortest time to time out a client           */
      /* If server or client time outs are zero, then servers/clients never  */
      /* time out (good for testing).                                        */
      if (currentClient == NULL)
	time_out = serverTimeOut ;
      else
	{
	  if (clientTimeOut == ACESERV_INFINITE_TIMEOUT)
	    time_out = clientTimeOut ;
	  else
	    time_out = shortest_timeout ;
	}

      /* Setting tv_ptr to NULL means "never time out", setting time_out to  */
      /* zero causes the select to return immediately.                       */
      if ((currentClient == NULL && serverTimeOut == 0)
	  || (currentClient != NULL && clientTimeOut == 0))
	tv_ptr = NULL ;
      else
	{
	  tv.tv_sec = time_out ;			    /* LINUX systems modify tv param. */
	  tv.tv_usec = 0 ;
	  tv_ptr = &tv ;
	}

      /* Block waiting for a socket to become ready, or we time out.         */
      nready = aceSocketSelect(maxfd+1, &rset, &wset, NULL, tv_ptr) ;
      if (nready == 0)
	{
	  /* If there are no clients then the server has timed out and we    */
	  /* shutdown, otherwise we fall through to time out clients.        */
	  if (currentClient == NULL)
	    {
	      /* Call server to alert to timeout                                 */
	      currTarget = ACESOCK_SERVER, currRequest = ACESOCK_TIMEDOUT ;
	      shutMessage = serverRequestFunc(serverData, client, 0, &currTarget, &currRequest,
					      NULL, NULL) ;
	      /* We could allow server to cancel timeout by returning NULL       */
	      /* request...currently timeout => shutdown.                        */
	      shuttingDown = TRUE ;
	      serverRequest = ACESOCK_TIMEDOUT ;
	      goto shutdownCode ;
	    }
	}



      /* Is there a new client connected to the listening descriptor ? Note  */
      /* that accept CAN fail even if select says there is a new client,     */
      /* maybe because client aborted, in this case we ignore the connect.   */
      /* We can only have up to (FD_SETSIZE - 1) sockets, if the new socket  */
      /* would exceed this, we accept the connection, log the problem and    */
      /* close the connection in an orderly way. Currently the unlucky user  */
      /* will get no message from us but their client should tell them.      */
      if (FD_ISSET(listenfd, &rset))
	{
	  clilen = sizeof(cliaddr);

	  /* connfd may get interrupted and so may not return a valid fd.    */
	  connfd = clientAccept(listenfd, (struct sockaddr *)&cliaddr, &clilen) ;
	  if (connfd >= 0)
	    {
	      /* uugghh, wrong place, we don't even need to accept the       */
	      /* above connection, I guess I put it here because I was going */
	      /* to send an error message.                                   */
	      if (shuttingDown)
		{
		  /* No message, we don't know what state server is in so no */
		  /* point in calling it.                                    */
		  messdump("Client attempting to connect but server is closing down "
			   "so connection refused.") ;
		  Shutdown(connfd, SHUT_RDWR) ;
		}
	      else if(connfd >= FD_SETSIZE)
		{
		  /* Oops, too many clients (FD_SETSIZE is usually 256). */
		  /* Currently, unlucky connection is just terminated.   */
		  messdump("Refused connection from client because it would have "
			   "exceeded the maximum number of files/sockets"
			   "allowed by the 'select' call (FD_SETSIZE), which is  %d",
			   FD_SETSIZE) ; 
		  Shutdown(connfd, SHUT_RDWR) ;
		}

#ifdef ED_G_NEVER_INCLUDE_THIS_CODE
	      else if (!setClientInetAddressInfo(&cliaddr, &dotted_name, &host_name))
		{
		  messdump("Refused connection from client because could not "
			   "resolve its hostname.") ;
		  Shutdown(connfd, SHUT_RDWR) ;
		}
#endif /* ED_G_NEVER_INCLUDE_THIS_CODE */

	      else
		{
		  /* Accept the connection and make a new client.        */
		  client = sClientCreate (&currentClient) ;
		  client->state = CLCON_CONN ;
		  client->fd = connfd ;
		  client->socket = cliaddr ;		    /* struct copy */
		  
		  /* Set up select call params to read from new clients sock.*/
		  if (connfd > maxfd)
		    maxfd = connfd;			    /* needed by select call */
		  FD_SET(client->fd, &allreads);	    /* new clients are always in read set. */
		}
 	    }
	}
      
      /* No clients and client connect failed so just select again.          */
      if (!currentClient)
	continue ;

      /* OK, process some clients...                                         */
      /* We always try to start our round with the next client in the        */
      /* circular list to attempt some fairness in dispatching.              */
      client = currentClient = currentClient->next ;

      /* We need to find out the client with the shortest time until it times*/
      /* out so we reset the curr client time out and check it against each  */
      /* client.                                                             */
      shortest_timeout = clientTimeOut ;

      /* Cycle through all the clients once, clients may move through several*/
      /* states in one "turn", e.g. read->error->close, both to ensure some  */
      /* events are processed in a timely manner _and_ to avoid getting stuck*/
      /* in a state not reliant on a read or write to a socket, e.g. request.*/
      do
	{

	lao:

	  /* If client is in one of the states associated with a socket then */
	  /* check for timeout (we don't check for the others as they happen */
	  /* synchronously). If client has not timed out AND its socket is   */
	  /* ready (for reading or writing) then reset its timer.            */
	  if (client->state == CLCON_CONN || client->state == CLCON_READWAIT
	      || client->state == CLCON_WRITEWAIT)
	    {
	      if (timeOut(client, clientTimeOut, &shortest_timeout))
		{
		  client->state = CLCON_CLOSE ;
		  client->close = CL_TIMEDOUT ; 
		}
	      else if (FD_ISSET(client->fd, &rset) || FD_ISSET(client->fd, &wset))
		{
		  resetTime(client) ;
		}
	    }


	  switch (client->state)
	    {
	    case CLCON_CONN:
	      if (FD_ISSET(client->fd, &rset))		    /* New data from client. */
		{
		  client->state = CLCON_READWAIT ;
		  client->message->isNew = TRUE ;
		  goto lao ;
		}
	      break ;
	    case CLCON_READWAIT:
	      if (FD_ISSET(client->fd, &rset))		    /* New/more data from client. */
		{
		  S_MSGState msgState ;

		  if (client->message->isNew)
		    sMessageSetRead(client->message) ;

		  msgState = sMessageSocketRead(client->fd, client->message) ;
		  if (msgState == SMSG_WAIT)
		    {
		      client->state = CLCON_READWAIT ;
		    }
		  else if (msgState == SMSG_DONE)
		    {
		      client->state = CLCON_REQUEST ;
		      goto lao ;			    /* Process request immediately. */
		    }
		  else if (msgState == SMSG_ERROR)
		    {
		      client->state = CLCON_CLOSE ;
		      client->close = CL_ERROR ;
		      goto lao ;			    /* Process error immediately. */
		    }
		}
	      break ;
	    case CLCON_REQUEST:
	      /* Send a request to the server code and get its reply to send */
	      /* back to the client.                                         */
	      stackDestroy (client->message->writeStack) ;
	      currTarget = ACESOCK_CLIENT, currRequest = ACESOCK_REQ ;

	      client->message->writeStack = serverRequestFunc(serverData,
							      client, client->fd,
							      &currTarget, &currRequest,
							      &(client->message->ah), 
							      client->message->readBuffer) ;
	      client->mustCallAce = TRUE ;		    /* Remember that we have called ace. */

	      /* Make sure that we get something sensible back, exit if not. */
	      checkServerReturn(currTarget, currRequest) ;

	      if (currTarget == ACESOCK_SERVER)
		{
		  /* Our interface shutdowns are generated at the request of */
		  /* this client.                                            */
		  if (currRequest == ACESOCK_SHUTFORCE)
		    {
		      /* Immediate kill of the server & all clients.         */
		      shuttingDown = TRUE ;
		      serverRequest = currRequest ;
		      if (client->message->writeStack)
			shutMessage = stackCopy(client->message->writeStack, 0) ;
		      goto shutdownCode ;		    /* N.B. the goto !! */
		    }
		  else if (currRequest == ACESOCK_SHUTQUIESCE)
		    {
		      /* Kill this client as it requested the shutdown and   */
		      /* then wait for all clients to quit.                  */
		      shuttingDown = TRUE ;
		      serverRequest = currRequest ;
		      client->mustCallAce = FALSE ;
		      client->state = CLCON_CLOSE ;
		      client->close = CL_NORMAL ;	
		      if (client->message->writeStack)
			  client->message->isNew = TRUE ;
		    }
		}
	      else
		{
		  if (currRequest == ACESOCK_REQ)
		    {
		      if (client->message->writeStack)
			{
			  /* a new message should be written */
			  client->state = CLCON_WRITEWAIT ;
			  client->message->isNew = TRUE ;
			  FD_CLR(client->fd, &allreads) ;
			  FD_SET(client->fd, &allwrites) ;
			}
		      else
			client->state = CLCON_CONN ;	    /* request with no reply from acedb. */
		    }
		  else if (currRequest == ACESOCK_KILL)
		    {
		      /* Server wishes to kill client so no need to call it  */
		      /* again.                                              */
		      client->mustCallAce = FALSE ;
		      client->state = CLCON_CLOSE ;
		      client->close = CL_NORMAL ;	
		      if (client->message->writeStack)
			  client->message->isNew = TRUE ;     /* a new message should be written */
		    }
		}
	      goto lao ;
	      break ;
	    case  CLCON_WRITEWAIT:
	      if (FD_ISSET(client->fd, &wset))
		{
		  S_MSGState msgState ;

		  if (client->message->isNew)
		    sMessageSetWrite(client->message) ;
		  msgState = sMessageSocketWrite(client->fd, client->message) ;
		  if (msgState == SMSG_WAIT)
		    {
		      client->state = CLCON_WRITEWAIT ;
		    }
		  else if (msgState == SMSG_DONE)
		    {
		      client->state = CLCON_CONN ;
		      FD_CLR(client->fd, &allwrites) ;	    /* Swop from write to read set. */
		      FD_SET(client->fd, &allreads) ;
		    }
		  else if (msgState == SMSG_ERROR)
		    {
		      client->state = CLCON_CLOSE ;
		      client->close = CL_ERROR ;  
		      goto lao ;
		    }
		}
	      break ;
	    case CLCON_CLOSE:
	      /* We do all closes immediately after the above operations other-  */
	      /* wise we may block on the select and not close the client for    */
	      /* a long time.                                                    */

	      /* Note that mustCallAce is TRUE only for timedout and error   */
	      /* closes. normal close => we have already called ace.         */
	      if (client->mustCallAce)
		{
		  /* If client is being closed by the socket layer, call ace */
		  /* server so it can free up its client resources.          */
		  currTarget = ACESOCK_CLIENT ;
		  if (client->close == CL_TIMEDOUT)
		      currRequest = ACESOCK_TIMEDOUT ;
		  else if (client->close == CL_ERROR)
		    currRequest = ACESOCK_KILL ;
		  client->message->writeStack =
		    serverRequestFunc(serverData,
				      client, client->fd,
				      &currTarget, &currRequest, &(client->message->ah), NULL) ;
		  if (client->message->writeStack && client->close == CL_TIMEDOUT)
		    client->message->isNew = TRUE ;
		}
	      
	      /* If there's a termination message send it, but only to a     */
	      /* client with no socket errors.                               */
	      if (client->close != CL_ERROR && client->message->isNew)
		{
		  /* Just try to write once, don't even bother to check, perhaps */
		  /* a bit hard on a normal quitter.                             */
		  /* Could insert my code to try writing a number of times */
		  /* and then just stop...                               */
		  sMessageSetWrite(client->message) ;
		  sMessageSocketWrite(client->fd, client->message) ;
		}

	      FD_CLR(client->fd, &allreads) ;
	      FD_CLR(client->fd, &allwrites) ;

	      sClientDestroy(&currentClient, client) ;	    /* Does the close as well. */
	      client = NULL ;

	      if (currentClient == NULL)		    /* Last client in list ? */
		client = NULL ;
	      else
		client = currentClient ;

	      break ;

	    default:
	      /* Bad news, somehow assigned a client state we don't recognise*/
	      messcrash ("Internal logic error, client connection FSM is in an unknown state.") ;
	      break ;
	    }

	  if (client != NULL)				    /* Any clients left ?? */
	    client = client->next ;

	} while (client != currentClient) ;		    /* Loop through clients once */

    }  /* loop back on Select */


  messcrash("We should never reach the end of this function (aceSocketListen), "
	    "logic error in code.") ;
  return ;
}


/************************************************************/

/* Only really called when we are run by inetd                               */
/* We try to clear any pending data on the listening socket, if we don't do  */
/* this then inetd will immediately start up the server again, which will    */
/* crash again only to be restarted again, etc.etc.                          */
/*                                                                           */
/* Note that the listening socket MUST NOT be closed here because this may   */
/* cause inetd to crash.                                                     */
/*                                                                           */
void aceSocketInetdCleanup(void)
{
  if (isInetdDaemon(NULL))
    {
      int nready, maxfd, connfd, listenfd ;
      fd_set rset ;
      socklen_t clilen;
      struct sockaddr_in cliaddr ;
      struct timeval tv ;

      /* BECAUSE its inetd, we _know_ the listening socket file desc. == 0   */
      listenfd = 0 ;

      maxfd = listenfd ;
      FD_ZERO(&rset) ;
      FD_SET(listenfd, &rset) ;
      tv.tv_sec = 0 ;					    /* Don't wait at all. */
      tv.tv_usec = 0 ;
      nready = aceSocketSelect(maxfd+1, &rset, NULL, NULL, &tv) ;
      if (nready == 0)
	{
	  messdump("Could not read listening socket even though there should be data"
		   "on it, suggests severe problem with inetd/listening socket") ;
	}
      else
	{
	  if (FD_ISSET(listenfd, &rset))		    /* redundant, must be true really. */
	    {
	      clilen = sizeof(cliaddr);

	      /* connfd may get interrupted and so may not return a valid fd.    */
	      connfd = clientAccept(listenfd, (struct sockaddr *)&cliaddr, &clilen) ;
	      if (connfd >= 0)
		{
		  Shutdown(connfd, SHUT_RDWR) ;
		}
	    }
	}

    }

  return ;
}



/************************************************************/
/***************  Internal routines  ************************/

/* Clients are organized as a circular doubly linked list
 * with currentClient pointing to zero (no clients in loop)
 * or to the currently active client in the list.
 */

static S_CLIENT sClientCreate (S_CLIENT* currentClientp)
{
  S_CLIENT client = (S_CLIENT) messalloc (sizeof(ClientConnectStruct)) ;

  if (!*currentClientp) /* Create a ring of a single client ! */
    {
      *currentClientp = client ;
      client->previous = client ;
      client->next = client ;
    }
  else   /* install the new client just under the current one */
    {
      client->next = (*currentClientp)->next ;
      client->next->previous = client ;
      (*currentClientp)->next = client ;
      client->previous = *currentClientp ;
    }

  /* Set client state to 'uninitialised'. */
  client->state = CLCON_EMPTY ;
  client->close = CL_NONE ;
  client->fd = -1 ;
  client->time = timeNow() ;
  client->message = sMessageCreate () ;
  client->mustCallAce = FALSE ;

  return client ;
}


/************************************************************/

static void sClientDestroy (S_CLIENT* currentClientp, S_CLIENT client) 
{
  client->previous->next = client->next ;
  client->next->previous = client->previous ;

  if (*currentClientp == client)
    *currentClientp = (client->next != client) ? client->next : NULL ;

  Shutdown(client->fd, SHUT_RDWR) ;
  sMessageDestroy (client->message) ;
  messfree (client) ;
}


/************************************************************/

/* This is basically a pretty crap interface, the first client is the one    */
/* that sent the shutdown request so we don't call ace to close that one.    */
/* Hardly transparent, there is much common code in the close client stuff   */
/* that could be unified...                                                  */
static void closeAllClients (S_CLIENT* currentClientp, S_CLIENT shut_client,
			     fd_set allreads, fd_set allwrites, Stack shutMessage,
			     AceSocketServerRequestRoutine serverRequestFunc, void *serverData)
{
  S_CLIENT client;
  AceSocketRequestTarget currTarget =  ACESOCK_CLIENT ;
  AceSocketRequestType currRequest = ACESOCK_KILL ;

  while ((client = *currentClientp))
    {
      /* Call server request func to close all clients....but not the   */
      /* one that sent the shutdown...                                       */
      if (client != shut_client)
	{
	  client->message->writeStack =
	    serverRequestFunc(serverData, client, client->fd,
			      &currTarget, &currRequest, &(client->message->ah), NULL) ;
	}
      
      /* Note that we ignore the message returned by the call to the server  */
      /* and just output our own.                                            */
      if (shutMessage)
	{
	  stackDestroy (client->message->writeStack) ;
	  client->message->writeStack = stackCopy(shutMessage, 0) ;
	  client->message->isNew = TRUE ;
	  sMessageSetWrite(client->message) ;
	  sMessageSocketWrite(client->fd, client->message) ;
	}
      FD_CLR(client->fd, &allreads) ;
      FD_CLR(client->fd, &allwrites) ;
      sClientDestroy (currentClientp,*currentClientp) ;
    }
}


/************************************************************/

/* UUUGGGHHH, I don't really like this interface but it will have to do for  */
/* now....don't go altering it without looking at the cmd line arguments and */
/* select call.                                                              */
/*                                                                           */
/* routine returns:                                                          */
/*                    FALSE if client has not timed out or if special "no    */
/*                          time out" value set for clients (i.e. zero).     */
/*                                                                           */
/*                     TRUE if client has timed out                          */
/*                                                                           */
/* shortest_timeout is up updated each time if the time left for this client */
/* is shorter than the shortest so far, in particular, if a client has timed */
/* out then shortest_timeout is set to zero to ensure the client is killed   */
/* immediately.                                                              */
/*                                                                           */
static BOOL timeOut(S_CLIENT client, int clientTimeOut, int *shortest_timeout)
{
  BOOL result = TRUE ;

  if (clientTimeOut == ACESERV_INFINITE_TIMEOUT)
    result = FALSE ;
  else
    {
      int curr_time, elapsed_time, time_left ;

      curr_time = timeNow() ;
      elapsed_time = curr_time - client->time ;

      /* Check to see if client has timed out, if not then reset             */
      if (elapsed_time >= clientTimeOut)
	{
	  result = TRUE ;
	  *shortest_timeout = ACESERV_MIN_TIMEOUT ;	    /* Ensures client killed immediately. */
	}
      else
	{
	  time_left = clientTimeOut - elapsed_time ;
	  if (time_left < *shortest_timeout)
	    *shortest_timeout = time_left ;

	  result = FALSE ;
	}
    }

  return result ;
}

/* Reset clients timer to "now", this is only done when there is actually    */
/* about to be some communication between the server/client.                 */
static void resetTime(S_CLIENT client)
{

  client->time = timeNow() ;

  return ;
}



/**************************************************************/
/* disassociate for terminals, */

/* THIS DEFINITELY NEEDS SOME CHECKING...                                    */
static void setCleanEnv(void)
{
  /* in dec alpha documentation i found the following recommendations */
  setsid () ;   /* reset neutral environment */
  chdir ("/") ;
  umask (0) ;
  
  return ;
}


/************************************************************/
/* OK, create a socket that uses streams and bind to it.      
 * return -1: failed
 *         0: inetd case, inetd duplicates socket on stdin
 *       n>0: foreground socket
 *
 * NOTE, for inetd we ignore 'port', the port is set by in /etc/inetd.conf
 */
static int aceSocketCreate (int port)
{
  int  aceSocket = -1 ;
  int listen_port ;


  /* NOTE that ports and sockets are not the same numbers, the socket is the */
  /* file descriptor number for inetd this is 0, for non-inetd it could be   */
  /* any value.                                                              */
  if (isInetdDaemon(&listen_port))			    /* inetd case */
    {

      /* Really the fact that this is 0 should be encapsulated somewhere,    */
      /* the inetd routine knows this but its not the correct place.         */
      aceSocket = 0 ;


      setCleanEnv () ;

      messout("Server listening socket inherited from inetd, port %d created", listen_port) ;
    }
  else
    { 
      int setSockOpt = 1 ;
      struct sockaddr_in sname ;
      
      aceSocket = socket (AF_INET, SOCK_STREAM, 0) ;
      if (aceSocket < 0)
	messcrash ("Cannot open server listening socket") ;

      /* Set SO_REUSEADDR option otherwise we get an "address already in use" error  */
      /* on the bind because there is some connection from a previous server incarnation */
      /* that has not died yet. (inetd also does this.)                      */
      if (setsockopt(aceSocket, SOL_SOCKET, SO_REUSEADDR,
		     &setSockOpt, sizeof(setSockOpt)) < 0)
	messcrash("Cannot set SO_REUSEADDR option for listening socket: %s", messSysErrorText()) ;

      /* create name with wildchar */
      bzero((char *)&sname, sizeof(sname));
      sname.sin_family = AF_INET ;
      sname.sin_addr.s_addr = htonl(INADDR_ANY) ;
      sname.sin_port = htons(port) ;
      
      if (bind(aceSocket, (struct sockaddr *) &sname, sizeof sname) < 0)
	messcrash ("Cannot bind to server listening socket %s", messSysErrorText()) ;

      if (listen(aceSocket, LISTENQ) < 0)
	messcrash ("Cannot listen to server listening socket %s",messSysErrorText()) ;

      messout("Server listening socket %d created", aceSocket) ;
    }

  return aceSocket ;
}


/**************************************************************/
/* Our version of accept where we catch the additional EWOULDBLOCK error     */
/* because we have set our listening socket non-blocking.                    */
/*                                                                           */
/* When Posix1.g appears then this function will need to be changed because  */
/* ECONNABORTED will mean the client has abborted and EPROTO will mean there */
/* is a probably unrecoverable stream error so we should terminate the server*/
/* But for now this is a real mish/mash on different systems so we will just */
/* do our best.                                                              */
static int clientAccept(int fd, struct sockaddr *sa, socklen_t *salenptr)
{
  int socket = -1 ;

  if ( (socket = accept(fd, sa, salenptr)) < 0)
    {
      if (errno != EINTR && errno != EWOULDBLOCK
#ifdef	EPROTO
	  && errno != EPROTO && errno != ECONNABORTED
#else
	  && errno != ECONNABORTED
#endif
	  )
	messcrash("Unrecoverable error in accept() for new client, %s", messSysErrorText());
    }

  return socket ;

}


/* Check validity of return request from acedb code. Add any new checks to   */
/* here rather then cluttering up main socket code.                          */
static void checkServerReturn(AceSocketRequestTarget target, AceSocketRequestType request)
{

  if (target != ACESOCK_SERVER && target != ACESOCK_CLIENT)
    messcrash("Internal logic error, socket code received bad target from acedb server code:\n"
	      "currTarget currTarget has invalid value: %d \n"
	      "(should be %d or %d).",
	      target, ACESOCK_SERVER, ACESOCK_CLIENT) ; 

  if (target == ACESOCK_SERVER
      && (request != ACESOCK_SHUTFORCE && request != ACESOCK_SHUTQUIESCE))
    messcrash("Internal logic error, socket code received bad request from acedb server code:\n"
	      "currTarget == ACESOCK_SERVER, but currRequest has invalid value: %d \n"
	      "(should be %d or %d).",
	      request, ACESOCK_SHUTFORCE, ACESOCK_SHUTQUIESCE) ; 

  if (target == ACESOCK_CLIENT
      && (request != ACESOCK_REQ && request != ACESOCK_KILL))
    messcrash("Internal logic error, socket code received bad request from acedb server code:\n"
	      "currTarget == ACESOCK_CLIENT, but currRequest has invalid value: %d \n"
	      "(should be %d or %d).",
	      request, ACESOCK_REQ, ACESOCK_KILL) ; 


  return ;
}



/**************************************************************/
/**************************************************************/
