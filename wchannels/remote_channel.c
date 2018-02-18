/*  File:  remote_channel.c
 *  Author:  D and J Thierry-Mieg (mieg@ncbi.nlm.nih.gov)
 *  Copyright (C) D and J Thierry-Mieg and Mark Sienkiewicz 2015
 * -------------------------------------------------------------------
 * AceView is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 * or see the on-line version at http://www.gnu.org/copyleft/gpl.txt
 * -------------------------------------------------------------------
 * This file is part of the AceView project developped by
 *	 D and J Thierry-Mieg (mieg@ncbi.nlm.nih.gov)
 *
 *  This should be linked with the acedb libraries
 *  available from http://www.aceview.org
 *
 *  Please direct all question about this code to
 *  mieg@ncbi.nlm.nih.gov
 *
 * File remote_channel.c, 
 *    Depends on libtsfree.a (channels), libwego.a (multi-trhreading), libtcpfn.a (tcp communcations), 
 * Purpose of this code
 *  Generalizes the channel class, inspired by the 'go' language channels
 *  to allow robust synchronisation of a multi-tasking programs accross a computer network
 *  with the provision that channels should be the only mean of communication between threads 
 *  and between tasks. No other type of synchronization calls should be needed. 
 *  
 ***  The public functions are, on the server side:
 *     TASK *t0 = taskServerInit () ;                       // server side, required initialization 
 *     TASK *t  = taskDispatch (program_name, parameters) ;    // server side launch n program(s) on unspecified computer(s)
 *** On the client side :
 *     TASK *t = taskClientInit (int *argcp, const char **argv, AC_HANDLE h) ; // client side, required, connect to the server tcp port
 *** On both sides
 *  The server and client must then open R/W channels with matching types and names using:
 *     w = taskCreateWriterChannel (t, name, type, max) ;  // will transmit data to that program 
 *     r = taskCreateReaderChannel (t, name, type, max) ;  // will receive data from that program 
 *  When all records have been put to the w channel, call
 *     channelClose (w) ;
 *     pending records will be received by some client(s), and all client will be notified that the channel is closed.
 *
 *  Then use channelGet/Put to transmit data exactly as in the multithreaded interface
 *  The difference is that memory is not shared accross processes, so the data
 *  are transmitted by value, not by address, and the latencies are longer. 
 *
 *  The user guide for this library is provided at the top of the file wh/remote_channel.h
 *
 * Implementationn stategy:
 *  On each side, the channel is used with the standard interface, 
                            channelPut(w,...) ;  channelGet(r, ...) ;
 *  so the renote_channel_class 'inherits' from channel, in C it means we just extend 
 *  the channel struct.
 *  Notice that only the long distance tranmission is unidirectional, it is legal
 *  to Get from the w channel, in parallel with the remote task which is probably
 *  also Getting from this channel, if it has opened as expected a corresponding 
 *  taskReaderChannel. It is also legal to put to the r channels, and these records 
 *  will be consumed locally in addition to the records transmitted from the probably 
 *  existing corresponding taskWriterChannel of the remote task.
 *
 *  Behind the scene, we have 3 running threads linked via channels, plus a control channel

           controller   <----------------------------------------------- ready
                  v                                                        ^
      [w]->echo_writer->{tcp-write-transmit-read}->[reader_channel]->echo_reader->[r]

 *  The 2 echo functions are blocking, amd automatically run in their own thread.
 *  Echo-writer gets a record from the w classic channel (when available), 
 *  and throws it down the tcp connection of the first available client reader (when ready).
 *  Echo-reader, reads a record from the tcp connection (when available), 
 *  puts it on the claasic r channel and sends back a ready signal to the same client.
 *  The tcp.fd file-descriptor is r/w and private to the channel/client/server triplet
 *  If N tasks, (possibly tunning on distinct machines), connect back to the same
 *  writer channel w, they will get records in a round-robin way
 *  each using its own tcp.fd connector. They can in parallel send records back to the server
 *  each using its own clieint-writer, to the same server-reader channel.
 */

#include "ac.h"
#include "wego_.h"
#include <pthread.h>
#include "remote_channel.h"
#include "remote_channel_.h"
#include "remote_task_.h"

extern int connect_socket (const char *host, int port, BOOL no_delay) ;
static int TASK_MAGIC = 12345 ;

static BOOL remote_channel_acknowledge (int fd, const char *channelName, char type, int size, int sig) ;
typedef struct fdmStruct {int fd ; char *buf; } FDM ;

/*********************************************************************/
/*********************************************************************/
/* these routines must be called called by task_init () */
/* find the ip of this machine */
static char *getHostIP (void)
{
  static char ip[64] ;
  static int firstPass = 1 ;
  const char *ccp ;

  if (firstPass)
    {
      AC_HANDLE h = ac_new_handle () ;
      ACEIN ai = aceInCreateFromPipe ("hostname -i", "r", 0, h) ;
      if (ai && aceInCard (ai) &&
	  (ccp = aceInWord (ai)) &&
	  strlen (ccp) < 63
	  )
	{
	  fprintf (stderr, "This host ip is %s\n", ccp) ; 
	  strncpy (ip, ccp, 63) ;
	}
      else
	messcrash ("Command \"hostname -i\" failed, I cannot find the ip of the host I am running on") ;
      ac_free (h) ;
    }
  return ip ;
} /* getHostIP */

/*********************************************************************/
/*********************************************************************/
/* echo functions
 * These functions must be started as separate threads, because they are blocking
 * They serve as transition between the standard channel (blocking) and
 * the shared tcp thread which should never be blocked
 */

typedef enum { ZERO=0, KILL_CLIENT } LIFE_TYPE ;

static void remote_channel_client_echo_writer (void *vp)
{
  CHAN *w = *(CHAN **) vp ;
  int n, n1, nn, nn2 ;
  char *cp, reply[2] ;
  const char *channelName = w->r.channelName ;
  int namLength = strlen (channelName) ;
  BOOL debug = FALSE ;
  int sig = w->r.task->signature ;

  remote_channel_acknowledge (w->r.fd, w->r.channelName, 'W', w->c.size, sig) ;

  nn = w->c.size + 19 + namLength ;
  cp = w->r.writeBuffer ;
  cp[0] = nn & 0xff ; cp[1] = (nn >> 8) & 0xff ; cp[2] = (nn >> 16) & 0xff ;  cp[3] = (nn >> 24) & 0xff ; cp[4] = ':' ;
  cp[5] = 'S' ; cp[6] = ':' ;
  n = w->c.size ;
  cp[7] = n & 0xff ; cp[8] = (n >> 8) & 0xff ; cp[9] = (n >> 16) & 0xff ;  cp[10] = (n >> 24) & 0xff ; cp[11] = ':' ;
  cp[12] = sig & 0xff ; cp[13] = (sig >> 8) & 0xff ; cp[14] = (sig >> 16) & 0xff ; cp[15] = (sig >> 24) & 0xff ; cp[16] = ':' ;
  memcpy (cp+17, channelName, namLength) ;
  cp[17 + namLength] = ':' ;

  cp[nn - 1] = ':' ; /* verif byte */
  
  while ((n1 =  uChannelMultiGet (w, w->r.writeBuffer + 18 + namLength, w->c.size, 1, TRUE)) >= 0)
    {
      BOOL isClosed = n1 > 0 ? FALSE : TRUE ;
      n = 0 ;
      if (isClosed)  
	{
	  cp = w->r.writeBuffer ; cp[16] = '*' ;
	  if (1) 
	    fprintf (stderr, "client_echo_writer snending a channelClose signal %s\n", timeShowNow()) ;
	}
      while (n < nn)
	{
	  int n1 = write (w->r.fd, w->r.writeBuffer, nn - n) ;
	  if (n1 > 0) n += n1 ;
	  else break ;
	}
      if (n != nn)
	{
	  fprintf (stderr, "error n = %d != %d + 6 while writing in remote_channel_echo_writer, reply = %c" 
		   ,n ,w->c.size, reply[0]
		   ) ;		    
	  break ;
	}
      if (debug)
	{
	  void *vp = w->r.writeBuffer + 18 + namLength ;
	  fprintf (stderr, "remote_channel_client_echo_writer exported %d bytes in channel %s  int:%d float:%f double:%g \n"
		   , n
		   , w->r.channelName 
		   , *(int*)vp, *(float*)vp, *(double*)vp
		   ) ;
	  
	}
      if (1)
	{
	  n = 0 ; nn2 = 1 ; memset (reply, 0, 2) ;
	  while (n < nn2)
	    {
	      int n2 = read (w->r.fd, reply + n, nn2 - n) ;
	      if (n2 > 0) n += n2 ;
	      if (n2 <= 0)
		break ;
	    }
	  if (n < nn2 || reply[0] != 'x')
	    {
	      fprintf (stderr, "error n = %d in remote_channel_echo_writer reply[0]=%c\n", n, reply[0]) ;
	      break ;
	    }     
	  if (debug)
	    fprintf (stderr, "--- remote_channel_echo_writer exported a record in channel %s, got x \n", w->r.channelName) ;
	}
      /* this record was successfully transmitted and acknowledged */
    }
  /* the classic w channel is closed, double close it  */
  pthread_mutex_lock (&(w->c.mutex)) ;
  w->r.fd = 0 ;  

  pthread_mutex_unlock (&(w->c.mutex)) ;
  channelClose (w) ;
  return ; /* closes this thread */
} /* remote_channel_client_echo_writer  */ 

/*********************************************************************/

static void remote_channel_client_echo_reader (void *vp)
{
  CHAN *r = *(CHAN **) vp ;
  int i, n, nn ;
  BOOL ok = TRUE;
  char *cp ;
  char prefix[10] ;
  BOOL debug = FALSE ;
  int sig =  r->r.task->signature ;

  remote_channel_acknowledge (r->r.fd, r->r.channelName, 'R', r->c.size, sig) ;
  
  nn = r->c.size + 11 ;
  cp = prefix ;
  cp[0] = nn & 0xff ; cp[1] = (nn >> 8) & 0xff ; cp[2] = (nn >> 16) & 0xff ;  cp[3] = (nn >> 24) & 0xff ; cp[4] = ':' ;
  cp[5] = sig & 0xff ; cp[6] = (sig >> 8) & 0xff ; cp[7] = (sig >> 16) & 0xff ;  cp[8] = (sig >> 24) & 0xff ; cp[9] = ':' ;
  if (debug) 
    fprintf (stderr, "########CR######## Client_echo_reader %s fd=%d\n",  r->r.channelName, r->r.fd) ;
  while (ok)
    {
      int n1 = 0 ;
      if (debug)
	fprintf (stderr, "Client_echo_reader %s waiting for a record\n", r->r.channelName) ;
      n = 0 ;
      while (n  < nn)
	{
	  n1 = read (r->r.fd,  r->r.readBuffer + n, nn - n) ;
	  if (n1 > 0) n += n1 ;
	  if (debug)
	    fprintf (stderr, "remote_channel_client_echo_reader %s got %d/%d bytes : %d\n", r->r.channelName, n, nn, r->r.readBuffer[0]) ;
	  if (n1 <= 0)
	    break ;
	}
      cp = r->r.readBuffer ;
      if (cp[9] == '*')  /* the server classic writer channel is closed */
	{
	  remote_channel_acknowledge (r->r.fd, r->r.channelName, 'R', r->c.size, r->r.task->signature) ;
	  break ;
	}
	  
      if (n1 <= 0 || n<= 0)
	break ;
      ok = TRUE ;
      for (i = 0 ; i < 9 ; i++)
	if (cp[i] != prefix[i])
	  {
	    if (debug) 
	      fprintf (stderr, "remote_channel_client_echo_reader error in %s received prefix char %d in remote_channel_echo_reader",  r->r.channelName, i) ;
	    ok = FALSE ;
	  } 
      if (cp[n - 1] != ':')
	  {
	    fprintf (stderr, "error in received suffix char %d in remote_channel_echo_reader", n - 1) ;
	    ok = FALSE ;
	  }
      if (debug)
	{
	  void *vp = r->r.readBuffer + 5 ;
	  fprintf (stderr, "Client_echo_reader %s received a record, int:%d float:%f double:%g\n"
		   , r->r.channelName
		   , *(int*)vp, *(float*)vp, *(double*)vp
		   ) ;
	}
      if (ok)  /* this record was successfully received, please acknowledge */
	{
	  uChannelMultiPut (r, r->r.readBuffer + 10, r->c.size, 1, TRUE) ;
	  remote_channel_acknowledge (r->r.fd, r->r.channelName, 'R', r->c.size, r->r.task->signature) ;
	}
    }
  /* close the classic r channel is close */
  if (r->c.debug)
    fprintf (stderr, " remote_channel_client_echo_reader closing channel %s %s\n",  r->r.channelName, timeShowNow ()) ;
  channelClose (r) ;
  return ;  /* closes this thread */
} /* remote_channel_client_echo_reader  */

/*********************************************************************/

static void remote_channel_server_echo_reader (void *vp)
{
  CHAN *r = *(CHAN **) vp ;
  FDM fdm ;
  char reply[2] = "x" ;
  BOOL debug = FALSE ;
  int n ;

  /* the tcp thread signals when the read buffer is ready */
  while (channelMultiGet (r->r.fdm, &fdm, 1, FDM) >= 1)
    {
      if (debug)
	{
	  void *vp = &fdm.buf ;
	  fprintf (stderr, " remote_channel_server_echo_reader got fdm fd=%d chan=%s  int:%d float:%f double:%g\n"
		   , fdm.fd, r->r.channelName
		   , *(int*)vp, *(float*)vp, *(double*)vp) ;
	}
      if (! uChannelMultiPut (r, fdm.buf, r->c.size, 1, TRUE))
	break ;
      n = 0 ;
      if (1)
	n = write (fdm.fd, reply, 1) ;

      free (fdm.buf) ; 
      if (debug)
	fprintf (stderr, " remote_channel_server_echo_reader exported %d bytes to fd=%d chan=%s\n"
		 , n, fdm.fd, r->r.channelName) ;
    }
  /* the fdm channel is closed */
  channelClose (r) ;
  return ; /* closes this thread */
} /* remote_channel_server_echo_reader  */

/*********************************************************************/

static void remote_channel_server_echo_writer (void *vp)
{
  CHAN *w = *(CHAN **) vp ;
  int n, n1, nn, fd = 0 ;
  char *cp ;
  BOOL debug = FALSE ;
  int sig = w->r.task->signature ;

  nn = w->c.size + 11 ;
  cp = w->r.writeBuffer ;
  cp[0] = nn & 0xff ; cp[1] = (nn >> 8) & 0xff ; cp[2] = (nn >> 16) & 0xff ;  cp[3] = (nn >> 24) & 0xff ; cp[4] = ':' ;
  cp[5] = sig & 0xff ; cp[6] = (sig >> 8) & 0xff ; cp[7] = (sig >> 16) & 0xff ;  cp[8] = (sig >> 24) & 0xff ; cp[9] = ':' ;
  cp[nn - 1] = ':' ; /* verif byte */
 
  while ((n1 =  uChannelMultiGet (w, w->r.writeBuffer + 10, w->c.size, 1, TRUE)) >= 0)
    {
      BOOL isClosed = n1 > 0 ? FALSE : TRUE ;
      if (isClosed)  
	{
	  cp = w->r.writeBuffer ; cp[9] = '*' ;
	  if (1) 
	    fprintf (stderr, " remote_channel_server_echo_writer sending a channelClose signal %s\n", timeShowNow()) ;
	}
          /* the tcp thread signals when the distant reader is ready */
      if (debug)
	{
	  void *vp =  w->r.writeBuffer + 5 ;
	  fprintf (stderr, "remote_channel_server_echo_writer waiting for r.signal in channel %s :: %s int:%d float:%f double:%g \n"
		   , w->r.channelName, timeShowNow()
		   , *(int*)vp, *(float*)vp, *(double*)vp
		   ) ;
	}
      if (channelMultiGet (w->r.signal, &fd, 1, int)) /* block untill the client channel is ready */
	{
	  n = write (fd, w->r.writeBuffer,  nn) ;
	  if (n != nn)
	    {
	      fprintf (stderr, "remote_channel_server_echo_writer error n = %d != %d + 6 while writing in remote_channel_echo_writer\n"
		       ,n ,w->c.size
		       ) ;		    
	    }
	}
      else
	break ;
     if (debug)
       fprintf (stderr, "remote_channel_server_echo_writer sent record to fd=%d : in channel %s ::%s\n", fd, w->r.channelName, timeShowNow()) ;
    }
  /* the classic w channel is closed */
  return ;  /* closes this thread */
} /* remote_channel_server_echo_writer  */

/*********************************************************************/
/* these calls are protected by the task mutex */
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netinet/tcp.h>
#include <arpa/inet.h>
#include <sys/time.h>
#include <sys/types.h>
#include <unistd.h>
#include <signal.h>
#include <ctype.h>
#include <sys/select.h> 

/*
 * the max fd and the fd_set used to select on client connections.
 * also includes the listen socket.
 */
static int tcp_server_fd_count;
static fd_set open_sockets;
static int listen_fd ;

extern int listen_socket(int) ;
/***********************************************************************
 * tracking clients
 *
 */

static int next_client_serial_number = 1;

static int current_active_clients = 0;

/*
 * There is a struct client for each connection
 */
struct client
{
  int	fd;
  /*
   * file descriptor of the socket to the client
   */

  int	serial_number;
  /*
   * what serial number of this connection
   */

  int	last_transaction_time;
  /*
   * time of the last client transaction
   */
  int	access_type;
  /*
   * access type granted to this client
   *	'r'	may read database
   *	'w'	may read or write database
   */

  char *challenge;
  /*
   * challenge used for authenticating requests
   * from this client.
   */

  /*  AceCommand	look;

   * the magic COMMAND_LOOK structure that a server needs.
   */

  int	saved_option;
  /*
   * the "option" value saved at the end of the last
   * transaction.  In the event of a "encore" command,
   * this option has some meaning about what the command
   * to be "encored" was.
   */

  struct incoming_transaction_rec *incoming_trans;
  struct outgoing_transaction_rec *outgoing_trans;
  /*
   * acetcp transaction state blocks
   */
};

/*
 * a pointer to the struct client for the connection is stored in the
 * array clients.  The index is the socket number.
 */
static Array clients;

static void accept_client (int listen_fd)
{
  struct sockaddr_in there;
  int s;
  struct client *c;
  int one;
  socklen_t n ;
  FILE *log_file = stderr ;

  n = sizeof(there);
  s = accept(listen_fd, (struct sockaddr *) &there, &n);
  if (s < 0)
    {
      perror("accept");
      return;
    }

  /* printf("accept client %d\n",s); */


  /*
   * it is important to set TCP_NODELAY on the socket.  if you
   * do not, transactions will be VERY slow.  Without TCP_NODELAY,
   * the system hopes you will write more data to the socket
   * very soon, so it does NOT transmit the data immediately -- it
   * wants to put the results of your next write in the same
   * packet.
   *
   * 6 is what getprotobyname returns for "tcp" - technically, I'm 
   * cheating but TCP is _always_ protocol number 6 because that
   * is the IP protocol number for TCP.
   */
  one = 1;
  if (setsockopt( s, 6, TCP_NODELAY , &one, sizeof(int)) < 0)
    {
      /*
       * this setsockopt can never fail, so it is ok to crash if it does.
       */
      messcrash("could not set TCP_NODELAY");
    }

  /*
   * make a client record with incoming and outgoing transaction
   * states and other basic initialization.
   */
  c = halloc(sizeof(struct client), 0);
  c->fd = s;

  /*
  c->incoming_trans = incoming_transaction_create();
  c->outgoing_trans = outgoing_transaction_create(s);
  */

  c->saved_option = 0;
  c->serial_number = next_client_serial_number++;
  array(clients, s, struct client *) = c;
  if (s >= tcp_server_fd_count)
    tcp_server_fd_count = s+1;
  FD_SET(s, &open_sockets);

  /*
   * check if their host has access.  We do it this late so that we 
   * have everything set up to send back a properly formatted 
   * error message.  From a computer security perspective, I would
   * just close the connection without saying anything but I think
   * that is likely to cause a lot of confusion when people don't
   * have their configuration correct.
   *
   * The type of access is stored in the client block, but that can
   * be changed later by the client if it responds to the challenge/response
   * authentication.  Otherwise the access remains at the default
   * granted to the host.
   */
  /*   c->access_type = check_tcp_access_list( (unsigned char *) & there.sin_addr ); */
  if (c->access_type < 0)
    {
      /*
       * no access - bummer for you
  
       static char msg[] = "your computer does not have access";
      */

      fprintf(log_file, "\n%s Refusing New Client %d - source IP %d.%d.%d.%d\n",
	      timeShowNow(),  c->serial_number
	      , ((unsigned char *) & there.sin_addr)[0]
	      , ((unsigned char *) & there.sin_addr)[1]
	      , ((unsigned char *) & there.sin_addr)[2]
	      , ((unsigned char *) & there.sin_addr)[3]
	      ) ;

      /* last param used to be  inet_ntoa(*(struct in_addr *) &there));
       * but that is not decoded on my linux 64
       */

      /*
       * We increment the client count because we will decrement it
       * when we close this socket.
       */
      ++current_active_clients;

      /*
       * send back a meaningful error message.
       */
      /*
      outgoing_transaction_send(c->outgoing_trans, 
				'E', msg, sizeof(msg) );
      */
      shutdown(c->fd, 2);
    }
  else
    {
      /*
       * some kind of valid access - send the challenge.
       */
      fprintf(log_file, "\n%s New Client %d, %d active clients\n", 
	      timeShowNow(), c->serial_number, ++current_active_clients);
      /*
      outgoing_transaction_send(c->outgoing_trans, 
				'C', c->challenge, strlen(c->challenge)+1);
      */
    }
}

static BOOL task_reply_ok (int fd, int ic)
{
  int n, nn ;
  char reply[13] ;
  reply[0] = 12 ;
  reply[1] = 0 ;
  reply[2] = 0 ;
  reply[3] = 0 ;
  reply[4] = ':' ;
  reply[5] = 1 ;
  reply[6] = ':' ;
  reply[7] = ic & 0xff ;
  reply[8] = (ic >> 8) & 0xff ;
  reply[9] = (ic >> 16) & 0xff ;
  reply[10] = (ic >> 24) & 0xff ;
  reply[11] = ':' ;
  reply[12] = 0 ;
  nn = 12 ; n = 0 ;
  while (n < nn)
    {
      int n1 = write (fd, reply, nn - n) ;
      if (n1 > 0) n += n1 ;
      if (n1 <= 0)
	return FALSE ;
    }
  return TRUE ;
}

static BOOL task_process_incoming (TASK *task, char *buffer, int size, int fd)
{ 
  int ic ;
  int type = buffer[0] ;
  int signature = 0 ;
  CHAN *c = 0 ;
  BOOL ok = FALSE ;
  BOOL debug = FALSE ;

  if (debug)
    fprintf (stderr, "task_process_incoming type=%c\n", type) ;
  buffer[size-1] = 0 ;  /* remove the terminal ':' */
  switch (type)
    {
    case 'R':  /* Read check */
    case 'U':  /* Write first check */
    case 'W':  /* Write check */
    case 'S':  /* Sent record */
      size = (buffer[2]  & 0xff) + ((buffer[3] & 0xff) << 8) + ((buffer[4] & 0xff) << 16) + ((buffer[5] & 0xff) << 24) ;
    case 'C':   /* close the channel */
    case 'A':  /* acknowledge reception of a record */
      signature = (buffer[7]  & 0xff) + ((buffer[8] & 0xff) << 8) + ((buffer[9] & 0xff) << 16) + ((buffer[10] & 0xff) << 24) ;
      if (type == 'S')
	{
	  char *cq = strchr (buffer+12, ':') ;
	  if (cq) *cq = 0 ;
	}
      dictFind (task->channelDict, buffer+12, &ic) ;
      if (ic >= 0 && task->channels &&
	  ic < arrayMax (task->channels) &&
	  (c = arr (task->channels, ic, CHAN*))
	  )
	{
	  if (type != 'A' && type != 'R' && c->c.isClosed)
	    fprintf (stderr, "......bizarre %s is closed type = %c\n", buffer + 12, type) ;
	  if (type == 'A' || (size == c->c.size && signature == c->r.task->signature))
	    {
	      if (c->r.readBuffer  && (type == 'U' || type == 'W'))
		{
		  /* we can now start the echo_reader */ 
		  c->r.fd = fd ;
		  if (task_reply_ok (fd, ic))
		    {
		      if (! c->r.wegoReader)  /* we are on server side */
			c->r.wegoReader = wego_group_go (task->group, remote_channel_server_echo_reader, &c, CHAN *) ;
		      ok = TRUE ;
		    }
		}
	      if (c->r.writeBuffer && type == 'R')
		{
		  /* we can now start the echo_writer */ 
		  c->r.fd = fd ;
		  
		  if (1)
		    {
		      pthread_mutex_lock (&(c->c.mutex)) ;
		      c->r.fd = fd ; 
		      pthread_mutex_unlock (&(c->c.mutex)) ;
		    }



		  if (task_reply_ok (fd, ic))
		    {
		      if (! c->r.wegoWriter) /* we are on server side */
			c->r.wegoWriter = wego_group_go (task->group, remote_channel_server_echo_writer, &c, CHAN *) ;
		      ok = TRUE ; 
		      if (debug) 
			fprintf (stderr, "#####SR######### echo_writer %s fd = %d type =%c\n", c->r.channelName, fd, type) ;
		    }
		}
	      if (type == 'S' && c->r.readBuffer)
		{ 
		  int namLength = strlen (buffer+12) ;
                  FDM fdm ;

		  fdm.fd = fd ;

		  if (debug) 
		    fprintf (stderr, "#####SR######### echo_writer %s fd = %d type =%c buffer=%s\n", c->r.channelName, fd, type, buffer+13) ;

		  fdm.buf = malloc (size) ;
		  /* write the record */
		  memcpy (fdm.buf, buffer + 13 + namLength, size) ;
		  /* signal */
		  channelPut (c->r.fdm, &fdm, FDM) ; 
		         /* this reply is now sent by  remote_channel_server_echo_reader after 
			  * copying fdm.buf to the classic channel
			  * we use one fdm.buf per fdm
			  * since several client could want to write in the same reader channel
			  */
		  	
		  ok = TRUE ;
		}	      
	      if (type == 'A')
		{
		  /* signal a new record may now be sent */
		  if (debug) 
		    fprintf (stderr, "#####SW######### task_process_incoming  %s fd = %d type =%c\n", c->r.channelName, fd, type) ;
		  if (c->r.writeBuffer)
		    channelPut (c->r.signal, &fd, int) ; 
		  /* no reply needed */
		  ok = TRUE ;
		  /* this record was successfully transmitted and acknowledged */
		}	      
	      if (type == 'C')
		{
		  /* the client does not have exclusive access, so nothing should be done */
		  /* no reply needed */
		  ok = TRUE ;
		  /* this record was successfully transmitted and acknowledged */
		}	      
	    }
	}
      break ;
    }
  return ok ;
}

static int taskTcpInit (TASK *task)
{
  int x ;
  
  /*
   * find which hosts have any access at all to the database
   */
  /* read_tcp_access_list(); */


  /*
   * IGNORE sigpipe - return values from the write tell us
   * when we hit a broken pipe.
   */
  signal (SIGPIPE, SIG_IGN);

  /*
   * set up to listen for incoming connections
   * initialize file descriptors
   */
  FD_ZERO (&open_sockets);
  clients = arrayHandleCreate (256, struct client *, task->h) ;

  /* starting from the port proposed on the command line, find a free port */
  for (x = 0 ; x < 1000 ; x++, task->port += 3)
    {
      listen_fd = listen_socket (task->port) ;
      if (listen_fd >= 0)  /* select this port, it is available */
	break ;
    }
  if (listen_fd < 0)
    messcrash("unable to listen on port %d", task->port);
  FD_SET(listen_fd, &open_sockets);
  tcp_server_fd_count = listen_fd + 1 ;

  /* clients = arrayCreate(tcp_server_fd_count+10, struct clients *) ; */

  fprintf (stderr, "waiting for client on port %d\n", task->port) ;
 
  return task->port ;
}


static void taskTcpListen (void *vp)
{
  TASK *task = *(TASK **) vp ;
  int x;
  int selecterror;
  static fd_set read_fds;
  FILE *log_file = stderr ;
  int safe = 1024 ;
  char *buffer = halloc (safe, 0) ;
  BOOL debug = FALSE ;

  /*
   * count select errors - we get one once in a while for reasons that
   * are not clear to me.  I think it may relate to alarm signals 
   * originating somewhere else in the database.  The server only
   * quits if it gets too many consecutive select errors.
   */
  selecterror = 0;

  for (;;)
    {
      read_fds = open_sockets ;

      /*
       * wait for activity from the listening socket or any client.
       *
       * If you need a timeout for server inactivity, do it in this select.
       */
      if (select (tcp_server_fd_count, &read_fds, NULL, NULL, NULL) < 0)
	{
	  perror("select");
	  if (selecterror++ > 5)
	    messcrash("too many select errors");
	}
      else
	selecterror = 0;

      /*
       * if there is a new connection, accept it
       */
      if (FD_ISSET(listen_fd, &read_fds))
	{
          if (1) FD_CLR (listen_fd, &read_fds) ;
	  accept_client (listen_fd) ;
	}

      /*
       * process transactions from any connection that has data in it.
       */
      for (x=0; x<tcp_server_fd_count; x++)
	{
	  struct client *c;

	  /*
	   * nothing on this fd
	   */
	  if (! FD_ISSET(x, &read_fds))
	    continue;

	  /*
	   * find the client
	   */
	  c = arr(clients, x, struct client *);
	  if (!c)
	    {
	      /*
	       * how odd - there is no client, but our select says
	       * we have data waiting on that socket.
	       */
	      close(x);
	      continue;
	    }

	  c->last_transaction_time = time(0);
	  /* 
	   * process client data
	   */
	  if (1) /* incoming_transaction(c->incoming_trans, x, FALSE, do_process_incoming , c)) */
	    {
	      char *cp, prefix[5] ;
	      int n, nn, n1 = 1 ; ;

	      nn = 5 ; n = 0 ; 
	      while (n1 > 0 && n < nn)
		{
		  n1 = read (c->fd, prefix, nn - n) ;
		  if (n1 > 0) n += n1 ;
		}
	      if (n == nn)
		{
		  cp = prefix ;
		  nn = (cp[0]  & 0xff) + ((cp[1] & 0xff) << 8) + ((cp[2] & 0xff) << 16) + ((cp[3] & 0xff) << 24) ;
		  if (debug)
		    fprintf (stderr, "taskTcpListen client %d message length %d bytes\n", x, nn) ;
		  n = 0 ; nn -= 5 ; 
		  if (nn > 0 && cp[4] == ':')
		    {
		      int nerr = 0 ;
		      if (nn > safe)
			{
			  ac_free (buffer) ;
			  safe = nn ;
			  buffer = halloc (safe, 0) ;
			}
		      while (n < nn)
			{
			  int n1 = read (c->fd, buffer + n, nn - n) ;
			  if (n1 > 0) n += n1 ;
			  if (n < 0 || n1 <= 0) break ;
			  if (n1 == 0 && nerr++ > 5) break ;
			}
		      
		      if (n == nn && nn >= 7 && buffer[1] == ':' && buffer[6] == ':' && buffer[11] == ':' && buffer[nn-1] == ':' &&
			  task_process_incoming (task, buffer, nn, c->fd)) 
			continue ;
		    }
		}

	      /* badly formatted query or 
	       * too many EOF - shut it down
	       */
 	      fprintf(log_file, "\n%s Closing Client %d, %d active clients\n",
		      timeShowNow(), c->serial_number, --current_active_clients);

	      /*
	      incoming_transaction_free(c->incoming_trans);
	      outgoing_transaction_free(c->outgoing_trans);
	      */
	      if (debug) printf("close client %d\n",x);
	      close(x);
	      FD_CLR(x, &open_sockets);
	      /*   bug: leaking struct client *c    */
	      messfree(c->challenge) ;
	      messfree(c);
	      arr(clients, x, struct client *) = 0;
	    }
	  else
	    {
	      /*
	       * it is a valid transaction - it was handled by the callback
	       * process_incoming()
	       */
	    }
	}
    } /* infinite loop end, for (;;) */
  ac_free (buffer) ;
  return ;
} /* taskTcpListen  */

/*********************************************************************/
/* check that the server side accepts the definition of this channel */
static BOOL remote_channel_check (TASK *task, int fd, const char *channelName, char type, int size)
{
  BOOL ok = TRUE ;
  int n, n1, nn = strlen (channelName) + 18 ;
  char *cp, buf[nn+1], reply[18] ;
  int sig = task->signature ;
  int namLength = strlen (channelName) ;

  cp = buf ;
  cp[0] = nn & 0xff ; cp[1] = (nn >> 8) & 0xff ; cp[2] = (nn >> 16) & 0xff ; cp[3] = (nn >> 24) & 0xff ; cp[4] = ':' ;
  cp[5] = type ; cp[6] = ':' ;
  cp[7] = size & 0xff ; cp[8] = (size >> 8) & 0xff ; cp[9] = (size >> 16) & 0xff ; cp[10] = (size >> 24) & 0xff ; cp[11] = ':' ;
  cp[12] = sig & 0xff ; cp[13] = (sig >> 8) & 0xff ; cp[14] = (sig >> 16) & 0xff ; cp[15] = (sig >> 24) & 0xff ; cp[16] = ':' ;
  memcpy (cp+17, channelName, namLength) ;

  cp[nn - 1] = ':' ; /* verif byte */
  cp[nn] = 0 ;
  n = write (fd, buf, nn) ;
  if (n != nn)
    messcrash ("remote_channel_check could only write %d/%d bytes during authentication", n, nn) ;
  
  n = 0 ; nn = 12 ; memset (reply, 0, 13) ;
  while (n  < nn)
    {
      n1 = read (fd,  reply + n, nn - n) ;
      if (n1 <= 0) break ;
      n += n1 ; 
    }
  if (n < nn || reply[0] != 12 || reply[4] != ':' || reply[5] != 1 || reply[6] != ':')
    {
      fprintf (stderr, "remote_channel_check received an invalid reply %d %d %d %d %d %d %d\n"
	       , reply[0]
	       , reply[1]
	       , reply[2]
	       , reply[3]
	       , reply[4]
	       , reply[5]  
	       , reply[6]
	       ) ;
      ok = FALSE ;
    }

  return ok ;
} /* remote_channel_check */

/*********************************************************************/

static BOOL remote_channel_acknowledge (int fd, const char *channelName, char type, int size, int sig)
{
  BOOL ok = TRUE ;
  int n, nn = strlen (channelName) + 18 ;
  char *cp, buf[nn+1] ;

  cp = buf ; 
  cp[0] = nn & 0xff ; cp[1] = (nn >> 8) & 0xff ; cp[2] = (nn >> 16) & 0xff ; cp[3] = (nn >> 24) & 0xff ; cp[4] = ':' ;
  cp[5] = 'A' ; cp[6] = ':' ;
  cp[7] = size & 0xff ; cp[8] = (size >> 8) & 0xff ; cp[9] = (size >> 16) & 0xff ; cp[10] = (size >> 24) & 0xff ; cp[11] = ':' ;
  cp[12] = sig & 0xff ; cp[13] = (sig >> 8) & 0xff ; cp[14] = (sig >> 16) & 0xff ; cp[15] = (sig >> 24) & 0xff ; cp[16] = ':' ;

  memcpy (cp+17, channelName, strlen (channelName)) ;
  cp[nn - 1] = ':' ; /* verif byte */
  cp[nn] = 0 ;
  n = write (fd, buf, nn) ;
  if (n != nn)
    messcrash ("remote_channel_acknowledge could only write %d/%d bytes during authentication", n, nn) ;

  return ok ;
} /* remote_channel_acknowledge */

/*********************************************************************/
/****************** Public interface *********************************/
/*********************************************************************/
static pthread_mutex_t taskMutex ;

static void taskDestroy (void *vp)
 { 
   TASK *task = (TASK *) vp ;

   if (task)
     {
       AC_HANDLE h = task->h ;
       int i ;
       task->h = 0 ; /* avoid looping */
       if (task->serverSide)
	 {
	   for (i = 0 ; i < arrayMax (task->clientTasks) ; i++)
	     {
	       TASK *tC = array (task->clientTasks, i, TASK *) ;
	       ac_free (tC) ;
	     }
	 }
       ac_free (h) ;
     } 
 } /* taskDestroy */

static void taskDispatchDestroy (void *vp)
 { 
   TASK *task = (TASK *) vp ;
   if (task && task->h)
     {
       AC_HANDLE h = task->h ;
       LIFE_TYPE t = KILL_CLIENT ;

       task->h = 0 ; /* avoid looping */
       if (task->serverSide &&  task->lifeChannel)
	   channelPut ( task->lifeChannel, &t, LIFE_TYPE) ; 
       ac_free (h) ;
     }
 } /* taskDispatchDestroy */

/*********************************************************************/
/* Client side: wait for a server request to kill the client */
static void taskClientLifeManager (void *vp)
{
  TASK *task = *(TASK **) vp ;
  LIFE_TYPE t = 0 ;
  int n = -1 ;

  n = channelGet (task->lifeChannel, &t, LIFE_TYPE) ;
    
  if (t == KILL_CLIENT)
    fprintf (stderr, "KILL_CLIENT requested on task->lifeChannel, I exit(102)\n") ;
  else if (t)
    fprintf (stderr, "Data received on task->lifeChannel, I exit(%d)\n", 100 + t) ;
  if (n == 0)  /* server is dead, we exit */
    fprintf (stderr, "lifeChannel is dead, I exit(101) %s\n", timeShowNow()) ;
  exit (101 + n) ;
} /* taskClientLifeManager */

/*************/
/*
        tcagggtccacacaaagctctcggatcccc
             aggacagcaaagccacaatgttc
   gaacattgtggctttgctgtcct
*/
/* create the tcp connection to the parent */
TASK *taskClientInit (int *argcp, const char **argv, AC_HANDLE h)
{
  static TASK *task = 0 ;
  int fd ;
  BOOL debug = FALSE ;
  int max_threads = 8 ; 
  const char *host = 0 ; int port = 0 ; int signature = 0 ;  int cycle = 0 ;

  if (task)
    return task ;
  task = (TASK *) halloc (sizeof (struct taskStruct), h) ;  

  pthread_mutex_lock (&(taskMutex)) ;

  messErrorInit (argv[0]) ;
  wego_max_threads (8) ; /* impose at least 8 */
  getCmdLineInt (argcp, argv, "--max_threads", &max_threads);
  wego_max_threads (max_threads) ;
  arrayReport (-2) ;  /* blocks static array counts */

  if (!getCmdLineOption (argcp, argv, "+rChan_h", &host))
    { fprintf (stderr, "Missing argument +rChan_h")  ; return 0 ; }
  if (!getCmdLineInt (argcp, argv, "+rChan_p", &port))
    { fprintf (stderr, "Missing argument +rChan_p")  ; return 0 ; }
  if (! getCmdLineInt (argcp, argv, "+rChan_s", &signature))
    { fprintf (stderr, "Missing argument +rChan_s")  ; return 0 ; }
  if (! getCmdLineInt (argcp, argv, "+rChan_c", &cycle))
    { fprintf (stderr, "Missing argument +rChan_c")  ; return 0 ; }
 
  task->h = ac_new_handle () ;
  task->magic = &TASK_MAGIC ;
  task->host = (host ? host  : "localhost") ;
  task->port = port ; 
  task->lifeCycle = cycle ;
  task->signature = signature ;
  task->timeout = 300 ;  /* unfortunately, we don't use the timeout in acetcp */
  
  task->clientSide = TRUE ;
  task->serverSide = FALSE ;
  task->group = wego_new_group () ;
  task->channelDict = dictHandleCreate (32, task->h) ;
  task->channels = arrayHandleCreate (32, CHAN*, task->h) ;
  blockSetFinalise (task, taskDestroy) ;
  fd = connect_socket (task->host, task->port, TRUE) ;
  if (fd <= 0)
    messprintf ("taskClientInit %s %d failed : %s\n", host, port, "connection refused") ;
  else
    {
      if (debug)
	fprintf (stderr, "taskClientInit %s %d  succeeded : %s\n", host, port, timeShowNow ()) ;
      close (fd) ;
    }
  task->client_fd = fd ;
  
  pthread_mutex_unlock (&(taskMutex)) ;
  if (debug) fprintf (stderr, "taskClientInit tries to taskCreateReaderChannel (task, __Life_channel_%d__)\n", task->lifeCycle) ;
  if (1)
    {
      char chNam[64] ;
      sprintf (chNam, "__Life_channel_%d__", task->lifeCycle) ;
      task->lifeChannel = taskCreateReaderChannel (task, chNam, 12, LIFE_TYPE) ;
      if (task->lifeChannel)
	{
	  if (debug) fprintf (stderr, "taskClientInit tries to taskCreateReaderChannel (task, __Life_channel_%d__) success\n", task->lifeCycle) ;
	  if (1) wego_group_go (0 /* task->group */, taskClientLifeManager, &task, TASK *) ; /* group==0: on close do not wait on this task */
	  if (0) channelDebug ( task->lifeChannel, TRUE, "LifeChannel") ;
	}
      else
	{
	  fprintf (stderr, "taceClientInit FAILED to connect to the server lifeChannel %s, exiting\n", chNam) ;
	  exit (0) ;
	}
    }
  return task ;
} /* taskClientInit  */

/*********************************************************************/

void taskSwitchContext (TASK *task)
{
  char buf[64] ;
  fprintf (stderr, "taskSwitchContext\n") ;  
  
  if (1) pthread_mutex_lock (&(taskMutex)) ;

  if (task &&  task->serverSide &&  task->lifeChannel)
    {  
      channelClose (task->lifeChannel) ;
      sleep (3) ;
      task->lifeCycle++ ;
      if (0)
	{
	  task->port = 0 ;
	  pthread_mutex_unlock (&(taskMutex)) ;
	  taskServerInit (0, 0, 0) ;  /* change port etc */
	  pthread_mutex_lock (&(taskMutex)) ;
	}
    }
  if (1) pthread_mutex_unlock (&(taskMutex)) ;
  sprintf (buf, "__Life_channel_%d__", task->lifeCycle) ;
  task->lifeChannel = taskCreateWriterChannel (task, buf,  12, LIFE_TYPE) ;
  sleep (1) ;
  return ;
} /* taskSwitchContext  */

/*********************************************************************/
/* launch a task server which will listen and communicate with the remote tasks */
/* Server side
 *
 * Open a listening port,  start a listener
 * This port will be called each time a remote task wishes to establish
 * a remote channel. After authentification, a read/write file descriptor 
 * private to this conection is established and an echo funtions
 * is launched as separate go-routines to transmit the records of the
 * appropriate channel down the connection 
 */
TASK *taskServerInit (const char *program_name, TASK_CONFIG config, AC_HANDLE h)
 {
   static TASK *task = 0 ;

   if (! task || ! task->port)
     {
       int k ;
       int cycle = task ? task->lifeCycle : 0 ;
       char buf[64] ;

        if (! task)
	 {
	   pthread_mutex_lock (&(taskMutex)) ;

	   messErrorInit (program_name) ;
	   wego_max_threads (8)  ;
	   arrayReport (-2) ;  /* blocks static array counts */

	   task = (TASK *) halloc (sizeof (TASK), h) ;
	   
	   task->program_name = program_name ;
	   task->config = config ;
	   task->host = getHostIP () ;
	   task->serverSide = TRUE ;
	   task->clientSide = FALSE ;
	   task->port = 12349 + task->lifeCycle ;
	 }
       else
	 {
	   cycle++ ;
	   wego_task_destroy (task->wegoListener) ;
	   if (0) ac_free (task->h) ;
	   task->port += cycle ;
	 }

       task->h = ac_new_handle () ;
       task->lifeCycle = cycle ;
       task->channelDict = dictHandleCreate (256, task->h) ;
       task->channels = arrayHandleCreate (256, CHAN *, task->h) ;
       task->clientTasks = arrayHandleCreate (256, TASK *, task->h) ;
       if (taskTcpInit (task))  /* we are on server side */
	 task->wegoListener = wego_group_go (task->group, taskTcpListen, &task, TASK *) ;
       k =  task->port ^ ((int)(0xffffffff & ((char*)(&task) - (char *)0))) ^ ((int)getpid ()) ;
       { int i = k & 0xf ; while (i--) randint () ; }
       task->signature = ((randint () ) ^ k) ;

       pthread_mutex_unlock (&(taskMutex)) ;
       /* used for shake hand */
       sprintf (buf, "__Life_channel_%d__", task->lifeCycle) ;
       task->lifeChannel = taskCreateWriterChannel (task, buf, 12, LIFE_TYPE) ;
     }
   sleep (1) ; 
   return task ;
 } /* taskServerInit */

/*********************************************************************/
/* taskDispatch defines the server side
 * 
 * launch a remote task which will call us on our taskServerPort 
*/
#define myPort 12347
static const char *taskProgramName (TASK *task0, const char *program_name, AC_HANDLE h)
{
  const char *pnam = 0 ;
  char *pnam2 = 0 ;
  if (! program_name)
    pnam = task0->program_name ;
  else if (program_name[0] == '/')
    pnam = program_name ;
  else if (task0->program_name)
    {
      /* if task0->program_name == ./bin.LINUX/foo, and program_name == bar, 
       * return ./bin.LINUX/bar, i.e assume server and client are in same directory
       */
      const char *ccp = task0->program_name + strlen (task0->program_name) - 1 ;
      while (*ccp != '/' && ccp > task0->program_name) ccp-- ;
      if (*ccp == '/')
	{
	  pnam2 = halloc (2 + (ccp - task0->program_name) + strlen (program_name), h) ;
	  memcpy (pnam2, task0->program_name, ccp - task0->program_name + 1) ;
	  strcpy (pnam2 + (ccp - task0->program_name + 1), program_name) ;
	}
      else /* no / in server name, assume ./ */
	{
	  pnam2 = halloc (3 + strlen (program_name), h) ;
	  memcpy (pnam2, "./", 2) ;
	  strcpy (pnam2 + 2, program_name) ;
	}
      pnam = pnam2 ;
    }
  if (! pnam)
    messcrash ("taskServerInit and taskProgramName both called with NULL program_name") ;
 #include <unistd.h>
  if (access (pnam, X_OK))  /* access returns the error type, 0 means ok, if error undo the file name modification */
    {
      pnam = program_name ; 
      if (access (pnam, X_OK))
	messcrash ("taskDispatch cannot locate the execulatable file %s, sorry\n", pnam) ;
    }
  return pnam ;  
} /* taskProgramName */

/*************************/

TASK *taskDispatch (const char *programName, char *params, AC_HANDLE h)
 {
   TASK *task0, *task = 0 ;
   char *cp ;
   
   if (1) pthread_mutex_lock (&(taskMutex)) ;

   task0 = taskServerInit (0, 0, h) ;  /* can be called several times */
   if (! task0->signature)
     {
       pthread_mutex_unlock (&(taskMutex)) ;
       return 0 ;
     }
   task = (TASK *) halloc (sizeof (TASK), h) ;
   array (task0->clientTasks, arrayMax (task0->clientTasks), TASK *) = task ;

   task->h = ac_new_handle () ;
   task->host = task0->host ;
   task->port = task0->port ;
   task->signature = task0->signature ;
   task->clientSide = FALSE ;
   task->serverSide = TRUE ;
   task->lifeChannel = task0->lifeChannel ;
   task->lifeCycle = task0->lifeCycle ;

   task->command = hprintf (task->h
			    , "%s +rChan_h %s +rChan_p %d +rChan_s %d +rChan_c %d %s"
			    , taskProgramName (task0, programName, h)
			    , task->host
			    , task->port 
			    , task->signature
			    , task->lifeCycle
			    , params ? params : ""
			    ) ;

   task->channelDict = task0->channelDict ;
   task->channels = task0->channels ;
   blockSetFinalise (task, taskDispatchDestroy) ;

   if (1) pthread_mutex_unlock (&(taskMutex)) ;

   switch (task0->config)
     {
     case TASK_SUBMIT:
       cp = hprintf (task->h, "submit taskDispatch.%d_%d \"%s\" \n"
		     , task0->lifeCycle
		     , arrayMax (task0->clientTasks)
		     , task->command
		     ) ;
       if (1) fprintf (stderr, "##taskDispatch:  %s\n", cp) ;
       system (cp) ;
       break ;
     case TASK_LOCAL: 
     default:
       cp = hprintf (task->h, "%s &"
		     , task->command
		     ) ;
       system (cp) ;
       break ;
     }

   return task ;
 } /* taskDispatch */

/*********************************************************************/ 

static void uRemoteChannelDestroy (void *vp)
{
  CHAN *c = (CHAN *) vp ;
  if (c && c->c.magic)
    {
      c->c.magic = NULL ;
      
      if (c->r.task && c->r.id >= 0  && 
	  arrayExists (c->r.task->channels) &&
	  c->r.id < arrayMax (c->r.task->channels)
	  )
	array (c->r.task->channels, c->r.id, CHAN*) = 0 ; /* prevent double deallocation */

      if (c->r.wegoReader)
	wego_task_destroy (c->r.wegoReader) ;
      if (c->r.wegoWriter)
	wego_task_destroy (c->r.wegoWriter) ;
      c->r.wegoReader = c->r.wegoWriter = 0 ;

      if (c->r.fd) { close (c->r.fd) ; c->r.fd = 0 ; }
      pthread_mutex_destroy (&(c->c.mutex)) ;
      pthread_cond_destroy (&(c->c.notFull)) ;
      pthread_cond_destroy (&(c->c.notEmpty)) ;

      ac_free (c->c.h) ;
    }
  return ;
} /* uRemoteChannelDestroy */

/*********************************************************************/ 
/* clientSide : the client wants to connect to the server
 *      construct immediatly the tcp connection
 * serverSide : the server hopes to be contacted by the future distant task
 *      the tcp connection will be created only when the client establishes the connection
 */

CHAN *uTaskCreateReaderChannel (TASK *task, const char *channelName, int cmax, int size)
{
  CHAN *c = 0 ;
  int id, fd = 0 ;
  char *chNam =  hprintf (task->h, "%s#%d", channelName, task->lifeCycle) ;

  if (1)   pthread_mutex_lock (&(taskMutex)) ;
  if (!task || task->clientSide)
    {
      fd = connect_socket (task->host, task->port, TRUE) ;
      if (fd < 0)
	{ c = 0 ; goto done ; }
      fprintf (stderr, "taskCreateWriterChannel %s successfully connected to %s %d\n"
	       , chNam, task->host, task->port) ;
      if (! remote_channel_check (task, fd, chNam, 'R', size))
	{  c = 0 ; goto done ; }
      fprintf (stderr, "taskCreateWriterChannel %s successfully talks to %s %d\n"
	       , chNam, task->host, task->port) ;
    }
  
  if (strlen(chNam) > 255)
    messcrash ("The name of a remote channel should not exceed 255 characters: %s", chNam) ;
  if (! dictAdd (task->channelDict, chNam, &id))
    messcrash ("double declaration of the same channel name in taskCreateReaderChannel: %s", chNam) ;
  c = uChannelCreate (cmax, size, task->h) ;
  blockSetFinalise (c, uRemoteChannelDestroy) ;
  c->r.readBuffer = halloc (size + 5, c->c.h) ;

  array (task->channels, id, CHAN*) = c ;
  pthread_mutex_lock (&(c->c.mutex)) ;
  c->r.fd = fd ; c->r.id = id ; 
  c->r.task = task ;

  strcpy (c->r.channelName, chNam) ;
  pthread_mutex_unlock (&(c->c.mutex)) ;

  if (task->clientSide)  /* the connection with the server is established, start wroking */
    /* group==0: on close do not wait on this task */
    c->r.wegoReader = wego_group_go (0 /* task->group */, remote_channel_client_echo_reader, &c, CHAN *) ;
  else
    {
      c->r.signal = uChannelCreate (1, sizeof(int), c->c.h) ;
      c->r.fdm = uChannelCreate (1024, sizeof(FDM), c->c.h) ;
    }  

 done:
  if (1)
    fprintf (stderr, "taskCreateReaderChannel created %s\n", chNam) ;
  if (1) pthread_mutex_unlock (&(taskMutex)) ;
   return c ;
} /* uTaskCreateReaderChannel */

/*********************************************************************/
/* clientSide : the client wants to connect to the server
 *      construct immediatly the tcp connection
 * serverSide : the server hopes to be contacted by the future distant task
 *      the tcp connection will be created only when the client establishes the connection
 */

CHAN *uTaskCreateWriterChannel (TASK *task, const char *channelName, int cmax, int size)
{
  CHAN *c = 0 ;
  int id, fd = 0 ;
  char *chNam =  hprintf (task->h, "%s#%d", channelName, task->lifeCycle) ;
  int namLength = strlen (chNam) ;
  char buf[64] ;

  if (1) pthread_mutex_lock (&(taskMutex)) ;
  if (! task || task->clientSide)
    {
      fd = connect_socket (task->host, task->port, TRUE) ;
      if (fd < 0)
	{  c = 0 ; goto done ; }
      fprintf (stderr, "taskCreateWriterChannel %s successfully connected to %s %d\n"
	       , chNam, task->host, task->port) ;
      if (! remote_channel_check (task, fd, chNam, 'U', size))
	{  c = 0 ; goto done ; }
       fprintf (stderr, "taskCreateWriterChannel %s successfully talks to %s %d\n"
	       , chNam, task->host, task->port) ;
   }
  
  sprintf (buf, "__Life_channel_%d__", task->lifeCycle) ;
  if (strlen(chNam) > 255)
    messcrash ("The name of a remote channel should not exceed 255 characters: %s", chNam) ;
  if (! dictAdd (task->channelDict, chNam, &id) && strcmp (chNam, buf))
    messcrash ("double declaration of the same channel name in taskCreateWriterChannel: %s", chNam) ;
  c = uChannelCreate (cmax, size, task->h) ;
  blockSetFinalise (c, uRemoteChannelDestroy) ;
  c->r.writeBuffer = halloc (size + 15 + namLength, c->c.h) ;
  strcpy (c->r.channelName, chNam) ;

  array (task->channels, id, CHAN*) = c ;
  c->r.fd = fd ; c->r.id = id ;
  c->r.task = task ;

  if ( task->clientSide)
    {
      /* group == task->group: on close do  wait on this task */
      c->r.wegoWriter = wego_group_go (task->group, remote_channel_client_echo_writer, &c, CHAN *) ;
    }
  else
    {
      c->r.signal = uChannelCreate (1024, sizeof(int), c->c.h) ;
    }

 done:
  if (1)
    fprintf (stderr, "taskCreateWriterChannel created %s\n", chNam) ;
  if (1) pthread_mutex_unlock (&(taskMutex)) ;
  return c ;
} /* uTaskCreateWriterChannel */

/*********************************************************************/
/*********************************************************************/
/*********************************************************************/


#ifdef JUNK
// example from the web
Example - Coin Toss

Lets start off with a Coin Toss example. Suppose we would like to flip a coin N times in parallel. Using par, we could write the following:

import (
  "code.google.com/p/go.net/context"
  "github.com/savaki/par"
  "math/rand"
)

func flipCoin(results chan CoinFlip) par.RequestFunc {
  return func(context.Context) error {
    if rand.Intn(2) == 0 {
      results <- Heads
    } else {
      results <- Tails
    }
    return nil
  }
}

func simple() {
  flips := 10

  // 1. create a channel to hold your results
  results := make(chan CoinFlip, flips)

  // 2. create a channel of requests
  requests := make(chan par.RequestFunc, flips)
  for flip := 0; flip < flips; flip = flip + 1 {
    requests <- flipCoin(results)
  }
  close(requests)

  // 3. execute the flips in parallel
  _ = par.Requests(requests).Do()

  // 4. results channel now has your results and can be ranged over
}

#endif
