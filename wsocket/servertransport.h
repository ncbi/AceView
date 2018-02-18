/*  File: transport.h
 *  Author: Ed Griffiths (edgrif@sanger.ac.uk)
 *  Copyright (c) J Thierry-Mieg and R Durbin, 2000
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.crbm.cnrs-mop.fr
 *
 * Description: This header defines the interface that a transport layer 
 *              must provide for the acedb server code. The transport
 *              layer might be sockets, RPC or whatever, but it must
 *              provide certain services for the server.
 * HISTORY:
 * Last edited: Jun  5 08:49 2000 (edgrif)
 * Created: Fri Feb 11 09:59:06 2000 (edgrif)
 * CVS info:   $Id: servertransport.h,v 1.1.1.1 2002/07/19 20:23:35 sienkiew Exp $
 *-------------------------------------------------------------------
 */
#ifndef DEF_TRANSPORT_H
#define DEF_TRANSPORT_H

/* It's assumed in the protocol that ints are 4 bytes, this will stop the    */
/* code from compiling, but we also do a run time check (although I can't    */
/* imagine that a 32-bit binary can really be run with ints that aren't 4    */
/* bytes...).                                                                */
/* n.b. solaris seems unable to cope with this test.                         */
/*                                                                           */
#ifndef SOLARIS
#if (UINT_MAX != 4294967295U)
#error "WARNING !! You cannot compile this program, it requires 4 byte integers and this machine does not have them."
#endif
#endif


/* MESSAGE DEFINITION:                                                       */
/*                                                                           */
/* Messages are passed between client and server, each message has a header  */
/* and a body. The header contains fields that show whether the header needs */
/* swopping, what sort of message it is and the length of the body of the    */
/* message.                                                                  */

/* serverVersion in the message header is set by the server to this number   */
/* to indicate the level of the protocol used. We may or may not use this    */
/* to modify the protocol at a later date.                                   */
enum {ACESERV_CURRENT_SERVER_VERSION = 1} ;


/* The client sets swapMagic to OK_MAGIC
   if on the server side it appears as SWAP_MAGIC
   then the server will swap bytes in the headers
   in both directions, so the client never needs to swap
   */
#define OK_MAGIC   0x12345678	
#define SWAP_MAGIC 0x78563412	


/* We strings rather than enums, easier for others to interface to...        */
/* No ACESERV_MSGnnn string should be longer than                            */
/*                             (ACESERV_MSGTYPE_BUFLEN - 1)                  */
/* you must not lightly change this length, this is the size clients will    */
/* be expecting, to change it means changing all clients.                    */
/*                                                                           */
enum {ACESERV_MSGTYPE_BUFLEN = 30} ;

/* Client only, its messages are either normal requests (commands) or are    */
/* ace data that needs to be parsed in.                                      */
#define ACESERV_MSGREQ    "ACESERV_MSGREQ"
#define ACESERV_MSGDATA   "ACESERV_MSGDATA"

/* Server only, it may just be sending or a reply or it may be sending an    */
/* instruction, such as "operation refused".                                 */
#define ACESERV_MSGOK     "ACESERV_MSGOK"
#define ACESERV_MSGENCORE "ACESERV_MSGENCORE"
#define ACESERV_MSGFAIL   "ACESERV_MSGFAIL"
#define ACESERV_MSGKILL   "ACESERV_MSGKILL"



/* The message header struct, this has to be packed into the the first       */
/* ACE_HEADER_BYTES of a message sent between client & server.               */
/* ATTENTION if you modify this struct, modify also headerPack, headerUnpack */
/* and the two constants below.                                              */

/* Number of fields in the header that may need to be byte swopped.          */
enum {ACE_HEADER_SWAPFIELDS = 5} ;

/* Because we want a CONSTANT size for the buffer, we define this very       */
/* tediously...we have 6 ints of 4 bytes each and a char buffer to hold the  */
/* message type which is ACESERV_MSGTYPE_MAXLEN long.                        */
enum {ACE_HEADER_BYTES = ((4 * ACE_HEADER_SWAPFIELDS) + ACESERV_MSGTYPE_BUFLEN)} ;

typedef struct _ACE_HEADER
{
  int swapMagic ;					    /* Shows if header needs swapping */
  int length ;						    /* Length of data being sent in bytes. */
  int serverVersion ;					    /* not used ?? */
  int clientId ;					    /* Unique client identifier. */
  int maxBytes ;					    /* Maximum request size. */
  char msgType[ACESERV_MSGTYPE_BUFLEN] ;		    /* See the msgs defined above. */
} *ACE_HEADER ;





/* SERVER <-> TRANSPORT LAYER DEFINITION:                                    */
/*                                                                           */
/* The acedb server code communicates with the tranport layer via this       */
/* interface.                                                                */
/*                                                                           */

/* Opaque pointer to client connection details, needed here because acedb    */
/* server layer needs sometimes to call the transport layer and pass in      */
/* details of client connection.                                             */
/*                                                                           */
typedef struct _ClientConnectStruct *ClientConnect ;


/* Type defines an id for a connection, typically this will be an int for a  */
/* file descriptor. It needs to be something that can be used in a simple    */
/* comparision. The acedb layers needs to check (for security reasons) that  */
/* this connection id remains the same for a particular client id throughout */
/* the connection. Otherwise someone could start faking client ids in order  */
/* to gain admin access.                                                     */
typedef int ClientConnectionID ; 


/* These types are to implement communication between the socket layer and   */
/* the server. This communication all passes back and forth between          */
/* the two via a single function call/return. What this all says is that     */
/* each request has a target (either a single client or the whole server),   */
/* and it also has an operation which needs to be performed.                 */
/* Some operations are specific to a client or a socket/server.              */
/*                                                                           */

/* A request can be for action on a client or server.                        */
typedef enum _AceSocketRequestTarget
{
  ACESOCK_NONE, ACESOCK_SERVER, ACESOCK_CLIENT
} AceSocketRequestTarget ;

typedef enum _AceSocketRequestType
{
  /* Initial invalid request.                                                */
  ACESOCK_NULL,

  /* Server/socket specific                                                  */
  ACESOCK_SHUTQUIESCE, ACESOCK_SHUTFORCE, ACESOCK_EXIT,

  /* Client specific                                                         */
  ACESOCK_REQ, ACESOCK_KILL,

  /* common to both                                                          */
  ACESOCK_TIMEDOUT

} AceSocketRequestType ;



/* The definition of the acedb server callback routine, this is the routine  */
/* that the transport layer calls each time it receives a request from the   */
/* network.                                                                  */
typedef Stack (*AceSocketServerRequestRoutine)(void *serverData,
					       ClientConnect client,
					       ClientConnectionID connection_id,
					       AceSocketRequestTarget *type_inout,
					       AceSocketRequestType *req_inout,
					       ACE_HEADER ah, char *request) ;



/* When the server is first started, it will initialise and then call this   */
/* routine which never returns (like the X Windows style of program).        */
/* This routine in the transport layer loops listening for requests which it */
/* then passes on to the acedb server. The server passes in time outs and    */
/* a callback function pointer and a pointer to some data that the server    */
/* wants passed to its callback routine, again all very X like.              */
void aceSocketListen (int port, int serverTimeOut, int clientTimeOut,
		      AceSocketServerRequestRoutine processRequestFunc, void *serverData) ;



/* Were we started by inetd ?  For lots of reasons acedb needs to know this. */
/* If port_out is non-NULL then the listening port number will be returned.  */
BOOL isInetdDaemon(int *port_out) ;


/* If we are controlled by inetd and crash on startup without clearing the   */
/* listening socket, inetd will immediately try to restart us, which is      */
/* very bad news, this routine must clear up the listening socket for acedb  */
/* so that this does not happen.                                             */
void aceSocketInetdCleanup(void) ;


/* Returns the dotted decimal and human readable forms of the clients host   */
/* name. This is just a service routine really, but its all network stuff    */
/* so it goes here.                                                          */
BOOL getClientInetAddressInfo(ClientConnect connect_struct,
			      char **dotted_name, char **host_name) ;



/* SERVER <-> TRANSPORT MESSAGE DEFINITIONS:                                 */
/*                                                                           */
/* The server expects various messages during logging on, not strictly       */
/* necessary, but it all helps.                                              */
/*                                                                           */

/* On connection to the server, the client must send this text.              */
#define ACESERV_CLIENT_HELLO "bonjour"

/* The server will reply with this string when password verification is      */
/* complete.                                                                 */
#define ACESERV_SERVER_HELLO "et bonjour a vous"



/* SERVER <-> TRANSPORT GENERAL DEFINITIONS:                                 */
/*                                                                           */

/* Timeouts, setting a time out of zero means that client and/or server will */
/* never be timedout by the server. Timeouts must be a positive integer      */
/* number of seconds.                                                        */
enum {ACESERV_MIN_TIMEOUT = 0, ACESERV_INFINITE_TIMEOUT = 0} ;
#define ACESERV_INFINITE_TIMEOUT_STR "infinite"


#endif /* DEF_TRANSPORT_H */ 

    
