/*  File: sclientlib.c
 *  Author: Ed Griffiths (edgrif@sanger.ac.uk)
 *  Copyright (c) J Thierry-Mieg and R Durbin, 2000
 *-------------------------------------------------------------------
 * Acedb is free software; you can redistribute it and/or
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
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.crbm.cnrs-mop.fr
 *
 * Description: 
 * Exported functions: See XXXXXXXXXXXXX.h
 * HISTORY:
 * Last edited: Apr 27 14:05 2000 (edgrif)
 * Created: Wed Apr 26 16:21:44 2000 (edgrif)
 * CVS info:   $Id: sclientlib.c,v 1.1.1.1 2002/07/19 20:23:35 sienkiew Exp $
 *-------------------------------------------------------------------
 */


#include <wh/regular.h>
#include <wh/aceversion.h>


/* Needed for S_MESSAGE, poor structuring because client is not encapsulated */
/* properly...my fault...used to be encapsulated by askServerBinary()        */
#include <wsocket/acesocket_.h>

#include <wsocket/serverclientutils.h>
#include <wsocket/sclientlib.h>



/************************************************************/

static void handShake(S_MESSAGE msg,
		      char *userid, char *passwd_hash, char **nonce_hash,
		      int time_out) ;

/************************************************************/

enum {WANT_ENCORE = -1, DROP_ENCORE = -2} ;



static int staticSock = -1 ;


/************************************************************/
/************************************************************/

void *openServer (char *host, int port, int timeOut, int *sock_ret, Client client)
{ 
  int sock ;
  struct sockaddr_in sname ;
  S_MESSAGE msg ;
  struct hostent *hp ;

  /* init message system details.                                            */
  messSetMsgInfo("saceclient", aceGetVersionString(), getLogin(TRUE), getSystemName()) ;


  /* create a socket for writing/reading */
  sock = socket(AF_INET, SOCK_STREAM, 0) ;
  if (sock < 0)
    messcrash("cannot open write socket, %s", messSysErrorText()) ;

  /* create absolute socket name */
  bzero((char *)&sname, sizeof(sname)) ;
  sname.sin_family = AF_INET ;
  sname.sin_port = htons(port) ;

  hp = gethostbyname (host) ;
  if (!hp)
    messcrash("unknown host %s", host) ;

  bcopy ((char*)hp->h_addr, (char*) &sname.sin_addr, hp->h_length) ;

  /* OK, connect to the socket and start sending.                            */
  Connect(sock, (struct sockaddr *)&sname, sizeof(sname));

  msg = sMessageCreate() ;
  sMessageSetMagic(msg) ;

  staticSock = sock ;
  handShake(msg, client->userid, client->passwd_hash, &(client->nonce_hash),
	    timeOut) ;

  *sock_ret = sock ;
  return msg ;
}

/************************************************************/

void closeServer (void *handle)
{
  S_MESSAGE msg = (S_MESSAGE) handle ;

  if (staticSock >= 0)
    Shutdown(staticSock, SHUT_RDWR) ;
  sMessageDestroy (msg) ;
}


/************************************************************/

/* All the error handling here is rubbish, it never worked....               */
/*                                                                           */
/* For now I have changed messerrors to messcrash's so at least we don't     */
/* just carry blithely on...                                                 */
/*                                                                           */
int askServerBinary (void *handle, char *request, char **answerp, 
		     int *lengthp, int *encorep, int chunkSize, int time_out) 
{
  S_MESSAGE msg = (S_MESSAGE) handle ;
  int retval = 0 ;     /* positive number -> error */
  char *answer ;
  int length ;
  S_MSGState state ;
  fd_set rset ;
  int maxfd ;
  struct timeval tv ;
  int sock = staticSock ;
  int encore = 0 ;


  /* This whole routine needs redoing, its a hangover from the rpc stuff,    */
  /* it doesn't set up the message header properly etc. etc.                 */
  /* it doesn't error check etc. etc. uuuggghhh....                          */


  FD_ZERO(&rset) ;
  
  if (!strncasecmp(request,"encore",6)) 
    {
      /* encore request */
      encore = WANT_ENCORE;

    } 
  else if (!strncasecmp(request,"noencore",8)) 
    {
      /* encore request */
      encore = DROP_ENCORE;

    } 

  if (encore == 3)
    encore = -3 ;

  length = strlen(request) ;
  directReply (msg, request) ;

  sMessageSetWrite(msg) ;

  setMsgType(msg->ah.msgType, ACESERV_MSGREQ) ;

  /* note that this particular construct may block forever if we can;t write */
  /* for some reason...                                                      */
  state = SMSG_WAIT ;
  while(state == SMSG_WAIT)
    state = sMessageSocketWrite(sock, msg) ;
      
  if (state == SMSG_ERROR)
    {
      retval = 1 ;
      messcrash("Error in writing to socket, closing down connection.") ;
    }


  /* Block waiting for server to send us stuff...this implies that no notice */
  /* will taken of stdin....like EOF etc....                                 */
  FD_SET(sock, &rset) ;
  maxfd = getMaxSocket(fileno(stdin), sock) + 1 ;
  tv.tv_sec = 300 ;					    /* LINUX systems modify tv param. */
  tv.tv_usec = 0 ;
  aceSocketSelect(maxfd, &rset, NULL, NULL, &tv) ;
  

  if (FD_ISSET(sock, &rset))
    {
      msg->isNew = TRUE ;
      sMessageSetRead(msg) ;
      state = SMSG_WAIT ;

      /* Note we may end up blocked here if the server freezes...should add  */
      /* some time out code to include the select above which would do the   */
      /* timeout.                                                            */
      while(state == SMSG_WAIT)
	{
	  state = sMessageSocketRead(sock, msg) ;
	}
      if (state == SMSG_ERROR)
	{
	  retval = 1 ;
	  messcrash("Error in reading from socket, closing down connection.") ;
	}
    }


  *lengthp = msg->ah.length ;
  
  /* i do this to have exactly the same interface we had in rpc */
  if (msg->ah.length > 0 && msg->readBuffer)
    {    
      answer = (char *) malloc ( msg->ah.length) ;
      memcpy (answer,  msg->readBuffer, msg->ah.length) ;
    }
  else
    answer = 0 ;

  *answerp = answer ;
  messfree (msg->readBuffer) ;
  
  return retval ;
}




/************************************************************/

/* Send an initial message to the client, we will then get a reply which     */
/* contains a 'nonce' or key with which we must encode the users userid and  */
/* password. Then get the users userid/password and do the encryption and    */
/* send it back to the server for verification. We should then get another   */
/* reply from the server to OK this.                                         */
/* If this routine fails then we exit.                                       */
/*                                                                           */
static void handShake(S_MESSAGE msg, char *userid, char *passwd_hash,
		      char **nonce_hash_out, int time_out)
{
  char *reply = NULL ;
  int replyLength = 0, encore = 0 ;
  int rc ;

#ifdef ED_G_NEVER_INCLUDE_THIS_CODE
  BOOL got_passwd ;
  char *passwd = NULL ;
#endif /* ED_G_NEVER_INCLUDE_THIS_CODE */

  enum {HASH_STRING_NUM = 2} ;
  char *hash_strings[HASH_STRING_NUM] ;
  char *hash_user, *hash_nonce ;
  char *request ;
  char *words[2] ;


  /* This is rubbish, handshake should be part of normal looping for client  */
  /* but with a 'select' in the loop according to where we are for processing*/

  /* First contact with server, server should reply with a single word which */
  /* is encryption nonce. We should check for a single word here...          */
  rc = askServerBinary(msg, ACESERV_CLIENT_HELLO, &reply, &replyLength, &encore, 10000,
		       time_out) ;
  if (rc != 0 || !reply)
    messcrash("Client aborted, server did not reply correctly during handshake") ;


  /* Removing passwd prompt from here...                                     */
  hash_user = passwd_hash ;

#ifdef ED_G_NEVER_INCLUDE_THIS_CODE
  /* User may have supplied the userid/passwd hash as an argument to the     */
  /* program so no need to prompt user for the passwd.                       */
  if (*userid_out != NULL)
    {
      userid = *userid_out ;
      hash_user = *passwd_hash_out ;
    }
  else
    {
      /* Get the userid and passwd from the user.                                */
      got_passwd = getUseridPasswd(&userid, &passwd) ;
      if (!got_passwd)
	{
	  messExit("Client terminating, userid/passwd not entered correctly") ;
	}

      /* hash the userid and password together and convert into a hex    */
      /* string.                                                         */
      hash_strings[0] = userid ;
      hash_strings[1] = passwd ;
      hash_user = hashAndHexStrings(hash_strings, HASH_STRING_NUM) ;
    }
#endif /* ED_G_NEVER_INCLUDE_THIS_CODE */



  /* Now hash the hex userid/passwd and the nonce together and       */
  /* convert to a hex string.                                        */
  hash_strings[0] = hash_user ;
  hash_strings[1] = reply ;
  hash_nonce = hashAndHexStrings(hash_strings, HASH_STRING_NUM) ;

  /* Now make a string containing the userid and hash to send back.  */
  words[0] = userid, words[1] = hash_nonce ;
  request = makePhrase(words, 2) ;

  /* Now send back the encrypted password to the server.                     */
  rc = askServerBinary(msg, request, &reply, &replyLength, &encore, 10000,
		       time_out) ;
  if (rc != 0 && !reply)
    messcrash("Client aborted, server did not reply correctly to handshake") ;

  if (strcmp(reply, ACESERV_SERVER_HELLO) != 0)
    messExit("Client terminated, server did not reply correctly during handshake,"
	     "\n// it should have said: %s"
	     "\n//   but actually said: %s", ACESERV_SERVER_HELLO, reply) ;

  /* Return the userid and passwd hash for later reuse.                      */
  *nonce_hash_out = hash_nonce ;

  return ;
}



/************************************************************/

int getMaxSocket(int sock1, int sock2)
{
  int result ;

  result = (sock1 > sock2 ? sock1 : sock2) ;

  return result ;
}




/* Takes an array of strings (which are each assumed to be a single word with*/
/* no blanks) and makes a single blank separated string of all the words.    */
/*                                                                           */
char *makePhrase(char *words[], int num_words)
{
  char *phrase = NULL ;
  int phrase_len ;
  int i ;
  
  for (i = 0, phrase_len = 0 ; i < num_words ; i++)
    {
      phrase_len += (strlen(words[i]) + 1) ;
    }

  phrase = messalloc(phrase_len) ;

  strcpy(phrase, words[0]) ;
  for (i = 1 ; i < num_words ; i++)
    {
      strcat(phrase, " ") ;
      strcat(phrase, words[i]) ;
    }

  return phrase ;
}

