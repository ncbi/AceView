/*  File: saceclient.c
 *  Author: Ed Griffiths (edgrif@sanger.ac.uk)
 *          & Jean Thierry-Mieg (mieg@kaa.cnrs-mop.fr)
 *  Copyright (c) J Thierry-Mieg and R Durbin, 2000
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.crbm.cnrs-mop.fr
 *
 * Description: 
 * Exported functions: main
 * HISTORY:
 * Last edited: Feb  5 16:17 2001 (edgrif)
 * Created: Thu Jan  6 21:16:11 2000 (edgrif)
 * CVS info:   $Id: saceclient.c,v 1.2 2011/10/09 03:13:28 mieg Exp $
 *-------------------------------------------------------------------
 */
#include <wh/regular.h>
#include <wh/array.h>
#include <wh/aceio.h>
#include <wh/version.h>
#include <wh/aceversion.h>
#include <w4/command_.h>
#include <wsocket/acesocket_.h>
#include <wsocket/serverclientutils.h>
#include <wsocket/sclientlib.h>
#include <wsocket/saceclient_.h>





/* Client timeouts (must be specified in seconds), after the given times,    */
/* the client will exit.                                                     */
/* CLIENT_WAITUSER_TIMEOUT    - time we wait for user to do something        */
/* CLIENT_WAITNORMAL_TIMEOUT  - time we wait for server to respond to a      */
/*                              normal request                               */
/* CLIENT_WAITEXIT_TIMEOUT    - time we wait for server to respond to a      */
/*                              request to exit                              */
/*                                                                           */
enum 
{
  CLIENT_WAITUSER_TIMEOUT =   60 * 20,			    /* 20 mins */
  CLIENT_WAITNORMAL_TIMEOUT = 60 * 10,			    /* 10 mins */
  CLIENT_WAITEXIT_TIMEOUT =   60 *  3			    /* 3 mins */
} ;


typedef enum _ClientState {CLIENT_INIT,
			   CLIENT_WAITUSER, CLIENT_READUSER, CLIENT_WRITEUSER,
			   CLIENT_WAITSERVER, CLIENT_READSERVER, CLIENT_WRITESERVER,
			   CLIENT_WAITEXIT, CLIENT_EXIT} ClientState ;

#define SOCKET_ERROR_MSGPREFIX "\n// Client terminated abnormally - "
#define SOCKET_ERROR_MSG  "// This may have been caused by one of the following:\n" \
		          "// -  the database being shutdown by the administrator\n" \
		          "// -  the acedb server program timing out through lack of client requests\n" \
		          "// -  the acedb server program crashing\n" \
		          "// -  or the server machine or network going down.\n"


static void setQuitMessage(ClientState *client_state, char **cp, BOOL silent) ;
static char* getAceLevel (int level, int limit, BOOL mask) ;
static char* getAceData () ;
static void printPrompt(char *host) ;
static void printMsgPrefix(char *start) ;
static char* getOrder(int level, BOOL silent, 
		      BOOL isAceIn, BOOL isAceOut, BOOL isReport, char *host) ;
static void checkStdin(BOOL *parseInput_out, BOOL *sendRequest_out, BOOL silent,
		       ClientState *client_state_out, char **cp_out) ;
static void askUserForPasswd(char **userid_out, char **passwd_hash_out) ;
static char *makeCmdFromFile(char *filename) ;

BOOL debug_G = FALSE ;					    /* debug on/off for whole client. */


static int nAceIn = 0 ;


/*************************************************************************/

int main(int argc, char *argv[])
{
  void *handle; /* JC do not need to know what it is */
  char *host;
  char *answer;
  int n, length, t ;
  int user_timeout, server_timeout ;
  int level ;
  unsigned long port = DEFAULT_PORT;
  unsigned long p ;
  BOOL silent = FALSE, doWrite = TRUE , isReport = FALSE, encoring = FALSE,
       isAceIn = FALSE, isAceOut = FALSE ;
  char *cp = 0 ;
  FILE *fil = 0 ;
  Stack s = 0 ;
  char *start_msg = "" ;

  fd_set rset ;
  int maxfd ;
  struct timeval tv ;
  int sock ;
  Client client = NULL ;
  ClientState client_state = CLIENT_INIT ;
  char *msg_type ;
  int exit_status ;
  BOOL got_passwd = FALSE ;
  char *in_file = NULL, *out_file = NULL ;


  /* Does not return if user typed "-version" option.                        */
  checkForCmdlineVersion(&argc, argv) ;


  /* Set up our client state data.                                           */
  /* (if userid != NULL then we prompt for passwd)                           */
  client = (Client)messalloc(sizeof(ClientStruct)) ;
  client->userid = client->nonce = client->passwd_hash = client->nonce_hash = NULL ;


  freeinit();

  if ( argc < 2 ) 
    goto usage ;
  host = argv[1] ;

  n = 2 ;
  if ( argc > n + 1 && !strcmp("-port", argv[n]))
    { if (sscanf(argv[n+1],"%lu",&p)== 1)
      { port = p ;
      n += 2 ;
      }
    else
      goto usage ;
    }


  user_timeout = CLIENT_WAITUSER_TIMEOUT ;
  server_timeout = CLIENT_WAITNORMAL_TIMEOUT ;
  if ( argc > n + 1 && !strcmp("-time_out", argv[n]))
    {
      if (sscanf(argv[n+1],"%d",&t)== 1)
	{
	  user_timeout = t ;

	  cp = argv[n+1] ;	  
	  if ((cp = strstr(argv[n+1],":")) && sscanf(++cp,"%d",&t) == 1)
	    {
	      server_timeout = t ;
	    }
	  n += 2 ;
	}
      else
	goto usage ;
    }


  if (argc > n && (!strcmp("-passwd", argv[n])))
    {
      client->userid = argv[n+1] ; n++ ;
      client->passwd_hash = argv[n+1] ; n++ ;
      n++ ;

      /* some checking needed here to make sure we have two strings...       */
      /*                                                                     */
      /* Also, James wanted two cmd line flags, one for userid, one for      */
      /* the hash. This is probably best...                                  */

      got_passwd = TRUE ;
    }

  if (argc > n && (!strcmp("-ace_out", argv[n]) || !strcmp("ace_out", argv[n])))
    {
      isAceOut = TRUE ;
      silent = TRUE ;
      n++ ;
    }

  if (argc > n && (!strcmp("-ace_in", argv[n]) || !strcmp("ace_in", argv[n])))
    {
      isAceIn = TRUE ;

      /* silent = TRUE ; */
      fil = stdin ;
      n++ ;
    }
 
  /* Is there an input file to be read ??                                    */
  if (argc > n)
    { if (!strcmp("-f", argv[n]))
      { n++ ;
      silent = TRUE ;
      isReport = TRUE ;
      fil = fopen(argv[n],"r") ;
      if (!fil)
	{
	  messerror("open file %s given on command line \n", argv[n]) ;
	  fflush(stdout);
	  goto usage ;
	}
      }
    else
      goto usage ;
    }


  s = stackCreate(50) ;
  pushText(s,"") ;

  if (argc > n)
    { int i ;
    n++ ;
    for (i= n ; i < argc ; i++)
      { catText(s, argv[i]) ;
      catText(s, " ") ;
      }
    }

  /* Note, either we read from stdin, or from the input file, not both.      */
  if (!fil)
    fil = stdin ;

  /* Set up the file to be read in.                                          */
  level = freesetfile(fil, stackText(s,0)) ;
  if (isReport)
    freespecial("\n\t\"\\@%$") ;				    /* no "//" (to echo empty lines), */
							    /* allow client-side sub shells */
  else
    freespecial("\n\t\"\\/@%$") ;				    /* allow client-side sub shells */



  /* If we haven't go the userid/passwd then get them....                    */
  /*                                                                         */
  /* NOTE I need to add stuff to get them as cmd line args, we can use       */
  /* glib stuff I guess to do all this now...                                */
  /*                                                                         */
  if (got_passwd == FALSE)
    {
      /* We exit if this fails.                                              */
      askUserForPasswd(&(client->userid), &(client->passwd_hash)) ;
    }


  /* Make contact with the server.                                           */
  if ((handle = openServer(host, port, CLIENT_WAITNORMAL_TIMEOUT, &sock, client)) == NULL)
    {
      messerror("cannot establish connection\n") ; 
      goto abort ;
    }
  else if (debug_G)
    messout("opened connection to %s on port %lu\n\n",host,port);


  /* Mucky having this here...but how else to get the first prompt...        */
  if (!silent && !isAceIn)
    {
      printPrompt(host) ;
      start_msg = "\n" ;
    }

  exit_status = EXIT_SUCCESS ;
  client_state = CLIENT_WAITUSER ;
  while(client_state != CLIENT_EXIT)
    {
      int sock_rc ;
      int timeout ;

      /* Set up the file descriptors.                                        */
      /* Always monitor socket, don't monitor stdin or input file if client  */
      /* is quitting.                                                        */
      FD_ZERO(&rset) ;
      FD_SET(sock, &rset) ;
      if (client_state != CLIENT_WAITEXIT)
	{
	  FD_SET(fileno(fil), &rset) ;
	  maxfd = getMaxSocket(fileno(fil), sock) + 1 ;
	}
      else
	maxfd = sock + 1 ;

      /* Set appropriate timeout.                                            */
      if (client_state == CLIENT_WAITUSER)
	timeout = user_timeout ;
      else if (client_state == CLIENT_WAITSERVER)
	timeout = server_timeout ;
      else if (client_state == CLIENT_WAITEXIT)
	timeout = CLIENT_WAITEXIT_TIMEOUT ;


      /* Block waiting for stdin or server input...                          */
      tv.tv_sec = timeout ;				    /* Reset every time, LINUX alters */
      tv.tv_usec = 0 ;					    /* them... */
      sock_rc = aceSocketSelect(maxfd, &rset, NULL, NULL, &tv) ;

#ifdef ED_G_NEVER_INCLUDE_THIS_CODE
      /* We need to sort out writing to socket like we have reading....      */
#endif /* ED_G_NEVER_INCLUDE_THIS_CODE */

      /* Either we've timed out or there is something to read from socket or */
      /* stdin/input file.                                                   */
      if (sock_rc == 0)
	{
	  /* OK, we have timed out, so exit. Currently we just bomb out but  */
	  /* we could send the server a "quit" if we timed out waiting for   */
	  /* the user to do something.                                       */
	  char *message ;

	  if (client_state  == CLIENT_WAITUSER)
	    message = "You have not entered a command for a long time." ;
	  else if (client_state == CLIENT_WAITSERVER)
	    message = "The server has not responded to your last request." ;
	  else if (client_state == CLIENT_WAITEXIT)
	    message = "The server has not responded to your request to quit." ;

	  printMsgPrefix(start_msg) ;
	  messerror(SOCKET_ERROR_MSGPREFIX "you have timed out\n"
		    "// because there has been no client activity for %d secs.\n"
		    "// Cause was:\n"
		    "//             %s\n", timeout, message) ;

	  client_state = CLIENT_EXIT ;
	}
      else if (FD_ISSET(fileno(fil), &rset) && client_state == CLIENT_WAITUSER)
        {
	  /* Either stdin needs to be read or the user specified an input    */
	  /* file and that is ready for reading.                             */
	  S_MSGState state ;
	  BOOL parseInput = TRUE ;			    /* Data to be parsed ? */
	  BOOL sendRequest = TRUE ;			    /* Anything to send to server ? */

	  if (debug_G)
	    messout("reading from stdin\n") ;

	  if (client_state == CLIENT_WAITSERVER)
	    {
	      /* Ignore user input if we are waiting for the server.         */
	      /* We could allow quitting and Cntl-D here at some time.       */
	    }
	  else
	    {
	      /* If reading from stdin, check what's there for specials,     */
	      /* e.g. Cntl-D etc.                                            */
	      cp = NULL ;
	      if (fil == stdin)
		{
		  checkStdin(&parseInput, &sendRequest, silent, &client_state, &cp) ;
		  start_msg = "" ;
		}


	      /* Anything to be parsed from the input ?                          */
	      if (parseInput)
		{
		  enum {MAX_REQUEST_WORDS = 10} ;
		  char *words[MAX_REQUEST_WORDS] ;	    /* for construction of requests. */

		  client_state = CLIENT_WAITSERVER ;
		  doWrite = TRUE ;


		  cp = getOrder(level, silent, isAceIn, isAceOut, isReport, host) ;

		  /* This lot should go in a subroutine...                   */
		  if (cp == NULL)			    /* EOF => quit */
		    {
		      setQuitMessage(&client_state, &cp, silent) ;
		    }
		  else if (strcmp(cp, "quit") == 0)	    /* quit command. */
		    {
		      setQuitMessage(&client_state, &cp, silent) ;
		    }
		  else if (!strncasecmp(cp, "shutdown", 8))
		    {
		      client_state = CLIENT_WAITEXIT ;
		      messout("Client sending shutdown to server") ;
		    }
		  else if (!strncasecmp(cp, "passwd", 6))
		    {
		      char *new_hash ;

		      new_hash = getNewPasswd(client->userid, client->passwd_hash) ;
		      if (!new_hash)
			{
			  sendRequest = FALSE ;
			  messout("Sorry, password change failed, please check current password\n"
				  "// and make sure you correctly enter new password.") ;
			  client_state = CLIENT_WAITUSER ;
			}
		      else
			{
			  words[0] = "passwd", words[1] = client->userid,
			    words[2] = client->nonce_hash, words[3] = new_hash ;

			  cp = makePhrase(words, 4) ;

			  /* This is dodgy, really we need to wait until the     */
			  /* command has completed before doing this update. We  */
			  /* should keep state and then assign....               */
			  /* We could now do this using the WAIT_SERVER state*/
			  /* etc. we would need a flag to say that we were   */
			  /* waiting to update the password.                 */
			  client->passwd_hash = new_hash ;
			}
		    }
		  else if (!strncasecmp(cp, "user", 4))
		    {
		      BOOL status ;

		      status = getUserUpdate(&cp) ;
		      if (!status)
			{
			  sendRequest = FALSE ;
			  messout("Sorry, incorrect arguments to \"user\" command.\n"
				  "// if you entered a password, "
				  "try again and make sure you re-enter it correctly.") ;
			  /* We should output some usage stuff here.         */
			  client_state = CLIENT_WAITUSER ;
			}
		    }
		  else if (!strncasecmp(cp, "domain", 4))
		    {
		      BOOL status ;

		      status = getDomainUpdate(&cp) ;
		      if (!status)
			{
			  sendRequest = FALSE ;
			  messout("Sorry, incorrect arguments to \"domain\" command.\n"
				  "// please check syntax and try again.") ;
			  client_state = CLIENT_WAITUSER ;
			  /* We should output some usage stuff here.         */
			}
		    }
		  else if (!strncasecmp(cp, "parse", 5) && strNumWords(cp) == 1)
		    {
		      /* If someone just enters "parse" then get the ace data from stdin. */
		      cp = getAceData("parse> ") ;
		      if (!cp)
			sendRequest = FALSE ;
		    }
		  else if (isAceOut && strncasecmp(cp,"Write",5) && strncasecmp (cp, "Table", 5))
		    {
		      doWrite = FALSE ;
		    }
		}

	      /* Request to send ??                                              */
	      if (sendRequest)
		{
		  CmdFilespecArg status ;
		  char *new_cp ;

#ifdef ED_G_NEVER_INCLUDE_THIS_CODE
		  /* Now check the cmd for file args...                      */
		  if (debug_G) messout("users command: %s", cp) ;

		  /* THERE IS A BUG IN THE BELOW ROUTINE, IT DOES NOT HANDLE */
		  /* ACE_IN PROPERLY...returns an error....                  */
		  status = aceCommandCheckForFiles(cp, &new_cp, &in_file, &out_file) ;
		  if (status == CMD_FILE_ERR)
		    {
		      messout("Sorry, you specified the command or one of its arguments incorrectly.\n"
			      "// Type \"?\" for help.") ;
		      client_state = CLIENT_WAITUSER ;
		    }
		  else if(status == CMD_FILE_FOUND)
		    {
		      char *file_text ;

		      if (debug_G) messout("  new command: %s,  in_file: %s  out_file: %s",
					   new_cp,
					   (in_file ? in_file : "no in_file"),
					   (out_file ? out_file : "no out_file")) ;
		      
		      /* if theres an infile we need to read it into a an   */
		      /* acein here...then we need to retrieve the text     */
		      /* and shove it on to the stack as "-c <text> "       */
		      /* then we can send it to the server which will grab  */
		      /* it.                                                */
		      if (in_file)
			file_text = makeCmdFromFile(in_file) ;


#ifdef ED_G_NEVER_INCLUDE_THIS_CODE
		      cp = new_cp ;
#endif /* ED_G_NEVER_INCLUDE_THIS_CODE */



		    }
		  else
		    {
		       if (debug_G) messout("Command does not use files...") ;
		    }
#endif /* ED_G_NEVER_INCLUDE_THIS_CODE */



		  /* For testing purposes...                                 */
		  status = CMD_FILE_NONE ;


		  if (status != CMD_FILE_ERR)
		    {
		      encoring = FALSE ;
		      if (isAceIn)
			msg_type = ACESERV_MSGDATA ;
		      else
			msg_type = ACESERV_MSGREQ ;

		      state = writeToSocket(sock, handle, msg_type, cp) ;
		      if (state == SMSG_ERROR)
			{
			  client_state = CLIENT_EXIT ;
			  printMsgPrefix(start_msg) ;
			  messerror(SOCKET_ERROR_MSGPREFIX "unable to write to socket") ;
			}


		      /* ACTUALLY I THINK THE STATES ARE ALL MESSED UP HERE, IF  */
		      /* ITS ACEIN WE MUST NOT READ ANYMORE FROM THE INPUT FILE  */
		      /* BECAUSE IT GETS CLOSED BY THE INITIAL READ. WE SHOULD   */
		      /* USE THIS FACT TO ADJUST OUR STATE I THINK THIS LOT IS   */
		      /* PROBABLY IN THE WRONG PLACE...AND THAT PERHAPS WE SHOULD*/
		      /* HAVE A SEPARATE FLAG FOR NOT READING THE INPUT FILE...  */
		      /*                                                         */
		      /* PERHAPS THIS WOULD HELP TO MAKE THINGS MORE NATURAL...  */
		      /* THINK ABOUT IT...                                       */
		      
		      /* If we were reading from an input file, we 'quit' as soon*/
		      /* as we've sent the data to the server.                   */
		      if (isAceIn)
			{
			  setQuitMessage(&client_state, &cp, silent) ;
			  
			  msg_type = ACESERV_MSGREQ ;
			  state = writeToSocket(sock, handle, msg_type, cp) ;
			  if (state == SMSG_ERROR)
			    {
			      client_state = CLIENT_EXIT ;
			      printMsgPrefix(start_msg) ;
			      messerror(SOCKET_ERROR_MSGPREFIX "unable to write to socket") ;
			    }
			}
		    }
		} /* sendRequest */

	      /* Don't try to read from stdin or any input file if its time to   */
	      /* quit.                                                           */
	      if (client_state == CLIENT_WAITEXIT || client_state == CLIENT_EXIT)
		{
		  FD_CLR(fileno(fil), &rset) ;
		}
 
	    }
	}
      else if (FD_ISSET(sock, &rset))
	{
	  /* There's something on the socket from the server.                */
	  S_MSGState state ;

	  if (debug_G)
	    messout("reading from socket\n") ;

	  state = readFromSocket(sock, handle, &msg_type, &answer, &length) ;
	  if (state == SMSG_ERROR)
	    {
	      /* Problem with the socket.                                    */
	      char *msg ;

	      if (client_state == CLIENT_WAITEXIT)
		msg = "error in reading from socket while waiting for quit from server." ;
	      else
		msg = "error in reading from socket." ;

	      client_state = CLIENT_EXIT ;
	      printMsgPrefix(start_msg) ;
	      messerror(SOCKET_ERROR_MSGPREFIX "%s\n" SOCKET_ERROR_MSG, msg) ;
	    }
	  else if (state == SMSG_WAIT)
	    {
	      /* Nothing to do except wait for rest of message to come.      */
	    }
	  else /* state == SMSG_DONE */
	    {
	      /* Must mean that socket replied OK.                           */
	      if (testMsgType(msg_type, ACESERV_MSGKILL))
		{
		  client_state = CLIENT_EXIT ;
		  if (!silent)
		      messout("Client sent termination signal by server.") ;
		}
	      else
		{
		  /* Because we sent a quit straight away to the server for  */
		  /* isAceIn we now need to make sure we quit here.          */
		  if (isAceIn)
		    client_state = CLIENT_WAITEXIT ;
		  else if (isReport && testMsgType(msg_type, ACESERV_MSGFAIL))
		    {
		      /* We exit here because we are processing an input file*/
		      /* and the user won't know about any error unless we   */
		      /* report it and exit.                                 */
		      client_state = CLIENT_EXIT ;
		      exit_status = EXIT_FAILURE ;
		      doWrite = TRUE ;
		      silent = FALSE ;

		      /* Output a message whatever happens so that user      */
		      /* gets to see something.                              */
		      messerror("server has returned an error probably caused by the input file.") ;
		    }
  		  else
		    client_state = CLIENT_WAITUSER ;
		}

	      if (length && !encoring && !silent)
		messout ("Response: %d bytes.", length) ;


	      s = stackReCreate(s,length) ;
	      stackTextForceFeed(s, length) ;
	      memcpy(stackText(s,0), answer, length) ;
	      free(answer);

	      stackCursor(s, 0) ;
		  
	      if (doWrite)
		{
		  while ((cp = stackNextText(s)))
		    {

		      if ((silent || testMsgType(msg_type, ACESERV_MSGENCORE))
			  && (!*cp || *cp == '/'))
			continue ;
			  
		      if (silent)
			{
			  /* Don't print comments (// lines)                 */
			  register char *cp1 = cp, *cp2, cc ;
			  while (*cp1)
			    {
			      cp2 = cp1 ;
			      while (*cp2 && *cp2 != '\n') cp2++ ;
			      cc = *cp2 ; *cp2 = 0 ;
			      if (!(*cp1 == '/' && *(cp1 + 1) == '/'))
				puts (cp1) ;
			      if (cc) *cp2++ = cc ;
			      cp1 = cp2 ;
			    }
			}
		      else
			{ 
			  /* print everything.                               */
			  if (*cp)
			    { register char *cp1 = cp + strlen(cp) - 1 ;
			    if (cp1 >= cp && *cp1 == '\n')
			      *cp1 = 0 ; /* because puts adds a terminal \n */
			    if (*cp) puts (cp) ; else putchar('\n') ;
			    }
			}
		    }
		}
	       

	      if (isAceIn)
		continue ;
	      else if (testMsgType(msg_type, ACESERV_MSGENCORE))
		{
		  encoring = TRUE ;

		  /* Try deleting this altogether...seems redundant...       */
		  /* but...good for debugging...                             */
		  cp = "encore" ;

		  /*if encoring we should write to the socket here...and then wait for more... */
		  state = writeToSocket(sock, handle, msg_type, cp) ;
		  /* error handling wally...                                 */

		  client_state = CLIENT_WAITSERVER ;
		}
	      else
		encoring = FALSE ;
	    }

	} /* socket ready */


      if (client_state != CLIENT_WAITSERVER && client_state != CLIENT_EXIT
	  && client_state != CLIENT_WAITEXIT && !silent && !encoring)
	{
	  printPrompt(host) ;
	  start_msg = "\n" ;
	}
    }

  if (isAceIn)
    messout ("Passed %d objects to the server.\n", nAceIn) ;

  closeServer(handle);

  if (!silent)
    messout("Please report problems to acedb@sanger.ac.uk\n// Bye\n") ;

  return(exit_status) ;

 usage:
  messerror("usage: saceclient host [-port number] [-time_out secs] "
	    "[-access_info] [-ace_out] [-ace_in] [-f reportfile parameters]\n") ;
 abort:
  return(EXIT_FAILURE) ;
}


/***************************************************************************/

/* Defines a routine to return the compile date of this file.                */
UT_MAKE_GETCOMPILEDATEROUTINE()


/***************************************************************************/

/* Set up state/message that we send to server to say that client wants to   */
/* quit.                                                                     */
static void setQuitMessage(ClientState *client_state, char **cp, BOOL silent)
{
  if (!silent)
    messout("Closing connection to server.") ;
  *client_state = CLIENT_WAITEXIT ;
  *cp = "quit" ;

  return ;
}


/*************************************************************************/

/* Check stdin for the few cases that are not some sort of command to acedb  */
/* but instead are the user just pressing 'enter' or whatever.               */
/* Note that is we get an error we do not try to exit cleanly by telling the */
/* server, its too risky, we don't know what state we are in.                */
/*                                                                           */
static void checkStdin(BOOL *parseInput_out, BOOL *sendRequest_out, BOOL silent,
		       ClientState *client_state_out, char **cp_out)
{
  int file_rc ;

  *parseInput_out = FALSE ;

  file_rc = getc(stdin) ;
  if (file_rc < 0)
    {
      /* There is some problem with stdin.                                   */
      if (ferror(stdin))
	  messcrash("unable to read stdin, client is terminating\n") ;
      else if (feof(stdin))				    /* user pressed 'EOF', (Cntl-D) */
	{
	  /* Prepare normal quit message to send to server.                  */
	  setQuitMessage(client_state_out, cp_out, silent) ;
	}
    }
  else if (file_rc == '\n')
    {
      /* User just pressed 'enter' with no data.                             */
      *sendRequest_out = FALSE ;
    }
  else
    {
      /* User entered some data.                                             */
      if (ungetc(file_rc, stdin) == EOF)
	  messcrash("unable to reset stdin, client is terminating\n") ;
      else
	{
	  *parseInput_out = TRUE ;
	  *sendRequest_out = TRUE ;
	}
    }

  return ;
}


/*************************************************************************/

static char* getAceLevel (int level, int limit, BOOL mask)
{ 
  static Stack s = 0 ;
  BOOL in = FALSE ;
  char *cp, cutter ;
  int nn = 1 ;
 
  s = stackReCreate (s, 30000) ;

  if (mask)
    pushText (s, "parse = ") ;

 while (freecard (level))
   {
     cp = freewordcut("\n" ,&cutter) ;
     if (cp && *cp)
       { if (!in)
	 { nn++ ;
	 nAceIn++ ; in = TRUE ;
	 }
       catText (s, cp) ;
       if (mask) catText (s, "\\n") ;
       else catText (s, "\n") ;
       }
     else
       {
	 if (mask)
	   catText (s, "\\n") ;
	 else
	   catText (s, "\n") ;
	 in = FALSE ;
	 if (limit && !(nn % limit))
	   break ;
       }
   }

 return nn > 1 ? stackText (s, 0) : 0 ;
}

/*************************************************************************/

static char* getAceData ()
{ 
  int level = freesetfile (stdin,"") ;
  char *cp ;

  freespecial ("\n/\\\"\t@") ;
  
  cp = getAceLevel (level, 0, TRUE) ;
  freeclose (level) ;
  return cp ;
}

/*************************************************************************/

/* A pair of noddy routines to attempt to make the prompt and messages       */
/* appear on the correct lines...                                            */
/*                                                                           */
/* In theory we should use readline/acein etc. but this doesn't work at the  */
/* moment....                                                                */
static void printPrompt(char *host)
{
  printf("acedb@%s> ", host) ; fflush(stdout);

  return ;
}

static void printMsgPrefix(char *start)
{
  printf(start) ;

  return ;
}

/*************************************************************************/

static char* getOrder(int level, BOOL silent, 
		      BOOL isAceIn, BOOL isAceOut, BOOL isReport, char *host)
{ char *cp = 0, cutter, *errtext;
  int np ;
  static Stack s = 0 ;

  /* It's a an ace file, so just read the lot in one go (and then it will be */
  /* sent in one go to the server).                                          */
  if (isAceIn)
    return getAceLevel (level, 500, FALSE) ;


  /* Ok, input must be a command or series of commands from a file or stdin  */
  /* (remember this could be piped input).                                   */
 lao:

  cutter = 0 ;   
  cp = freewordcut(isReport ? "#" : "" ,&cutter) ;

  if ((!cp && !cutter) || (cp && *cp == '/' && *(cp+1) == '/'))
    {
      if (!level || !freecard (level))
	return 0 ;                       /* just get a card */

      if (isReport && (!cp || *cp != '/'))
	{
	  putchar('\n'), fflush(stdout) ;
	}

      goto lao ;
    }

  if (!silent || isAceOut)
    return cp ;   /* the whole line will be processed by the server */
  
    /* Now silent = TRUE and we are in report mode*/

  if (!cutter)    /* just print out the line */
    {
      if (cp)
	{ messout("%s",cp) ;fflush(stdout); }
      goto lao ;
    }
  if (cp)
    { messout("%s",cp) ; /* print the beginning */
      fflush(stdout);
    }
  s = stackReCreate(s, 50) ;
  pushText(s,"") ;
  
  if (!freestep('('))
    { putchar (cutter) ; goto lao ; }

  np = 1 ;
  while (np && (cp = freewordcut("()" ,&cutter)))
    { catText(s, cp) ;
      switch (cutter)
	{ 
	case '(':
	  catText(s,"(") ;
	  np++ ;
	  break ;
	case ')':
	  np-- ;
	  if (np)
	    catText(s,")") ; 
	  break ;
	case 0:
	  if (np)
	    { errtext = "Missing right parenthese" ;
	      goto error ;
	    }
	}
    }
  cp = stackText(s,0) ;

  return cp ;  /* let the server process this command */


 error:
  messerror("%s after %% in %s\n", errtext, stackText(s,0)) ;
  cp = freepos() ;
  if (cp)
    messerror("%s",cp) ;
  fflush(stdout);
  goto lao ;

}

/*************************************************************************/

/* Either works or we exit.                                                  */
/*                                                                           */
static void askUserForPasswd(char **userid_out, char **passwd_hash_out)
{
  char *userid = NULL ;
  char *passwd = NULL ;
  char *hash_user = NULL ;
  enum {HASH_STRING_NUM = 2} ;
  char *hash_strings[HASH_STRING_NUM] ;

  if (!getUseridPasswd(&userid, &passwd))
    messExit("Client terminating, userid/passwd not entered correctly") ;
  
  /* hash the userid and password together and convert into a hex    */
  /* string.                                                         */
  hash_strings[0] = userid ;
  hash_strings[1] = passwd ;
  hash_user = hashAndHexStrings(hash_strings, HASH_STRING_NUM) ;

  *userid_out = userid ;
  *passwd_hash_out = hash_user ;

  return ;
}


/*****************************************************************************/

/* Get the contents of a text file and put them in a C string which has the  */
/* form " -c <text from file> ", i.e. we are creating a command line arg     */
/* from a file.                                                              */
/*                                                                           */
static char *makeCmdFromFile(char *filename)
{
  char *result = NULL ;
  ACEIN infile ;
  Stack file_text ;
  char *next_line ;

  infile = aceInCreateFromFile(filename, "r", "", 0) ;

  file_text = stackCreate(2000) ;
  stackTextOnly(file_text) ;

  pushText(file_text, "") ;
  while ((next_line = aceInCard(infile)) != NULL)
    {
      catText(file_text, next_line) ;
    }

  result = strnew(stackText(file_text, 0), 0) ;
  
  aceInDestroy(infile) ;
  stackDestroy(file_text) ;

  return result ;
}



/*************************************************************************/
/*************************************************************************/
 
 
