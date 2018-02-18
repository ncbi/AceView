/*  File: sclient.c
 *  Author: Ed Griffiths (edgrif@sanger.ac.uk)
 *  Copyright (c) J Thierry-Mieg and R Durbin, 1999
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.crbm.cnrs-mop.fr
 *
 * Description: A dummy client for testing the socket server.
 * HISTORY:
 * Last edited: Jan  4 12:04 2001 (edgrif)
 * Created: Thu Jul 22 14:12:15 1999 (edgrif)
 * CVS info:   $Id: sclient.c,v 1.1.1.1 2002/07/19 20:23:35 sienkiew Exp $
 *-------------------------------------------------------------------
 */
#include <wh/mystdlib.h>
#include <termios.h>					    /* for non-echoing terminal code. */
#include <readline/readline.h>				    /* for getting user input. */

#include <wh/regular.h>


#ifdef ED_G_NEVER_INCLUDE_THIS_CODE
#include <wh/array.h>
#include <wsocket/acesocket_.h>
#endif /* ED_G_NEVER_INCLUDE_THIS_CODE */

#include <wsocket/serverclientutils.h>
#include <wsocket/saceclient_.h>


/************************************************************/

static void promptUser(int num_prompts, char *prompts[], char *answers[], int passwd_prompt) ;



/************************************************************/
/************************************************************/
static char *myreadline(char *prompt)
{
  int BUFSIZE = 2000 ;
  static char *s = 0 ;

  if (!s) s = messalloc(BUFSIZE) ;

  printf ("%s", prompt) ;
  s = fgets(s,BUFSIZE,stdin) ;
  return s ;
}

/* Get the userid and password from the user.                                */
/*                                                                           */
/* Try as I might I can't see how to get the readline package _not_ to echo  */
/* the password back to the screen, very annoying, none of the documentation */
/* mentions this very basic action. So instead I fiddle with the terminal    */
/* characteristics directly.                                                 */
/*                                                                           */
BOOL getUseridPasswd(char **userid, char **passwd)
{
  BOOL got_passwd = FALSE ;
  char *useridPrompt = "Please enter userid: " ;
  char *passwdPrompt = "Please enter passwd: " ;
  struct termios term ;

  /* once I know what to do with readline to do non-echoing etc. I will try  */
  /* doing this using acein                                                 */


  /* OK, get the password.                                                   */
  *userid = myreadline(useridPrompt) ;


  /* Now set the terminal non-echoing, get the password and reset the        */
  /* terminal.                                                               */
  if (tcgetattr(fileno(stdin), &term) < 0)
    messcrash("Unable to get terminal attributes for input of password.");

  term.c_lflag &= ~(ECHO) ;				    /* Turn echo off. */
  if (tcsetattr(fileno(stdin), TCSADRAIN, &term) < 0)
    messcrash("Unable to set terminal attributes for input of password.");

  *passwd = myreadline(passwdPrompt) ;

  term.c_lflag |= ECHO ;				    /* Turn on again. */
  if (tcsetattr(fileno(stdin), TCSADRAIN, &term) < 0)
    messcrash("Unable to reset terminal attributes after input of password.");

  printf("\n") ;					    /* newline after password is not
							       echoed to terminal so send one. */

  if (*userid && *passwd)
    got_passwd = TRUE ;

  return got_passwd ;
}




/************************************************************/

char *getNewPasswd(char *userid, char *curr_passwd_hash)
{
  char *passwd_str = NULL ;
  enum {PROMPT_NUM = 3, PASSWD_PROMPT = 0} ;
  char *prompts[PROMPT_NUM] = {"Please enter current passwd: ",
			       "Please enter new passwd: ", 
			       "Please re-enter new passwd: "} ;
  char *answers[PROMPT_NUM] ;
  BOOL status = TRUE ;
  enum {HASH_STRING_NUM = 2} ;
  char *hash_strings[HASH_STRING_NUM] ;
  char *passwd_hash ;


  /* Now get the current password and the new password from the user doing   */
  /* all the non-echoing stuff...                                            */
  /*                                                                         */
  messout("Changing password for %s", userid) ;
  promptUser(PROMPT_NUM, prompts, answers, PASSWD_PROMPT) ;

  /* Check that current password matches the one person logged on with.      */
  hash_strings[0] = userid ;
  hash_strings[1] = answers[0] ;
  passwd_hash = hashAndHexStrings(hash_strings, HASH_STRING_NUM) ;
  if (strcmp(curr_passwd_hash, passwd_hash) != 0)
    status = FALSE ;
  messfree(passwd_hash) ;				    /* no longer needed. */

  /* Check that new password was reentered as the same string                */
  if (status)
    {
      if (strcmp(answers[1], answers[2]) != 0)
	status = FALSE ;
    }

  /* hash the new password                                                   */
  if (status)
    {
      hash_strings[0] = userid ;
      hash_strings[1] = answers[1] ;
      passwd_str = hashAndHexStrings(hash_strings, HASH_STRING_NUM) ;
    }

  return passwd_str ;
}


/************************************************************/

BOOL getUserUpdate(char **cp)
{
  BOOL status = TRUE ;
  char *passwd_str = NULL ;
  enum {MAX_PROMPTS = 4, USERID_PROMPT = 0} ;
  char *prompts[MAX_PROMPTS] = {"Please enter userid to be updated: "} ;
  char *groupnew = "Please enter name of users new group: " ;
  char *group = "Please enter name of users group: " ;
  char *passwd1 = "Please enter users new passwd: " ;
  char *passwd2 = "Please re-enter users new passwd: " ;
  char *answers[MAX_PROMPTS] ;
  int num_prompts, passwd_prompt = -1 ;
  enum {HASH_STRING_NUM = 2} ;
  char *hash_strings[HASH_STRING_NUM] ;


  /* NOTE THAT THIS IS NOT SECURE AT THE MOMENT, REALLY WE SHOULD ASK THE    */
  /* ADMIN USER FOR THEIR OWN PASSWORD EVERY TIME THEY WANT TO UPDATE ANY    */
  /* ASPECT OF THE PASSWORD FILE AT ALL...OTHERWISE ALL SORTS OF TRICKS      */
  /* CAN BE PLAYED, SUCH AS MOVING A DIFFERENT USER INTO THE ADMIN GROUP OR  */
  /* CHANGING THE ADMIN PASSWORD OR SETTING THE DATABASE TO WORLD WRITEABLE. */
  /*                                                                         */
  /* BUT, if we make it secure, this will become a pain for the admin person */
  /* because they will have to enter their password  all the time and        */
  /* actually I think root users on UNIX are not asked their password when   */
  /* they change someone elses...                                            */


  /* Now get the current password and the new password from the user doing   */
  /* all the non-echoing stuff...                                            */
  /*                                                                         */
  messout("Updating user...") ;

  if (strstr(*cp, "passwd"))
    {
      prompts[1] = passwd1 ;
      prompts[2] = passwd2 ;
      num_prompts = 3 ;
      passwd_prompt = 1 ;
    }
  else if (strstr(*cp, "group"))
    {
      prompts[1] = groupnew ;
      num_prompts = 2 ;
    }
  else if (strstr(*cp, "new"))
    {
      prompts[1] = group ;
      prompts[2] = passwd1 ;
      prompts[3] = passwd2 ;
      num_prompts = 4 ;
      passwd_prompt = 2 ;
    }
  else if (strstr(*cp, "delete"))
    {
      /* No extra prompts.                                                   */
      num_prompts = 1 ;
    }
  else
    {
      status = FALSE ;
    }

  promptUser(num_prompts, prompts, answers, passwd_prompt) ;

#ifdef ED_G_NEVER_INCLUDE_THIS_CODE
  /* THIS IS DISABLED CURRENTLY, SEE COMMENTS ABOVE....                      */

  /* Check that current password matches the one person logged on with.      */
  hash_strings[0] = userid ;
  hash_strings[1] = answers[0] ;
  passwd_hash = hashAndHexStrings(hash_strings, HASH_STRING_NUM) ;
  if (strcmp(curr_passwd_hash, passwd_hash) != 0)
    status = FALSE ;
  free(passwd_hash) ;					    /* no longer needed. */
#endif /* ED_G_NEVER_INCLUDE_THIS_CODE */


  /* For some operations we need to hash a userid & password.                */
  if (status && passwd_prompt >= 0)
    {
      /* Check that new password was reentered as the same string                */
      if (status)
	{
	  if (strcmp(answers[passwd_prompt], answers[passwd_prompt + 1]) != 0)
	    status = FALSE ;
	}

      /* hash the new password                                                   */
      if (status)
	{
	  hash_strings[0] = answers[USERID_PROMPT] ;
	  hash_strings[1] = answers[passwd_prompt] ;
	  passwd_str = hashAndHexStrings(hash_strings, HASH_STRING_NUM) ;
	}
    }

  /* Now construct the request to the server.                                */
  if (status)
    {
      int i, length ;
      char *tmp ;
      int p = (passwd_prompt < 0) ? num_prompts : passwd_prompt ;

      length = strlen(*cp) ;
      for (i = 0 ; i < num_prompts ; i++)
	{
	  length += (1 + strlen(answers[i])) ;
	}
      length += 1 ;
      tmp = messalloc(length) ;
      
      if (!strcpy(tmp, *cp))
	messcrash("string operation failed while copying initial user request.") ;
      for (i = 0 ; i < p ; i++)
	{
	  if (!strcat(tmp, " ") || !strcat(tmp, answers[i]))
	    messcrash("string operation failed while concatenating non-passwd part of user request.") ;
	  free(answers[i]) ;
	}
      if (passwd_prompt >= 0)
	{
	  if (!strcat(tmp, " ") || !strcat(tmp, passwd_str))
	    messcrash("string operation failed while concatenating passwd part of user request.") ;
	}

      *cp = tmp ;					    /* Return the result. */
    }

  return status ;
}


/************************************************************/

BOOL getDomainUpdate(char **cp)
{
  BOOL status = TRUE ;
  enum {MAX_PROMPTS = 2} ;
  char *prompts[MAX_PROMPTS] = {"Please enter domain to be updated: "} ;
  char *groupnew = "Please enter name of domains new group: " ;
  char *group = "Please enter name of domains group: " ;
  char *answers[MAX_PROMPTS] ;
  int num_prompts, passwd_prompt = -1 ;

  /* NOTE THAT THIS IS NOT SECURE AT THE MOMENT, REALLY WE SHOULD ASK THE    */
  /* ADMIN USER FOR THEIR OWN PASSWORD EVERY TIME THEY WANT TO UPDATE ANY    */
  /* ASPECT OF THE PASSWORD FILE AT ALL...OTHERWISE ALL SORTS OF TRICKS      */
  /* CAN BE PLAYED, SUCH AS MOVING A DIFFERENT USER INTO THE ADMIN GROUP OR  */
  /* CHANGING THE ADMIN PASSWORD OR SETTING THE DATABASE TO WORLD WRITEABLE. */
  /*                                                                         */
  /* BUT, if we make it secure, this will become a pain for the admin person */
  /* because they will have to enter their password  all the time and        */
  /* actually I think root users on UNIX are not asked their password when   */
  /* they change someone elses...                                            */


  /* Now get the current password and the new password from the user doing   */
  /* all the non-echoing stuff...                                            */
  /*                                                                         */
  messout("Updating domain...") ;

  if (strstr(*cp, "group"))
    {
      prompts[1] = groupnew ;
      num_prompts = 2 ;
    }
  else if (strstr(*cp, "new"))
    {
      prompts[1] = group ;
      num_prompts = 2 ;
    }
  else if (strstr(*cp, "delete"))
    {
      /* No extra prompts.                                                   */
      num_prompts = 1 ;
    }
  else
    {
      status = FALSE ;
    }

  promptUser(num_prompts, prompts, answers, passwd_prompt) ;

  /* THIS IS DISABLED CURRENTLY, SEE COMMENTS ABOVE....                      */
#ifdef ED_G_NEVER_INCLUDE_THIS_CODE
  /* Check that current password matches the one person logged on with.      */
  hash_strings[0] = userid ;
  hash_strings[1] = answers[0] ;
  passwd_hash = hashAndHexStrings(hash_strings, HASH_STRING_NUM) ;
  if (strcmp(curr_passwd_hash, passwd_hash) != 0)
    status = FALSE ;
  free(passwd_hash) ;					    /* no longer needed. */
#endif /* ED_G_NEVER_INCLUDE_THIS_CODE */

  /* Now construct the request to the server.                                */
  if (status)
    {
      int i, length ;
      char *tmp ;

      length = strlen(*cp) ;
      for (i = 0 ; i < num_prompts ; i++)
	{
	  length += (1 + strlen(answers[i])) ;
	}
      length += 1 ;
      tmp = messalloc(length) ;
      
      if (!strcpy(tmp, *cp))
	messcrash("string operation failed while copying initial user request.") ;
      for (i = 0 ; i < num_prompts ; i++)
	{
	  if (!strcat(tmp, " ") || !strcat(tmp, answers[i]))
	    messcrash("string operation failed while concatenating non-passwd part of user request.") ;
	  free(answers[i]) ;
	}

      *cp = tmp ;					    /* Return the result. */
    }

  return status ;
}



/************************************************************/
/************************************************************/

/* uugghh, horrible interface, we assume that password stuff will be last,   */
/* we only turn on non-echoing for password stuff.                           */
static void promptUser(int num_prompts, char *prompts[], char *answers[], int passwd_prompt)
{
  struct termios term ;
  int i ;

  /* Get anything that doesn't need echo'ing to be turned off.               */
  if (passwd_prompt != 0)
    {
      int p = (passwd_prompt < 0) ? num_prompts : passwd_prompt ;
	
      for (i = 0 ; i < p ; i++)
	{
	  answers[i] = myreadline(prompts[i]) ;
	}
    }

  /* Now get the password stuff where we set the terminal non-echo'ing.      */
  if (passwd_prompt >= 0)
    {
      /* Get the terminal attributes to manipulate echo'ing.                     */
      if (tcgetattr(fileno(stdin), &term) < 0)
	messcrash("Unable to get terminal attributes for input of password.");

      term.c_lflag &= ~(ECHO) ;				    /* Turn echo off. */
      if (tcsetattr(fileno(stdin), TCSADRAIN, &term) < 0)
	messcrash("Unable to set terminal attributes for input of password.");
      
      for (i = passwd_prompt ; i < num_prompts ; i++)
	{
	  answers[i] = myreadline(prompts[i]) ;

	  printf("\n") ;				    /* newline after password is not
							       echoed to terminal so send one. */
	}

      /* Really some basic password checking should be done here or somewhere*/
      /* to make sure the user does not just enter a NULL string etc.        */


      term.c_lflag |= ECHO ;				    /* Turn echo on again. */
      if (tcsetattr(fileno(stdin), TCSADRAIN, &term) < 0)
	messcrash("Unable to reset terminal attributes after input of password.");
    }

  return ;
}



/* THESE SHOULDN'T BE HERE....                                               */
/* these two routines should be in a socket lib....                          */
/*                                                                           */
/* They are only accessed from saceclient.c at the moment...                 */
/*                                                                           */
S_MSGState writeToSocket(int sock, void *handle, char *msg_type, char *request)
{
  S_MESSAGE msg = (S_MESSAGE) handle ;
  S_MSGState state ;

  directReply (msg, request) ;
  sMessageSetWrite(msg) ;

  setMsgType(msg->ah.msgType, msg_type) ;

  /* note that this particular construct may block forever if we can;t write */
  /* for some reason other than as error.                                    */
  state = SMSG_WAIT ;
  while(state == SMSG_WAIT)
    state = sMessageSocketWrite(sock, msg) ;
      
  return state ;
}


S_MSGState readFromSocket(int sock, void *handle, char **msg_type_out, 
			  char **answerp, int *lengthp)
{
  S_MESSAGE msg = (S_MESSAGE) handle ;
  char *answer ;
  static S_MSGState state = SMSG_DONE ;

  if (state == SMSG_DONE)
    {
      msg->isNew = TRUE ;
      sMessageSetRead(msg) ;
    }
  else
    msg->isNew = FALSE ;

  if (debug_G)
    messout("about to read from socket, hbytes pending: %d, mbytespending: %d\n",
	    msg->hBytesPending, msg->mBytesPending) ;

  state = sMessageSocketRead(sock, msg) ;

  if (debug_G)
    messout("done read from socket, hbytes pending: %d, mbytespending: %d\n",
	    msg->hBytesPending, msg->mBytesPending) ;
  
  if (state == SMSG_DONE)
    {
      *msg_type_out = msg->ah.msgType ;

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
    }

  return state ;
}



/************************************************************/
/************************************************************/





