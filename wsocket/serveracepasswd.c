/*  File: serveracepasswd.c
 *  Author: Ed Griffiths (edgrif@sanger.ac.uk)
 *  Copyright (c) J Thierry-Mieg and R Durbin, 2000
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.crbm.cnrs-mop.fr
 *
 * Description: Implements the new password/group handling functions for
 *              the server.
 * Exported functions: See serverace_.h
 * HISTORY:
 * Last edited: Feb 28 14:32 2001 (edgrif)
 * Created: Mon Mar  6 21:40:26 2000 (edgrif)
 * CVS info:   $Id: serveracepasswd.c,v 1.2 2006/06/22 20:21:34 mieg Exp $
 *-------------------------------------------------------------------
 */
#include <wh/regular.h>
#include <wh/aceio.h>
#include <wh/command.h>
#include <wh/session.h>
#include <wh/dbpath.h>
#include <wsocket/serverace_.h>
#include <wsocket/serverclientutils.h>

/************************************************************/

/* Updates to the password file are either for users or for domain name      */
/* updates. There are different actions for the two types of update.         */
typedef enum _PasswdGroupEntryType {USER_ENTRY, DOMAIN_ENTRY} PasswdGroupEntryType ;



/* PUT THIS INTO REGULAR.H or a new UTILS.H ??                               */
/* Noddy macro to produce an array large enough to hold the printed form of  */
/* an integer. The relationship between number of bytes for the integer and  */
/* its maximum number of digits (for 8 bit bytes !) is:                      */
/*                                                                           */
/*             digits = 0.49 + (2.42 * bytes)                                */
/*                                                                           */
/* a safe version of this is:  (int)(0.5 + (2.5 * bytes))                    */
/*                                                                           */
/* but then you need the "+ 2" for the string terminator and the sign.       */
/*                                                                           */
/* Actually you may need to remember a lot more if you choose the exotic     */
/* forms of sprintf to print your numbers (e.g. with commas: 1,000,242 etc). */
/* In other words, this macro is for integers in their simplest form only:   */
/* 222343243 or -4440584054 etc.                                             */
/*                                                                           */
/* you supply the type, e.g.    MAKE_NUM_STRING(unsigned long)               */
/*                                                                           */
#define MAKE_INT_STRING(type) \
(((int)((sizeof(type) * 2.5) + 0.5)) + 2)

/************************************************************/

static AceUserGroup getUserDBGroup(AceServ server, PasswdGroupEntryType update_type, char *userid) ;
static ACEIN findAndOpenServerFile(AceServ server, char **path_out) ;
static char *getFileBackupName(char *passwd_path) ;
static BOOL addUser(Array passwd_entries, char *userid, char *hash) ;
static BOOL updateUserEntry(Array passwd_file, char *userid, char *hash) ;
static BOOL removeUser(Array passwd_entries, char *userid) ;
static BOOL addEntryToGroup(Array passwd_entries, char *userid, char *group) ;
static BOOL removeEntryFromGroup(Array passwd_entries, char *userid) ;
static BOOL updatePasswdFile(AceServ server, PasswdGroupEntryType update_type,
			     char *userid, char *group, char *passwd_hash) ;


/* Add code to retrieve error messages, so calling code can just call a      */
/* routine to find out what the rc was, perhaps we should return something   */
/* other than BOOL, DON'T use a static errno....this code could end up being */
/* threaded....                                                              */


/************************************************************/
/************************************************************/



/* OK THIS COPY ROUTINE NEEDS RATIONALISING WITH THE OTHER ONE IN THIS       */
/* CODE...SHOULD BE FAIRLY STRAIGHT FORWARD...                               */

/* Make an update to the server file to change the database global read/write*/
/* permissions.                                                              */
/*                                                                           */
/* Change the database global read/write permissions in the server file. We  */
/* assume that the line containing the users entry has this format:          */
/*             perm_type perm_level                                          */
/*             |                                                             */
/*              --- beginning of line.                                       */
/*                                                                           */
/* We want any altering of the server file to be as failsafe as possible     */
/* so we:                                                                    */
/*      - make a backup copy of server file in the wspec directory.          */
/*      - read the contents of the server file into memory                   */
/*      - alter the in-memory contents                                       */
/*      - write the in-memory contents out to a temporary file               */
/*      - overwrite the original passwd file with the new one (N.B. rename() */
/*        cannot be used because it does not span filesystems)               */
/*      - if the overwrite worked, remove the backup server file and the     */
/*        temporary file                                                     */
/*                                                                           */
/* This seems long-winded but remember just how much NFS is used for         */
/* database disks !!                                                         */
/*                                                                           */
BOOL  updateGlobalDBPerms(AceServ server, char *perm_type, char *perm_level)
{
  BOOL status = TRUE ;
  ACEIN server_file = NULL ;
  char *server_path = NULL ;
  char *backup_name = NULL ;
  Array contents ;
  int i ;
  ACETMP tmp_copy = NULL ;
  char *tmp_name = NULL ;


  /* Check arguments.                                                        */
  /* Currently we insist on upper case in the server file which makes this a */
  /* little long winded...                                                   */
  if (strcasecmp(perm_type, ACESERV_READ) == 0)
    perm_type = ACESERV_READ ;
  else if (strcasecmp(perm_type, ACESERV_WRITE) == 0)
    perm_type = ACESERV_WRITE ;
  else
    status = FALSE ;

  if (status)
    {
      if (strcasecmp(perm_level, ACESERV_NONE) == 0)
	perm_level = ACESERV_NONE ;
      else if (strcasecmp(perm_level, ACESERV_PASSWD) == 0)
	perm_level = ACESERV_PASSWD ;
      else if (strcasecmp(perm_level, ACESERV_WORLD) == 0)
	perm_level = ACESERV_WORLD ;
      else
	status = FALSE ;
    }

  /* Try and open the server file.                                         */
  if (status)
    {
      server_file = findAndOpenServerFile(server, &server_path) ;
      if (!server_file)
	messcrash("Cannot continue without creating temporary copy of server file.") ;
      server_path = strnew(server_path, 0) ;		    /* Need a copy because filXXX routines */
							    /* are not reentrant. */
      aceInSpecial(server_file, "\n") ;
    }

  /* Make a backup copy of the server file in the server file directory. */
  if (status)
    {
      backup_name = getFileBackupName(server_path) ;
      status = filCopyFile(server_path, backup_name) ;
      if (!status)
	messfree(backup_name) ;
    }

  /* Read the server server file into an array.                            */
  if (status)
    {
      contents = arrayCreate(64, char*) ;

      while (aceInCard(server_file))
	{
	  char *tmp = strnew(aceInPos(server_file), 0) ;
	  array(contents, arrayMax(contents), char*) = tmp ;
	}
    }

  /* Change the global permissions.                                          */
  if (status)
    {
      int i ;

      /* Find the global entry, delete and add the new one.                      */
      for (i = 0, status = FALSE ; i < arrayMax(contents) && status == FALSE ; i++)
	{
	  if (strncasecmp(array(contents, i, char*), perm_type, strlen(perm_type)) == 0)
	    {
	      messfree(array(contents, i, char*)) ;
	      array(contents, i, char*) = messalloc(strlen(perm_type) + 1 + strlen(perm_level) + 1) ;
	      if (strcpy(array(contents, i, char*), perm_type) == NULL
		  || strcat(array(contents, i, char *), " ") == NULL
		  || strcat(array(contents, i, char*), perm_level) == NULL)
		messcrash("Standard str func failed while trying to construct "
			  "global database permissions entry for server file.") ;
	      status = TRUE ;
	    }
	}
    }


  /* Put the updated server file contents into a temporary file.             */
  if (status)
    {
      int file_rc ;
      ACEOUT tmp_out ;

      tmp_copy = aceTmpCreate("w", 0) ;
      if (!tmp_copy)
	messcrash("Cannot continue without creating temporary copy of server file.") ;

      tmp_name = aceTmpGetFileName(tmp_copy) ;

      /* There is a bug here where, we end up writing an extra newline to the*/
      /* file, I reckon this is a problem with reading stuff in probably and */
      /* getting an extra line for the end of file...check this out...       */
      tmp_out = aceTmpGetOutput(tmp_copy) ;
      for (i = 0, file_rc = 0 ; i < arrayMax(contents) && file_rc == 0 ; i++)
	{
	  /* We may have set some elements to NULL because of not being able */
	  /* to use arrayRemove (see removeUser() ), so don't print these.   */
	  if (array(contents, i, char*) != NULL)
	    file_rc = aceOutf(tmp_out, "%s\n", array(contents, i, char*)) ;
	}

      if (file_rc != 0)
	{
	  status = FALSE ;
	  
	  messcrash("Could not copy server server file updates to temporary file \"%s\", "
		    "reason was: %s", tmp_name, messSysErrorText()) ;
	}

      aceTmpClose(tmp_copy) ;
    }

  /* Now replace the server file with the new one.                           */
  if (status)
    {
      status = filCopyFile(tmp_name, server_path) ;
      if (!status)
	  messcrash("Could not replace server server file with new server file \"%s\", "
		    "reason was: %s. "
		    "There is a backup copy in \"%s\".",
		    tmp_name, messSysErrorText(), backup_name) ;
    }

  /* clear up                                                                */
  if (server_file)
      aceInDestroy(server_file) ;

  if (server_path)
    messfree(server_path) ;

  if (contents)
    arrayDestroy(contents) ;

  if (tmp_copy)
    aceTmpDestroy(tmp_copy) ;				    /* n.b. removes tmp_copy as well. */


  /* Remove the server backup file if the update worked.                     */
  if (status)
    {
      if (unlink(backup_name) != 0)
	messcrash("Cannot remove backup server file \"%s\".", backup_name) ;
    }

  return status ;
}





/************************************************************/

/* Routine to change users password.                                         */
/* Arguments are, userid, nonce for that                                     */
/*                                                                           */
/*                                                                           */
BOOL changeUserPasswd(AceServ server, char *userid, char *nonce, char *old_hash, char *new_hash,
		      ClientAccessStatus *status_out)
{
  BOOL status = TRUE ;

  /* Is the user in a group in the passwd file ?                             */
  if (getUserDBGroup(server, USER_ENTRY, userid) == ACEUSER_NONE)
    {
      status = FALSE ;
      *status_out = CLIENTACCESS_UNKNOWN ;
    }

  /* Is the old password correct ?                                           */
  if (status == TRUE)
    {
      if (!checkUserHash(server, userid, old_hash, nonce))
	{
	  status = FALSE ;
	  *status_out = CLIENTACCESS_BADHASH ;
	}
    }

  /* ok, change the password.                                                */
  if (status == TRUE)
    {
      status = updatePasswdFile(server, USER_ENTRY, userid, NULL, new_hash) ;
    }

  return status ;
}


/************************************************************/

/* Make changes to the server password file for individual users.            */
/*                                                                           */
/* words[0] = "user"  words[1] = "passwd|group|new|delete" words[2] = userid */
/* words[3} & words[4] depend on what words[1] is, but will be either a      */
/*   group name ("admin|write|read") or a password hash.                     */
/*                                                                           */
/* status_out is where we return the type of error if there was one.         */
/*                                                                           */
BOOL updateUser(AceServ server, Array words, ClientAccessStatus *status_out)
{
  BOOL status = TRUE ;
  int num_words = arrayMax(words) ;
  char *request, *userid ;

  request = array(words, 1, char*) ;
  userid = array(words, 2, char*) ;

  /* Do some verification:                                                   */
  /* Must not allow any "." chars in userid because this could be used to    */
  /* add users with names that are the same as a domain name, this would then*/
  /* allow access to users of that whole domain because we have both userids */
  /* and domain names in the groups.                                         */
  /* We also don't allow ":" in userid because it is the separator we use    */
  /* for groups. And we don't userids that are the same as the three group   */
  /* names.                                                                  */
  if (strstr(userid, ".") != NULL
      || strstr(userid, ":") != NULL
      || strcmp(userid, ACEGROUP_ADMIN_TEXT) == 0
      || strcmp(userid, ACEGROUP_WRITE_TEXT) == 0
      || strcmp(userid, ACEGROUP_READ_TEXT) == 0)
    {
      status = FALSE ;
      *status_out = CLIENTACCESS_BADHASH ;
    }

  /* We don't check the userid for permissions because we already know that  */
  /* this is an 'admin' id. so we just get on and update the user.           */
  if (status)
    {
      if ((strcmp(request, "passwd") == 0)
	  && num_words == 4)
	{
	  status = updatePasswdFile(server, USER_ENTRY, userid, NULL, array(words, 3, char*)) ;
	}
      else if ((strcmp(request, "group") == 0)
	       && num_words == 4)
	{
	  status = updatePasswdFile(server, USER_ENTRY, userid, array(words, 3, char*), NULL) ;
	}
      else if ((strcmp(request, "new") == 0)
	       && num_words == 5)
	{
	  status = updatePasswdFile(server, USER_ENTRY, userid, array(words, 3, char*),
				    array(words, 4, char*)) ;

	}
      else if ((strcmp(request, "delete") == 0)
	       && num_words == 3)
	{
	  status = updatePasswdFile(server, USER_ENTRY, userid, NULL, NULL) ;
	}
      else
	{
	  status = FALSE ;
	  *status_out = CLIENTACCESS_BADUSER_REQUEST ;
	}
    }

  return status ;
}


/************************************************************/

/* Change a domain name entry in the password file.                          */
/* We are expecting: "domain group|new|delete name .."                       */
/*                      (where ".." is dependent on the 2nd word)            */
/*                                                                           */
BOOL updateDomain(AceServ server, Array words, ClientAccessStatus *status_out)
{
  BOOL status = FALSE ;
  int num_words = arrayMax(words) ;
  char *str_pos, *domain_name, *request ;

  request = array(words, 1, char*) ;
  domain_name = array(words, 2, char*) ;

  /* Do some verification:                                                   */
  /* Domain names MUST HAVE at least one "." to distinguish them from user   */
  /* ids and to ensure that they are sensible, i.e. if they had no dot this  */
  /* would mean a domain along the lines of the whole of the "uk" or "com"   */
  /* was being allowed....not good.                                          */
  /* The dot must be somewhere inside the string, not at beginning or end.   */
  str_pos = strstr(domain_name, ".") ;
  if (str_pos != NULL
      && str_pos != domain_name
      && str_pos != (domain_name + (strlen(domain_name) - 1)))
    status = TRUE ;

  if (status)
    {
      if ((strcmp(request, "group") == 0)
	  && num_words == 4)
	{
	  /* Change the group of an existing domain name.                    */
	  status = updatePasswdFile(server, DOMAIN_ENTRY, domain_name, array(words, 3, char*), NULL) ;
	}
      else if ((strcmp(request, "new") == 0)
	       && num_words == 4)
	{
	  /* Add a new domain name.                                          */
	  status = updatePasswdFile(server, DOMAIN_ENTRY, domain_name, array(words, 3, char*), NULL) ;

	}
      else if ((strcmp(request, "delete") == 0)
	       && num_words == 3)
	{
	  /* Delete an existing domain name.                                 */
	  status = updatePasswdFile(server, DOMAIN_ENTRY, domain_name, NULL, NULL) ;
	}
      else
	{
	  /* should set a reason here to be returned via acc_status.              */
	  status = FALSE ;
	  *status_out = CLIENTACCESS_BADDOMAIN_REQUEST ;
	}
    }

  return status ;
}


/************************************************************/

/* Routine to construct a "nonce", this is a 'random' value passed to the    */
/* client which the client must hash with its hash of the userid/passwd of   */
/* the user before passing it back.                                          */
/*                                                                           */
char *makeNonce(void *client_addr)
{
  char *nonce = NULL ;
  unsigned long addr ;
  int random ;
  mytime_t time ;
  enum {NONCE_ITEMS = 3} ;
  char addr_str[MAKE_INT_STRING(addr)] ;
  char random_str[MAKE_INT_STRING(random)] ;
  char time_str[MAKE_INT_STRING(time)] ;
  char *hash_strings[NONCE_ITEMS] ;

  hash_strings[0] = &addr_str[0], hash_strings[1] = &random_str[0],
    hash_strings[2] = &time_str[0] ;


  /* Items to make up our nonce, we want values that are not related to the  */
  /* client but also values that are unique per client.                      */
  /*                                                                         */
  /*   - addr is the address of the memory malloc'd to hold the client which */
  /*     should be unique to a client and not predictable.                   */
  /*   - random is as good a random integer as the system can provide        */
  /*   - time will also be unique per client because it is ever increasing   */
  /*                                                                         */
  addr = (unsigned long)client_addr ;
  random = randint() ;
  time = timeNow() ;
  
  /* Put the values in to string form as plain integers.                     */
  if (sprintf(addr_str, "%lu", addr) < 1
      || sprintf(random_str, "%d", random) < 1
      || sprintf(time_str, "%u", time) < 1)
    messcrash("a sprintf to create strings for nonce value failed, reason was: %s",
	      messSysErrorText()) ;

  /* make an md5 hash of the values as a hex string.                         */  
  nonce = hashAndHexStrings(hash_strings, NONCE_ITEMS) ;

  return nonce ;
}


/**************************************************/

/* Read the server configuration file wspec/serverconfig.wrm and find out what     */
/* the global read/write permissions for the database are.                   */
/* It is a serious error if the file cannot be found/opened and we return    */
/* FALSE, otherwise we return TRUE.                                          */
/* Note that permissions will default to 'passwd' if the read or write       */
/* keywords cannot be found or are found but either no option or an incorrect*/
/* option is specified.                                                      */
/*                                                                           */
/* READ & WRITE options are independent except if WRITE is set to WORLD, in  */
/* which case READ is also set to WORLD. This is because anyone with WRITE   */
/* access also has READ access so if anyone can write then anyone can also   */
/* read.                                                                     */
/*                                                                           */
BOOL checkGlobalDBAccess(AceServPermOptions *read_acc, AceServPermOptions *write_acc)
{
  BOOL status = TRUE ;
  char *filename ;
  ACEIN acefi ;

  /* Can we find the file ?                                                  */
  if (status == TRUE)
    {
      filename = dbPathStrictFilName("wspec", ACESERV_CONFIG_NAME, ACESERV_CONFIG_EXT, "r", 0) ;
      if (filename == NULL)
	{
	  messdump("Server cannot find the "
		   "\"wspec/"ACESERV_CONFIG_NAME"."ACESERV_CONFIG_EXT"\" configuration file, "
		   "the server will not run on this database until it can find the file.") ;
	  status = FALSE ;
	}
    }

  /* Can we open the file ?                                                  */
  if (status == TRUE)
    {
      acefi = aceInCreateFromFile(filename, "r", "", 0) ;
      if (acefi == NULL)
	{
	  messdump("Server cannot open the "
		   "\"wspec/"ACESERV_CONFIG_NAME"."ACESERV_CONFIG_EXT"\" configuration file, "
		   "the server will not run on this database until it can open this file.") ;
	  status = FALSE ;
	}
    }

  /* Read the permissions.                                                   */
  if (status == TRUE)
    {
      /* Set defaults to be safe, i.e. password protected.                   */
      *write_acc = *read_acc = PERM_PASSWD ;

      while (aceInCard(acefi))
	{
	  char *word ;
	  
	  if ((word = aceInWord(acefi)))
	    {
	      if (strcmp(word, ACESERV_WRITE) == 0)
		{
		  if ((word = aceInWord(acefi)))
		    {
		      if (strcmp(word, ACESERV_NONE) == 0)
			*write_acc = PERM_NONE ;
		      else if (strcmp(word, ACESERV_PASSWD) == 0)
			*write_acc = PERM_PASSWD ;
		      else if (strcmp(word, ACESERV_WORLD) == 0)
			*write_acc = PERM_WORLD ;
		    }
		}
	      else if (strcmp (word, ACESERV_READ) == 0)
		{
		  if ((word = aceInWord(acefi)))
		    {
		      if (strcmp(word, ACESERV_NONE) == 0)
			*read_acc = PERM_NONE ;
		      else if (strcmp(word, ACESERV_PASSWD) == 0)
			*read_acc = PERM_PASSWD ;
		      else if (strcmp(word, ACESERV_WORLD) == 0)
			*read_acc = PERM_WORLD ;
		    }
		}
	    }
	}

      /* write_acc == PERM_WORLD implies _anyone_ can read or write so the   */
      /* read permission should also be world.                               */
      if (*write_acc == PERM_WORLD)
	*read_acc = PERM_WORLD ;
    }

  if (filename)
    messfree(filename) ;

  if (acefi)
    aceInDestroy(acefi) ;

  return status ;
}
 

/**************************************************/ 

/* This routine encapsulates the rules for giving access to a user given     */
/* the global access permissions for the database and the users possible     */
/* entry in the password file. The intention is to keep all this logic in    */
/* one routine and not scatter it around the code.                           */
/* On failure ACEUSER_NONE is returned, the status_out indicates the reason  */
/* for the failure.                                                          */
/*                                                                           */
AceUserGroup getUserDBAccess(AceServ server,
			     char *userid, char *hash, char *nonce,
			     AceServPermOptions globalReadAcc,
			     AceServPermOptions globalWriteAcc,
			     ClientConnect connection,
			     ClientAccessStatus *status_out)
{
  AceUserGroup group = ACEUSER_NONE ;

  /* Is the user in a group in the passwd file ?                             */
  group = getUserDBGroup(server, USER_ENTRY, userid) ;

  if (group == ACEUSER_ADMIN)
    {
      /* admin users always have read/write permissions, but must have a     */
      /* password.                                                           */
      if (!checkUserHash(server, userid, hash, nonce))
	{
	  group = ACEUSER_NONE ;
	  *status_out = CLIENTACCESS_BADHASH ;
	}
    }
  else if (group == ACEUSER_WRITE)
    {
      if (globalWriteAcc == PERM_WORLD)
	{
	  /* Nothing to do                                                   */
	}
      else if (globalWriteAcc == PERM_PASSWD)
	{
	  if (checkUserHash(server, userid, hash, nonce))
	    group = ACEUSER_WRITE ;
	  else
	    {
	      group = ACEUSER_NONE ;
	      *status_out = CLIENTACCESS_BADHASH ;
	    }
	}
      else
	{
	  group = ACEUSER_NONE ;
	  *status_out = CLIENTACCESS_BLOCKED ;
	}
    }
  else if (group == ACEUSER_READ)
    {
      if (globalReadAcc == PERM_WORLD)
	{
	  /* Nothing to do                                                   */
	}
      else if (globalReadAcc == PERM_PASSWD)
	{
	  if (checkUserHash(server, userid, hash, nonce))
	    group = ACEUSER_READ ;
	  else
	    {
	      group = ACEUSER_NONE ;
	      *status_out = CLIENTACCESS_BADHASH ;
	    }
	}
      else
	{
	  group = ACEUSER_NONE ;
	  *status_out = CLIENTACCESS_BLOCKED ;
	}
    }
  else /* group == ACEUSER_NONE */
    {
      /* We can only give access to an unknown user if the database has      */
      /* 'world' permissions for either read or write, or the user is        */
      /* calling from one of an allowed list of host machines/domains.       */
      if (globalWriteAcc == PERM_WORLD)
	group = ACEUSER_WRITE ;
      else if (globalReadAcc == PERM_WORLD)
	group = ACEUSER_READ ;
      else
	{
	  /* Can only do this if the NO_HOSTNAME_RESOLUTION option is _not_  */
	  /* set.                                                            */
	  if (checkForServerConfigValue(server, ACESERV_NODOMAINRESOLVE) == FALSE)
	    {
	      char *dotted, *plain ;

	      if (getClientInetAddressInfo(connection, &dotted, &plain))
		{
		  group = getUserDBGroup(server, DOMAIN_ENTRY, plain) ;
		  if (group == ACEUSER_ADMIN)
		    {
		      /* admin group not allowed to have domain names in it. */
		      group = ACEUSER_NONE ;
		      *status_out = CLIENTACCESS_ADMINHOSTNAME ;
		    }
		  else if (group == ACEUSER_NONE)
		    {
		      *status_out = CLIENTACCESS_BLOCKED ; 
		    }
		  /* Otherwise they must in the read or write groups.        */
		}
	    }
	  else
	    {
	      group = ACEUSER_NONE ;
	      *status_out = CLIENTACCESS_UNKNOWN ;
	    }
	}
    }

  return group ;
}


/************************************************************/

/* Reads the file <database>/wspec/serverpasswd.wrm to see if the supplied   */
/* userid is any of the read, write or admin groups in the file.             */
/* Returns ACEUSER_NONE if it cannot find the userid in the groups.          */
static AceUserGroup getUserDBGroup(AceServ server, PasswdGroupEntryType group_type, char *userid)
{
  AceUserGroup group = ACEUSER_NONE ;
  BOOL status = TRUE ;
  ACEIN acefi = NULL ;

  /* Can we open the server password file ?                                  */
  if (status == TRUE)
    {
      acefi = findAndOpenPasswdFile(server, NULL) ;
      if (!acefi)
	status = FALSE ;
    }

  /* Look for the userid in the file.                                        */
  if (status == TRUE)
    {
      BOOL found = FALSE ;

      while (aceInCard(acefi) && found == FALSE)
	{
	  char *word ;
	  
	  if ((word = aceInWord(acefi)))
	    {
	      /* Does this line begin with one of the groups ?               */
	      if (strcmp(word, ACEGROUP_ADMIN) == 0)
		group = ACEUSER_ADMIN ;
	      else if (strcmp(word, ACEGROUP_WRITE) == 0)
		group = ACEUSER_WRITE ;
	      else if (strcmp(word, ACEGROUP_READ) == 0)
		group = ACEUSER_READ ;
	      else
		group = ACEUSER_NONE ;

	      /* If we found a group then look along the line for the userid.*/
	      /* or domain.                                                  */
	      if (group != ACEUSER_NONE)
		{
		  while ((word = aceInWord(acefi)) && found == FALSE)
		    {
		      if (group_type == USER_ENTRY)
			{
			  if (strcmp(word, userid) == 0)
			    found = TRUE ;
			}
		      else /* DOMAIN_ENTRY */
			{
			  if (strstr(word, ".") != NULL)
			    {
			      /* Is the end of the users domain the same as  */
			      /* the domain in the password file ?           */
			      if (strstr(userid, word) == userid + strlen(userid) - strlen(word))
				found = TRUE ;
			    }
			}
		    }
		}
	    }
	}

      if (!found)
	group = ACEUSER_NONE ;
    }

  if (acefi)
    aceInDestroy(acefi) ;

  return group ;
}


/************************************************************/

/* This routine checks out whether the correct password for a user has been  */
/* sent by the client.                                                       */
/* It does this by getting an encrypted form of the password from the file   */
/* wspec/serverpasswd.wrm and comparing it to the encrypted password sent    */
/* by the client.                                                            */
/* The routine will return FALSE if it cannot open the password file or if   */
/* it can't find the clients userid in the file or if it found the userid    */
/* but there was no hash following the id.                                   */
/*                                                                           */
BOOL checkUserHash(AceServ server, char *userid, char *client_hash, char *nonce)
{
  BOOL hash_ok = FALSE ;
  BOOL status = TRUE ;
  ACEIN acefi = NULL ;
  char *server_hash = NULL ;
  enum {HASH_ITEMS = 2} ;
  char *hash_strings[HASH_ITEMS]  ;
  char *digest_str ;

  /* Can we open the server password file ?                                  */
  if (status == TRUE)
    {
      acefi = findAndOpenPasswdFile(server, NULL) ;
      if (!acefi)
	status = FALSE ;
    }

  /* Look for the userid in the file and get the passwd hash.                */
  if (status == TRUE)
    {
      BOOL found = FALSE ;

      while (aceInCard(acefi) && found == FALSE)
	{
	  char *word ;

	  /* Try and find the userid and hash in the file.                   */
	  if ((word = aceInWord(acefi)))
	    {
	      if (strcmp(word, userid) == 0)
		{
		  found = TRUE ;			    /* Stop now we've found the userid. */
		  if ((word = aceInWord(acefi)))
		    server_hash = strnew(word, 0) ;
		}
	    }
	}
    }

  if (acefi)
    aceInDestroy(acefi) ;

  /* Did the client send the correct password ?                              */
  /* The comparison involves:                                                */
  /*       1) MD5 the nonce with server copy of the userid/passwd hash       */
  /*       2) compare this server hash with the one the client sent          */
  /*                                                                         */
  if (server_hash != NULL)
    {
      hash_strings[0] = server_hash ;
      hash_strings[1] = nonce ;
      digest_str = hashAndHexStrings(hash_strings, HASH_ITEMS) ;
      
      if (strcmp(digest_str, client_hash) == 0)
	hash_ok = TRUE ;
    }

  return hash_ok ;
}


/************************************************************/

/* Make an update to the password file, note that this function can be used  */
/* in four ways depending on which parameters are NULL:                      */
/*                                                                           */
/*      1) To add a new user, pass userid, group and passwd                  */
/*      2) To change an existing users password, pass userid & passwd, set   */
/*         group to null                                                     */
/*      3) To change an existing users group, pass userid and group, set     */
/*         passwd to null                                                    */
/*      4) To delete a user, pass just the userid.                           */
/*                                                                           */
/* We want any altering of the password file to be as failsafe as possible   */
/* so we:                                                                    */
/*      - make a backup copy of password file in the wspec directory.        */
/*      - read the contents of the password file into memory                 */
/*      - alter the in-memory contents                                       */
/*      - write the in-memory contents out to a temporary file               */
/*      - overwrite the original passwd file with the new one (N.B. rename() */
/*        cannot be used because it does not span filesystems)               */
/*      - if the overwrite worked, remove the backup password file and the   */
/*        temporary file                                                     */
/*                                                                           */
/* This seems long-winded but remember just how much NFS is used for         */
/* database disks !!                                                         */
/*                                                                           */
static BOOL updatePasswdFile(AceServ server, PasswdGroupEntryType group_type,
			     char *userid, char *group, char *passwd_hash)
{
  BOOL status = TRUE ;
  ACEIN passwd_file = NULL ;
  char *passwd_path = NULL ;
  char *backup_name = NULL ;
  Array contents ;
  int i ;
  ACETMP tmp_copy = NULL ;
  char *tmp_name = NULL ;


  /* As this function has several different uses, check the input.           */
  /* (currently trivial)                                                     */
  if (userid == NULL)
    messcrash("updatePasswdFile() received bad call, arguments were  "
	      "userid: %s   group: %s  passwd_hash: %s",
	      (userid ? userid : "NULL"),
	      (group ? group : "NULL"),
	      (passwd_hash ? passwd_hash : "NULL")) ;

  /* Try and open the password file.                                         */
  passwd_file = findAndOpenPasswdFile(server, &passwd_path) ;
  if (!passwd_file)
    messcrash("Cannot continue without creating temporary copy of password file.") ;
  passwd_path = strnew(passwd_path, 0) ;		    /* Need a copy because filXXX routines */
							    /* are not reentrant. */
  aceInSpecial(passwd_file, "\n") ;

  /* Make a backup copy of the password file in the password file directory. */
  backup_name = getFileBackupName(passwd_path) ;
  status = filCopyFile(passwd_path, backup_name) ;
  if (!status)
    messfree(backup_name) ;

  /* Read the server password file into an array.                            */
  if (status)
    {
      contents = arrayCreate(64, char*) ;

      while (aceInCard(passwd_file))
	{
	  char *tmp = strnew(aceInPos(passwd_file), 0) ;
	  array(contents, arrayMax(contents), char*) = tmp ;
	}
    }

  /* Add new data to the array.                                              */
  if (status)
    {
      if (group_type == USER_ENTRY)
	{
	  if (userid != NULL && passwd_hash != NULL && group != NULL)
	    {
	      /* new user                                                            */
	      status = addUser(contents, userid, passwd_hash) ;
	      if (status)
		status = addEntryToGroup(contents, userid, group) ;
	    }
	  else if (passwd_hash != NULL)
	    {
	      /* Change user password.                                               */
	      status = updateUserEntry(contents, userid, passwd_hash) ;
	    }
	  else if (group != NULL)
	    {
	      /* Change user group.                                                  */
	      status = removeEntryFromGroup(contents, userid) ;
	      status = addEntryToGroup(contents, userid, group) ;
	    }
	  else
	    {
	      /* Delete a user.                                                      */
	      status = removeUser(contents, userid) ;
	      status = removeEntryFromGroup(contents, userid) ;
	    }
	}
      else /* group_type == DOMAIN_ENTRY */
	{
	  if (userid != NULL && group != NULL)
	    {
	      /* N.B. this is slightly hokey, if this is a new domain name   */
	      /* the call to remove it will fail, we allow this and then     */
	      /* do the add which must be done whether the domain name is a  */
	      /* new one or not.                                             */
	      status = removeEntryFromGroup(contents, userid) ;
	      status = addEntryToGroup(contents, userid, group) ;
	    }
	  else
	    {
	      /* Delete a domain name                                                */
	      status = removeEntryFromGroup(contents, userid) ;
	    }
	}
    }

  /* Put the updated passwd file contents into a temporary file.             */
  if (status)
    {
      int file_rc ;
      ACEOUT tmp_out ;

      tmp_copy = aceTmpCreate("w", 0) ;
      if (!tmp_copy)
	messcrash("Cannot continue without creating temporary copy of password file.") ;

      tmp_name = aceTmpGetFileName(tmp_copy) ;

      /* There is a bug here where, we end up writing an extra newline to the*/
      /* file, I reckon this is a problem with reading stuff in probably and */
      /* getting an extra line for the end of file...check this out...       */
      tmp_out = aceTmpGetOutput(tmp_copy) ;
      for (i = 0, file_rc = 0 ; i < arrayMax(contents) && file_rc == 0 ; i++)
	{
	  /* We may have set some elements to NULL because of not being able */
	  /* to use arrayRemove (see removeUser() ), so don't print these.   */
	  if (array(contents, i, char*) != NULL)
	    file_rc = aceOutf(tmp_out, "%s\n", array(contents, i, char*)) ;
	}

      if (file_rc != 0)
	{
	  status = FALSE ;
	  
	  messcrash("Could not copy server password file updates to temporary file \"%s\", "
		    "reason was: %s", tmp_name, messSysErrorText()) ;
	}

      aceTmpClose(tmp_copy) ;
    }

  /* Now replace the passwd file with the new one.                           */
  if (status)
    {
      status = filCopyFile(tmp_name, passwd_path) ;
      if (!status)
	  messcrash("Could not replace server password file with new password file \"%s\", "
		    "reason was: %s. "
		    "There is a backup copy in \"%s\".",
		    tmp_name, messSysErrorText(), backup_name) ;
    }

  /* clear up                                                                */
  if (passwd_file)
      aceInDestroy(passwd_file) ;

  if (passwd_path)
    messfree(passwd_path) ;

  if (contents)
    arrayDestroy(contents) ;

  if (tmp_copy)
    aceTmpDestroy(tmp_copy) ;				    /* n.b. removes tmp_copy as well. */


  /* Remove the passwd backup file if the update worked.                     */
  if (status)
    {
      if (unlink(backup_name) != 0)
	messcrash("Cannot remove backup passwd file \"%s\".", backup_name) ;
    }

  return status ;
}


/************************************************************/

/* Two noddy functions that try to find and open either the password file    */
/* or the server file, this just saves having the strings embedded all over  */
/* place. The actual work is done by the third function (findAndOpenFile).   */
/*                                                                           */
/* They return NULL on failure or a valid ACEIN on success. Reason for       */
/* failure is logged.                                                        */
/*                                                                           */

ACEIN findAndOpenPasswdFile(AceServ server, char **path_out)
{
  ACEIN passwd_file = NULL ;
  char *filename = ACESERV_PASSWD_NAME ;
  char *filepath = NULL ;

  passwd_file = findAndOpenFile(server, filename, &filepath) ;
  if (passwd_file && path_out != NULL)
    *path_out = filepath ;

  return passwd_file ;
}

static ACEIN findAndOpenServerFile(AceServ server, char **path_out)
{
  ACEIN server_file = NULL ;
  char *filename = ACESERV_CONFIG_NAME ;
  char *filepath = NULL ;

  server_file = findAndOpenFile(server, filename, &filepath) ;
  if (server_file && path_out != NULL)
    *path_out = filepath ;

  return server_file ;
}

/************************************************************/

/* Make the name of the backup password file.                                */ 
/*                                                                           */ 
static char *getFileBackupName(char *file_path)
{
  char *backup_name = NULL ;
  char *backup_ext = ".BACKUP" ;

  backup_name = messalloc(strlen(file_path) + strlen(backup_ext) + 1) ;
  
  if (strcpy(backup_name, file_path) == NULL
      || strcat(backup_name, backup_ext) == NULL)
    messcrash("Standard str func failed while trying to construct \"backup\" file name.") ;

  return backup_name ;
}


/************************************************************/

/* Add a new user to the password file. The format of the entry added is:    */
/*             userid   passwd_hash                                          */
/*             |                                                             */
/*              --- beginning of line.                                       */
/*                                                                           */
/* This entry is simply added to the bottom of the password file.            */
/*                                                                           */
static BOOL addUser(Array passwd_entries, char *userid, char *hash)
{
  BOOL status = TRUE, found ;
  int i ;

  /* Is the entry already there ?, its an error if it is...                  */
  if (status)
    {
      for (i = 0, found = FALSE ; i < arrayMax(passwd_entries) && found == FALSE ; i++)
	{
	  /* We may have set some elements to NULL because of not being able */
	  /* to use arrayRemove (see removeUser() ), so don't process these. */
	  if (array(passwd_entries, i, char*) != NULL)
	    {
	      if (strstr(array(passwd_entries, i, char*), userid) == array(passwd_entries, i, char*))
		{
		  found = TRUE ;
		  status = FALSE ;
		}
	    }
	}
    }

  if (status)
    {
      i = arrayMax(passwd_entries) ;

      array(passwd_entries, i, char*)
	= messalloc(strlen(userid) + 1 + strlen(hash) + 1) ;

      if (strcpy(array(passwd_entries, i, char*), userid) == NULL
	  || strcat(array(passwd_entries, i, char *), " ") == NULL
	  || strcat(array(passwd_entries, i, char*), hash) == NULL)
	messcrash("Standard str func failed while trying to construct "
		  "userid entry for passwd file.") ;
    }

  return status ;
}


/************************************************************/

/* Change the users password in the password file. We assume that the line   */
/* containing the users entry has this format:                               */
/*             userid   passwd_hash                                          */
/*             |                                                             */
/*              --- beginning of line.                                       */
/*                                                                           */
/* If there is white space in front of the userid then the code will fail.   */
/* (note that strtok is no good, it mangles the strings..)                   */
static BOOL updateUserEntry(Array passwd_entries, char *userid, char *hash)
{
  BOOL status = FALSE ;
  int i ;

  /* Find the userid entry, delete and add the new one.                      */
  for (i = 0, status = FALSE ; i < arrayMax(passwd_entries) && status == FALSE ; i++)
    {
      /* We may have set some elements to NULL because of not being able */
      /* to use arrayRemove (see removeUser() ), so don't print these.   */
      if (array(passwd_entries, i, char*) != NULL)
	{
	  if (strstr(array(passwd_entries, i, char*), userid) == array(passwd_entries, i, char*))
	    {
	      messfree(array(passwd_entries, i, char*)) ;
	      array(passwd_entries, i, char*) = messalloc(strlen(userid) + 1 + strlen(hash) + 1) ;
	      if (strcpy(array(passwd_entries, i, char*), userid) == NULL
		  || strcat(array(passwd_entries, i, char *), " ") == NULL
		  || strcat(array(passwd_entries, i, char*), hash) == NULL)
		messcrash("Standard str func failed while trying to construct "
			  "userid entry for passwd file.") ;
	      status = TRUE ;
	    }
	}
    }

  /* Should set an error message here....                                    */

  return status ;
}


/************************************************************/

/* Remove a user from the password file. We assume that the line             */
/* containing the users entry has this format:                               */
/*             userid   passwd_hash                                          */
/*             |                                                             */
/*              --- beginning of line.                                       */
/*                                                                           */
/* If there is white space in front of the userid then the code will fail.   */
static BOOL removeUser(Array passwd_entries, char *userid)
{
  BOOL status = FALSE ;
  int i ;

  /* Find the userid entry and delete it.                                    */
  for (i = 0, status = FALSE ; i < arrayMax(passwd_entries) && status == FALSE ; i++)
    {
      /* We may have set some elements to NULL because of not being able */
      /* to use arrayRemove (see removeUser() ), so don't print these.   */
      if (array(passwd_entries, i, char*) != NULL)
	{
	  if (strstr(array(passwd_entries, i, char*), userid) == array(passwd_entries, i, char*))
	    {
	      messfree(array(passwd_entries, i, char*)) ;
	      array(passwd_entries, i, char*) = NULL ;
	      status = TRUE ;
	      /* I would like to user the following arrayRemove but this seems   */
	      /* only to work on sorted arrays and you have to provide a         */
	      /* function to say if two elements are <, > or == to each other    */
	      /* to help with the sorting...useless to me, I just want to say    */
	      /* "remove this element"...sometime I will stick such a call in to */
	      /* array...now I don't have the time...                            */
#ifdef ED_G_NEVER_INCLUDE_THIS_CODE
	      status = arrayRemove(passwd_entries, array(passwd_entries, i, char*), NULL) ;
	      if (!status)
		messcrash("Could not remove element from array holding password entries while "
			  "trying to remove user userid entry from passwd file.") ;
#endif /* ED_G_NEVER_INCLUDE_THIS_CODE */
	    }
	}
    }

  /* Should set an error message here....                                    */

  return status ;
}


/************************************************************/

/* Remove a user from a group in the password file. We assume that the groups*/
/* containing the users entry have this format:                              */
/*             <groupname>: user1 user2 user3                                */
/*             |                                                             */
/*              --- beginning of line.                                       */
/*                                                                           */
/* where <groupname> is one of  read|write|admin                             */
/*                                                                           */
/* If there is white space in front of the group then the code will fail.    */
/*                                                                           */
static BOOL removeEntryFromGroup(Array passwd_entries, char *userid)
{
  BOOL status = FALSE, group_found, userid_found ;
  int i ;

  /* Find the group entry with the user and delete the user from it.         */
  for (i = 0, group_found = FALSE, userid_found = FALSE ;
       i < arrayMax(passwd_entries) && userid_found == FALSE ; i++)
    {
      /* We may have set some elements to NULL because of not being able */
      /* to use arrayRemove (see removeUser() ), so don't print these.   */
      if (array(passwd_entries, i, char*) != NULL)
	{
	  /* Find a group and look for the userid in it, its an error if we don't    */
	  /* find it.                                                                */
	  if (strstr(array(passwd_entries, i, char*), ACEGROUP_ADMIN) == array(passwd_entries, i, char*)
	      || strstr(array(passwd_entries, i, char*), ACEGROUP_WRITE) == array(passwd_entries, i, char*)
	      || strstr(array(passwd_entries, i, char*), ACEGROUP_READ) == array(passwd_entries, i, char*))
	    {
	      int num_words ;
	      Array words ;
	      
	      group_found = TRUE ;
	      
	      /* Check to see if the userid is in this group.                    */
	      num_words = parseWords(array(passwd_entries, i, char*), &words) ;
	      if (num_words == 1 || stringInWords(words, userid) < 0)
		{
		  userid_found = FALSE ;
		}
	      else
		{
		  char *tmp = array(passwd_entries, i, char*) ;
		  int j ;

		  userid_found = TRUE ;
		  
		  array(passwd_entries, i, char*) = messalloc(strlen(tmp) + 1 - strlen(userid)) ;
							    /* Don't need space for removed userid */
		  
		  /* copy the group entry missing out userid to be deleted.         */
		  j = 0 ;
		  if (strcpy(array(passwd_entries, i, char*), array(words, j, char*)) == NULL)
		    messcrash("strcpy() failed while trying to construct "
			      "group entry for passwd file.") ;
		  for (j = 1 ; j < num_words ; j++)
		    {
		      if (strcmp(array(words, j, char*), userid) != 0)
			{
			  if (strcat(array(passwd_entries, i, char *), " ") == NULL
			      || strcat(array(passwd_entries, i, char*), array(words, j, char*)) == NULL)
			    messcrash("strcat() failed while trying to construct "
				      "group entry for passwd file.") ;
			}
		    }
		  
		  messfree(tmp) ;
		}
	      if (num_words > 0)
		freeWords(words) ;
	    }
	}
    }

  if (userid_found)
    status = TRUE ;
  else if (group_found == FALSE || userid_found == FALSE)
    status = FALSE ;
  /* Should set an error message here....about not finding group..           */
  /* Actually this should be a messcrash really....                          */


  return status ;
}


/************************************************************/

/* Change the users group in the password file. We assume that the groups    */
/* containing the users entry have this format:                              */
/*             <groupname>: user1 user2 user3                                */
/*             |                                                             */
/*              --- beginning of line.                                       */
/*                                                                           */
/* where <groupname> is one of  read|write|admin                             */
/*                                                                           */
/* If there is white space in front of the group then the code will fail.    */
/*                                                                           */
static BOOL addEntryToGroup(Array passwd_entries, char *userid, char *group)
{
  BOOL status, group_found, userid_found ;
  int i ;

  /* Find the group entry, delete and add the new one.                       */
  for (i = 0, status = FALSE, group_found = FALSE, userid_found = FALSE ;
       i < arrayMax(passwd_entries) && status == FALSE ; i++)
    {
      /* We may have set some elements to NULL because of not being able */
      /* to use arrayRemove (see removeUser() ), so don't print these.   */
      if (array(passwd_entries, i, char*) != NULL)
	{
	  /* Find the group.                                                     */
	  if (strstr(array(passwd_entries, i, char*), group) == array(passwd_entries, i, char*))
	    {
	      int num_words ;
	      Array words ;

	      group_found = TRUE ;
	  
	      /* If there are any userids in the group, check that the new       */
	      /* userid is not the same as any of them. If its not, then add it  */
	      /* to the group.                                                   */
	      num_words = parseWords(array(passwd_entries, i, char*), &words) ;
	      if (num_words > 1 && stringInWords(words, userid) >= 0)
		{
		  userid_found = TRUE ;
		  status = TRUE ;
		  /* should set an error return here....                         */
		}
	      else
		{
		  char *tmp = array(passwd_entries, i, char*) ;

		  array(passwd_entries, i, char*) = messalloc(strlen(tmp) + 1 + strlen(userid) + 1) ;
		  if (strcpy(array(passwd_entries, i, char*), tmp) == NULL
		      || strcat(array(passwd_entries, i, char *), " ") == NULL
		      || strcat(array(passwd_entries, i, char*), userid) == NULL)
		    messcrash("Standard str func failed while trying to construct "
			      "userid entry for passwd file.") ;
		  messfree(tmp) ;
		  status = TRUE ;
		}
	      freeWords(words) ;
	    }
	}
    }

  if (group_found == FALSE || userid_found == TRUE)
    status = FALSE ;

  return status ;
}


/************************************************************/

/* This mechanism could be made more general for lots of "packages" so we    */
/* had message catalogs but we should try and use stuff like from the        */
/* "Windows Error Messages" book.                                            */
/*                                                                           */
/* Why a message function ? because passwords are likely to cause trouble    */
/* and we need to log problems without necessarily reporting the exact       */
/* problem to the user.                                                      */
/*                                                                           */
char *getPasswdErrorMessage(ClientAccessStatus msg_num)
{
  typedef struct _PasswdMessages {ClientAccessStatus msg_ID ; char *msg_text ;} PasswdMessages ;
  const PasswdMessages msg_catalog[] = 
  {
    {CLIENTACCESS_FIRST, NULL},
    {CLIENTACCESS_UNKNOWN, "Client userid not found in database password file and database has no \"world\" permissions"},
    {CLIENTACCESS_BADHASH, "Client userid was found but client password hash was wrong"},
    {CLIENTACCESS_BLOCKED, "Client userid was found but database global access permissions have blocked clients access"},
    {CLIENTACCESS_ADMINHOSTNAME, "Error in password file, machine host names are not allowed for the admin group"},
    {CLIENTACCESS_BADUSER_REQUEST, "Request to update user is incorrect, request was: %s"},
    {CLIENTACCESS_BADDOMAIN_REQUEST, "Request to update domain is incorrect, request was: %s"},
    {CLIENTACCESS_LAST, NULL}
  } ;
  enum {MSG_TOTAL = UtArraySize(msg_catalog)} ;
  BOOL found ;
  int i, index ;

  if (msg_num <= CLIENTACCESS_FIRST || msg_num >= CLIENTACCESS_LAST)
    messcrash("bad logic in call to error message finding routine, "
	      "message number was: %d", msg_num) ;

  for (i = 0, found = FALSE ; i < MSG_TOTAL && found == FALSE ; i++)
    {
      if (msg_catalog[i].msg_ID == msg_num)
	{
	  found = TRUE ;
	  index = i ;
	}
    }

  if (found == FALSE)
    messcrash("bad logic in call to error message finding routine, "
	      "message number was: %d", msg_num) ;

  return (msg_catalog[index].msg_text) ;
}


/************************************************************/
/************************************************************/


