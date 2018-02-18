/*  File: socketace_.h
 *  Author: Ed Griffiths (edgrif@sanger.ac.uk)
 *  Copyright (c) J Thierry-Mieg and R Durbin, 1999
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.crbm.cnrs-mop.fr
 *
 * Description: Internal header for the ACEDB server code. This code
 *              has been written for the sockets based version of the
 *              client/server but in theory could be made to work with
 *              other communication protocols.
 * HISTORY:
 * Last edited: Apr  5 19:16 2001 (edgrif)
 * Created: Fri Nov 12 10:53:40 1999 (edgrif)
 * CVS info:   $Id: serverace_.h,v 1.1.1.1 2002/07/19 20:23:35 sienkiew Exp $
 *-------------------------------------------------------------------
 */
#ifndef DEF_ACESERVER_PRIV_H
#define DEF_ACESERVER_PRIV_H

#include <wsocket/servertransport.h>			    /* Defines server <-> transport layer. */



/* wspec/serverconfig.wrm keywords, configure the server.                    */
/* ----------------------                                                    */
/*                                                                           */


#define ACESERV_CONFIG_NAME      "serverconfig"
#define ACESERV_CONFIG_EXT       "wrm"

/* switch on debugging output from server.                                   */
#define ACESERV_DEBUG     "DEBUG"

/* database is read only, so don't try to use log files etc. in database.    */
#define ACESERV_READONLY_DATABASE  "READONLY_DATABASE"

/* working directory of server, if not set then defaults to database dir.    */
#define ACESERV_WORKING_DIRECTORY  "WORKING_DIRECTORY"

/* switch off host domain name resolution to allow any user from a domain    */
/* specified in the passwd file to log in. This is switchable because        */
/* it can take a long time to execute, it could well be a security hole and  */
/* it often fails because the tables for IP address to domain resolution are */
/* not up to date.                                                           */
#define ACESERV_NODOMAINRESOLVE  "NO_HOSTNAME_RESOLUTION"


/* The database can be locked so that the server cannot be started on it for */
/* a number of reasons, e.g. the administrator shutdown the database to      */
/* prevent possible corruption. This locking is implemented via "lock" files */
/* which are written into the /database directory and if present prevent the */
/* acedb server from starting up. The name of the lock file shows the reason */
/* the database is locked. The database cannot be restarted until the        */
/* lockfile is manually removed.                                             */
/* N.B. these file names will have the usual ".wrm" extension added.         */

#define ACESERV_NORESTART     "NO_RESTART"
#define ACESERV_NORESTART_OPT "norestart"

#define ACESERV_LOCKFILE_DIR      "database"
#define ACESERV_LOCKFILE_EXT      "wrm"
#define ACESERV_SHUTDOWN_LOCKFILE       "SERVER_SHUTDOWN_LOCKFILE"
#define ACESERV_CRASH_LOCKFILE          "SERVER_CRASH_LOCKFILE"
#define ACESERV_UNDETECTABLECRASH_LOCKFILE    "SERVER_UNDETECTABLECRASH_LOCKFILE"


#define ACESERV_SHUTDOWN_LOCKFILE_MSG "The database has previously been shutdown " \
                                      "with the \""ACESERV_NORESTART"\" option, "            \
                                      "it cannot be restarted until you remove the lockfile:  %s"

#define ACESERV_CRASH_LOCKFILE_MSG "The database crashed last time it was run, "     \
                                   "but the code detected the crash and will have tried "  \
                                   "to save the database. "  \
                                   "You should check the integrity of the database with " \
                                   "a stand alone program such as tace before restarting " \
                                   "the server. Make sure it is safe to restart "  \
				   "and then remove the lockfile: %s"

#define ACESERV_UNDETECTABLECRASH_LOCKFILE_MSG "The \"Undetectable crash\" lockfile has been " \
                                         "found in the database sub-directory which indicates " \
                                         "that the server may have crashed uncontrollably. " \
                                         "You should check the integrity of the database with " \
                                         "a stand alone program such as tace before restarting " \
                                         "the server. Make sure it is safe to restart "  \
				         "and then remove the lockfile: %s"

#define ACESERV_CRASH_INETD_MSG "Note that the server was inetd started and may " \
                                "repeatedly restart and crash because inetd detects new clients " \
                                "trying to connect and automatically attempts to restart the server."

/* GLOBAL READ/WRITE database access permissions.                            */
/*                                                                           */
/* WRITE and READ are keywords which may both have one of the NONE, PASSWD   */
/* or WORLD options.                                                         */
#define ACESERV_WRITE  "WRITE"
#define ACESERV_READ   "READ"

#define ACESERV_NONE   "NONE"
#define ACESERV_PASSWD "PASSWD"
#define ACESERV_WORLD  "WORLD"
typedef enum _AceServPermOptions {PERM_NONE, PERM_PASSWD, PERM_WORLD} AceServPermOptions ;




/* wspec/serverpasswd.wrm keywords, user passwords/groups.                   */
/* ----------------------                                                    */
/*                                                                           */

#define ACESERV_PASSWD_NAME      "serverpasswd"
#define ACESERV_PASSWD_EXT       "wrm"


/* Group membership of users.                                                */
/* The below nonsense with having to add the ":" separately is because DECs  */
/* strncmp interface is rubbish, it returns zero (i.e. OK) if a comparison   */
/* ends prematurely...sigh. So we have to use strcmp..see updateUser.        */
#define ACEGROUP_ADMIN_TEXT "admin"
#define ACEGROUP_WRITE_TEXT "write"
#define ACEGROUP_READ_TEXT  "read"
#define ACEGROUP_ADMIN ACEGROUP_ADMIN_TEXT":"
#define ACEGROUP_WRITE ACEGROUP_WRITE_TEXT":"
#define ACEGROUP_READ  ACEGROUP_READ_TEXT":"



/* CLIENT data                                                               */
/* The client structure, describes access permissions etc. for a single      */
/* client.                                                                   */
/*                                                                           */

/* Describes where client is in its life cycle, e.g. server is waiting for   */
/* password from client (n.b. init and killed states are not needed because  */
/* client is allocated/deallocated immediately).                             */
/* ACECLIENT_WAITPASSWD => Server waiting for client to send password.       */
/* ACECLIENT_ACTIVE     => Client has registered correctly.                  */
typedef enum _AceClientState {ACECLIENT_WAITPASSWD, ACECLIENT_ACTIVE} AceClientState ;

/* The level of the clients access to the database.                          */
typedef enum _AceUserGroup {ACEUSER_NONE, ACEUSER_READ, ACEUSER_WRITE, ACEUSER_ADMIN} AceUserGroup ;

/* Error message stuff...                                                    */
/* Reasons why a client can be refused access to the database.               */
/* (keep enums and strings in step because enums index the message array.    */
/* SEE the function  getPasswdErrorMessage(), this is rubbish and will be    */
/* improved.                                                                 */
/*                                                                           */
typedef enum _ClientAccessStatus {CLIENTACCESS_FIRST = 0,
				  CLIENTACCESS_UNKNOWN, CLIENTACCESS_BADHASH,
				  CLIENTACCESS_BLOCKED, CLIENTACCESS_ADMINHOSTNAME,
				  CLIENTACCESS_BADUSER_REQUEST, CLIENTACCESS_BADDOMAIN_REQUEST,
				  CLIENTACCESS_LAST} ClientAccessStatus ;





/* Definition of an ace client.                                              */
typedef struct _AceClientStruct
{
  int clientId ;
  AceClientState state ;
  mytime_t lastAccess;
  AceCommand aceCommand ;
  KEY lastCommand ;
  AceServPermOptions globalReadAcc, globalWriteAcc ;
  char *userid ;
  char *nonce ;
  AceUserGroup group ;
  ClientConnect connection ;
  ClientConnectionID connection_id ;
  BOOL compress ;
} AceClientStruct, *AceClient ;


/*                                                                           */
/* Server state is kept in this control block which is passed around rather  */
/* than keeping a load of statics.                                           */
/* Also include here some defaults for various bits of state, if you alter   */
/* any of these you should change the corresponding bit in the usage string  */
/* defined below.                                                            */
/*                                                                           */
/* Notes:                                                                    */
/*                                                                           */
/* We hold the name of the crash lock file statically because we want to     */
/* construct it at database startup when memory should be plentiful, not     */
/* when we are crashing and don't know what state we are in.                 */
/*                                                                           */
/* OK to have port here, it is not specific to sockets, rpc, udp etc. also   */
/* use port nos.                                                             */

enum {ACESERV_PORT             = 23100,			    /* should be > 49151 */
      ACESERV_CLIENTTIMEOUT    = 600,			    /* secs, infinite if zero */
      ACESERV_SERVERTIMEOUT    = 600,			    /* secs, infinite if zero */
      ACESERV_AUTOSAVEINTERVAL = 600			    /* secs */
} ;

typedef struct _AceServStruct
{
  BOOL debug ;
  ACEOUT server_log ;
  char *dbpath ;
  char *working_dir ;
  BOOL inetd_started ;
  BOOL readonly_DB ;
  int port ;
  int clientTimeOut ;
  int serverTimeOut ;
  int maxKbytes ;
  int autoSaveInterval ;
  int nActiveClients; 
  BOOL shuttingDown ;
  Associator clientAss ;
  char *crashLockFile ;
  char *undetectable_lockfile ;
  BOOL crashfile_found_at_startup ;
  Array wspecTar ;
  Array wspecTarEncoded ;
} AceServStruct, *AceServ ;


/* Usage string for help with wrong args to server prog., see defines above  */
/* for defaults for port etc.                                                */
/*                                                                           */
/* reason for lack of -option notation is that inetd only allows 5 params in some 
   implementations 
*/
#define USAGE_STR "Usage: saceserver database_dir [port_number [params]]\n"                 \
		  "  params are clientTimeout:serverTimeout:maxKbytes:autoSaveInterval\n"   \
		  "  defaults are port             = 23100 \n"                              \
		  "               clientTimeout    = 600 seconds\n"                         \
		  "               serverTimeout    = 600 seconds\n"                         \
		  "               maxKbytes        = 0 [no limit]\n"                        \
		  "               autoSaveInterval = 600 seconds\n"                         \
		  "  autoSave checks are only certain at serverTimeout intervals\n"         \
		  "  example:  aceserver /local1/worm  23100  1200:1200:100\n"



/* Functions private to the acedb server layer but shared between server     */
/* layer files.                                                              */
int parseWords(char *request, Array *word_array_out) ;
int stringInWords(Array words, char *string) ;
void freeWords(Array word_array) ;
ACEIN findAndOpenFile(AceServ server, char *filename, char **filepath_out) ;
ACEIN findAndOpenPasswdFile(AceServ server, char **path_out) ;
BOOL checkUserHash(AceServ server, char *userid, char *hash, char *nonce) ;
BOOL checkGlobalDBAccess(AceServPermOptions *read_acc, AceServPermOptions *write_acc) ;
char *makeNonce(void *client_addr) ;
BOOL parseUseridHash(char *request, char **userid_out, char **hash_out) ;
AceUserGroup getUserDBAccess(AceServ server,
			     char *userid, char *hash, char *nonce,
			     AceServPermOptions globalReadAcc,
			     AceServPermOptions globalWriteAcc,
			     ClientConnect connection,
			     ClientAccessStatus *status_out) ;
char *getPasswdErrorMessage(ClientAccessStatus msg_num) ;
BOOL changeUserPasswd(AceServ server, char *userid, char *nonce, char *old_hash, char *new_hash,
		      ClientAccessStatus *status_out) ;
BOOL updateUser(AceServ server, Array user_text, ClientAccessStatus *status_out) ;
BOOL updateDomain(AceServ server, Array words, ClientAccessStatus *status_out) ;
BOOL updateGlobalDBPerms(AceServ server, char *perm_type, char *perm_level) ;
BOOL checkForServerConfigValue(AceServ server, char *keyword) ;
BOOL getUserUpdate(char **cp_out) ;

#endif /* DEF_ACESERVER_PRIV_H */
