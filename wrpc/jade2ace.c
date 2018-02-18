/*  File: jade2ace.c
 *  Author: Lincoln Stein (Whitehead) and Jean Thierry-Mieg
 *  Copyright (C) L Stein and J Thierry-Mieg, 1996
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
 * Modified from aceclient to handle communication
 * with a java client
 * 
 * Exported functions:
 * HISTORY:
 * Last edited: Dec  7 12:56 1998 (fw)
 * * Dec  2 16:33 1998 (edgrif): Corrected decl. of main, added code to
 *              record build time of this module.
 * * Oct 16 16:20 1998 (fw): reverted strerror back to messSysErrorText()
 *                           because it resolves system differences
 *                           at that level (SunOS has no strerror)
 * * Jul  7 10:26 1998 (edgrif): Replace references to sys_errlist with
 *                               strerror function.
 * Created: Wed Nov 25 20:02:45 1992 (mieg)
 *-------------------------------------------------------------------
 */

 /* $Id: jade2ace.c,v 1.2 2017/03/18 15:31:36 mieg Exp $ */

#if !defined(WIN32)
#include <unistd.h>
#endif

#include <string.h>
#include <sys/types.h>
#include <ctype.h>
#include <errno.h>
#include <malloc.h>
#include "regular.h"
#include "array.h"
#include "../wh/aceclient.h"
#include "version.h"

#define OK 200
#define GOODBYE 201
#define DEBUG 301
#define COMMENT 302
#define REDIRECT 303
#define PARAMERROR 401
#define SYNTAXERROR 402
#define UNIMPLEMENTED 403
#define TIMEOUTERROR 501
#define COMMERROR 502
#define MEMERROR 503
#define FILEERROR 504

#define FIND "find"
#define SHOW "show"
#define QUERY "query"
#define UNDO "undo"
#define FOLLOW "follow"
#define TABLE "table"
#define NTABLE "ntable"
#define AQL "AQL"
#define MODEL "model"
#define GIFACE "giface"
#define IS "is"
#define LIST "list"
#define QUIT "quit"
#define HELP "help"

#define COMMENTSTRING "//"
#define ERRORSTRING "//!! "
#define KEYSETSTRING "KeySet"
#define EOT "."

/* I promise never to make my commands larger than this! */
#define COMMAND_MAX 128
/* I promise that Ace class names can't get larger than this! */
#define MAX_CLASSNAME 512
/* I promise never to send a message back to the user larger than this! */
#define MAX_MESSAGE 512

#ifndef FILENAME_MAX
#define FILENAME_MAX 125
#endif

static int chunkSize = 10;
/* Opaque database handle JC do not need to know what it is */
static void *handle; 

void tStatus (void) { return ; }

void writeStatus (int status,char* message) {
  fprintf(stdout,"%d %s\r\n",status,message);
  fflush(stdout);
}

void writeln(char* line) {
  fprintf(stdout,"%s\r\n",line);
  fflush(stdout);
}

void die (int status, char* message) {
  writeStatus(status,message);
  exit(-1);
}

void reallocate (void** buffer, int old_size, int new_size) {
  void* newbuffer;
  int bytes_to_copy;

  bytes_to_copy = (old_size < new_size) ? old_size : new_size;

  newbuffer = malloc(new_size);
  if (newbuffer == NULL)
    die(MEMERROR,"ran out of memory allocating input buffer");
  if (*buffer != NULL) {
    memcpy(newbuffer,*buffer,bytes_to_copy);
    free(*buffer);
  }
  *buffer = newbuffer;
}

char* readln() {
  static char* linebuffer = NULL;
  static int linebufsize = 512;
  int length;
  char* result;
  int oldsize = 0;
  int offset = 0;

  if (linebuffer == NULL) {
    reallocate((void**)&linebuffer,oldsize,linebufsize);
  }

  while (1) {

    result = fgets(linebuffer+offset,linebufsize-offset,stdin);
    if (result == NULL)
      return NULL;

    length = strlen(linebuffer);

    /* empty string is OK */
    if (length == 0)
      return linebuffer;

    /* huge buffers present a security risk */
    if (length > (1<<20))
      { free (linebuffer) ; return NULL ; }


    /* if the last character is a \n, then we got the end of the line
       so we're done */
    if (linebuffer[length-1] == '\n') {
      if ( (length >= 2) && (linebuffer[length-2] == '\r') )
	linebuffer[length-2] = '\0';
      else
	linebuffer[length-1] = '\0';
      return linebuffer;
    }

    /* otherwise we got a truncated, and we need to make the buffer bigger
     and try again */
    oldsize = linebufsize;
    linebufsize *= 2;
    reallocate((void**)&linebuffer,oldsize,linebufsize);
    offset = length;
  }
       
}

char* readCommand(char** extraData) {
  static char command[COMMAND_MAX]; /* command can't be longer than 512 characters */
  char *line;
  char nextchar;
  int i=0,j=0;

  if ( (line = readln()) == NULL)
    return NULL;

  /* Skip over whitespace and other non-alphanumeric chars at the beginning of the line.  */
  while ((i < COMMAND_MAX-1) && line[i] && !isalnum((int)line[i]))
    i++;


  while (i < COMMAND_MAX-1 && line[i] && !isspace((int)(nextchar = line[i]))) {
    command[j++] = nextchar;
    i++;
  }

  command[j]='\0';
  *extraData = isspace((int)line[i]) ? line + i + 1 : line + i;

  return command;
}

BOOL streq (char* s1, char* s2) {
  int length;
  length = strlen(s2);
  if (strncasecmp(s1,s2,length) == 0)
    return TRUE;
  return FALSE;
}

/*************************************************************************/
/* Find the "%d Active Objects" line and return the count */
/*************************************************************************/
int countActive (char* text) {
  int objectCount;
  while (text) {
    if (sscanf(text,"// %d Active Objects",&objectCount))
      return objectCount;
    /* if we didn't find it, skip to the next line and try again */
    text = strchr(text,'\n');
    if (!text) break;
    if (!*(++text)) break;
  }
  /* if we get here, we never found the string.  Return -1 */
  return -1;
}

/*************************************************************************/
/* replace newlines with CRLF pairs.  Strip // comments entirely. */
/*************************************************************************/
void stripAndPrint (char* text,int length,
		    BOOL stripBlanks,BOOL stripLastBlank, int encore) {
  char* end_of_line;
  char* current = text;
  
  if (encore && *text) /* mieg: because writeln adds an extra \n */
    {
      end_of_line = text + strlen(text) - 1 ;
      if (*end_of_line == '\n') *end_of_line = 0 ;
    }
  
  if (current[0] == '\n') current++;

  while (current < text + length) {

    if ((strncmp(current,COMMENTSTRING, 2) == 0) ||
	(strncmp(current,KEYSETSTRING, 6) == 0) ) 
      {

       if ((end_of_line = strchr(current,'\n')) == NULL)
	 return;

       current = end_of_line + 1;
       continue;
    }

    end_of_line = strchr(current,'\n');

    /* sorry, this is a little weird because the number of lines that come back
       from the ace server is not entirely predictable */
    if (end_of_line != NULL) {

      *end_of_line = '\0';
      
      if (!stripBlanks || strlen(current) > 0)
	writeln(current);

      current = end_of_line + 1;

    } else {

      if ((!stripLastBlank && !stripBlanks) || (strlen(current) > 0) )
	writeln(current);

      return;

    }

  }

}

/*************************************************************************/
/* All-purpose query handler.  If show is TRUE, then  */
/*************************************************************************/
void doQuery(char* query,BOOL show) {
  char message[MAX_MESSAGE];
  char subclass[MAX_CLASSNAME],superclass[MAX_CLASSNAME];
  unsigned char* answer;
  int retval, length, active = 0, encore = 0;
  BOOL encoring = FALSE;

  do {
    retval = askServerBinary(handle,query,&answer,&length,&encore,chunkSize);

    if (retval > 0) {
      sprintf(message,"ACE: error code %d",retval);
      writeStatus(TIMEOUTERROR,message);
      if (answer != NULL) free(answer);
      break;
    }

    if (show && !encoring) {

      /* Some error detection and reporting.  This has to be much expanded.*/
      *message = '\0'; /* paranoia */
      if (2 == sscanf((char*)answer,"// %s is a sub class of %[a-zA-Z_0-9]",subclass,superclass)) {
	sprintf(message,"%s is the superclass for %s",superclass,subclass);
	writeStatus(REDIRECT,message);
	free(answer);
	return;
      } else if (1 == sscanf((char*)answer,"// Sorry, %[a-z A-Z_0-9:]",message)) {
	writeStatus(SYNTAXERROR,message);
	free(answer);
	return;
      } else if (!strncmp((char*)answer,"// This table has a problem",22)) {
	writeStatus(SYNTAXERROR,"table syntax error");
	free(answer);
	return;
      }

      writeStatus(OK,"Results follow, terminated by \".\"");
    }

    encoring = encore;

    if (show)
      stripAndPrint((char*)answer,length,FALSE,FALSE, encore);
    else
      if ((active = countActive((char*)answer)) >= 0) {
	sprintf(message,"%d objects",active);
	writeStatus(OK,message);
	encoring = FALSE;
      }

    free(answer);

    if (!encore)
      break;
    else
      query = "encore";

  } while (encoring);

  if (show) writeln(EOT);
}

/*************************************************************************/
/* Handling of externally-defined tables:
   1. opens up a temporary file using tempnam() function
   2. writes the user's table into that file
   3. asks the server to process the file
   4. prints the result
   5. unlinks the table
   */
#define TMP "/tmp"
void doTable() {
  char* tempName;
  char message[FILENAME_MAX + COMMAND_MAX + 1];
  char* line;
  FILE* temp;

  if (!(tempName = tempnam("/var/tmp", ".acenetcl"))) {
    sprintf(message,"Error getting temp file name: %s", messSysErrorText()) ;
    writeStatus(FILEERROR,message);
    return;
  }
  unlink(tempName); /* just to be safe */
  temp = fopen(tempName,"w");

  if (temp == NULL) {
    sprintf(message,"Error creating temp table file \"%s\": %s", tempName, messSysErrorText()) ;
    writeStatus(FILEERROR,message);
    return;
  }

  /* Here is where we accept lines from the remote and write them into the
     temp file.  This is very, very hacky */
  writeStatus(OK,"Enter table data.  End with \".\" on a line by itself.");
  while ((line = readln()) && (line[0] != '.') ) {
    fprintf(temp,"%s\n",line);
  }
  fclose(temp);

  /* it's possible that the input suddenly died.  In this case we just
     clean up and get out of here */
  if (line != NULL) {
    /* Perform a "table" style query using the current file for input
     This is the obviously gross part. */
    sprintf(message,"table -j %s",tempName);
    doQuery(message,TRUE);
  }

  unlink(tempName);
  free((void*) tempName);
}

/*************************************************************************/
void printHelp () {
  char* helpText[] = {
    "- Commands:",
    "-    FIND   SHOW   FOLLOW    TABLE   NTABLE MODEL   GIF",
    "-    LIST     IS   QUERY     AQL     UNDO    QUIT   HELP",
    "- You need to know the Ace query language to use this",
    "- service.  Sorry.",
    "End of HELP info",
    NULL
  };
  char** m = helpText;
  while (*m != NULL)
    writeStatus(COMMENT,*m++);
}

/*************************************************************************/
void unimplemented (char* command) {
  char message[MAX_MESSAGE];
  sprintf(message,"%s is currently unimplemented",command);
  writeStatus(UNIMPLEMENTED,message);
}

/*************************************************************************/
BOOL processCommand(char* command,char* data) {

  /* the following stuff is a slightly awkward dynamic
     string copy.  It works and is not particularly inefficient
     since the commands are really quite small.*/
  static char* message = NULL;
  static int currentsize = 0;
  BOOL show = FALSE;
  BOOL pCheck = TRUE;

  int sizeneeded = strlen(data) + COMMAND_MAX + 1; 
  if ((message == NULL) || (currentsize < sizeneeded)) {
    reallocate((void**)&message,currentsize,sizeneeded);
    currentsize = sizeneeded;
  }
  *message = '\0';

  if (streq(command,QUIT))
    return TRUE;
  if (streq(command,HELP))
    printHelp();
  else if (streq(command,FIND))
    strcat(message,"find ");
  else if (streq(command,FOLLOW))
    strcat(message,"follow ");
  else if (streq(command,IS))
    strcat(message,"is ");
  else if (streq(command,QUERY))
    strcat(message,"query ");
  else if (streq(command,UNDO)) {
    strcat(message,"undo ");
    pCheck = FALSE;
  }
  else if (streq(command,SHOW)) {
    strcat(message,"show -j ");
    pCheck = FALSE;
    show = TRUE;
  }
  else if (streq(command,NTABLE)) {
    strcat(message,"table -j -n ");
    show = TRUE;
  }
  else if (streq(command,AQL)) {
    strcat(message,"bql -j ");
    show = TRUE;
  }
  else if (streq(command,MODEL)) {
    strcat(message,"model ");
    show = TRUE;
  }
  else if (streq(command,GIFACE)) {
    strcat(message,"gif ");
    show = TRUE;
  }
  else if (streq(command,LIST)) {
    strcat(message,"list -j ");
    pCheck = FALSE;
    show = TRUE;
  }
  else if (streq(command,TABLE))
    doTable();
  else if (!command[0]) /* do nothing for empty commands */
      ;
  else
    writeStatus(SYNTAXERROR,"unknown command");

  if (*message != '\0') {
    if (pCheck && !*data) { 
      writeStatus(PARAMERROR,"missing parameter(s)");
      return FALSE;
    }
    strcat(message,data);
    doQuery(message,show);
  }

  return FALSE;
}

/* Defines a routine to return the compile date of this file.                */
UT_MAKE_GETCOMPILEDATEROUTINE()


/*************************************************************************/
int main(int argc, char *argv[])
{
  char *host = "localhost";
  unsigned long port = DEFAULT_PORT;
  int timeOut = 300;
  BOOL isQuit = FALSE;
  char *command,*extraData ;

  /* Read command line parameters */
  while (argc > 1) {
    argv++; argc--;
    if ( (argc > 1) && !strcmp("-host",*argv) ) {
      argv++; argc--;
      host = *argv;
    } 
    else if ( (argc > 1) && !strcmp("-port",*argv) ) {
      argv++; argc--;
      port = atoi(*argv);
    }
    else if ( (argc > 1) && !strcmp("-time_out",*argv) ) {
      argv++; argc--;
      timeOut = atoi(*argv);
    }
    else {
      fprintf(stderr,"Usage: netclient [-host host] [-port port_num] [-time_out nn_in_seconds]\n");
      return (EXIT_FAILURE);
    }
  }

  if ((handle = openServer(host, port, timeOut)) == NULL) {
    writeStatus(COMMERROR,"cannot establish connection");
    return (EXIT_FAILURE);
  }

  writeStatus(OK,"Ace server ready and waiting");
  while (!isQuit && (command = readCommand(&extraData)) )
    isQuit = processCommand(command,extraData);

  closeServer(handle);
  writeStatus(GOODBYE,"A bientot");
  return(EXIT_SUCCESS);
}

/*************************************************************************/
/*************************************************************************/
 
