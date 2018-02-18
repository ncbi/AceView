/*  File: jadeserver.c
 *  Author: J D Barnett, Lincoln Stein  and Jean Thierry-Mieg
 *  Copyright (C) L Stein and J Thierry-Mieg, 1996
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
 * Runs tace in the back with a netclient dialog on stdin/out
 * Effectivelly allowing Jade to talk to a standalone tace
 * 
 * Exported functions:
 * HISTORY:
 * Last edited: Dec  2 16:34 1998 (edgrif)
 * * Dec  2 16:34 1998 (edgrif): Corrected decl. of main, added code to
 *              record build time of this module.
 * * Oct 16 17:27 1998 (fw): use messSysErrorText instead of strerror
 * * Jul  7 10:24 1998 (edgrif): Replace references to sys_errlist with
 *                               strerror function.
 * Created: Wed Nov 25 20:02:45 1992 (j barnett)
 *-------------------------------------------------------------------
 */

/* $Id: jadeserver.c,v 1.2 2003/07/07 15:50:47 mieg Exp $ % */

#include <unistd.h>
#include <string.h>
#include <sys/types.h>
#include <ctype.h>
#include <errno.h>
#include <malloc.h>
#include "regular.h"
#include "array.h"
#include "aceclient.h"
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
#define MODEL "model"
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

/*******************************************************************/
/********** Missing routines ***************************************/

void graphCleanUp(void)
{ return ;
}

void graphText (char *text, float x, float y)
{ printf(text) ;
}

void graphHelp(void)
{ return ;
}
void graphClear(void)
{ printf("\n\n") ;
  return ;
}

int graphExists(void)
{ return 0 ;
}

void graphRedraw(void)
{ printf("\n") ;
  return ;
}

void graphTextBounds(int w, int h)
{ return ;
}

void graphColor(int color)
{ return ;
}

void * graphCreate()
{ return NULL ;
}

BOOL graphInterruptCalled (void)
{ return FALSE ;
}

BOOL xClientGetKey (KEY key)
{ return TRUE ;
}

/*  end of graph stubs */


static chunkSize = 10;
/* Opaque database handle JC do not need to know what it is */
static Stack mainStack;
/* static KEY option = 0; */

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
  while ((i < COMMAND_MAX-1) && line[i] && !isalnum(line[i]))
    i++;


  while (i < COMMAND_MAX-1 && line[i] && !isspace(nextchar = line[i])) {
    command[j++] = nextchar;
    i++;
  }

  command[j]='\0';
  *extraData = isspace(line[i]) ? line + i + 1 : line + i;

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
void doQuery(void* commandLook, char* query,BOOL show) {
  char message[MAX_MESSAGE];
  char subclass[MAX_CLASSNAME],superclass[MAX_CLASSNAME];
  char* answer;
  int retval, length, active = 0, level;
  BOOL encoring = FALSE;

    level = freesettext (query,"");
    mainStack = stackReCreate(mainStack, 4000);
    aceCommandDoExecute (commandLook, level, FALSE, 0, 0);
    answer = stackText (mainStack, 0);

    /*  What to do about retval??? */
    /*      retval = askServerBinary(handle,query,&answer,&length,&encore,chunkSize);

    if (retval > 0) {
      sprintf(message,"ACE: error code %d",retval);
      writeStatus(TIMEOUTERROR,message);
      if (answer != NULL) free(answer);
      break;
    }
    */
    if (show) {

      /* Some error detection and reporting.  This has to be much expanded.*/
      *message = '\0'; /* paranoia */
      if (2 == sscanf((char*)answer,"// %s is a sub class of %[a-zA-Z_0-9]",subclass,superclass)) {
	sprintf(message,"%s is the superclass for %s",superclass,subclass);
	writeStatus(REDIRECT,message);
	return;
      } else if (1 == sscanf((char*)answer,"// Sorry, %[a-z A-Z_0-9:]",message)) {
	writeStatus(SYNTAXERROR,message);
	return;
      } else if (!strncmp((char*)answer,"// This table has a problem",22)) {
	writeStatus(SYNTAXERROR,"table syntax error");
	return;
      }

      writeStatus(OK,"Results follow, terminated by \".\"");
      stripAndPrint((char*)answer,strlen(answer),FALSE,FALSE, 0);
    }
    else
      if ((active = countActive((char*)answer)) >= 0) {
	sprintf(message,"%d objects",active);
	writeStatus(OK,message);
      }
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
void doTable( void* commandLook)
{
  char* tempName;
  char message[FILENAME_MAX + COMMAND_MAX + 1];
  char* line;
  FILE* temp;

  if (!(tempName = tempnam("/var/tmp", ".acenetcl"))) {
    sprintf(message,"Error getting temp file name: %s", 
	    messSysErrorText()) ;
    writeStatus(FILEERROR,message);
    return;
  }
  unlink(tempName); /* just to be safe */
  temp = fopen(tempName,"w");

  if (temp == NULL) {
    sprintf(message,"Error creating temp table file \"%s\":%s", 
	    tempName, messSysErrorText()) ;
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
    doQuery(commandLook, message,TRUE);
  }

  unlink(tempName);
  free((void*) tempName);
}

/*************************************************************************/
void printHelp () {
  char* helpText[] = {
    "- Commands:",
    "-    FIND   SHOW   FOLLOW   TABLE  NTABLE MODEL",
    "-    LIST     IS   QUERY     UNDO    QUIT  HELP",
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
BOOL processCommand(void* commandLook, char* command,char* data) {

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
  else if (streq(command,MODEL)) {
    strcat(message,"model ");
    show = TRUE;
  }
  else if (streq(command,LIST)) {
    strcat(message,"list -j ");
    pCheck = FALSE;
    show = TRUE;
  }
  else if (streq(command,TABLE))
    doTable(commandLook);
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
    doQuery(commandLook, message,show);
  }

  return FALSE;
}


/* Defines a routine to return the compile date of this file.                */
UT_MAKE_GETCOMPILEDATEROUTINE()


/*************************************************************************/
int main(int argc, char *argv[] )
{
  char *host = "localhost";
  unsigned long port = DEFAULT_PORT;
  int timeOut = 300;

  BOOL isQuit = FALSE;
  char *command,*extraData ;

  void *commandLook;

  sessionInit (argc>1 ? argv[1]: 0) ;

  mainStack = stackCreate(4000);
  commandLook =( void*)aceCommandDoCreate (1, NULL, mainStack);

  writeStatus(OK,"Ace server ready and waiting");
  while (!isQuit && (command = readCommand(&extraData)) )
    isQuit = processCommand(commandLook, command,extraData);

  writeStatus(GOODBYE,"A bientot");

  return(EXIT_SUCCESS) ;
}

/*************************************************************************/
/*************************************************************************/
