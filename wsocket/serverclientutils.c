/*  File: acecliservutils.c
 *  Author: Ed Griffiths (edgrif@sanger.ac.uk)
 *  Copyright (c) J Thierry-Mieg and R Durbin, 1999
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.crbm.cnrs-mop.fr
 *
 * Description: Provides utility routines required by both the client
 *              and the server.
 * Exported functions: See clientservutils.h
 * HISTORY:
 * Last edited: Dec  6 10:04 2000 (edgrif)
 * Created: Wed Nov 17 15:00:47 1999 (edgrif)
 * CVS info:   $Id: serverclientutils.c,v 1.1.1.1 2002/07/19 20:23:35 sienkiew Exp $
 *-------------------------------------------------------------------
 */

#include <wh/regular.h>
#include <wmd5/global.h>
#include <wmd5/md5.h>
#include <wsocket/serverclientutils.h>
#include <wsocket/servertransport.h>

static char *convertMD5toHexStr(unsigned char digest[]) ;



/*****************************************************************************/
/* Encryption and hashing code.                                              */
/*                                                                           */
/* We use MD5 for encryption and this can be found in the wmd5 directory.    */
/*                                                                           */
/* Note that MD5 produces as its output a 128 bit value that uniquely        */
/* represents the input strings. We convert this output value into a hex     */
/* string version of the 128 bit value for several reasons:                  */
/*                                                                           */
/*     - the md5 algorithm requires strings as input and we need to use      */
/*       some of the md5 output as input into a new md5 hash.                */
/*     - it makes all handling of the encrypted data much simpler            */
/*     - using hex means that the string we produce consists entirely of the */
/*       digits 0-9 and the letters a-f (i.e. no unprintable chars)          */
/*     - the passwd hash can be kept in a plain text file in the database    */
/*                                                                           */


/* Takes an array of strings, hashes then together using MD5 and then        */
/* produces a hexstring translation of the MD5 output.                       */
/*                                                                           */
char *hashAndHexStrings(char *strings[], int num_strings)
{
  char *hex_string = NULL ;
  MD5_CTX Md5Ctx ;
  unsigned char digest[MD5_HASHLEN] ;
  int i ;

  MD5Init(&Md5Ctx) ;

  for (i = 0 ; i < num_strings ; i++)
    {
      MD5Update(&Md5Ctx, (unsigned char *)strings[i], strlen(strings[i])) ;
    }

  MD5Final(digest, &Md5Ctx) ;

  hex_string = convertMD5toHexStr(&digest[0]) ;

  return hex_string ;
}


/************************************************************/

/* Takes the array of unsigned char output by the MD5 routines and makes a   */
/* hexadecimal version of it as a null-terminated string.                    */
/*                                                                           */
static char *convertMD5toHexStr(unsigned char digest[])
{
  char *digest_str ;
  int i ;
  char *hex_ptr ;

  digest_str = halloc(MD5_HEX_HASHLEN, 0) ;
  for (i = 0, hex_ptr = digest_str ; i < MD5_HASHLEN ; i++, hex_ptr+=2)
    {
      sprintf(hex_ptr, "%02x", digest[i]) ;
    }

  return digest_str ;
}



/* Instead of using an enum for the message type we use a string, this will  */
/* hopefully make the protocol more robust for clients written in other      */
/* languages, e.g. perl, java etc.                                           */
/*                                                                           */
/* These two routines set the type and test it, they could be macros but     */
/* performance is not the problem here.                                      */
/*                                                                           */
/* They cope with caller supplying msgType which itself points to buffer.    */
/*                                                                           */
void setMsgType(char buffer[], char *msgType)
{
  char *msg_copy = NULL ;

  /* Ugly bug here...what if msgType points to buffer ? We'll do belt and    */
  /* braces and clean the buffer anyway.                                     */
  if (msgType == &(buffer[0]))
    {
      msg_copy = strnew(msgType, 0) ;
    }

  /* Reset the message type section of the header to be zeroed so there is   */
  /* no extraneous bumpf if the previous message was longer than this one.   */
  /* This is important for perl and other non-C languages that have to parse */
  /* this bit out of the message buffer.                                     */
  memset(buffer, 0, ACESERV_MSGTYPE_BUFLEN) ;

  if (msg_copy != NULL)
    msgType = msg_copy ;

  if (strcpy(buffer, msgType) == NULL)
    messcrash("copy of message type failed, message was: %s",  msgType) ;

  if (msg_copy != NULL)
    messfree(msg_copy) ;

  return ;
}


BOOL testMsgType(char buffer[], char *msgType)
{
  BOOL result = FALSE ;

  if (msgType == &(buffer[0]))
    result = TRUE ;
  else if (strcmp(buffer, msgType) == 0)
    result = TRUE ;

  return result ;
}


/************************************************************/
/************************************************************/
