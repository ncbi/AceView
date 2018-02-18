/*  File: acecliservutils.h
 *  Author: Ed Griffiths (edgrif@sanger.ac.uk)
 *  Copyright (c) J Thierry-Mieg and R Durbin, 1999
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.crbm.cnrs-mop.fr
 *
 * Description: Interface for client/server utilities.
 * HISTORY:
 * Last edited: Feb 11 12:06 2000 (edgrif)
 * Created: Wed Nov 17 15:02:22 1999 (edgrif)
 * CVS info:   $Id: serverclientutils.h,v 1.1.1.1 2002/07/19 20:23:35 sienkiew Exp $
 *-------------------------------------------------------------------
 */
#ifndef DEF_ACECLIENTSERV_UTILS_H
#define DEF_ACECLIENTSERV_UTILS_H



/* The MD5 code returns an array of unsigned char of size 16, the value 16   */
/* has no symbolic constant in md5.h so I define one, plus the size of the   */
/* string required to hold the hexadecimal string version of the array.      */
/*                                                                           */
enum {MD5_HASHLEN = 16, MD5_HEX_HASHLEN = ((16 * 2) + 1)} ;

char *hashAndHexStrings(char *strings[], int num_strings) ;


/* Noddy functions to set/test the string forms of the message types.        */
void setMsgType(char buffer[], char *msgType) ;
BOOL testMsgType(char buffer[], char *msgType) ;




#endif /* DEF_ACECLIENTSERV_UTILS_H */
