/*  File: sclientlib.h
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
 * Description: Library used by all clients for the ace server.
 *              
 * HISTORY:
 * Last edited: Feb  5 16:15 2001 (edgrif)
 * Created: Wed Apr 26 16:23:47 2000 (edgrif)
 * CVS info:   $Id: sclientlib.h,v 1.1.1.1 2002/07/19 20:23:35 sienkiew Exp $
 *-------------------------------------------------------------------
 */
#ifndef DEF_CLIENTLIB_H
#define DEF_CLIENTLIB_H


/* Client information we want to hold onto between transactions.             */
typedef struct _ClientStruct
{
  char *userid ;
  char *nonce ;
  char *passwd_hash ;
  char *nonce_hash ;
} ClientStruct, *Client ;


/* Handling connection and communication with the server.                    */
/*                                                                           */
/* Open the connection to the server, userid/passwd passed in via client.    */
void *openServer(char *host, int port, int timeOut, int *sock_ret, Client client) ;

/* Send a request to the server, wait for and return the servers reply.      */
int askServerBinary (void *handle, char *request, char **answerp, 
		     int *lengthp, int *encorep, int chunkSize, int time_out) ;

/* Terminate connection to the server.                                       */
void closeServer (void *handle) ;


/* Utilities.                                                                */

/* Return the higher numbered socket, avoid mucky macro.                     */
int getMaxSocket(int sock1, int sock2) ;

/* Glue together an array of words in to a single string.                    */
char *makePhrase(char *words[], int num_words) ;



#define DEFAULT_PORT 23100



#endif /* DEF_CLIENTLIB_H */
