/*  File: get_url.c
 *  Author: Jean Thierry-Mieg (mieg@ncbi.nlm.nih.gov)
 *  Copyright (c) J Thierry-Mieg and R Durbin, 1999
 * -------------------------------------------------------------------
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
 * -------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.crbm.cnrs-mop.fr
 *
 * Description:  Performs one transaction with remote console <strHostName>:<nHostPort>.
	Writes  <writeBuf> content and then waits for transaction until <timeout> 
	or until gets <maxLen> bytes.
	Returns the received buffer. 
 * Exported functions: see wh/vtxt.h
 * HISTORY:
 * Created: Thu May 8 2003: mieg
 *-------------------------------------------------------------------
 */

#include <wh/ac.h>
#include <wh/vtxt.h>
#include <wtcp/tcp_connect.h>
#include <sys/ioctl.h>

/* Adapted from accclinet_actcp.c, originaly by Mark */
static BOOL vtxtUrlTransaction (vTXT txt, char *urlHost, int port, char *urlData, int timeOut, char **errpp) 
{
    int fd = -1 ;
    int i ;
    int bufSize = 64000 ;
    int availBytes, ln = 1 ; /* number of bytes announced/obtained */
    char *buf = (char *) messalloc (bufSize) ;
    BOOL ok = FALSE ;
    char *waitStr1="</HTML>", *waitStr2="</html>" ;

    if (timeOut < 1) 
      timeOut = 1; 
    
    /* make a socket */
    fd = connect_socket (urlHost, port, TRUE) ; 
    if (fd < 0)
      {
	if (errpp)
	  *errpp = "cannot create a socket" ;
	goto abort ;
      }

    /* write and read from socket into buffer  */
    if (urlData && *urlData)
      write (fd, urlData, strlen(urlData)+1);

    while (ln) 
      { 
	/* int xxx = FIONREAD ;  for debugging */
        ioctl (fd, FIONREAD, &availBytes) ;
        for (i = 0; i < timeOut && availBytes < 1; i++)
	  {
	    sleep (1) ;
	    /* ioctl is called ioctlsocket on win32 */
	    ioctl (fd, FIONREAD, &availBytes) ;
	  }
       
        if (availBytes < 1)
	  break;
	
	if (availBytes >= bufSize)
	  {
	    bufSize = availBytes + 1 ;
	    messfree (buf) ;
	    buf = (char *) messalloc (bufSize) ;
	  }
        ln = read (fd, buf, availBytes) ;
	buf [ln] = 0 ; /* insure a final 0 */
	vtxtPrintf (txt, "%s", buf) ;
	if (strstr (buf, waitStr1) || strstr (buf, waitStr2))
	  break ;
      }
        
    ok = TRUE ;

 abort:
    if (fd > 0)
      close (fd) ;
    messfree (buf) ;

  return ok ;
}

/*********************************************************************/
/* Parse the url into its host/port/data components
*/

static BOOL vtxtParseUrl (char *url, char **hostp, int *portp, char **datap)
{
  int x = 0 ;
  char cc = 0, *cp, *cq ;

  cp = strstr (url, "http://") ;
  if (!cp)
    return FALSE ;

  cp += 7 ;
  *hostp = cq = cp ;
  
  /* clip host */
  while (*cp && *cp != ':' && *cp != '/') cp++ ; 
  cc = *cp ; *cp++ = 0 ;

  *portp = 80 ; /* standard HTTP port */
  if (cc == ':') /* get port number */
    {
      x = 0 ;
      cq = cp ; while (*cq && *cq != '/') cq++ ; cc = *cq ; *cq++ = 0 ;
      if (sscanf (cp, "%d", &x) == 1  && x > 0)
	*portp = x ;
      cp = cq ;
    }
  if (cc == '/') /* get local details */
    {
      *datap = (char *) messalloc (strlen(cp) + 1000) ;
      sprintf (*datap
	       , "GET /%s HTTP/1.0 \n"
	       "Accept: text/html, text/plain\n"
	       "User-Agent: xFiche\n"
	       "Host: %s:%d\n"
	       "Connection: Keep-Alive\n\n"
	       , cp, *hostp, *portp) ;
      return TRUE ;
    }
  return FALSE ;
}

/*********************************************************************/
/* 
	downloads the whole ccontent of the url using CGI/POST <data> 
	until <timeOut> seconds (default 60) 
	or </HTML>" or "</html>" is reached
*/


vTXT vtxtHandleGetURL (const char *url, int timeOut, AC_HANDLE h)
{
  vTXT txt = 0 ;
  char *urlCopy = 0, *host, *data, *errp ;
  int port = 80 ;

  if(timeOut <= 0) timeOut = 60 ; 

  if (url && *url)
    {
      urlCopy = strnew (url, 0) ;
      txt = vtxtHandleCreate (h) ;
      
      if (!vtxtParseUrl (urlCopy, &host, &port, &data) ||
	  !vtxtUrlTransaction (txt, host, port, data, timeOut, &errp)) 
	vtxtDestroy (txt) ; 
      messfree (urlCopy) ;
    }

  return txt ;
} /* vtxtHandleGetURL */

/*********************************************************************/

vTXT vtxtGetURL (const char *url, int timeOut)
{
  return 
    vtxtHandleGetURL (url, timeOut, 0) ;
}/* vtxtGetURL */

/*********************************************************************/
/*********************************************************************/
/*********************************************************************/
