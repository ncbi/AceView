/*  File: htmltest.c
 *  Author:  D and J Thierry-Mieg (mieg@ncbi.nlm.nih.gov)
 *  Copyright (C) D and J Thierry-Mieg 2004
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
 * This file is part of the COMPARATOR genome database package, written by
 *	 D and J Thierry-Mieg (mieg@ncbi.nlm.nih.gov)
 *
 *  This code works as a client against an acedb server
 *  and should be linked with the acedb libraries
 *  available from http://www.acedb.org
 *
 *  The present code and the acedb schema for this project
 *  are available from http://www.aceview.org/aligold
 */

#include "vtxt.h"

/***************************************************************************************/
/***************************** Public interface **************************************/
/*************************************************************************************/

static void htmlTest (const char *url)
{   
  vTXT txt ;
  char *buf2 ;
 
  txt =  vtxtGetURL (url, 6) ; /* get the content of the url page */

  buf2 = txt ? vtxtPtr (txt) : "vtxt: page not found" ; 
  printf ("%s\n", buf2) ;

  return  ;
}

/*************************************************************************************/

static void usage (void)
{
  fprintf (stderr, "// Usage: htmltest -url url\n") ;
  exit (1) ;
}

/*************************************************************************************/

int main (int argc, const char **argv)
{
  const char *url = 0 ;

  freeinit () ; /* needed init messalloc and even to link in LINUX_4_OPT */

  /* consume optional args */
  url = getArgV (&argc, argv, "-url") ;

  if (url)
    htmlTest (url) ;
  else
    usage () ;

  return 0 ;
}

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/

