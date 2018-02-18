/*  File: strsubs.c
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
 * Description: String utilities, we do not have some basic stuff for
 *              strings. Some of these maybe simple covers to the
 *              utility libraries (e.g. glib) in the end.
 *              This is thread safe so far, please try to keep it so.
 *              
 * Exported functions: See strsubs.h
 * HISTORY:
 * Last edited: Jul 11 15:52 2000 (edgrif)
 * Created: Tue Jul 11 14:07:45 2000 (edgrif)
 * CVS info:   $Id: strsubs.c,v 1.2 2010/11/22 00:49:13 mieg Exp $
 *-------------------------------------------------------------------
 */

#include <wh/regular.h>
#include <wh/aceio.h>
#include <wh/strsubs.h>



/* Return the number of "space" delimited words in a string, currently this  */
/* means space as determined by aceInWord(). This could be made more         */
/* efficient by using strtok or even cruising through the string one char    */
/* at a time, but this will do for now.                                      */
/*                                                                           */
int strNumWords(char *user_string)
{
  int num_words = 0 ;
  ACEIN string = NULL ;
  char *next_word ;

  string = aceInCreateFromText(user_string, "", 0) ;

  if (string == NULL)
    messcrash("Could not create aceIn from text: %s", user_string) ;

  if (aceInCard(string))
    {
      while ((next_word = aceInWord(string)))
	{
	  num_words++ ;
	}
    }

  aceInDestroy(string) ;

  return num_words ;
}





