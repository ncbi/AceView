/*  File: liste.c
 *  Author: Richard Durbin (rd@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1995
 * -------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 * -------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: 
     You can add and remove from a liste in a loop
     without making the list grow larger then the max current number
     this is more economic than hashing, where removing is inneficient
     and faster than an ordered set

     This library will not check for doubles,
      i.e. it maintains a list, not a set.
 * Exported functions:
     Liste listeCreate (AC_HANDLE hh) ;
     listeDestroy(_ll) 
     listeMax(_ll) 
     int listeFind (Liste liste, void *vp)  ;
     int listeAdd (Liste liste, void *vp) ;
     void listeRemove (Liste liste, void *vp, int i)  ;

 * HISTORY:
 * Last edited: Nov 23 11:21 1998 (fw)
 * Created: oct 97 (mieg)
 *-------------------------------------------------------------------
 */

/* $Id: liste.c,v 1.7 2009/06/22 19:33:13 mieg Exp $ */

#include "liste.h"

/* A liste is a list of void*
 * It is less costly than an associator if
 * the things in it are transient
 * because in an associator, you cannot easilly destroy
 * The liste does NOT check for doubles.
 */

static int LISTEMAGIC = 0 ;

static void listeFinalize(void *vp)
{
  Liste liste = (Liste) vp ;
  /* we need to destroy in case liste was allocated on hh = 0 
   */
  if (liste && liste->magic == &LISTEMAGIC)
    {
      arrayDestroy(liste->a) ;
      liste->magic = 0 ;
    }
}

Liste listeCreate (AC_HANDLE hh)
{
  Liste liste = (Liste)  halloc(sizeof(struct listeStruct), hh) ;
  blockSetFinalise(liste,listeFinalize) ;
  
  liste->magic = &LISTEMAGIC ;
  liste->i = 1 ;

  /* NOTE: the array is a freefloating allocation and must not be 
   * allocated upon the given handle - 
   * it is destroyed explicitly by finalisation */
  liste->a = arrayCreate (32, void*) ;


  array (liste->a,0,void*) = 0 ; /* avoid zero */

  return liste ;
} /* listeCreate */

void listeRemove (Liste liste, void *vp, int i) 
{
  Array a = liste->a ;
  void *wp = array(a, i, void*) ;

  if (vp != wp) messcrash ("Confusion in listeRemove") ;
  array (a,i,void*) = 0 ;
  if (i < liste->i)
    liste->i = i ;
}
  
int listeAdd (Liste liste, void *vp) 
{
  Array a = liste->a ;
  int i = liste->i ;
  void **vpp = arrayp(a, i, void*) ;
  int n = arrayMax(a) ;  /* comes after arrayp of above line ! */
  while (*vpp && i++ < n) vpp++ ;
  array (a,i,void*) = vp ;
  return i ;
} /* listeAdd */
  
int listeFind (Liste liste, void *vp) 
{
  Array a = liste->a ;
  int i = arrayMax (liste->a) ;
  void **vpp = arrayp(a, i - 1, void*) + 1 ;
  
  while (vpp--, i--)
    if (vp == *vpp) return i ;
  return 0 ;
} /* listeFind */
  
/******************** end of file **********************/

 
 
