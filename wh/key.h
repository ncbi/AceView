/*  File: key.h
 *  Author: Fred Wobus (fw@sanger.ac.uk)
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
 * Description: defines values that represent various keys on keyboard,
 *              these are independent of the underlying windowing system.
 *              The windowing system code must remap its key values to
 *              these values that the application will use.
 * HISTORY:
 * Last edited: May  9 11:44 2003 (rnc)
 * Created: Thu Aug 26 17:30:10 1999 (fw)
 * CVS info:   $Id: key.h,v 1.4 2008/03/14 23:11:57 mieg Exp $
 *-------------------------------------------------------------------
 */

#ifndef ACEDB_KEY_H
#define ACEDB_KEY_H


/* Modifier keys, i.e. Cntl, Alt, Shift, can be OR'd together. */
typedef enum 
{
  NULL_MODKEY,
  SHIFT_MODKEY = 1 ,
  CAPSLOCK_MODKEY = 2 ,
  CNTL_MODKEY  = 4 ,
  ALT_MODKEY   = 8 ,
  ALTGRAPH_MODKEY   = 16 ,
  NUMLOCK_MODKEY   = 32 ,
  META_MODKEY   = 64 ,
} GraphKeyModifier ;


/* Non-modifier keys. */

/* these were probably meant to match up to the usual ascii values for these keys
 * but this is a non-portable thing to do, sometime this should be fixed up
 * and done properly. */
#define      SPACE_KEY           32
#define     ESCAPE_KEY           27
#define     RETURN_KEY           13
#define  BACKSPACE_KEY            8
#define        TAB_KEY          '\t'

/* The numbers seem to be just the alphabetic position of the character.
** Setting up all alpha characters although not many used so far.      */
#define          A_KEY            1
#define          B_KEY            2  /* treedisp: find backwards */
#define          C_KEY            3
#define          D_KEY            4
#define          E_KEY            5
#define          F_KEY            6  /* treedisp: find forwards */
#define          G_KEY            7
#define          H_KEY            8
#define          I_KEY            9
#define          J_KEY           10
#define          K_KEY           11
#define          L_KEY           12
#define          M_KEY           13
#define          N_KEY           14
#define          O_KEY           15
#define          P_KEY           16
#define          Q_KEY           17
#define          R_KEY           18
#define          S_KEY           19
/* #define          T_KEY           20 alpha complains about this definition */
#define          U_KEY           21  /* treedisp: Update Mode */
#define          V_KEY           22
#define          W_KEY           23
#define          X_KEY           24
#define          Y_KEY           25
#define          Z_KEY           26
#define      COMMA_KEY           44  /* treedisp: find prev == CTL/< */
#define        DOT_KEY           46  /* treedisp: find next == CTL/> */

/* remainder are arbitrary here since their values vary depending
   on the graphics system used -- they must be remapped in 
   the event procedure that calls the function registered by 
   graphRegister (KEYBOARD, f) ;
*/

#define         UP_KEY          257
#define       DOWN_KEY          258
#define       LEFT_KEY          259
#define      RIGHT_KEY          260
#define       HOME_KEY          261
#define        END_KEY          262
#define    PAGE_UP_KEY          263
#define  PAGE_DOWN_KEY          264
#define     INSERT_KEY          265
#define     DELETE_KEY          127

#endif /* !ACEDB_KEY_H */
