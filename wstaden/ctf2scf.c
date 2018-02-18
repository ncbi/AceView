/*  File: ctf2scf.c
 *  Author: mieg@crbm.cnrs-mop.fr
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
 * Description: 
 * Exported functions: <standalone>
 * HISTORY:
 * Last edited: Aug 26 17:09 1999 (fw)
 * Created: Thu Aug 26 17:08:53 1999 (fw)
 * CVS info:   $Id: ctf2scf.c,v 1.1.1.1 2002/07/19 20:23:37 sienkiew Exp $
 *-------------------------------------------------------------------
 */

#include "regular.h"
#include "Read.h"

int main(int argc, char **argv)
{ 
  Read *rr = 0 ;

  freeinit() ;

  if (argc>1)
    goto usage ;

  rr = fread_reading(stdin, "stdin", 0) ; 
  fwrite_reading (stdout, rr, TT_SCF) ; 

  return 0 ;

usage:
  printf("Usage:  ctf2scf < ff.ctf > ff.scf\n") ;
  return 1 ;
}

/***************************************************************/
/***************************************************************/
