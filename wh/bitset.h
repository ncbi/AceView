/*  File: bitset.h
 *  Author: Richard Durbin (rd@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1997
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
 * 	Richard Durbin (MRC LMB, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * SCCS: $Id: bitset.h,v 1.14 2016/10/15 17:03:45 mieg Exp $
 * Description: This file should be included in codes taht explicitely access bitsets
 *   the actual descriotion of the bitset interface is now in wh/bitset_.h
 * Created: Sat Feb  1 11:07:23 1997 (rd)2015_09_25 (mieg)
 *-------------------------------------------------------------------
 */

#ifndef BITSET_H
#define BITSET_H

#include "bitset_.h"
static unsigned int bitField[32] = {0x1, 0x2, 0x4, 0x8,                              
				    0x10, 0x20, 0x40, 0x80,                          
				    0x100, 0x200, 0x400, 0x800,                      
				    0x1000, 0x2000, 0x4000, 0x8000,                  
				    0x10000, 0x20000, 0x40000, 0x80000,              
				    0x100000, 0x200000, 0x400000, 0x800000,          
				    0x1000000, 0x2000000, 0x4000000, 0x8000000,      
				    0x10000000, 0x20000000, 0x40000000, 0x80000000 } ;


#ifdef JUNK
static void showBit (BitSet bb, int ii)
{
  int i ;
  char buf[11] ;

  if (bb)
    {
      for (i=0; i<10;i++)
	{
	  if (bitt(bb,ii+i)) 
	    buf[i]='+' ; 
	  else 
	    buf[i]='-' ;
	}
      buf[i]=0 ;
      fprintf(stderr, "BitSet %d:: %s\n", ii, buf) ;
      if (FALSE)    /* never ! */
	showBit (0, 0) ; /* to please the compiler */
    }
}
#endif


#endif

/***************** end of file *****************/
 
 
