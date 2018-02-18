/*  File: aceBinary.c
 * $Header: /home/mieg/aa/CVS_ACEVIEW/ace/w6/acebinary.c,v 1.1 2011/05/10 16:40:58 mieg Exp $
 *  Author: Richard Durbin (rd@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1994
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
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
 * Exported functions:
 * HISTORY:
 * Created: Tue May 10 23:59:13 1994 (rd)
 *-------------------------------------------------------------------
Thanks Jean,

I got your file.  
I also update several other tables in the aceBinary.c file:
the molecular weight of Sec(U) is 168.05, so the molecularWeight[] is 
changed
to 168 for U, I also update Z, B, and X molecular weight.
PKa of NH and COOH of U is update according to that of C
I did not look for the exact number but U and C shuld be close.

The molecular weight of aceBinary can be simply calculated by
sum(Mw of amino acids) - (length-1)*Mw of Water
if you want to be more accurate, then molecularWeight[] should be
declared as float or double.
The aceBinary.c file is calculating the number of atoms (kind of
 unnecessary and slower).
Kemin Zhu

 */

/* $Id: acebinary.c,v 1.1 2011/05/10 16:40:58 mieg Exp $ */

#include <ctype.h>
#include <wh/acedb.h>
#include <wh/a.h>
#include <wh/aceBinary.h>
#include <wh/lex.h>
#include <whooks/sysclass.h>

#include "freeout.h"
#include <wh/pick.h>   /* pickWord2Class */
#include <wh/dump.h>
#define  SIZEOFINT (8 * sizeof (int))

