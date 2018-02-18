/*  File: makemrna.h
 *  Author: mieg
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
 * Exported functions:
 * HISTORY:
 * Created: Mars 2001 (mieg)
 * CVS info:   $Id: makemrna.h,v 1.4 2017/02/15 20:39:28 mieg Exp $
 *-------------------------------------------------------------------
 */
/* %W% %G% */

#ifndef MAKEMRNA_H_DEF
#define MAKEMRNA_H_DEF
#include "cdnatr.h"

KEY  makeMrnaGene (S2M *s2m, SC* sc, SMRNA *gmrna, Array smrnas, KEYSET clipTops, KEYSET clipEnds, Array linkPos) ;

void mrnaAnalyseClusterAllPg (char *cp) ;
void mrnaAnalyseAllNeighbours (void) ;
void mrnaAnalyseNeighbours (KEYSET genes) ;
int  mrnaQualityEvaluate (int ln, int ali, int  nerr, BOOL isMrna, int *ixp, int *iyp) ;
KEY mrnaIntronBoundary (Array dnaD, Array dnaR,  int a1, int a2) ;
BOOL mrnaDesignUsingCompositeStrategy (S2M *s2m, SC* sc, SMRNA *gmrna, Array smrnas) ;
#endif
