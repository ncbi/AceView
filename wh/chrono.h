/*  File: chrono.h
 *  Author: mieg (mieg@sanger.ac.uk)
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
 * Last edited: Aug 26 17:05 1999 (fw)
 * Created: Thu Aug 26 17:05:18 1999 (fw)
 * CVS info:   $Id: chrono.h,v 1.2 2003/06/30 03:33:43 mieg Exp $
 *-------------------------------------------------------------------
 */

/* $Id: chrono.h,v 1.2 2003/06/30 03:33:43 mieg Exp $ */


#ifndef ACEDB_CHRONO_H
#define ACEDB_CHRONO_H

void chronoDoReturn(void);
void chronoSwitch(char *procedureName) ;
BOOL chronoStart(void) ;
void chronoStop(void) ;
void chronoReport (void) ;

#ifdef CHRONO
#define chrono(cp) chronoSwitch(cp)
#define chronoReturn() chronoDoReturn() 
#else	  /* suppress any call to chrono */
#define chrono(cp)
#define chronoReturn() 
#endif

#endif /* !ACEDB_CHRONO_H */
 
 


