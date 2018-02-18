
/*  File: win32thread.h
 *  Author: Richard Bruskiewich (rbrusk@octogene.medgen.ubc.ca)
 *  Copyright (C) R. Bruskiewich, J Thierry-Mieg and R Durbin, 1997
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
 *	WIN32 multithreading functionality for winace et al.
 *
 * Exported functions:
 * HISTORY:
 * Last edited:Mar 30 14:30 1997 (rbrusk)
 * Created: Mar 30 14:30 1997 (rbrusk)
 *-------------------------------------------------------------------
 */
/* $Id: win32thread.h,v 1.1.1.1 2002/07/19 20:23:27 sienkiew Exp $ */
/* WIN32 multithreading support */


typedef void *hACETHREAD ;

#if defined(WIN32THREAD_CPP)
extern "C" void parseLine (int isSkip, KEYSET ks) ;

#else

extern void parseGraphBoxDraw (int k, int fcol, int bcol) ;
extern hACETHREAD initParseThread(FILE *fil) ;
extern void win32ParseLine(int isSkip, KEYSET ks) ;
extern void getParseThreadLine(hACETHREAD ptHandle, int isSkip, KEYSET ks) ;
extern void closeParseThread(hACETHREAD handle) ;

#endif

 
