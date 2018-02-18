/*  File: graphcolour.c
 *  Author: Simon Kelley (srk@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1991
 * -------------------------------------------------------------------
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
 * Description: maps graph package colours to RGB / grey triplets/
 *-------------------------------------------------------------------
 */

#include <colours.h>
#include <w2/graphcolour.h>


/* NB grey values are arranged so that (eg Red |= Green != Blue to make them
   distinguishable on greyscale display */
  /* These colors match in the same order those declared 
   * whooks/systags.h
   * whooks/sysclass.c twice
   * wh/colours.h
   * wh/colorRGB.h
   * w2/graphcolour.c  
   * w4/palette.c
   * waceview/swfc.c
   */

static int colorTable[4*NUM_TRUECOLORS]= {
255,	255,	255,   255, 	   /* WHITE           */
0,	0,	0,     0,	   /* BLACK           */
200,	200,	200,   178,	   /* LIGHTGRAY       */
100,	100,	100,   76, 	   /* DARKGRAY        */
255,	0,	0,     102,        /* RED             */
0,	255,	0,     128,        /* GREEN           */
0,	0,	255,   142,	   /* BLUE            */
255,	255,	0,     229,	   /* YELLOW          */
0,	255,	255,   178,	   /* CYAN            */
255,	0,	255,   153, 	   /* MAGENTA         */

255,	113,	113,   191,	   /* LIGHTRED        */ /* mieg 22 sept 2004 */
0,	201,	196,   204,	   /* LIGHTGREEN      */ /* mieg 22 sept 2004 */
133,	159,	255,   216,        /* LIGHTBLUE       */ /* mieg 22 sept 2004 */
204,	0 ,	51,   38,	   /* DARKRED         */ /* mieg 22 sept 2004 */
0,	153,	153,   51,	   /* DARKGREEN       */ /* mieg 22 sept 2004 */
22,	0,	200,   63,	   /* DARKBLUE        */ /* mieg 22 sept 2004 */
255,	209,	209,   229,	   /* PALERED         */ /* mieg 22 sept 2004 */
201,	255,	209,   229,	   /* PALEGREEN       */ /* mieg 22 sept 2004 */
205,	236,	255,   242,	   /* PALEBLUE        */ /* mieg 22 sept 2004 */
255,	255,	195,   242,	   /* PALEYELLOW      */ /* mieg 22 sept 2004 */
189,	255,	255,   242,	   /* PALECYAN        */ /* mieg 22 sept 2004 */
255,	200,	255,   242,	   /* PALEMAGENTA     */ /* mieg 22 sept 2004 */


168,     40,    0,     102,        /* BROWN           */ /* mieg 22 sept 2004 */
255,	61,	0,     153,	   /* ORANGE          */ /* mieg 22 sept 2004 */
255,	214,	201,   204,	   /* PALEORANGE      */ /* mieg 22 sept 2004 */
164,	0,	164,   127,	   /* PURPLE          */ /* mieg 22 sept 2004 */
161,	94,	228,   178,	   /* VIOLET          */ /* mieg 22 sept 2004 */
228,	209,	247,   229,	   /* PALEVIOLET      */ /* mieg 22 sept 2004 */
150,	150,	150,   140,	   /* GRAY            */
235,	235,	235,   216,	   /* PALEGRAY        */
255,	0,	92,    127,	   /* CERISE          */ /* mieg 22 sept 2004 */
86,	178,	222,   127,	   /* MIDBLUE         */ /* mieg 22 sept 2004 */
255,	125 ,	255,   200,	   /* LIGHTMAGENTA    */ /* mieg 22 sept 2004 */
127,	255 ,	255,   200,	   /* LIGHTCYAN       */ /* mieg 22 sept 2004 */
123,	35,	211,    38,	   /* DARKVIOLET      */ /* mieg 22 sept 2004 */
65,     106,    255,   140,        /* LAVANDER        */

  /* these colour scales have a dew independant grey value */
236,238,254,255, /* BLUE1 */
201,207,251,224, /* BLUE2 */
165,175,249,192, /* BLUE3 */
130,144,246,128, /* BLUE4 */
94,112,244,238, /* BLUE5 */
59,81,241,172, /* BLUE6 */
22,48,238,118, /* BLUE7 */
14,35,190,32,  /* BLUE8 */

232,254,247,255, /* GREEN1 */
197,251,233,224, /* GREEN2 */
161,249,220,192, /* GREEN3 */
106,246,199,160, /* GREEN4 */
35,241,172,128, /* GREEN5 */
13,205,141,96, /* GREEN6 */
10,152,105,64, /* GREEN7 */
7,115,79,32, /* GREEN8 */


255,217,217,255, /* RED1 */
255,179,179,224, /* RED2 */
255,139,139,192, /* RED3 */
255,101,101,160, /* RED4 */
255,63,63,128, /* RED5 */
255,5,5,96, /* RED6 */
222,0,0,64, /* RED7 */
184,0,0,32, /* RED8 */

199,	159,	239,   178,	   /* LIGHTVIOLET      */ /* mieg 2008 */
0,	215,    210,   178,	   /* DARKCYAN      */ /* mieg 2008 */
255,    170,    143,   178,	   /* LIGHTORANGE      */ /* mieg 2008 */
} ;


void graphGetColourAsFloat(int colour, float *r, float *g, float *b, float *gr)
{
  if (r)
    *r = ((float)colorTable[4*colour])/255.0;
  
  if (g)
    *g = ((float)colorTable[4*colour+1])/255.0;

  if (b) 
    *b = ((float)colorTable[4*colour+2])/255.0;

  if (gr)
    *gr = ((float)colorTable[4*colour+3])/255.0;


}

void graphGetColourAsInt(int colour, int *r, int *g, int *b, int *gr)
{
  if (r)
    *r = colorTable[4*colour];
  
  if (g)
    *g = colorTable[4*colour+1];

  if (b)
    *b = colorTable[4*colour+2];
  
  if (gr)
    *gr = colorTable[4*colour+3];
}

