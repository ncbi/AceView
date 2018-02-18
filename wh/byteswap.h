/*  File: byteswap.h
 *  Author: Simon Kelley (srk@sanger.ac.uk)
 * Copyright (C) J Thierry-Mieg and R Durbin, 1994
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: internal to external data representation conversion
 * HISTORY:
    mieg: suppressed acedb4 protection, 
          added reincude protection 
          made swapData a BOOL
 * Last edited: Nov 26 10:47 1998 (fw)
 * Created: Wed 1 Mar 1995 (srk)
 *-------------------------------------------------------------------
 */

/*  $Id: byteswap.h,v 1.3 2003/06/03 22:14:29 mieg Exp $  */

/* NB we assume that any type of BSdata can be swapped as if it were an int.
   this means that they must all be exactly 32 bits. True for the new lex disk.
 */

#ifndef BYTESWP_H_DEF
#define BYTESWP_H_DEF

#include "session.h"		/* for extern BOOL swapData */


#define swap_word(a) ( ((a) << 24) | \
                      (((a) << 8) & 0x00ff0000) | \
                      (((a) >> 8) & 0x0000ff00) | \
        ((unsigned int)(a) >>24) )

#define swap_half(a) ( ((a & 0xff) << 8) | ((unsigned short)(a) >> 8) )


#define swapInt(x) swap_word(x)
#define maybeSwapInt(x) (swapData ? swapInt(x) : (x))
#define swapKEY(x) swap_word(x)
#define maybeSwapKEY(x) (swapData ? swapKEY(x) : (x))
#define swapDISK(x) swap_word(x) /* only valid for ACEDB4 */
#define maybeSwapDISK(x) (swapData ? swapDISK(x) : (x))
#define swapNODE(x) swap_half(x)
#define maybeSwapNODE(x) (swapData ? swapNODE(x) : (x))


#endif
