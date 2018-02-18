/*  File: model.h
 *  Author: Fred Wobus (fw@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1998
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@crbm.cnrs-mop.fr
 *
 * Description: public header for model.c
 * Exported functions:
 * HISTORY:
 * Last edited: Dec  2 17:29 1998 (fw)
 * Created: Mon Nov 23 15:07:46 1998 (fw)
 *-------------------------------------------------------------------
 */

#ifndef _MODEL_H
#define _MODEL_H

#include "acedb.h"

BOOL readModels (void);
BOOL isModel (KEY key);

BOOL getNewModels (void);    /* used by update.c:updateDoAction() */

VoidRoutine modelChangeRegister (VoidRoutine func);

#endif /* _MODEL_H */
