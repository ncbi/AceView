/*  Last edited: Dec  4 14:50 1998 (fw) */

/* @(#)freeout.h	1.2 12/13/95 */

#ifndef FREEOUT_H_DEF
#define FREEOUT_H_DEF

#include "regular.h"

int freeOutSetFile (FILE *fil) ;
int freeOutSetStack (Stack s) ;

void freeOutInit (void) ;
void freeOut (const char *text) ;
void freeOutf (const char *format,...) ;
void freeOutxy (const char *text, int x, int y) ; 
void freeOutBinary (const char *data, int size) ;

void freeOutClose (int level) ;

int freeOutLine (void) ;
int freeOutByte (void) ;
int freeOutPos (void) ;

void freeOutShutDown (void) ;

#endif
