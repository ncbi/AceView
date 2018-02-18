/*  Last edited: Nov 18 18:38 1991 (rd) */

/* $Id: ace.h,v 1.1.1.1 2002/07/19 20:23:18 sienkiew Exp $ */


#ifndef DEFINE_ace_h
#define DEFINE_ace_h 1

char *nextField (char *pcutter) ;
void outItem (FILE *fil, char *tag, char *text) ;
void outList (FILE *fil, char *tag, char **items, int n) ;
void outStack (FILE *fil, char *tag, Stack s);

extern FILE *output ;

#endif
/******************************************************************/


