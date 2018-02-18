/* $Id: pref.h,v 1.1.1.1 2002/07/19 20:23:20 sienkiew Exp $ */

#ifndef DEF_PREF_H
#define DEF_PREF_H

void prefInit() ;
BOOL prefValue(char * name) ;
float prefFloat(char * name) ;
int prefInt(char * name) ;
int prefColour (char * name) ;
void editPreferences(void) ;

extern Array prefList;  /* bizare to have that public */


#endif
 
 
