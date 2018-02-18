#ifndef ANNOT_H
#define ANNOT_H

typedef struct noteStruct 
 { 
   int x; char *title ; char *tagName; KEY tag ; KEY type ; int length ; BSunit u ;
 } ANNOTATION ;

void annotate (KEY key, ANNOTATION *an) ;

#endif
