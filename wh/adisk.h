
/* $Id: adisk.h,v 1.3 2015/08/17 16:12:51 mieg Exp $ */

/* prototypes for the adisk routines
 *  new way to store data in test for acedb5
 */


#ifndef DEFINE_DISK_H
#define DISK_H

#ifndef  DEF_DISK       
typedef KEY  DISK; 
#define DEF_DISK           
#endif
int aCacheSaveAll(void) ;


void aDiskArrayStore (KEY key, KEY parent, Array a, Array b) ;
BOOL aDiskArrayGet (KEY key, KEY *pp, Array *ap, Array *bp) ;

void aDiskKill (KEY key) ;

void aDiskReadHeader (DISK d, KEY key, 
		     KEY *parent, 
		     int *n1, int *s1, 
		     int *n2, int *s2) ;
void aDiskReadData (DISK d,
		   char *p1, int n1, int s1,
		   char *p2, int n2, int s2) ;
DISK aDiskWrite (KEY key, KEY parent, 
		char* p1, int n1, int s1,
		char* p2, int n2, int s2) ;

int  aDiskGetGlobalAdress(void)  ;

DISK aDiskAssign (KEY key, int n1, int s1, int n2, int s2) ;
void aDiskAccessInit(BP bp) ;

#endif
