

#ifndef DEFINE_BU_H
#define DEFINE_BU_H

    /* A buTree is a bsTree plus surrounding ferns managing paths and Xref */

/* $Id: bu.h,v 1.1.1.1 2002/07/19 20:23:18 sienkiew Exp $ */

      /* Functions to transfer a buTree  in/out the block cache */      
OBJ_HANDLE    	buGet     (KEY key) ;
void  	buStore     (KEY key, OBJ_HANDLE  x) ;
                              /* does not imply destructiony */
void buStoreNoCopy  (KEY key, OBJ_HANDLE  x) ;
                              /* implies bsTreePrune(x->root) */

     /* Functions for the manipulations of whole buTrees */
void 		buDestroy(OBJ_HANDLE  x) ;
OBJ_HANDLE   	buCopy   (OBJ_HANDLE  x) ;


#endif
/******************************************************************/








