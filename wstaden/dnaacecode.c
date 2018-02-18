#include "regular.h"
#include "keyset.h"
#include "dna.h"

/* this is the mapping use to print out */
char dnaDecodeChar[] =
 { '-','a','t','w','g','r','k','d','c','m','y','h','s','v','b','n',
     0,  0,  0,  0, '.'} ;
char rnaDecodeChar[] =
 { '-','a','u','w','g','r','k','d','c','m','y','h','s','v','b','n',
     0,  0,  0,  0, '.'} ;
				/* RD added '.' to help alignment code */
char complementBase[] =	
 { 0, T_,A_,W_,C_,Y_,M_,H_,G_,K_,R_,D_,S_,B_,V_,N_ } ;
/* this is the mapping use to parse in */
char dnaEncodeChar[] =
{  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,

   0,  A_,  B_,  C_,  D_,   0,   0,  G_,  H_,   0,   0,  K_,   0,  M_,  N_,   0,
   0,   0,  R_,  S_,  T_,  U_,  V_,  W_,   0,  Y_,   0,   0,   0,   0,   0,   0,
   0,  A_,  B_,  C_,  D_,   0,   0,  G_,  H_,   0,   0,  K_,   0,  M_,  N_,   0,
   0,   0,  R_,  S_,  T_,  U_,  V_,  W_,   0,  Y_,   0,   0,   0,   0,   0,   0,
} ;


