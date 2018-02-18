
#ifndef VDEF_H_DEFINED /* vdef.h */

/*
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
_/
_/ In this include file the 
_/ basic types, constants, 
_/ variables, and methods 
_/ are defined
_/
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/

/** BASICS **/

    /* primary original constants declaration */
    #define vNULL    0             /* the absolute null integer number */
    #define vTRUE    1             /* true for logical variables */
    #define vFALSE   0             /* false for logical variables */

    #define vMAXINUM (0x7FFFFFFFL) 			/* maximum integer number module */
    #define vMININUM (0x80000000L)     		/* minimum integer number module */
    #define vMAXRNUM (1.7976931348623158e+308)   /* maximum real number module */
    #define vMINRNUM (2.2250738585072014e-308)   /* minimum real number module */

    #define vINVADR  vNULL         /* memory allocation fault result */
    #define vINVIVAL vMAXINUM      /* reserve as invalid integer number */
    #define vINVRVAL vMAXRNUM      /* reserve as invalid real number */

    /* numerical operations declaration */
    #define vDim(vrn_vn)           ( sizeof( (vrn_vn) ) / sizeof( (vrn_vn)[0] ) )    	/* dimesion of matrixes of any types */
    #define vMax(vrn_vn1,vrn_vn2)  ( (vrn_vn1) > (vrn_vn2) ? (vrn_vn1) : (vrn_vn2) )   /* maximum value from two variables */
    #define vMin(vrn_vn1,vrn_vn2)  ( (vrn_vn1) < (vrn_vn2) ? (vrn_vn1) : (vrn_vn2) )   /* minimum value from two variables */
    #define vAbs(vrn_vn)           ( (vrn_vn) > 0 ? (vrn_vn) : -(vrn_vn) )             /* absolute value of variable */
    #define vSi0(vrn_vn)           ( (vrn_vn) > 0 ? 1 : ( (vrn_vn) < 0 ? -1 : 0 ) )  	/* the sign of variable with 0 */
    #define vSig(vrn_vn)           ( (vrn_vn) >= 0 ? 1 : -1 )  						/* the sign of variable */
    #define vSif(vrn_vn1,vrn_vn2)  ( vAbs((vrn_vn1))*vSig((vrn_vn2)))				/* the firs variable with signt of the second one */
    #define vToStr(vrn_vn)         ( #vrn_vn)

	#define	vSet0(vrn_vn)			(memset(&vrn_vn,0,sizeof((vrn_vn))))
	#define	vToVoid(vrn_vn)			(void *)((char *)(0)+(vrn_vn))
	#define	vToInt(vrn_vn)			(int)((char *)(vrn_vn)-(char *)(0))
	#define vVarDef(vrn_vn,vrn_fmt) vToVoid(&(vrn_vn)),vToVoid(sizeof((vrn_vn))),vToVoid(#vrn_vn),vToVoid((vrn_fmt))

	#define	vCallVarg(v_res,v_fun,v_frmt) {va_list ap;va_start(ap,v_frmt);(v_res)=(v_fun)((v_frmt),ap);va_end(ap);}
	#define	vCallParaVarg(v_res,v_fun,v_para,v_frmt) {va_list ap;va_start(ap,v_frmt);(v_res)=(v_fun)((v_para),(v_frmt),ap);va_end(ap);}
	typedef void * (* vCallbackUniversal)(void * param,...);
	typedef void * (* vCallbackLogical)(void * param,...);

    /* type definitions */
    #ifdef  __cplusplus
        #define vIN  inline
    #else
        #define vIN
    #endif

    #define vclrRGB(r,g,b)           ( (((unsigned int)(b))&0xFF) | ((((unsigned int)(g))&0xFF)<<8) | ((((unsigned int)(r))&0xFF)<<16) )
    #define vclrB(rgb)               ( (((unsigned int)(rgb)) )&0xFF )
    #define vclrG(rgb)               ( (((unsigned int)(rgb)) >> 8)&0xFF )
    #define vclrR(rgb)               ( (((unsigned int)(rgb)) >> 16)&0xFF )

#define VDEF_H_DEFINED
#endif

