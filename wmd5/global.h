/*  Last edited: Jan  6 13:06 2000 (edgrif) */
/* GLOBAL.H - RSAREF types and constants
 */

/* PROTOTYPES should be set to one if and only if the compiler supports
  function argument prototyping.
The following makes PROTOTYPES default to 0 if it has not already
  been defined with C compiler flags.
 */
#ifdef ED_G_NEVER_INCLUDE_THIS_CODE
#ifndef PROTOTYPES
#define PROTOTYPES 0
#endif
#endif /* ED_G_NEVER_INCLUDE_THIS_CODE */

/* Default to prototypes.                                                    */
#ifndef PROTOTYPES
#define PROTOTYPES 1
#endif


/* POINTER defines a generic pointer type */
typedef unsigned char *POINTER;

/* UINT2 defines a two byte word */
typedef unsigned short int UINT2;

/* UINT4 defines a four byte word */
#ifdef ED_G_NEVER_INCLUDE_THIS_CODE
typedef unsigned long int UINT4;
#endif /* ED_G_NEVER_INCLUDE_THIS_CODE */
/* "long" is incorrect for 64 bit architectures, should just use int.        */
typedef unsigned int UINT4;


/* PROTO_LIST is defined depending on how PROTOTYPES is defined above.
If using PROTOTYPES, then PROTO_LIST returns the list, otherwise it
  returns an empty list.
 */
#if PROTOTYPES
#define PROTO_LIST(list) list
#else
#define PROTO_LIST(list) ()
#endif


