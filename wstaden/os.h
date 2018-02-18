/*
 * File: os.h
 *
 * Author: 
 *         MRC Laboratory of Molecular Biology
 *	   Hills Road
 *	   Cambridge CB2 2QH
 *	   United Kingdom
 *
 * Description: operating system specific type definitions
 *
 */

#ifndef _OS_H_
#define _OS_H_

#include <limits.h>


/* ACEDB modifications, jean thierry-mieg */
#if defined (ALPHA) || defined (LINUX) || defined (INTEL_SOLARIS) || defined (OPTERON) || defined (CYGWIN)
#if !defined (LITTLE_ENDIAN)
#define LITTLE_ENDIAN
#endif
#if defined (BIG_ENDIAN)
#undef BIG_ENDIAN
#endif
#else   /* SGI SUN SOLARIS IBM (I think) */
#if !defined (BIG_ENDIAN)
#define BIG_ENDIAN
#endif
#if defined (LITTLE_ENDIAN)
#undef LITTLE_ENDIAN
#endif
#endif

/*
 *-----------------------------------------------------------------------------
 * Typedefs for data sizes. Note there's umpteen versions of typedefs here
 * due to old code being supported. The ones that should be used everywhere
 * are {u,}int[124].
 *-----------------------------------------------------------------------------
 */

/*
 * One byte integers
 */ 
/*
typedef unsigned char	uint1;
typedef signed char	uint1;
*/
typedef unsigned char	int1;

/*
 * Two byte integers
 */
typedef signed short	int2;
typedef unsigned short	uint2;

/*
 * Four byte integers
 */
typedef signed int	int4;
typedef unsigned int	uint4;


/*
 * Backwards compatibility
 */
typedef signed char	int_1;
typedef unsigned char	uint_1;
typedef signed short	int_2;
typedef unsigned short	uint_2;
typedef signed int	int_4;
typedef unsigned int	uint_4;


/*
 *-----------------------------------------------------------------------------
 * The FORTRAN interface.
 *-----------------------------------------------------------------------------
 */

typedef int4 f_int;
typedef int4 f_implicit;
typedef void f_proc_ret;	/* procedure return value */

/* James Bonfield compatability mode */
typedef int4 int_f;		/* f_int */
typedef int4 int_fl;		/* f_implicit */

#define f_proc_return() return /* (f_proc_ret) 0 */

/*
 * Use when calling/defining a Fortran function from C.
 */
#ifdef VMS
#    define FORT(symbol) (symbol)
#else
#    define FORT(symbol) (_symbol)
#endif


/*
 *-----------------------------------------------------------------------------
 * Some handy definitions.
 *-----------------------------------------------------------------------------
 */

#define MAXINT4 (INT_MAX)
#define MAXINT2 (SHRT_MAX)


/*
 *-----------------------------------------------------------------------------
 * Machine specifics
 *-----------------------------------------------------------------------------
 */
/*
 * SunOS 4.x
 * Even though we use the ANSI gcc, we make use the the standard SunOS 4.x
 * libraries and include files, which are non-ansi
 */
#if defined(__sun__) && !defined(__svr4__)
#define SEEK_SET 0
#define SEEK_CUR 1
#define SEEK_END 2
#define FOPEN_MAX 64
#define FILENAME_MAX 1024
#endif

/* Microsoft Visual C++ */
#if defined(_MSC_VER)
#define popen _popen
#define pclose _pclose
typedef int mode_t;
#endif

/*
 * Linux defines both BIG_ENDIAN and LITTLE_ENDIAN in its include files.
 * We ought to change things to (eg) define ENDIAN to be big or little and
 * then check lower down on this.
 */
#if defined(__linux__) && defined(BIG_ENDIAN)
#undef BIG_ENDIAN
#endif

/*
 *=============================================================================
 * Anything below here should not be changed.
 *=============================================================================
 */

#define False 0
#define True 1


#ifdef OLD_SWAP
/* See below for reasons not to revert back. We ought to remove these. */

/* copy INT4 from src to dst byteswapping on the way */
#define swap_int4(src, dst) \
    do {\
	((char *)&(dst))[0] = ((char *) &(src))[3];\
	((char *)&(dst))[1] = ((char *) &(src))[2];\
        ((char *)&(dst))[2] = ((char *) &(src))[1];\
        ((char *)&(dst))[3] = ((char *) &(src))[0];\
    } while (0)

/* copy INT2 from src to dst byteswapping on the way */
#define swap_int2(src, dst) \
    do {\
        ((char *) &(dst))[0] = ((char *) &(src))[1];\
        ((char *) &(dst))[1] = ((char *) &(src))[0];\
    } while (0)

#else

/*
 * Our new swap runs at the same speed on Ultrix, but substantially faster
 * (300% for swap_int4, ~50% for swap_int2) on an Alpha (due to the lack of
 * decent 'char' support).
 *
 * They also have the ability to swap in situ (src == dst). Newer code now
 * relies on this so don't change back!
 */
#define swap_int4(src, dst) \
    dst = ((src & 0x000000ff) << 24) + \
          ((src & 0x0000ff00) <<  8) + \
          ((src & 0x00ff0000) >>  8) + \
          ((src & 0xff000000) >> 24)

#define swap_int2(src, dst) \
    dst = ((src & 0x00ff) << 8) + \
          ((src & 0xff00) >> 8)
#endif


/*
 * Slightly updated swap_int? routines that return results rather than
 * swapping from source to destination.
 */
#define iswap_int4(x) \
    (((x & 0x000000ff) << 24) + \
     ((x & 0x0000ff00) <<  8) + \
     ((x & 0x00ff0000) >>  8) + \
     ((x & 0xff000000) >> 24))

#define iswap_int2(x) \
    (((x & 0x00ff) << 8) + \
     ((x & 0xff00) >> 8))

/*
 * Macros to specify that data read in is of a particular endianness.
 * The macros here swap to the appropriate order for the particular machine
 * running the macro and return the new answer. These may also be used when
 * writing to a file to specify that we wish to write in (eg) big endian
 * format.
 *
 * This leads to efficient code as most of the time these macros are
 * trivial.
 */
#ifdef BIG_ENDIAN
#define be_int4(x) (x)
#define be_int2(x) (x)
#define be_int1(x) (x)

#define le_int4(x) iswap_int4((x))
#define le_int2(x) iswap_int2((x))
#define le_int1(x) (x)
#endif

#ifdef LITTLE_ENDIAN
#define be_int4(x) iswap_int4((x))
#define be_int2(x) iswap_int2((x))
#define be_int1(x) (x)

#define le_int4(x) (x)
#define le_int2(x) (x)
#define le_int1(x) (x)
#endif

#ifndef BIG_ENDIAN
#ifndef LITTLE_ENDIAN
#error Must define BIG_ENDIAN or LITTLE_ENDIAN in Makefile
#endif
#endif

#ifdef BIG_ENDIAN
#ifdef LITTLE_ENDIAN
#error Must only define one of BIG_ENDIAN and LITTLE_ENDIAN in Makefile
#endif
#endif

#endif /*_OS_H_*/
