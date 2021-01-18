/*  Last edited: Oct 13 13:00 1993 (ulrich) */

  /* topology.h */

/* $Id: topology.h,v 1.3 2020/03/12 00:18:40 mieg Exp $ */

typedef struct { KEY a , b ; int group , type ; } LINK ;
#define linkFormat "kkii"

typedef struct { KEY a ; int group, pos ; } VERTEX ;
#define vertexFormat "kii"

int topoVertexOrder (const void *v1, const  void *v2) ;
int topoLinkOrder (const void *v1, const  void *v2) ;
int topoConnectedComponents(Array links, Array vx) ;

/* sort a 2D matrix recursivelly by the total weights of its lines and columns
 * aa must be an array of float of size iMax*jMax
 * if iMax = keySetMax(lines)
 *    jMax = keySetMax(cols)
 * aa must be an array of float with arrayMax(aa) = iMax*jMax
 *
 * On exit, aa has the same shape but is fully reordered
 * the order is returned as 'lines' and 'cols'
 *  so that the value is aa[i*jMax + j]
 *  but     newLineTitle[i] = originalTitle[lines[i]]
 *          newColTitle[j]  = originalTitle[cols[j]]
 */
BOOL topoChronoOrder (Array aa, KEYSET lines, KEYSET cols) ;

