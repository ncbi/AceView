/*  File: pmap_.h
 *  Author: Neil Laister (nl1@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1992
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *      Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: objects and definitions internal to pmap module but shared
 * between its files
 * Exported functions: NONE
 * HISTORY:
 * Last edited: Dec 17 17:10 1998 (fw)
 *     - Added the Flag FLAG_DISPLAY_LOCI for toggling the display of loci.
 * * Feb  1 15:15 1997 (rd)
 * Created: Fri Jan 17 15:43:30 1992 (nl1)
 *-------------------------------------------------------------------
 */

/* $Id: pmap_.h,v 1.2 2014/11/30 03:21:16 mieg Exp $ */

/*..... ..... ..... .....
*/
#define MINUSINF (- ACEDB_MAXINT - 1)
#define FMINUSINF ((float)(-1000000000.0))

/*..... ..... ..... .....macros for handling flags:
*/

#define ClearReg(__obj) ((__obj)->flag=(unsigned int) 0)
#define Set(__obj, __m) ((__obj)->flag|=(__m))
#define Clear(__obj, __m) ((__obj)->flag&= ~(__m))
#define Toggle(__obj, __m) ((__obj)->flag^=(__m))
#define IsSet(__obj, __m) ((__obj)->flag & (__m))
#define IsntSet(__obj, __m) (!IsSet((__obj), (__m)))

/*for positioning operations..*/

/*..... ..... ..... .....displayed entities are represented as SEGs
*/
typedef struct
  { KEY key;
    KEY parent; /* clone associated with this item */
    float x0, x1;
      /*
      better store band coords of the ends for clones, otherwise I get
      a cumulative rounding error which quickly becomes significant
      under a lot of graphical manipulations; for YACs, the ends are
      computed arbitrarily on the midpoints of cosmids (so are subject
      to rounding errors); for some entities, I only store a midpoint
      for location on the map: since this is just for display purposes
      I keep it as a float (in x0): in these cases x1 is FMINUSINF 
      */
    unsigned int flag;
  } SEG ;

#define segFormat "kkffi"

#define SegMidPt(__seg) \
  ((float) (__seg)->x0==FMINUSINF? FMINUSINF: \
   (__seg)->x1==FMINUSINF? (__seg)->x0: 0.5*((__seg)->x0+(__seg)->x1))
  /*
  first case - the seg is not located on the map;
  second case - the seg is located by midpoint only;
  general case - centre of seg is computed;
  */
#define SegLen(__seg) ((__seg)->x1-(__seg)->x0)
#define SegRightEnd(__seg) ((__seg)->x1==FMINUSINF? (__seg)->x0: (__seg)->x1)

typedef struct PhysMapStruct
{
  magic_t *magic;        /* == &PhysMap_MAGIC*/
  KEY key; /*<<--the pmap, derived from contig or clone being displayed*/
  int activebox, min, max;
  int winHeight, winWidth, scrollBox;
  float centre;
  float cursor; /*position of the vertical Xor bar*/
  unsigned int flag; /*for switches dealing with the display entity, eg when "Show Buried Clones" gets selected*/
  Array /*of SEG*/ segs; /*content of all possible boxes*/
  Array /*of (SEG *)*/ box2seg; /*displayed boxes*/
  FREEOPT *menu;
  KEY currentSelected;
  Graph mainGraph; 
} *PhysMap;

/*..... ..... ..... .....
*/

typedef enum
{
  NORMAL, SELECTED /*currently-picked box*/, SISTER, HIGHLIT, BURIED_CLONE, DEBUG
} DrawType;

/*..... ..... ..... .....masks for LOOK flags:
*/

#define FLAG_SHOW_BURIED_CLONES     0x00000001
#define FLAG_BURIED_CLONES_ATTACHED 0x00000002
#define LOOK_WORK_REMARKS           0x00000004
#define FLAG_REVERT_DISPLAY         0x00000008
#define LOOK_NO_REMARKS             0x00000010
#define LOOK_NO_LOCUS               0x00000020
#define LOOK_NO_SCROLL_BAR          0x00000040
/*..... ..... ..... .....masks for SEG flags:
*/

#define FLAG_FINGERPRINT         0x00000001
#define FLAG_CANONICAL           0x00000002
#define FLAG_MORE_TEXT           0x00000004
#define FLAG_HIDE                0x00000008
#define FLAG_COSMID_GRID         0x00000010
#define FLAG_YAC_GRID            0x00000020
#define FLAG_WORK_REMARK         0x00000040
#define FLAG_PROBE               0x00000080
#define FLAG_HIGHLIGHT           0x00000100
#define FLAG_IS_BURIED_CLONE     0x00000200
#define FLAG_ERROR               0x00000400
#define FLAG_DISPLAY_BURIED_CLONES 0x00000800
  /*for marking (on the canonical) for display the buried clones of individual clones*/
#define FLAG_DISPLAY             0x00001000
  /*
  when adding clones, I may have to guess the band coords if I haven't enough data - I assume
  a "length" of 10 bands, centred on the current mark position, if there is one, else on
  the midpoint of the contig
  */
#define FLAG_SEQUENCED            0x00004000
#define FLAG_BURIED_IS_POSITIONED 0x00004000
  /*
  because order along contig may be affected (buried clones are positioned above their canonicals)
  */
#define FLAG_IS_COSMID 0x00010000
#define FLAG_IS_YAC    0x00020000
#define FLAG_IS_CDNA   0x00040000
#define FLAG_IS_FOSMID 0x00080000

#define FLAG_REMARK    0x00100000
/*..... ..... ..... .....masks for pickEntity_t flags:
*/
#define FLAG_SELECTED 0x00000001
#define FLAG_RELATED_SELECTION 0x00000002
  /*an apparently related match*/


/* pmapconvert.c */
int pMapCompareSeg (const void *a, const void *b);
Array pMapConvert(PhysMap look, KEY contig, KEY pmap);
void pMapRecalculateContig(PhysMap look, BOOL fromScratch) ;

/* pmapdisp.c */
BOOL pMapLiftBuriedClones(Array segs, SEG *canonical, BOOL setDisplayFlag);


#define CanonicalFromBuriedClone(__Objp, __canonical) \
  (bsGetKey((__Objp), _Exact_Match_to, &(__canonical)) || \
   bsGetKey((__Objp), _Approximate_Match_to, &(__canonical)) || \
   bsGetKey((__Objp), _Funny_Match_to, &(__canonical)))

#define CloneIsBuriedClone(__Objp) \
  (bsFindTag((__Objp), _Exact_Match_to) || \
   bsFindTag((__Objp), _Approximate_Match_to) || \
   bsFindTag((__Objp), _Funny_Match_to))

#define CloneIsGridded(__Objp, __key) bsGetKey((__Objp), _Hybridizes_to, &(__key))

#define CloneIsCosmid(__seg) (IsSet((__seg), FLAG_FINGERPRINT) || IsSet((__seg), FLAG_PROBE))


#define IsYAC(__seg) (class((__seg)->key)==_VClone && CloneIsYAC(__seg))

  /*are all these predicates right?*/

#define Max(__u, __v) ((__u)<(__v)? (__v): (__u))
#define Min(__u, __v) ((__v)<(__u)? (__v): (__u))
#define Trunc(__l, __v, __h) (((__v)<(__l))? (__l): ((__h)<(__v))? (__h): (__v))

/************************ end of file **************************/
 
