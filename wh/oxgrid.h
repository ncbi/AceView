/* $Id: oxgrid.h,v 1.2 2015/08/12 14:12:11 mieg Exp $ */
/*  File: oxgrid.h
 *  Written by Jo Dicks (jld@bioch.ox.ac.uk)
 *-------------------------------------------------------------------
 * This file is part of the ACeDB comparative mapping package, 
 * written by Jo Dicks and Michelle Kirby
 * The original ACeDB system was written by :
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmba.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *      Copyright (C) J Thierry-Mieg and R Durbin, 1991
 *
 * Description: header file for all comparative mapping functions
 * Exported functions: oxgridCreate (), oxgridDisplay (), pairMapCreate (),  
 *                     o2mCreate (), oxhomlist ()  
 * HISTORY:
 * Created: July 1993
 * Last edited: Dec  4 14:48 1998 (fw)
 *-------------------------------------------------------------------
 */

#include "regular.h"

typedef struct OXSTUFF
  { 
    int   flag, specflag;                      /* Oxford Grid */
    void  *magic ;
    int   searchBox, sp1Box, sp2Box, mapBox ;
    int   gridBox, h, high, wide, map_labels, mlb ;
    const char  *searchp, *mapset1p, *mapset2p, *specmapp ;
    char  findhom[20], findloc[20] ;
    BOOL  modified, hideFish, showRandom ;
    KEY   s1, s2, search, mapset1, mapset2, specmap ;
    Array segs, boxes, axis1, axis2, size1, size2, extent1, extent2; 
    Graph graph, chooseGraph ;

    int   listBox, list, stayBox ;                     /* Homology List */
    Graph hlGraph ;

    int   pairBox, pairhigh, pairwide, pair_labels, messageBox ; 
    char  pairhom[20], pairloc[20], messageText[80] ;      
    float c1min, c1max, c2min, c2max ;                 /* Pairmap */
    KEY   s3, s4, c1, c2 ;
    Array points, c1segs, c2segs ;
    Graph pairGraph ;

    int   manyBox, o2mhigh, o2mwide, o2m_labels, message2Box ;
    int   segsize ;
    char  o2mhom[20], o2mloc[20], message2Text[80] ;
    float c3min, c3max ;                               /* One-to-many Map */
    KEY   s5, c3 ;
    Array many, csegs, c3segs, axis3, size3, ohomols ;
    Graph manyGraph ;

    int   specBox, s, spechigh, specwide, spec_labels, message3Box ;
    char  spechom[20], specloc[20], message3Text[80] ; 
    float c4min, c4max ;                               /* Species Grid */
    KEY   c4 ;
    Array specsegs, specs, c4segs, homols ;                   
    Graph specGraph ;

  } *OX ;

typedef struct                /* standard seg for homologies */
  { KEY   key, spec2 ;
    KEY   gene1, gene2 ;
    KEY   map1, map2 ;
    float x1, dx1 ;
    float z1 , z2, z3, z4 ;  /* position on map1 map2 */
    int   box, col, hbox ;
    int   x, y ;
    int   flag ;
    BOOL isFish ;
  } SEG ;

typedef struct                /* contains data on chromosomes */
  { KEY key ;         
    float x, dx ;
    char* sym ;        
    unsigned int flag ; 
  } CHROMS ;

typedef struct                /* contains data on chromsome bands */
  { int chrom ;
    float xmin, xmax ;
    KEY key ;
    Array bands ;
  } BANDS ;

typedef struct                /* contains data on homologies */
  { int pbox, x, y ;
    KEY hom, key1, key2 ;
    KEY map1, map2 ;
    float x1, dx1, x2, dx2 ;
    unsigned int flag ;
  } HOMOL ;

typedef struct                /* contains data on axes of boxes in */
  { int box, x, y ;           /* Oxford Grid - handy when cell is empty */
    KEY key, map1, map2 ;
  } BOX ;

typedef struct                /* contains data on map_sets in a Species Grid */
  { int homol ;
    KEY key ;
    Array maps ;
  } SPEC ;

typedef struct                /* contains data on the location of a point */
  { KEY map ;                 /* boolean TRUE if a position on the map is known */
    BOOL pos ;
    float x1, dx1 ;
  } POS ;

typedef struct               /* used for sorting on the location of points on a */
{ KEY map ;                  /* Species Grid so that segment lines can be drawn */
  float y ;
} POINT ;

typedef struct               /* used for drawing homologies in different colours */
{ int box ;                  /* on a Species Grid */
  BOOL hid ;
} HIDBOX ;

extern int  OX_MAGIC ;

	/* per OX flags */
#define FLAG_PAIR_CHROMS        0x00001
#define FLAG_XBANDS             0x00002
#define FLAG_YBANDS             0x00004
#define FLAG_SEGS               0x00008
#define FLAG_XCEN               0x00010
#define FLAG_YCEN               0x00020
#define FLAG_OX_FLIP            0x00040
#define FLAG_PAIR_FLIP          0x00080
#define FLAG_BARS		0x00100
#define FLAG_SHADE		0x00200
#define FLAG_EQUAL              0x00400
#define FLAG_NOSIZE             0x00800
#define FLAG_O2M_BANDS          0x01000
#define FLAG_O2M_SIZE           0x02000
#define FLAG_O2M_CHROM          0x04000
#define FLAG_O2M_EQUAL          0x08000
#define FLAG_SELF               0x10000
#define FLAG_SPECBANDS          0x20000
#define FLAG_SPEC_CHROM         0x40000
#define FLAG_PARALOGY           0x80000

       /* per SEG flags */
#define FLAG_HIDE               0x00001
#define FLAG_LIGHT              0x00002
#define FLAG_SELECT             0x00004
#define FLAG_FISH               0x00008

        /* per CHROMS flags */
#define FLAG_CENTROMERE		0x001
#define FLAG_DARK_BAND		0x002
#define FLAG_NOR		0x004
#define FLAG_SYM                0x008

       /* per HOMOL flags */
#define FLAG_XPOS               0x00010  /* does an homology have an x */
#define FLAG_YPOS               0x00020  /* and a y position ? */
#define FLAG_DUPLICATE		0x00040
#define FLAG_HOMHIDE	        0x00080
#define FLAG_HIGHLIGHT		0x00100


void oxgridCreate (void) ;
void oxgridDisplay (void) ;
void oxhomlist (void) ;
void pairMapCreate (void) ;
void o2mCreate (void) ;
void japanese (char *text, float x, float y) ;
BOOL chkSpecies (KEY key1, KEY key2) ;
void specDisplay (void) ;

#define OXGET(name)      OX ox ; \
                         if (!graphAssFind (&OX_MAGIC, &ox)) \
		           messcrash ("graph not found in %s", name) ; \
			 if (!ox) \
                           messcrash ("%s received a null pointer", name) ; \
                         if (ox->magic != &OX_MAGIC) \
                           messcrash ("%s received a wrong pointer", name)


extern KEY
 _Ox_Grid,
 _Homology,
 _Labelling,
 _Chrom_Pair,
 _History,
 _Prefix,
 _MIM,
 _Autosomal_Chrom,
 _Band_data,
 _Long_name,
 _Loci,
 _Doc,
 _Protein,
 _Group,
 _Pairwise,
 _X_Chrom,
 _Y_Chrom, 
 _Symbol,
 _Chr1,
 _Chr2,
 _Breakpoints,
 _Acquired,
 _Constitutional,
 _Within, 
 _Chrom_Anomaly ;

extern int
 _VHomology_group,
 _VSpecies,
 _VPrefix,
 _VDoc ;




/********* end of file ********/
 
 
 
