/*  File:geneintronwalls.c
 *  Author:  D and J Thierry-Mieg (mieg@ncbi.nlm.nih.gov)
 *  Copyright (C) D and J Thierry-Mieg 2005
 * -------------------------------------------------------------------
 * AceView is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 * or see the on-line version at http://www.gnu.org/copyleft/gpl.txt
 * -------------------------------------------------------------------
 * This file is part of the AceView project developped by
 *	 D and J Thierry-Mieg (mieg@ncbi.nlm.nih.gov)
 *
 *  This should be linked with the acedb libraries
 *  available from http://www.acedb.org
 *
 *  Please direct all question about this code to
 *  mieg@ncbi.nlm.nih.gov
 */

#include "ac.h"
#include "cdna.h"
#include "vtxt.h"

typedef struct giwStruct { KEY key, map, tg, pg, type ; char subtype[9] ; BOOL isDown ; int a1, a2, nEst, nMrna, nR, nAV, nNM, nXM, x1, x2, nja, njb ; } GIW ;

/******************************************************************/
/******************************************************************/

static AC_DB giwDb (void)
{
  static AC_DB db = 0 ;

  if (!db)  
    db = ac_open_db ("local", 0) ;
  return db ;
}

/******************************************************************/

static int giwOrder (const void *va, const void *vb)
{
  const GIW *a= ( const GIW *)va, *b = ( const GIW *)vb ;
  int n ;
  n = a->map - b->map ; if (n) return n ;
  n = a->isDown - b->isDown ; if (n) return n ;
  n = a->a1 - b->a1 ; if (n) return n ;
  n = a->a2 - b->a2 ; if (n) return n ;
  return 0 ;
}

/******************************************************************/
/******************************************************************/
/* 2008_04_11
 * this code takes on the jensen database versus the next code
 * an unbelievable time, around 400 times too slow !
 * but the outputs seem identical
 *   Y 0s 74s reproduced 12 times
 *  21 1s  2h
 *   9 4s 10h
 *   1 7s >15h killed
 *
 * i need to see  if it comes from bql or from aceC
 */
/******************************************************************/
#include "bs.h"
#include "query.h"
#include "../whooks/systags.h"
#include "../whooks/tags.h"

static int giwAddIntronClass2 (AC_DB db, KEY gene)
{
   AC_HANDLE h = ac_new_handle () ;
   int ii, jj, nii, nn = 0, npgs, ntgs, pg1, pg2, tg1, tg2, a1, a2, x1, x2, oldx1, oldx2 ;
   KEYSET pgs = 0, tgs = 0 ;
   Array aa = 0, introns = 0, tmp = 0 ;
   OBJ Pg = 0, Tg = 0 ;
   KEY map, est, tg, pg, nmId ;
   BSunit *bb ;
   GIW *up, *vp ;
   const char *ccp ;
   char *mapNam, *cp ;
   KEY _NM_id = str2tag("NM_id") ;

   tgs = queryKey (gene, "{CLASS Transcribed_gene} SETOR {>Transcribed_gene}; Intron_boundaries && ! shedded_from") ;
   ntgs = keySetMax (tgs) ;
   pgs = queryKey (gene, "{CLASS Predicted_gene} SETOR {>Genefinder};{source_exons} SETOR {>Genefinder ; source_exons}") ;
   npgs = keySetMax (pgs) ;
   if (!ntgs && !npgs)
     goto done ;
   aa = arrayHandleCreate (1000, BSunit, h) ;

   /* register all the predicted introns */
   tmp = arrayHandleCreate (60, GIW, h) ;
   introns = arrayHandleCreate (60, GIW, h) ;
   nii = 0 ;
   if (npgs)
     {
       for (up = 0, ii = 0 ; ii < npgs ; ii++)
	 {
	   pg = keySet (pgs, ii) ;
	   Pg = bsCreate (pg) ;
	   if (!bsGetKey (Pg, _IntMap, &map) ||
	       !bsGetData (Pg, _bsRight, _Int, &pg1) ||
	       !bsGetData (Pg, _bsRight, _Int, &pg2)
	       )
	     { bsDestroy (Pg) ; continue ; }
	   bsGetArray (Pg, _Source_Exons, aa, 2) ;
	   /* attention, they exons may be in disorder */
	   tmp = arrayReCreate (tmp, 200, GIW) ;
	   for (jj = 0 ; jj < arrayMax(aa) ; jj += 2)
	     {
	       bb = arrp (aa, jj, BSunit) ;
	       vp = arrayp (tmp, jj/2, GIW) ;
	       vp->x1 = bb[0].i ; vp->x2 = bb[1].i ; /* i mean x */
	       vp->a1 = bb[0].i ; vp->a2 = bb[1].i ; /* set a for giwOrder */
	     }
	   arraySort (tmp, giwOrder) ;
	   arrayCompress (tmp) ;
	   for (jj = 0 ; jj + 1 < arrayMax(tmp) ; jj++)
	     {
	       vp = arrp (tmp, jj, GIW) ; /* tmp is an array of exons */
	       x1 = vp->x2 + 1 ; x2 = (vp+1)->x1 - 1 ;
	       if (x1 > x2) continue ;
	       bb = arrp (aa, jj, BSunit) ;
	       up = arrayp (introns, nii++, GIW) ;
	       up->map = map ;
	       up->pg = pg ;
	       if (bsGetKey (Pg, _NM_id, &nmId))
		 {
		   if (nmId && *name(nmId) == 'N') up->nNM++ ;
		   if (nmId && *name(nmId) == 'X') up->nXM++ ;
		   if (nmId && *name(nmId) == 'A') up->nAV++ ;
		 }
	       up->x1 = x1 ; up->x2 = x2 ;
	       if (pg1 < pg2)
		 {
		   up->a1 = pg1 + x1 - 1 ;
		   up->a2 = pg1 + x2 - 1 ;
		 }
	       else
		 {
		   up->a1 = pg1 - x1 + 1 ;
		   up->a2 = pg1 - x2 + 1 ;
		 }
	     }
	   bsDestroy (Pg) ;
	 }
     }

   /* register all the observed introns */
   if (ntgs)
     {
       KEY type ;
       int dx = 0 ;
       for (ii = 0 ; ii < ntgs ; ii++)
	 {
	   tg = keySet (tgs, ii) ;
	   Tg = bsCreate (tg) ;
	   if (!bsGetKey (Tg, _IntMap, &map) ||
	       !bsGetData (Tg, _bsRight, _Int, &tg1) ||
	       !bsGetData (Tg, _bsRight, _Int, &tg2)
	       )
	     { bsDestroy (Tg) ; continue ; }
	   up = 0 ; oldx1 = oldx2 = 0 ;
	   bsGetArray (Tg, _Intron_boundaries, aa, 6) ;
	   for (jj = 0 ; jj < arrayMax(aa) ; jj += 6)
	     {
	       bb = arrp (aa, jj, BSunit) ;
	       type = bb[0].k ;
	       if (type == _Other) dx = 1 ;
	       else dx = 0 ;
	       x1 = bb[2+dx].i ; x2 = bb[3+dx].i ;
	       if (!x1 || !x2) continue ;
	       if (x1 != oldx1 || x2 != oldx2)
		 {
		   oldx1 = x1 ; oldx2 = x2 ;
		   up = arrayp (introns, nii++, GIW) ;
		   up->tg = tg ;
		   up->map = map ;
		   up->type = type ;
		   if (type == _Other) strncpy (up->subtype, bb[1].s, 8) ;
		   up->x1 = x1 ; up->x2 = x2 ;
		   if (tg1 < tg2)
		     {
		       up->a1 = a1 = tg1 + x1 - 1 ;
		       up->a2 = a2 = tg1 + x2 - 1 ;
		     }
		   else
		     {
		       up->a1 = tg1 - x1 + 1 ;
		       up->a2 = tg1 - x2 + 1 ;
		     }
		 }
	       /* increment the observed support */
	       est = bb[4+dx].k ;
	       if (est)
		 {
		   if (!strncmp (name(est), "NM_", 3) ||
		       !strncmp (name(est), "XM_", 3) ||
		       !strncmp (name(est), "NR_", 3) ||
		       !strncmp (name(est), "XR_", 3)
		       )
		     up->nR++ ;
		   else 
		     up->nEst++ ;
		   /* Jensen MAQC 454 sequencing 2008_04_08 */
		   ccp = name(est) ;
		   if (ccp[0] == 'J' && ccp[2] == 'A' && ccp[4] == '_') up->nja++ ;
		   if (ccp[0] == 'J' && ccp[2] == 'B' && ccp[4] == '_') up->njb++ ;
		 }
	     }
	   bsDestroy (Tg) ;
	 }
     }

   if (nii)
     {
       vTXT bfr = vtxtHandleCreate (h) ;

       arraySort (introns, giwOrder) ;
       arrayCompress (introns) ;
       nii = arrayMax (introns) ;

       for (ii = 0 ; ii < arrayMax (introns) ; ii++)
	 {
	   up = arrp (introns, ii, GIW) ;
	   mapNam = ac_protect(name(up->map), h) ;
	   cp = mapNam + strlen (mapNam) - 1 ; *cp = 0 ;
	   vtxtPrintf (bfr, "Intron %s__%d_%d\"\n", mapNam, up->a1, up->a2) ;
	   if (up->tg)
	     {
	       cp = ac_protect(name(up->tg), h) ;
	       vtxtPrintf (bfr, "From_gene %s %d %d\n",  cp, up->x1, up->x2) ;
	     }
	   if (up->pg)
	     {
	       cp = ac_protect(name(up->pg), h) ;
	       vtxtPrintf (bfr, "From_genefinder %s %d %d\n",  cp, up->x1, up->x2) ;
	     }
	   if (up->nAV) vtxtPrint (bfr, "AV\n") ;
	   if (up->nNM) vtxtPrint (bfr, "NM\n") ;
	   if (up->nXM) vtxtPrint (bfr, "XM\n") ;
	   if (up->nR) vtxtPrintf (bfr, "RefSeq %d\n", up->nR) ; 
	   if (up->nMrna) vtxtPrintf (bfr, "mRNA %d\n", up->nMrna) ; 
	   if (up->nEst) vtxtPrintf (bfr, "EST %d\n", up->nEst) ; 
	   if (up->nja) vtxtPrintf (bfr, "Brain %d\n", up->nja) ; 
	   if (up->njb) vtxtPrintf (bfr, "Stratagene %d\n", up->njb) ; 
	   
	   ccp = name(up->map) ;
	   vtxtPrintf (bfr, "IntMap %s %d %d\n", ac_protect(ccp,h), up->a1, up->a2) ; 
	   
	   if (up->type == _Other && *up->subtype) vtxtPrintf (bfr, "Type %s %s\n", name(up->type), up->subtype) ;
	   else if (up->type) vtxtPrintf (bfr, "Type %s\n", name(up->type)) ;
	   
	   /*       ccp = ac_table_printable (tbl, ir, 10, 0) ;
		    if (ccp && strcmp (ccp, "FALSE") && strcmp (ccp, "NULL")) vtxtPrintf (bfr, "Tissue %s\n", ccp) ;
	   */
	   vtxtPrint (bfr, "\n") ;
	 }
      nn++ ;
      ac_parse (db, vtxtPtr (bfr), 0, 0, h) ;
     }
 done:
   keySetDestroy (pgs) ;
   keySetDestroy (tgs) ;
   ac_free (h) ;
   return nn ;
} /* giwAddIntronClass2  */

/******************************************************************/

static int giwAddProbeWalls (AC_DB db, KEY key)
{
   AC_HANDLE h = ac_new_handle () ;
   int nn = 0, i, j ;
   AC_TABLE tbl = 0 ; /* tblv3, tblv5 */
   GIW *up, *vp ;
   Array introns = 0, exons = 0 ;
   KEYSET walls = 0 ;
   vTXT bfr = vtxtHandleCreate (h) ;
   vTXT qry = vtxtHandleCreate (h) ;
   vTXT qryV3 = vtxtHandleCreate (h) ;
   vTXT qryV5 = vtxtHandleCreate (h) ;
   const char *type ;
   int ir, g1, g2, m1, x1, x2, dx, ni, nx ;


   vtxtPrintf (qry, " select g, g1, g2, m1, m2, x1, x2, type") ;
   vtxtPrintf (qry, " from g in class gene where g like %s", ac_protect(name(key), h)) ;
   vtxtPrintf (qry, " , gm in g->intmap, g1 in gm[1], g2 in g1[1] ") ;
   vtxtPrintf (qry, " , tg in g->transcribed_gene ") ;
   vtxtPrintf (qry, " , m in tg->mrna where m ") ;
   vtxtPrintf (qry, " , mm in m->intmap, m1 in mm[1], m2 in m1[1] ") ;
   vtxtPrintf (qry, " , x1 in m->splicing, x2 in x1[1], type in x2[3]" ) ;
  
   vtxtPrintf (qryV3, " select g, g1, g2, m, m1, m2, x1") ;
   vtxtPrintf (qryV3, " from g in class gene where g like %s", ac_protect(name(key), h)) ;
   vtxtPrintf (qryV3, " , gm in g->intmap, g1 in gm[1], g2 in g1[1] ") ;
   vtxtPrintf (qryV3, " , tg in g->transcribed_gene ") ;
   vtxtPrintf (qryV3, " , m in tg->mrna where m#Valid3p") ;
   vtxtPrintf (qryV3, " , mm in m->intmap, m1 in mm[1], m2 in m1[1] ") ;
   vtxtPrintf (qryV3, " , x1 in m->valid3p") ;
   
   vtxtPrintf (qryV5, " select g, g1, g2, m1, m2, x1") ;
   vtxtPrintf (qryV5, " from g in class gene where g like %s", ac_protect(name(key), h)) ;
   vtxtPrintf (qryV5, " , gm in g->intmap, g1 in gm[1], g2 in g1[1] ") ;
   vtxtPrintf (qryV5, " , tg in g->transcribed_gene ") ;
   vtxtPrintf (qryV5, " , m in tg->mrna where m#Valid5p") ;
   vtxtPrintf (qryV5, " , mm in m->intmap, m1 in mm[1], m2 in m1[1] ") ;
   vtxtPrintf (qryV5, " , x1 in m->valid5p") ;
   
   tbl = ac_bql_table (db, vtxtPtr (qry), 0, 0, 0, h) ;
   /*
     tblV3 = ac_bql_table (db, vtxtPtr (qryV3), 0, 0, 0, h) ;
     tblV5 = ac_bql_table (db, vtxtPtr (qryV5), 0, 0, 0, h) ;
   */

   /* un donneur qui a 2 accepteurs ou qui est dans un exon compte
    * un accepteur qui a 2 donneurs ou qui est dans un exon compte
    * un flag valid5p ou valid3p compte a 50bp pres
    * les extremites hors tout du gene comptent
    * les unspliced ne comptent pas
    * basta
    */

   /* cree un array des donneurs
    * compte ses accepteurs, si > 1, flag
    *                        sinon cheche si dans un exon
    * idem pour les accepteurs
    * rajouter les flags et les extremites hors tout
    */
   /* construire l'array des introns et des exons */
   ni = nx = 0 ;
   exons = arrayHandleCreate (12, GIW, h) ;
   introns = arrayHandleCreate (12, GIW, h) ;
   for (ir = 0 ; ir < tbl->rows ; ir++)
     {
       g1 = ac_table_int (tbl, ir, 1, 0) ;
       g2 = ac_table_int (tbl, ir, 2, 0) ;
       m1 = ac_table_int (tbl, ir, 3, 0) ;
       /*  m2 = ac_table_int (tbl, ir, 4, 0) ; */
       x1 = ac_table_int (tbl, ir, 5, 0) ;
       x2 = ac_table_int (tbl, ir, 6, 0) ;
       type = ac_table_printable (tbl, ir, 7, 0) ;
       if (g1 < g2) dx = m1 - g1 ;
       else dx = g1 - m1 ;
       if (type && strstr (type, "tron"))
	 {
	   up = arrayp (introns, ni++, GIW) ;
	   up->x1 = dx + x1 ; up->x2 = dx + x2 ;
	 }
       if (type && strstr (type, "xon"))
	 {
	   up = arrayp (exons, nx++, GIW) ;
	   up->x1 = dx + x1 ; up->x2 = dx + x2 ;
	 }
     }
   if (nx) { arraySort (exons, giwOrder) ; arrayCompress (exons) ; }
   if (ni) { arraySort (introns, giwOrder) ; arrayCompress (introns) ; }
   nx = arrayMax (exons) ;
   ni = arrayMax (introns) ;
   /* remplir l'array des walls avec les donneurs qui ont 2 accepteurs */
   for (i = 0 ; i < ni ; i++)
     {
       up = arrp (introns, i, GIW) ;
       x1 = up->x1 ; x2 = up->x2 ;
       if (!keySetFind (walls, x1, 0))
	 for (j = i+1, vp = up + 1 ; j < ni && vp->x1 == x1 ; j++, vp++)
	   if (vp->x2 != x2)
	     {
	       if (j > i+1) { keySetInsert (walls, x1) ;keySetInsert (walls, x2) ; }
	       keySetInsert (walls, vp->x2) ;
	     }
     }
   if (nn)
     {
       vtxtPrint (bfr, "\n") ;
       ac_parse (db, vtxtPtr (bfr), 0, 0, h) ;
     }

   ac_free (h) ;
   return nn ;
} /* giwAddProbeWalls  */

/******************************************************************/
/******************************************************************/
/* called with ks from the acembly menu, with key from makemrna.c */
static int giwAddStuff (KEYSET ks, KEY key, int stuff)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_DB db = giwDb () ;
  int ii, nn = 0 ;
  KEY *kp ;

  if (!ks && key)
    { 
      ks = keySetHandleCreate (h) ;
      keySet (ks, 0) = key ;
    }
  if (ks && keySetMax(ks))
    {
      kp = arrp (ks, 0, KEY) - 1 ;
      ii = keySetMax(ks) ;
      while (kp++, ii--)
	switch (stuff)
	  {
	  case 1: nn += giwAddIntronClass2 (db, *kp) ; break ;
	  case 2: nn += giwAddProbeWalls (db, *kp) ; break ;
	  }
    }

  ac_free (h) ;
  
  return nn ;  
} /* geneAddKeysetIntronClass */

int giwAddKeysetProbeWalls (KEYSET ks, KEY key)
{
  return 0 ;
  return giwAddStuff (ks, key, 2) ;
}

int giwAddKeysetIntronClass (KEYSET ks, KEY key)
{
  return giwAddStuff (ks, key, 1) ;
}

/* look for cassette introns */
int giwAddIntronHierarchy (void)
{
  AC_HANDLE h = ac_new_handle () ;
  int ii, jj, nn = 0, nii = 0 ; 
  GIW *up, *vp ;
  KEY key ;
  OBJ obj = 0 ;
  BSunit *bb ;
  Array aa = 0, introns = 0 ;
  KEYSET ks = query (0, "Find intron") ;
  vTXT bfr = vtxtHandleCreate (h) ;

  introns = arrayHandleCreate (100000, GIW, h) ;
  aa = arrayHandleCreate (1000, BSunit, h) ;

  for (ii = jj = 0 ; ii < keySetMax (ks) ; ii++)
    {
      key = keySet (ks, ii) ;
      if ((obj = bsCreate (key)))
	{
	  bsGetArray (obj, _IntMap, aa, 3) ;
	  if (arrayMax(aa) >= 3)
	    {
	      up = arrayp (introns, jj++, GIW) ;
	      bb = arrp (aa, 0, BSunit) ;
	      up->key = key ;
	      up->map = bb[0].k ; up->a1 = bb[1].i ; up->a2 = bb[2].i ;
	      if (up->a1 < up->a2) 
		up->isDown = TRUE ;
	      else
		{
		  int tmp = up->a1 ; up->a1 = up->a2 ; up->a2 = tmp ;
		  up->isDown = FALSE ;
		}
	    }
	  bsDestroy (obj) ;
	}
    }
  
  arraySort (introns, giwOrder) ;
  arrayCompress (introns) ;
  nii = arrayMax (introns) ;

  for (ii = 0 ; ii < nii ; ii++)
    {
      int n1 = 0 ;
      up = arrp (introns, ii, GIW) ;
      for (jj = ii + 1 ; jj < nii ; jj++)
	{
	  vp = arrp (introns, jj, GIW) ;
	  if (vp->a2 > up->a2 || vp->isDown != up->isDown || vp->map != up->map) break ;
	  if (!n1++)
	    {
	      vtxtPrintf (bfr, "Intron %s\n", ac_protect (name(up->key), h)) ;
	      vtxtPrintf (bfr, "Includes %s\n", ac_protect (name(vp->key), h)) ;
	    }
	}
      if (n1)
	{
	  nn++ ;
	  vtxtPrint (bfr, "\n") ;
	}
    }

  if (nn)
    ac_parse (giwDb (), vtxtPtr (bfr), 0, 0, h) ;
    
   keySetDestroy (ks) ;
   ac_free (h) ;
   return nn ;
} /* giwAddIntronHierarchy */

/*********************************************************************/
/*********************************************************************/

