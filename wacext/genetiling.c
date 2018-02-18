/*
 * Export from aceview a separation of the genome
 * into : intergenic/intronic/exonic regions
 *
 * for each region defined on the genome
 * associate a keyset of genes
 */


#include "ac.h"
#include "freeout.h"

typedef enum {ZERO, INTERGENIC, INTRON, EXON, CODINGEXON } zTYPE ;
char *zType[] = {"ZERO", "Intergenic", "Intronic", "Exonic non-coding", "Exonic coding" }  ;

typedef struct zoneStruct { AC_DB db ; Array aa ; const char *gene ; AC_HANDLE h ; BOOL mixStrands, noCloud, codingIgnored, remap, pg, probeType, probe2doubleTile;} ZZ ;
typedef struct zzoneStruct { int a1, a2 ; BOOL isDown ; zTYPE type ; KEY map, gId[12], gene[12] ; int nGid, nGene ; int touchedBy ; } ZZZ ;
typedef struct zzP2T { KEY key, map ; char type ; int a1, a2 ; } P2T ;

static void zzExport (ZZ *zz) ;

/*************************************************************************************/

static void zzRemap (ZZ *zz)
{
  AC_HANDLE h = ac_new_handle () ;
  int ir, jr, a1, a2, x1, x2, u1, u2, ne, ln ;
  const char *tile, *probe ;
  AC_TABLE tiles, probes ;
  AC_KEYSET aks ;
  KEYSET ks = keySetHandleCreate (h) ;
  const char *mapP, *mapT ;

  aks = ac_dbquery_keyset (zz->db, "Find sequence Is_Cosmid && IntMap", h) ;
  tiles = ac_bql_table (zz->db, "select s,m,a1,a2 from s in @, m in s->intmap, a1 in m[1], a2 in m[2] "
		      , aks, 0, 0, h) ;
  aks = ac_dbquery_keyset (zz->db, "Find probe IS R_* && Genome_approximate_hit", h) ;
  probes = ac_bql_table (zz->db, "select s,m,x1,x2, ne, ln from s in @, m in s->Genome_approximate_hit, x1 in m[1], x2 in x1[1], ne in x1[2], ln in x1[3] "
		      , aks, 0, 0, h) ;

  for (ir = 0 ; tiles && ir < tiles->rows ; ir++) 
    {
      tile = ac_table_printable (tiles, ir, 0, 0) ;
      mapT = ac_table_printable (tiles, ir, 1, "") ;
      a1 = ac_table_int (tiles, ir, 2, 0) ;
      a2 = ac_table_int (tiles, ir, 3, 0) ;
      for (jr = 0 ; probes && jr < probes->rows ; jr++) 
	{
	  if (keySet (ks, jr))
	    continue ; 
	  mapP = ac_table_printable (probes, jr, 1, "") ;
	  if (strcasecmp (mapP,mapT)) /* they are not in the same class, but have the  same name */
	    continue ;
	  x1 = ac_table_int (probes, jr, 2, 0) ;
	  x2 = ac_table_int (probes, jr, 3, 0) ;
	  if (x1 >= a1 && x2 >= a1 && x1 <= a2 && x2 <= a2)
	    {
	      ne = ac_table_int (probes, jr, 4, 0) ;
	      ln = ac_table_int (probes, jr, 5, 0) ;
	      if (a1 < a2)
		{ u1 = x1 - a1 + 1 ; u2 = x2 - a1 + 1 ; }
	      else
		{ u1 = a1 - x1 + 1 ; u2 = a1 - x2 + 1 ; }
	      probe = ac_table_printable (probes, jr, 0, 0) ;
	      keySet (ks, jr) = 1 ;
	      freeOutf ("Probe %s\nTiling_hit %s %d %d %d %d\n\n", probe, tile, u1, u2, ne, ln) ;
	      freeOutf ("Sequence %s\nProbe_hit %s %d %d %d %d\n\n", tile, probe, u1, u2, ne, ln) ;
	    }
	}
    }

  ac_free (h) ;
} /* zzRemap */

/*************************************************************************************/

static void zzProbeTypeExons (AC_OBJ probe, AC_OBJ gene, int p1, int p2, AC_HANDLE h)
{
  int ir, m1, m2, a1, a2, b1, b2, p, x1, x2, u1, u2, pass, cds1, cds2 ;
  AC_OBJ mrna = 0, product = 0 ;
  AC_TABLE tbl ;
  AC_ITER iter ;
  const char *ccp ;
  BOOL pDown, mDown ;

  pDown =  p1 < p2 ? TRUE : FALSE ;
  iter = ac_objquery_iter (gene, ">transcribed_gene; >mrna", h) ;
  while (ac_free (mrna), (mrna = ac_iter_obj (iter)))
    {
      tbl = ac_tag_table (mrna, "IntMap", 0) ;
      m1 = ac_table_int (tbl, 0, 1, 0) ;
      m2 = ac_table_int (tbl, 0, 2, 0) ;
      ac_free (tbl) ;

      mDown =  m1 < m2 ? TRUE : FALSE ;

      if (mDown && m1 <= p1 && m2 >= p1) ;
      else if (mDown && m1 <= p2 && m2 >= p2) ;
      else if (!mDown && m2 <= p1 && m1 >= p1) ;
      else if (!mDown && m2 <= p2 && m1 >= p2) ;
      else continue ;

      tbl = ac_tag_table (mrna, "Product", 0) ;
      cds1 = cds2 = 0 ;
      for (ir = 0 ; tbl && !cds1 && ir < tbl->rows ; ir++)
	{
	  product = ac_table_obj (tbl, ir, 0, 0) ;
	  if (ac_has_tag (product, "Best_product") &&
	      ac_has_tag (product, "Good_product"))
	    {
	      cds1 = ac_table_int (tbl, ir, 1, 0) ;
	      cds2 = ac_table_int (tbl, ir, 2, 0) ;
	    }
	  ac_free (product) ;
	}
      ac_free (tbl) ;
      tbl = ac_tag_table (mrna, "Splicing", 0) ;
      for (ir = 0 ; tbl && ir < tbl->rows ; ir++)
	{
	  ccp = ac_table_printable (tbl, ir, 4, "") ;
	  if (!strstr (ccp, "xon"))
	    continue ;
 	  a1 = b1 = ac_table_int (tbl, ir, 0, 0) ;
	  a2 = b2 = ac_table_int (tbl, ir, 1, 0) ;
	  x1 = ac_table_int (tbl, ir, 2, 0) ;
	  x2 = ac_table_int (tbl, ir, 3, 0) ;
	  if (x2 < cds1) { b1 = b2 = 0 ; }
	  else if (x1 > cds2) { b1 = b2 = 0 ; }
	  else 
	    { 
	      if (x1 < cds1 && x2 >= cds1) b1 += cds1 - x1 ;
	      if (x1 <= cds2 && x2 > cds2) b2 += cds2 - x2 ;
	    }

	  if (mDown) { u1 = m1 + a1 - 1 ;  u2 = m1 + a2 - 1 ; }
	  else { u1 = m1 - a2 + 1 ;  u2 = m1 - a1 + 1 ; }
	  for (pass = 1 ; pass <= 2 ; pass++)
	    { 
	      p = pass == 1 ? p1 : p2 ;
	      if (p >= u1 && p <= u2) 
		{
		  if (pDown == mDown)
		    freeOutf ("Exon_f%d\n", pass) ;
		  else
		    freeOutf ("Exon_r%d\n", pass) ;
		}
	    }
	  if (b1 < b2)
	    {	    
	      if (mDown) { u1 = m1 + b1 - 1 ;  u2 = m1 + b2 - 1 ; }
	      else { u1 = m1 - b2 + 1 ;  u2 = m1 - b1 + 1 ; }
	      if (u1 < u2)
		for (pass = 1 ; pass <= 2 ; pass++)
		  { 
		    p = pass == 1 ? p1 : p2 ;
		    if (p >= u1 && p <= u2) 
		      {
			if (pDown == mDown)
			  freeOutf ("Coding_f%d\n", pass) ;
			else
			  freeOutf ("Coding_r%d\n", pass) ;
		      }
		  }
	    }
	}
      ac_free (tbl) ;
    }
} /* zzProbeTypeExons */

/*************************************************************************************/

static void zzProbeType (ZZ *zz)
{
  AC_HANDLE h = 0, h0 = ac_new_handle () ;
  int ir, nn, a1, a2, p, p1, p2, g1, g2, pass ;
  AC_OBJ gene = 0, probe = 0 ;
  AC_TABLE tiles, tbl ;
  AC_ITER genes ;
  AC_KEYSET aks ;
  KEY mapKey ;
  BOOL pDown, gDown, ok ;

  /* it is more optimal to work by tiles */
  aks = ac_dbquery_keyset (zz->db, "Find sequence Is_Cosmid && IntMap && probe_hit", h0) ;
  tiles = ac_bql_table (zz->db, "select s, sm, a1, a2, p, p1,p2 from s in @active:1, sm in s->intmap, a1 in sm[1], a2 in a1[1], p in s->probe_hit, p1 in p[1], p2 in p1[1] "
		      , aks, 0, 0, h) ;

  for (ir = nn = 0 ; tiles && ir < tiles->rows ; ac_free (h), ir++) 
    { 
      h = ac_new_handle () ;
      nn = 0 ;
      mapKey = ac_table_key (tiles, ir, 1, 0) ;
      a1 = ac_table_int (tiles, ir, 2, 0) ;
      a2 = ac_table_int (tiles, ir, 3, 0) ;
      
      probe = ac_table_obj (tiles, ir, 4, h) ;
      if (0 && strcmp (ac_name(probe), "BBB_1560_754_1008"))
	continue ;
      p1 = ac_table_int (tiles, ir, 5, 0) ;
      p2 = ac_table_int (tiles, ir, 6, 0) ;
      /* we are in cosmid coords, move to intmap coords */
      if (a1 < a2)
	{ p1 = a1 + p1 - 1 ; p2 = a1 + p2 - 1 ; }
      else
	{ p1 = a1 - p1 + 1 ; p2 = a1 - p2 + 1 ; }
      pDown = p1 < p2 ? TRUE : FALSE ;
      /* place the exploring point 10bp inside the probe */
      if (pDown) { p1 += 5 ; p2 -= 5 ; }
      else  { p1 -= 5 ; p2 += 5 ; }

      genes = ac_objquery_iter (probe,
				">Tiling_hit; {Is_cosmid} $| {>in_junction} $| {>source;>subsequence;>Is_gene_tile};>genes"
				, h
				) ;
      while (ac_free (gene) , (gene = ac_iter_obj (genes)))
	{
	  if (0 && strcmp (ac_name(gene), "GDI1"))
	    continue ;
	  tbl = ac_tag_table (gene, "IntMap", h) ;
	  if (!tbl || mapKey != ac_table_key (tbl, 0, 0, 0))
	    continue ;
	  g1 = ac_table_int (tbl, 0, 1, 0) ;
	  g2 = ac_table_int (tbl, 0, 2, 0 );
	  gDown = g1 < g2 ? TRUE : FALSE ;
	  if (!gDown) { int dummy = g1 ; g1 = g2 ; g2 = dummy ; }
	  for (ok = 0, pass = 1 ; pass <= 2 ; pass++)
	    { 
	      p = pass == 1 ? p1 : p2 ;
	      if (p >= g1 && p <= g2) 
		{
		  if (!nn++)
		    freeOutf ("Probe %s\n", ac_protect (ac_name(probe), h)) ;
		  if (ac_has_tag (gene, "Cloud_gene"))
		    {
		      if (pDown == gDown)
			{
			  freeOutf ("Cloud_f%d\n", pass) ;
			  freeOutf ("Cloud_gene %s\n", ac_protect (ac_name(gene), h)) ;
			}
		      else
			{
			  freeOutf ("Cloud_r%d\n", pass) ;
			  freeOutf ("Cloud_antigene %s\n", ac_protect (ac_name(gene), h)) ;
			}
		    }
		  else
		    {
		      ok++ ;
		      if (pDown == gDown)
			{
			  freeOutf ("Gene_f%d\n", pass) ;
			  freeOutf ("Gene %s\n", ac_protect (ac_name(gene), h)) ;
			  {
			    int igid ;
			    AC_TABLE gids = ac_tag_table (gene, "GeneId", h) ;
			    for (igid = 0 ; gids && igid < gids->rows ; igid++)
			      freeOutf ("GeneId %s\n", ac_table_printable (gids, igid, 0, "")) ;
			  }
			}
		      else
			{
			  freeOutf ("Gene_r%d\n", pass) ; 
			  freeOutf ("Antigene %s\n", ac_protect (ac_name(gene), h)) ;
			  {
			    int igid ;
			    AC_TABLE gids = ac_tag_table (gene, "GeneId", h) ;
			    for (igid = 0 ; gids && igid < gids->rows ; igid++)
			      freeOutf ("AntigeneId %s\n", ac_table_printable (gids, igid, 0, "")) ;
			  }
			}
		    }
		}
	    }
	  if (ok) /* inside a gene, look for exons/cds/utr */
	    zzProbeTypeExons (probe, gene, p1, p2, h) ;
	}
      if (nn) freeOutf ("\n") ;
    }
  
  ac_free (h0) ;
} /* zzProbeType */

/*************************************************************************************/

static void zzGetTgData (ZZ *zz)
{
  int nn = 0, i, ir, g1, g2, a1, a2, ngid ;
  const char *ccp ;
  ZZZ *up ;
  BOOL isCloud = FALSE ;
  AC_OBJ Gene = 0, tg = 0 ;
  AC_ITER iter ;
  AC_TABLE tbl = 0 ;
  KEY gene, map, gids[12] ;
  
  if ((iter = ac_dbquery_iter (zz->db, messprintf ("Find tg %s", zz->gene),  0)))
    while (ac_free (tg) , tg = ac_iter_obj (iter))
      {
	gene = ac_tag_key (tg, "Gene", 0) ;

	ngid = 0 ;
	if ((Gene = ac_tag_obj (tg, "Gene", 0)))
	  {
	    isCloud = ac_has_tag (Gene, "Cloud_gene") ;
	    tbl = ac_tag_table (Gene, "GeneId", 0) ;
	    for (ir = 0 ; tbl && ir < tbl->rows  && ir < 12; ir++)
	      {
		gids[ngid++] = ac_table_key (tbl, ir, 0, 0) ;	    
	      }
	    ac_free (tbl) ;
	    ac_free (Gene) ;
	  }
	if (isCloud && zz->noCloud)
	  continue ;
	gene = ac_tag_key (tg, "Gene", 0) ;
	map = 0 ; g1 = g2 = -1 ;
	tbl = ac_tag_table (tg, "IntMap", 0) ;
	for (ir = 0 ; tbl && ir < tbl->rows ; ir++)
	  {
	    map = ac_table_key (tbl, ir, 0, 0) ;
	    g1 = ac_table_int (tbl, ir, 1, 0) ;
	    g2 = ac_table_int (tbl, ir, 2, 0) ;
	    
	  }
	ac_free (tbl) ;
	if (g2 == -1)
	  continue ;
	up = arrayp (zz->aa, nn++, ZZZ) ;
	up->gene[0] = gene ;
	up->nGene = 1 ;
	for (i=0 ; i < ngid ; i++)
	  up->gId[i] = gids[i] ;
	up->nGid = ngid ;
	up->map = map ;
	up->type = INTRON ; /*  could be 'gap' so all premessenger not covered by a confirmed intron would be a gap */
	if (g1 < g2)
	  {
	    up->isDown = TRUE ;
	    up->a1 = g1 ;
	    up->a2 = g2 ;
	  }
	else
	  {
	    up->isDown = FALSE ;
	    up->a2 = g1 ;
	    up->a1 = g2 ;
	  }

	tbl = ac_tag_table (tg, "Splicing", 0) ;
	for (ir = 0 ; tbl && ir < tbl->rows ; ir++)
	  {
	    ccp = ac_table_printable (tbl, ir, 2, "") ; 
	    if (! strstr (ccp, "xon"))
	      continue ; /* remove this if you want to label the gaps separatelly */
	    if (strstr (ccp, "xon"))
	    up = arrayp (zz->aa, nn++, ZZZ) ;
	    up->gene[0] = gene ;
	    up->nGene = 1 ;
	    for (i=0 ; i < ngid ; i++)
	      up->gId[i] = gids[i] ;
	    up->nGid = ngid ;
	    up->map = map ;
	    a1 = ac_table_int (tbl, ir, 0, 0) ;
	    a2 = ac_table_int (tbl, ir, 1, 0) ;
	    if (g1 < g2)
	      {
		up->isDown = TRUE ;
		up->a1 = g1 + a1 - 1 ;
		up->a2 = g1 + a2 - 1 ;
	      }
	    else
	      {
		up->isDown = FALSE ;
		up->a2 = g1 - a1 + 1 ;
		up->a1 = g1 - a2 + 1 ;
	      }
	    ccp = ac_table_printable (tbl, ir, 2, "") ;
	    if (strstr (ccp, "xon"))
	      up->type = EXON ;
	    else if (strstr (ccp, "tron"))
	      up->type = INTRON ;
	  }
	ac_free (tbl) ;
      }
  ac_free (iter) ;
} /* zzGetTgData */

/*************************************************************************************/

static void zzGetPgData (ZZ *zz)
{
  int nn = 0, ir, g1, g2, a1, a2 ;
  ZZZ *up ;
  AC_OBJ pg = 0 ;
  AC_ITER iter ;
  AC_TABLE tbl = 0 ;
  KEY map ;
  
  if ((iter = ac_dbquery_iter (zz->db, messprintf ("Find Predicted_gene %s", zz->gene),  0)))
    while (ac_free (pg) , pg = ac_iter_obj (iter))
      {
	map = 0 ; g1 = g2 = -1 ;
	tbl = ac_tag_table (pg, "IntMap", 0) ;
	for (ir = 0 ; tbl && ir < tbl->rows ; ir++)
	  {
	    map = ac_table_key (tbl, ir, 0, 0) ;
	    g1 = ac_table_int (tbl, ir, 1, 0) ;
	    g2 = ac_table_int (tbl, ir, 2, 0) ;
	    
	  }
	ac_free (tbl) ;
	if (g2 == -1)
	  continue ;
	up = arrayp (zz->aa, nn++, ZZZ) ;
	up->map = map ;
	up->type = INTRON ; /*  could be 'gap' so all premessenger not covered by a confirmed intron would be a gap */
	if (g1 < g2)
	  {
	    up->isDown = TRUE ;
	    up->a1 = g1 ;
	    up->a2 = g2 ;
	  }
	else
	  {
	    up->isDown = FALSE ;
	    up->a2 = g1 ;
	    up->a1 = g2 ;
	  }

	tbl = ac_tag_table (pg, "Source_exons", 0) ;
	for (ir = 1 ; tbl && ir < tbl->rows ; ir++)
	  {
	    a1 = ac_table_int (tbl, ir-1, 1, 0) + 1 ;
	    a2 = ac_table_int (tbl, ir, 0, 0) - 1 ;
	    if (a2 < a1+10) continue ;
	    up = arrayp (zz->aa, nn++, ZZZ) ;
	    up->map = map ;
	    if (g1 < g2)
	      {
		up->isDown = TRUE ;
		up->a1 = g1 + a1 - 1 ;
		up->a2 = g1 + a2 - 1 ;
	      }
	    else
	      {
		up->isDown = FALSE ;
		up->a2 = g1 - a1 + 1 ;
		up->a1 = g1 - a2 + 1 ;
	      }
	    /* ccp = ac_table_printable (tbl, ir, 2, "") ; */
	    up->type = INTRON ;
	  }
	for (ir = 0 ; tbl && ir < tbl->rows ; ir++)
	  {
	    a1 = ac_table_int (tbl, ir, 0, 0) ;
	    a2 = ac_table_int (tbl, ir, 1, 0) ;
	    
	    if (a2 < a1+10) continue ;
	    up = arrayp (zz->aa, nn++, ZZZ) ;
	    up->map = map ;
	    if (g1 < g2)
	      {
		up->isDown = TRUE ;
		up->a1 = g1 + a1 - 1 ;
		up->a2 = g1 + a2 - 1 ;
	      }
	    else
	      {
		up->isDown = FALSE ;
		up->a2 = g1 - a1 + 1 ;
		up->a1 = g1 - a2 + 1 ;
	      }
	    /* ccp = ac_table_printable (tbl, ir, 2, "") ; */
	    up->type = EXON ;
	  }
	ac_free (tbl) ;
      }
  ac_free (iter) ;
} /* zzGetPgData */

/*************************************************************************************/

static void zzGetCodingData (ZZ *zz)
{
  int nn = arrayMax (zz->aa), i, ir, jj, g1, g2, a1, a2, ngid ;
  const char *ccp; char cc ;
  ZZZ *up ;
  BOOL isCloud = FALSE ;
  AC_OBJ tg = 0, mrna = 0, product = 0, Gene = 0 ;
  AC_ITER iter ;
  AC_TABLE tbl = 0 ;
  KEY gene, map, gids[12] ;
  char buf1a[8], buf3a[8], buf5a[8] ;
  char buf1b[8], buf3b[8], buf5b[8] ;
  
  buf1a[1] = 0 ; buf3a[0] = '3' ; buf3a[2] = 0 ; buf5a[0] = '5' ; buf5a[2] = 0 ; 
  buf1b[1] = 0 ; buf3b[0] = '3' ; buf3b[2] = 0 ; buf5b[0] = '5' ; buf5b[2] = 0 ; 

  if ((iter = ac_dbquery_iter (zz->db, messprintf ("Find Product IS %s* && Good_product ; > mrna", zz->gene),  0)))
    while (ac_free (mrna) , mrna = ac_iter_obj (iter))
      {
	map = 0 ; g1 = g2 = -1 ;
	tbl = ac_tag_table (mrna, "IntMap", 0) ;
	for (ir = 0 ; tbl && ir < tbl->rows ; ir++)
	  {
	    map = ac_table_key (tbl, ir, 0, 0) ;
	    g1 = ac_table_int (tbl, ir, 1, 0) ;
	    g2 = ac_table_int (tbl, ir, 2, 0) ;
	    
	  }
	ac_free (tbl) ;
	if (g2 == -1)
	  continue ;

	tg = ac_tag_obj (mrna, "From_gene", 0) ;
	gene = ac_tag_key (tg, "Gene", 0) ;
	ngid = 0 ;
	if ((Gene = ac_tag_obj (tg, "Gene", 0)))
	  {
	    isCloud = ac_has_tag (Gene, "Cloud_gene") ;
	    tbl = ac_tag_table (Gene, "GeneId", 0) ;
	    for (ir = 0 ; ngid < 12 && tbl && ir < tbl->rows ; ir++)
	      {
		gids[ngid++] = ac_table_key (tbl, ir, 0, 0) ;	    
	      }
	    ac_free (tbl) ;
	    ac_free (Gene) ;
	  }
	ac_free (tg) ;
	if (isCloud && zz->noCloud)
	  continue ;


	tbl = ac_tag_table (mrna, "Product", 0) ;
	for (ir = jj = 0, cc = 'A' ; tbl && ir < tbl->rows && ir < 25 && jj < 2 ; cc++, ir++)
	  {
	    product = ac_table_obj (tbl, ir, 0, 0) ;
	    if (ac_has_tag (product, "Best_product"))
	      { 
		jj++ ;
		buf1a[0] = cc ;
		buf3a[1] = cc ;
		buf5a[1] = cc ;
	      }
	    else if (ac_has_tag (product, "Good_product"))
	      {
		jj++ ;
		buf1b[0] = cc ;
		buf3b[1] = cc ;
		buf5b[1] = cc ;
	      }
	    ac_free (product) ;	    
	  }
	ac_free (tbl) ;
	if (!jj)
	  continue ;
	tbl = ac_tag_table (mrna, "Coding", 0) ;
	for (ir = 0 ; tbl && ir < tbl->rows ; ir++)
	  {
	    ccp = ac_table_printable (tbl, ir, 4, "") ;
	    {
	      BOOL found = FALSE ;
	      
	      if (strstr (ccp, buf1a) && !strstr (ccp, buf5a)  && !strstr (ccp, buf3a))
		found = TRUE ;
	      if (jj > 1 && strstr (ccp, buf1b) && !strstr (ccp, buf5b)  && !strstr (ccp, buf3b))
		found = TRUE ;
	      if (!found)
		continue ;
	    }
	    up = arrayp (zz->aa, nn++, ZZZ) ;
	    up->gene[0] = gene ;
	    up->nGene = 1 ;
	    a1 = ac_table_int (tbl, ir, 0, 0) ;
	    a2 = ac_table_int (tbl, ir, 1, 0) ;
		if (g1 < g2)
	      {
		up->isDown = TRUE ;
		up->a1 = g1 + a1 - 1 ;
		up->a2 = g1 + a2 - 1 ;
	      }
	    else
	      {
		up->isDown = FALSE ;
		up->a2 = g1 - a1 + 1 ;
		up->a1 = g1 - a2 + 1 ;
	      }

	    for (i=0 ; i < ngid ; i++)
	      up->gId[i] = gids[i] ;
	    up->nGid = ngid ;
	    up->map = map ;
	    up->type = CODINGEXON ;
	  }
	ac_free (tbl) ;
      }
  ac_free (iter) ;
} /* zzGetCodingData */

/*************************************************************************************/

static int zzOrder (const void *a, const void *b)
{
  const ZZZ *up = (const ZZZ*)a, *vp = (const ZZZ *)b ;
  int n ;

  if (up->type && ! vp->type)
    return -1 ;
  if (vp->type && ! up->type)
    return 1 ;

  if (up->map < vp->map)
    return -1 ;
  if (up->map > vp->map)
    return 1 ;
  if (up->isDown && ! vp->isDown) 
    return -1 ;
  if (vp->isDown && ! up->isDown) 
    return 1 ;
  n = up->a1 - vp->a1 ;
  if (n)
    return n ;
  n = up->a2 - vp->a2 ;
  if (n)
    return n ;
  n = up->type - vp->type ;
  if (n)
    return n ;
  return up->gene - vp->gene ;
} /* zzOrder */

/*************************************************************************************/
/* the 2 strands are exported intermingled */
static int zzOrderMixed (const void *a, const void *b)
{
  const ZZZ *up = (const ZZZ*)a, *vp = (const ZZZ *)b ;
  int n ;

  if (up->map < vp->map)
    return -1 ;
  if (up->map > vp->map)
    return 1 ;
  n = up->a1 - vp->a1 ;
  if (n)
    return n ;
  n = up->a2 - vp->a2 ;
  if (n)
    return n ;
  n = up->type - vp->type ;
  if (n)
    return n ;
  return up->gene - vp->gene ;
} /* zzOrderMixed */

/*************************************************************************************/
/* merge the 2 strands, now called strand m */
static void zzMixStrands (ZZ *zz)
{
  int ii = zz->aa ? arrayMax (zz->aa) : 0 ;
  ZZZ *up ;

  if (ii)
    {
      up = arrp (zz->aa, 0, ZZZ) ;
      while (up++, ii--)
	up->isDown = TRUE ;
    }

  return ;
} /*  zzMixStrands */
      
 /*************************************************************************************/

static void zzMerge (ZZ *zz)
{
  int i, j, ii, jj, kk, a1, a2,  b2 ;
  BOOL oldstrand ;
  KEY gene, oldmap ;
  ZZZ *up, *vp, *wp ;
  KEYSET ks = keySetCreate () ;
  Array aa = zz->aa, bb ;

  arraySort (aa, zzOrder) ;
  
  /* merge contiguous coding segments and equivalent coding/exonic segments 
   *  for the moment there is a single gene per segment so nothing to merge
   *
   * first pass merge coding and exonic annotations 
   */
  for (ii = 0, up = arrp (aa, 0, ZZZ) ; ii < arrayMax (aa) ; up++, ii++)
    {
      if (!up->type)
	continue ;
      /* look forward for equivalent segment */
      for (kk = ii + 1, vp = up+1 ; kk < arrayMax (aa) ; kk++, vp++)
	{
	  if (!up->type)
	    continue ;
	  if (up->a1 == vp->a1 && up->a2 == vp->a2 &&
	      up->gene[0] == vp->gene[0])
	    {
	      /* merge the types */
	      if (vp->type > up->type)
		up->type = vp->type ;
	      vp->type = ZERO ;
	      continue ; 
	    }
	}
    }

  arraySort (aa, zzOrder) ;  
  /* clip the ZERO rejected at the bottom by zzOrder */
  for (ii = arrayMax (aa), up = arrp (aa, ii - 1, ZZZ) ; ii > 0 && up->type == ZERO ; ii--, up--) ;
  arrayMax (aa) = ii ;
    

  bb = arrayCopy (zz->aa) ;
  arrayMax (aa) = 0 ;
  /* split and duplicate the segments */
  for (ii = jj = 0, up = arrp (bb, 0, ZZZ) ; ii < arrayMax (bb) ; up++, ii++)
    {
      if (!up->type)
	continue ;
      /* copy as if complete */
      vp = arrayp (aa, jj++, ZZZ) ;	
      *vp = *up ;
      oldmap = up->map ; 
      b2 = up->a2 ;
      oldstrand = up->isDown ;
      /* now look if we must split */
      i = 0 ;
      /* look forward on possibly overlapping segments */
      for (kk = ii + 1, wp = up+1 ; kk < arrayMax (bb) ; kk++, wp++)
	{
	  if (wp->map != oldmap || wp->isDown != oldstrand || wp->a1 > b2)
	    break ;
	  if (!wp->touchedBy)
	    wp->touchedBy = ii + 1 ; /* +1 becuase first exon is ii==0 */
	  if (wp->a1 > up->a1 && wp->a1 <= up->a2)
	    keySet (ks, i++) = wp->a1 ;
	  if (wp->a2 >= up->a1 && wp->a2 < up->a2)
	    keySet (ks, i++) = wp->a2 + 1 ;
	}
      /* look backwards untill first segement that had touched me */
      for (kk = ii -1, wp = up-1 ; up->touchedBy && kk >= up->touchedBy - 1 ; kk--, wp--)
	{
	  if (wp->map != oldmap || wp->isDown != oldstrand)
	    break ;
	  if (wp->a1 > up->a1 && wp->a1 <= up->a2)
	    keySet (ks, i++) = wp->a1 ;
	  if (wp->a2 >= up->a1 && wp->a2 < up->a2)
	    keySet (ks, i++) = wp->a2 + 1 ;
	}
      if (!i)
	continue ;
      /* sort the collection of break points */
      keySetMax (ks) = i ;
      keySetSort (ks) ;
      keySetCompress (ks) ;
      /* duplicate the segment */
      a1 = up->a1 ; a2 = up->a2 ;
      for (i = 0 ; i < keySetMax (ks) ; i++)
	{ 
	  a2 = keySet (ks, i) ;
	  if (a2 > a1)
	    {
	      wp = arrayp (aa, jj++, ZZZ) ;	
	      *wp = *up ; /* do not use vp which may have been reallocated */
	      (wp-1)->a2 = a2 - 1 ;
	      a1 = wp->a1 = a2 ;
	    }
	}
      }
    
  /* sort again, then merge segments with identical coordinates */
  arraySort (aa, zzOrder) ;
  for (ii = 0, up = arrp (aa, 0, ZZZ) ; ii < arrayMax (aa) ; up++, ii++)
    {
      if (!up->type)
	continue ;
      /* look forward for equivalent segment */
      for (kk = ii + 1, vp = up+1 ; kk < arrayMax (aa) ; kk++, vp++)
	{
	  if (up->a1 != vp->a1 || up->a2 != vp->a2 ||
	      up->map != vp->map || up->isDown != vp->isDown )
	    break ;
	  if (vp->type < up->type) /* vp loses */
	    vp->type = ZERO ;
	  if (vp->type > up->type) /* vp wins */
	    {
	      up->type = vp->type ;
	      for (j = 0 ; j < vp->nGene ; j++)
		up->gene[j] = vp->gene[j] ;
	      for (j = 0 ; j < vp->nGid ; j++)
		up->gId[j] = vp->gId[j] ;
	      up->nGene = vp->nGene ;
	      up->nGid = vp->nGid ;
	      vp->type = ZERO ;
	    }
	  if (!vp->type) /* done with vp */
	    continue ;

	  /* equal types: merge the genes */ 
	  vp->type = ZERO ;
	  for (i = 0 ; i < vp->nGene ; i++)
	    {
	      gene = vp->gene[i] ;
	      if (up->nGene >= 12)
		gene = 0 ;
	      for (j = 0 ; gene && j < up->nGene ; j++)
		if (up->gene[j] == gene)
		  gene = 0 ;
	      if (gene) /* missing */
		{
		  up->gene[(up->nGene)++] = gene ;
		}
	    }
	  /* merge the gIds */
	  for (i = 0 ; i < vp->nGid ; i++)
	    {
	      gene = vp->gId[i] ;
	      if (up->nGid >= 12)
		gene = 0 ;
	      for (j = 0 ; gene && j < up->nGid ; j++)
		if (up->gId[j] == gene)
		  gene = 0 ;
	      if (gene) /* missing */
		{
		  up->gId[(up->nGid)++] = gene ;
		}
	    }
	}
    }
  arraySort (aa, zzOrder) ;  
  /* clip the ZERO rejected at the bottom by zzOrder */
  for (ii = arrayMax (aa), up = arrp (aa, ii - 1, ZZZ) ; ii > 0 && up->type == ZERO ; ii--, up--) ;
  arrayMax (aa) = ii ;

  /* we have to do this second pass after the completion of the annot merging */
  for (ii = 0, up = arrp (aa, 0, ZZZ) ; ii < arrayMax (aa) ; up++, ii++)
    {
      if (!up->type)
	continue ;
      /* look forward for contiguous segment */
      for (kk = ii + 1, vp = up+1 ; kk < arrayMax (aa) ; kk++, vp++)
	{
	  if (!up->type)
	    continue ;
	  if (up->map != vp->map)
	    break ;
	  if (up->a2+1 == vp->a1 && up->type == vp->type)
	    {
	      i = 12 ;
	      while (i--)
		if (up->gene[i] != vp->gene[i])
		  break ;
	      if (i == -1)
		{
		  /* merge the coords */
		  up->a2 = vp->a2 ;
		  vp->type = ZERO ;
		  continue ; 
		}
	    }
	}
    }

  keySetDestroy (ks) ;
  return ;
} /* zzMerge */

/*************************************************************************************/

static void zzAddIntergenic (ZZ *zz)
{
  int ii, jj, b2 ;
  ZZZ *up, *vp ;
  Array aa = zz->aa, bb = 0 ;
  BOOL oldstrand = FALSE ;
  KEY oldmap = 1 ; /* impossible value */

  bb = arrayCopy (aa) ;
  arraySort (bb, zzOrder) ; 

  for (ii = jj = 0, b2 = 0, up = arrp (bb, 0, ZZZ) ; ii < arrayMax (bb) ; up++, ii++)
    {
      if (!up->type)
	continue ;

      if (jj && b2 < up->a1 && oldstrand == up->isDown && oldmap == up->map)
	{
	  vp = arrayp (aa, jj++, ZZZ) ;	
	  vp->a1 = b2 ; vp->a2 = up->a1 - 1 ;
	  vp->isDown = oldstrand ;
	  vp->type = INTERGENIC ;
	  vp->nGene = vp->nGid = 0 ;
	  vp->map = up->map ;
	}
      vp = arrayp (aa, jj++, ZZZ) ;	
      *vp = *up ;
      oldstrand = up->isDown ;
      oldmap = up->map ;
      b2 = up->a2 + 1 ;
    }

  arrayMax (aa) = jj ;
  arrayDestroy (bb) ;
  return ;
} /* zzAddIntergenic */

/*************************************************************************************/

static void zzExport (ZZ *zz)
{
  int ii, i, j, nclo ;
  ZZZ *up ;
  char *cc ;
  const char *ccp ; 
  AC_KEYSET clones = 0 ;
  AC_OBJ Gene = 0 ;
  KEY gene ;
  KEYSET ks = keySetCreate () ; /* to export gene titles just once */

  arraySort (zz->aa, zzOrderMixed) ;

  freeOutf ("# Chromosome\tStrand\ta1\ta2\tType\tGeneIds\tGenes\n") ;
  for (ii = 0, up = arrp (zz->aa, 0, ZZZ) ; ii < arrayMax (zz->aa) ; up++, ii++)
    {
      if (!up->type)
	continue ;
      freeOutf ("%s\t%s\t%d\t%d\t%s\t"
		, ac_key_name(up->map) 
		, zz->mixStrands ? "0" : (up->isDown ? "+" : "-")
		, up->a1 
		, up->a2 
		, zType[up->type]
		) ;
      if (up->nGid)
	for (i = 0, cc="" ; i < up->nGid ; cc = "|", i++)
	  freeOutf ("%s%s", cc, ac_key_name(up->gId[i])) ;
      else
	if (up->type > 1) freeOut ("novel") ;
      freeOut ("\t") ;
      for (i = 0, cc="" ; i < up->nGene ; cc = "|", i++)
	freeOutf ("%s%s", cc, ac_key_name(up->gene[i])) ;
      freeOut ("\t") ;
      /* export the title of the new gene */
      nclo = -1 ;
      for (i = j = 0; !j && i < up->nGene ; i++)
	{
	  gene = up->gene[i] ;
	  if (! keySetFind (ks, gene, 0))
	    {
	      keySetInsert (ks, gene) ;
	      if (( Gene = ac_get_obj (zz->db, "gene", ac_key_name(gene), 0)))
		{
		  clones =  ac_objquery_keyset (Gene, ">transcribed_gene ; >read ; {IS *} SETOR {>buries;>buried_est;}; > cdna_clone", 0) ;
		  nclo = clones ? ac_keyset_count (clones) : 0 ;
		  ccp = ac_tag_printable (Gene, "Title", 0) ;
		  if (ccp)
		    {
		      j = 1 ;
		      freeOutf ("%d cDNA clone%s\t%s", nclo, nclo > 1 ? "s" : "", ccp) ;
		    }
		  ac_free (clones) ;
		  ac_free (Gene) ;
		}
	    }
	}
      if (nclo == -1)
	freeOut ("-\t-") ;
      freeOut ("\n") ;
    }
} /* zzExport */

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/ 

static int probesP2Torder (const void *a, const void *b)
{
  const P2T *up = (const P2T*)a, *vp = (const P2T*)b ;
  int n ;

  if (up->map < vp->map)
    return -1 ;
  if (up->map > vp->map)
    return 1 ;
  n = up->a1 - vp->a1 ;
  if (n) return n ;
  n = up->a2 - vp->a2 ;
  return n ;
} /* probesP2Torder */

/*************************************************************************************/ 

static void zzProbe2doubleTile (ZZ *zz)
{
  AC_HANDLE h = ac_new_handle() ;
  Array probes, tiles ;
  P2T *up, *vp ;
  int ii, jj, nn, a1, a2 ;
  const char *error_messages = 0 ;
  KEY probe, tile, map ;
  AC_KEYSET ks ;
  AC_TABLE tbl ;

  probes = arrayHandleCreate (400000, P2T, h) ;
  nn = 0 ;
  ks = ac_dbquery_keyset (zz->db, "Find probe intmap", h) ;
  tbl = ac_bql_table (zz->db, "select p, m, a1, a2 from p in @, m in p->intmap, a1 in m[1], a2 in m[2]", ks, 0, &error_messages, h) ;
  for (ii = 0 ; ii < tbl->rows ; ii++)
    {
      probe = ac_table_key (tbl, ii, 0, 0) ;
      for (jj = ii + 1 ;  jj < tbl->rows && probe == ac_table_key (tbl, jj, 0, 0) ; jj++) ;
      if (jj > ii + 1) { ii = jj - 1 ; continue ; } /* jump ambiguous probes */
      
      up = arrayp (probes, nn++, P2T) ;
      up->key = ac_table_key (tbl, ii, 0, 0) ;
      up->map = ac_table_key (tbl, ii, 1, 0) ;
      up->a1 = ac_table_int (tbl, ii, 2, 0) ;
      up->a2 = ac_table_int (tbl, ii, 3, 0) ;  
      up->type = 'p' ;
      if (up->a1 > up->a2) { int a0 = up->a1 ;  up->a1 =  up->a2 ;  up->a2 = a0 ; }
    }
  arraySort (probes, probesP2Torder) ;

  tiles = arrayHandleCreate (4000, P2T, h) ;
  nn = 0 ;
  ks = ac_dbquery_keyset (zz->db, "find sequence junction ", h) ;
  tbl = ac_bql_table (zz->db, "select t, m, a1, a2 from t in @, m in t->intmap, a1 in m[1], a2 in m[2]", ks, 0, &error_messages, h) ;
  for (ii = 0 ; ii < tbl->rows ; ii++)
    {
      tile = ac_table_key (tbl, ii, 0, 0) ;
      for (jj = ii + 1 ;  jj < tbl->rows && tile == ac_table_key (tbl, jj, 0, 0) ; jj++) ;
      if (jj > ii + 1) { ii = jj - 1 ; continue ; } /* jump ambiguous probes */
      
      up = arrayp (tiles, nn++, P2T) ;
      up->type = 't' ;
      up->key = ac_table_key (tbl, ii, 0, 0) ;
      up->map = ac_table_key (tbl, ii, 1, 0) ;
      up->a1 = ac_table_int (tbl, ii, 2, 0) ;
      up->a2 = ac_table_int (tbl, ii, 3, 0) ;
      if (up->a1 > up->a2) { int a0 = up->a1 ;  up->a1 =  up->a2 ;  up->a2 = a0 ; }
    }
  arraySort (tiles, probesP2Torder) ;

  for (ii = 0 ; ii < arrayMax (tiles) ; ii++)
    {
      nn = 0 ;
      up = arrayp (tiles, ii, P2T) ;
      a1 = up->a1 ; a2 = up->a2 ; map = up->map ;
       for (jj = 0 ; jj < arrayMax (probes) ; jj++)
	 {
	   vp = arrayp (probes, jj, P2T) ;
	   if (vp->map < map) continue ;
	   if (vp->map > map) break ;
	   if (vp->a1 < a1) continue ;
	   if (vp->a2 > a2) break ;
	   if (! nn++)
	     {
	       printf ("\nSequence %s\n"
		       , freeprotect (ac_key_name (up->key))
		       ) ;
	     }
	   printf ("Probe_hit %s\n"
		       , freeprotect (ac_key_name (vp->key))
		   ) ;
	 }
    }
  printf ("\n") ;
  ac_free (h) ;
} /* zzProbe2doubleTile */

/*************************************************************************************/
/***************************** Public interface **************************************/
/*************************************************************************************/

static void usage (void)
{
  fprintf (stderr, "// Usage: genetiling [-db ACEDB] [-noCloud] [-gene template] [-mixStrands] \n") ;
  fprintf (stderr, "// Example:  genetiling -db ZZ -chrom 4* \n") ;  
  fprintf (stderr, "//   -db ACEDB database\n") ;
  fprintf (stderr, "//   -gene template : just the genes with matching name (use ? and * as wild chars\n") ;
  fprintf (stderr, "//   -noCloud : ignore all the cloud genes\n") ;
  fprintf (stderr, "//   -pg : split according to the predicted_genes\n") ;
  fprintf (stderr, "//   -codingIgnored : do not split exons into coding/non-coding\n") ;
  fprintf (stderr, "//   -mixStrands : fuse the 2 strands\n") ;
  fprintf (stderr, "//   -remap : attribute the probes to their tiles based on their intmap\n") ;
  fprintf (stderr, "//   -probeType : classify the probes as exonic intergenic etc\n") ;
  fprintf (stderr, "//   -probe2doubleTile: associate the probes to a genomic double tile\n") ;
  exit (1) ;
} /* usage */

/*************************************************************************************/

int main (int argc, const char **argv)
{
  FILE *f = 0 ;
  int outlevel = 0 ;
  ZZ zz ;
  AC_HANDLE h = 0 ;
  const char *ici ;
  const char *s = "ok" ;

  freeinit () ; 
  h = ac_new_handle () ;
  memset (&zz, 0, sizeof (ZZ)) ;
  zz.h = h ;

  /* optional temple argument */
  zz.gene = "*" ;
  getCmdLineOption (&argc, argv, "-gene", &(zz.gene)) ;
  zz.noCloud = getCmdLineOption (&argc, argv, "-noCloud", 0) ;
  zz.mixStrands = getCmdLineOption (&argc, argv, "-mixStrands", 0) ;
  zz.codingIgnored = getCmdLineOption (&argc, argv, "-codingIgnored", 0) ;
  zz.remap = getCmdLineOption (&argc, argv, "-remap", 0) ;
  zz.pg = getCmdLineOption (&argc, argv, "-pg", 0) ;
  zz.probeType = getCmdLineOption (&argc, argv, "-probeType", 0) ;

  /* july 2007, analysis for boubou */
  zz.probe2doubleTile = getCmdLineOption (&argc, argv, "-probe2doubleTile", 0) ;
  outlevel = freeOutSetFile (stdout) ;	

 /* mandatory database descriptor */
  if (getCmdLineOption (&argc, argv, "-db", &ici))
    {
      zz.db = ac_open_db (messprintf("%s",ici), &s);
      if (!zz.db)
	messcrash ("Failed to open db %s, error %s", ici, s) ;
    }
  else 
     usage () ;

  zz.aa = arrayHandleCreate (20000, ZZZ, zz.h) ;

  fprintf (stderr, "// start: %s\n", timeShowNow()) ;

  if (zz.remap)
    zzRemap (&zz) ;
  else if (zz.probeType)
    zzProbeType (&zz) ;
  else if (zz.probe2doubleTile)
    zzProbe2doubleTile (&zz) ;
  else
    {
      if (zz.pg) zzGetPgData (&zz) ;
      else zzGetTgData (&zz) ;
      if (! zz.codingIgnored) zzGetCodingData (&zz) ;
      if (zz.mixStrands) zzMixStrands (&zz) ;
      
      fprintf (stderr, "// data parsed, start merging: %s\n", timeShowNow()) ;
      zzMerge (&zz) ;
      zzAddIntergenic (&zz) ;
      fprintf (stderr, "// merging done: %s\n", timeShowNow()) ;
      zzExport (&zz) ;
    }

  fprintf (stderr, "// done: %s\n", timeShowNow()) ;

  if (outlevel)
    freeOutClose (outlevel) ;
  if (f) filclose (f) ;
  ac_db_close (zz.db) ;
  ac_free (zz.h) ;
  if (0) sleep (1) ; /* to prevent a mix up between stdout and stderr */
  return 0 ;
}

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/

