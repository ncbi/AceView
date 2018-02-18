/*
 Exportation des statistiques acembly pour la conference japon 2003
*/
/* 
   beautiful gawk code
   gawk '/^Variant/{aa[$2]++;}
        END{for (v in aa) printf("%s\t%d\n",v,aa[v]);}'
      variant.ace | wc


anatomy of a gene
tbly . << EOF
  query find gene structure
  show -a -f gene_structure.ace structure
  quit
EOF

#count the mrna
cat ZH* /d14.structure.*.ace | 
gawk '/^Gene /{ngene++;next;}/^$/{next;}{aa[$1]+=$2;}END{printf("using %d genes\n", ngene); for (v in aa) printf("%s\t%d\n",v,aa[v]);}'

#count by cumul
cat ZH* /d14.structure.*.ace | 
gawk '/^Gene /{ngene++;next;}/^$/{next;}{aa[$1]+=1;}END{printf("using %d genes\n", ngene); for (v in aa) printf("%s\t%d\n",v,aa[v]);}' gene_structure.ace
#count by cumul truque pour les promoteurs
gawk '/^Gene /{ngene++;next;}/^$/{next;}{if($2>1)aa[$1]+=1;}END{printf("using %d genes\n", ngene); for (v in aa) printf("%s\t%d\n",v,aa[v]);}' 
*/

#include "acedb.h"
#include <errno.h>
#include "cdna.h"
#include "dict.h"
#include "freeout.h"
#include "utils.h"

static BOOL isWorm = FALSE ;

/*********************************************************************/

static int orderByA1 (const void *va, const void *vb)
{
  const HIT *up = (const HIT *)va, *vp = (const HIT *)vb ;

  if (up->cDNA_clone != vp->cDNA_clone)
    return up->cDNA_clone - vp->cDNA_clone ;  /* the cDNA_clone */
  else if (up->est != vp->est)
    return up->est - vp->est ;  /* the cDNA_clone */
  else if (up->a1 != vp->a1)
    return up->a1 - vp->a1 ;
  else if (up->a2 != vp->a2)
    return up->a2 - vp->a2 ;      
  else  
    return up->x1 - vp->x1 ;    
}

/***************************************************************************/

static void sshowHits(Array hits)
{
  int j ;
  HIT *vp ;

  if (hits && arrayExists(hits)) 
    for (j = 0  ; j < arrayMax(hits) ; j++)
      {
	vp = arrp(hits, j, HIT) ;
	printf("%2d: a1=%d a2 = %d x1 = %d x2 = %d", j, vp->a1, vp->a2, vp->x1, vp->x2) ;
	printf ("\n") ;
      }
} /* sshowHits */

 /* analyse cloud of small genes to see if they are inside real genes */
typedef struct { int ii, map, a1, a2 ; BOOL  hasIntron, isCoding, isDown, inSameStrand, inOtherStrand ; } CLOUD ;

/* we ant to plot for genes with a givem number of clones the rate ofgenes
 *with at least 1 2 3 4 5 10 20 forms
*/
static void sCountAlternatifs (AC_DB db)
{
  int i, j, k, nClones, nMrna, nTg = 0 ;
  int kk [] = { 1,2,3,4,5, 10, 20, 50, 100, 9999999 } ;
  int MAXCLO = 100 ;
  AC_KEYSET tgs, ks ;
  AC_ITER iter ;
  AC_OBJ tg ;
  AC_HANDLE h = handleCreate () ;
  Array aa = arrayHandleCreate (10000, int, h) ;
  Array bb = arrayHandleCreate (100, int, h) ;

  printf ("############################################################\n") ;
  printf ("## phase 40: sCountAlternatifs   Are we staturated in alternative forms[2]\n") ;

  tgs = ac_dbquery_keyset (db, "find tg  mrna &&  ! shedded_from && (gt_ag || gc_ag)", h) ;

  iter = ac_keyset_iter (tgs, TRUE, h) ;
  while ((tg = ac_next_obj (iter)))
    {
      ks = ac_objquery_keyset (tg, ">cdna_clone", 0) ;
      nClones = ac_keyset_count (ks) ;
      ac_free (ks) ;

      ks = ac_objquery_keyset (tg, ">mrna", 0) ;
      nMrna = ac_keyset_count (ks) ;
      ac_free (ks) ;

      /* j is log2 of nClones */
      for (i = 1, j = 0 ; j < MAXCLO && i < nClones ; j++, i <<= 1 ) ;
      (array (aa, MAXCLO * nMrna + j, int)) ++ ;
      (array (bb, j, int) )++ ;

      nTg++ ;
      ac_free (tg) ;
    }
  printf ("\n\nFor genes with at #n clones, What fraction have 1 2 3  n alternatives\n") ;
  printf ("%12s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\n",
	  "#clo", "#genes", "1", "2", "3", "4", "5", "6-10","11-20", "21-50") ;
  printf ("%12s\t%8d\n", "Total", nTg) ;
  for (i = 1, j = 0 ; j < arrayMax(bb) ; i <<= 1, j++)
    {
      int n, i1, j1, ng = 0 ;

      for (j1 = j; j1 <= j ; j1++)
	ng += arr (bb, j1, int) ;
      
      printf ("\t%s", messprintf("%d-%d", (i/2+1), i)) ;
      printf ("\t%d", ng) ; /* nb of gene in this bin */
      
      if (!ng) ng = 1 ;
      
      for (k = 0 ; k < 8 ; k++)
	{
	  n = 0 ;
	  for (j1 = j; j1 <= j ; j1++)
	    for (i1 = j1 + MAXCLO * kk[k] ; i1 < arrayMax(aa) && i1 <  j1 + MAXCLO * kk[k+1] ; i1 += MAXCLO)
	      n += arr (aa, i1, int) ;
	  printf ("\t%8.2f", 100.0 * n / ng) ;
	}
      printf ("\n") ;
    } 
  printf ("// done \n\n") ;
} /* sCountAlternatifs */

/*************************************************************************/
/*************************************************************************/

static int sAliDoStats (char *title, AC_ITER iter) /* added feb 25 */
{
  /*
    lg tot
    lg tot alignee
    nb ali
    nb totala moyen d erreurs
  */
  AC_HANDLE h = handleCreate () ;
  int ir, jr ;
  long int nnest = 0, nnestali, nbpali, nln = 0, nerr, bestAli, ali, err ;
  AC_OBJ est = 0 ;
  AC_TABLE mfg ;

  ir = jr = nnest = nnestali = nbpali = nln = nerr = bestAli = ali = err = 0 ;
  while (ac_free (est), (est = ac_next_obj (iter)))
    {
      nnest++ ;
      if (!ac_has_tag (est, "from_gene"))
	continue ;
      nnestali++ ;
      mfg = ac_tag_table (est, "from_gene", h) ;
      for (ir = 0, jr = -1, bestAli = 0 ; mfg && ir < mfg->rows ; ir++)
	{
	  ali = ac_table_int (mfg, ir, 3, 0) ;
	  if (bestAli < ali) { jr = ir ; bestAli = ali ; }
	}
      if (jr >= 0)
	{
	  nln += ac_table_int (mfg, jr, 1, 0) ;
	  nbpali += ac_table_int (mfg, jr, 3, 0) ;
	  nerr += ac_table_int (mfg, jr, 4, 0) ;
	} 
      ac_free (mfg) ;
    }

  freeOutf ("\n%s\t%8ld\t%8ld", title, nnest, nnestali) ;
  if (nnest)
    {
      freeOutf ("\t(%6.2f%%)\n", (100.0 * nnestali)/nnest) ;
      freeOutf ("\tlength to align =\t%12ld bp,\t%d bp/sequence\n"
		, nln, nln/nnestali) ;
      if (nln)
	freeOutf ("\tlength aligned =\t%12ld bp\t(%6.2f%%)\t%ld bp/sequence\n"
		  ,nbpali, (100.0 * nbpali)/nln, nbpali/nnestali) ;
      if (nbpali)
	freeOutf ("\tcumulated errors =\t%12ld bp\t(%6.2f%%)\t%ld err/sequence\n"
		  , nerr , (100.0 * nerr)/nbpali, nerr/nnestali) ;
    }
  ac_free (h) ;
  return 1 ;
}

static int sAliStats (AC_DB db)  /* added feb 25 */
{
  AC_HANDLE h = handleCreate () ;
  AC_ITER NMs, mrnas, ests ;
  
  printf ("############################################################\n") ;
  printf ("## phase 21: sAliStats\n") ;
  printf ("## Values are computed on the affectivelly aligned reads\n") ;
  if (! isWorm)
    {
      NMs = ac_dbquery_iter (db, "Find EST Ref_Seq", h) ;
      sAliDoStats ("NMs", NMs) ;
      mrnas = ac_dbquery_iter (db, "Find EST (IS_mrna || ref_mrna) && ! Ref_Seq", h) ;
      sAliDoStats ("mRNAs", mrnas) ;
      ests = ac_dbquery_iter (db, "Find EST !IS_mrna && !ref_mrna && ! Ref_Seq", h) ;
    }
  else
    ests = ac_dbquery_iter (db, "Find EST DNA && IS_read", h) ;
  sAliDoStats ("ESTs", ests) ;
   
  ac_free (h) ;
  return 1 ;
}

static int sAlterStats (AC_DB db) /* added feb 25 */
{
  AC_HANDLE h = handleCreate () ;
  AC_KEYSET tgs, ks1, ks12, ks13, ks2, ks22, ks23 ;
  int ng, nn1, nn12, nn13, nn2, nn22, nn23 ;

  printf ("############################################################\n") ;
  printf ("## phase 52: sAlterStats Are we staturated in alternative forms\n") ;


  tgs = ac_dbquery_keyset (db, "Find tg  mrna &&  ! shedded_from && (gt_ag || gc_ag)", h) ;
  ng = ac_keyset_count (tgs) ;
  ks1 = ac_ksquery_keyset (tgs, "COUNT cdna_clone > 1", h) ;
  nn1 = ac_keyset_count (ks1) ;
  ks12 = ac_ksquery_keyset (ks1, "COUNT mrna >= 2", h) ;
  nn12 = ac_keyset_count (ks12) ;
  ks13 = ac_ksquery_keyset (ks1, "COUNT {>mrna; gt_ag || gc_ag} >= 2", h) ;
  nn13 = ac_keyset_count (ks13) ;
  
  ks2 = ac_ksquery_keyset (tgs, "COUNT cdna_clone > 2", h) ;
  nn2 = ac_keyset_count (ks2) ;
  ks22 = ac_ksquery_keyset (ks2, "COUNT mrna >= 2", h) ;
  nn22 = ac_keyset_count (ks22) ;
  ks23 = ac_ksquery_keyset (ks2, "COUNT {>mrna; gt_ag || gc_ag} >= 2", h) ;
  nn23 = ac_keyset_count (ks23) ;
  
  freeOutf ("%d non-shed genes with gt-ag||gc_ag introns are defined in the database\n", ng) ;
  if (ng)
    {
      freeOutf ("\t%d\t%6.2f%% are supported by > 1 clones\n", nn1, (100.0 * nn1)/ng) ;
      if (nn1)
	freeOutf ("\t\t%d\t%6.2f%% of them have >= 2 mRNAs\n", nn12, (100.0 * nn12)/nn1) ;
      if (nn1)
	freeOutf ("\t\t%d\t%6.2f%% of them have >= 2 mRNAs with introns\n", nn13, (100.0 * nn13)/nn1) ;
      freeOutf ("\t%d\t%6.2f%% are supported by > 2 clones\n", nn2, (100.0 * nn2)/ng) ;
      if (nn2)
	freeOutf ("\t\t%d\t%6.2f%% of them have >= 2 mRNAs\n", nn22, (100.0 * nn22)/nn2) ;
      if (nn2)
	freeOutf ("\t\t%d\t%6.2f%% of them have >= 2 mRNAs with introns\n", nn23, (100.0 * nn23)/nn2) ;
    }

  ac_free (h) ;
  return 1 ;
}

static int sLocusIdStats  (AC_DB db) /* added 2005_08_30 */
{
  AC_HANDLE h = handleCreate () ;
  AC_KEYSET aks, aks1 ;
  int n, nn , nn1, ii, total ;
  /*
    AC_ITER iter ;
    AC_OBJ gene = 0, pg = 0 ;
    KEYSET ks = 0 ;
    KEYSET gene2tg2pg = keySetHandleCreate (h) ;
    KEYSET pg2tg2gene = keySetHandleCreate (h) ; 
  */

  printf ("############################################################\n") ;
  printf ("## phase 53:GeneID+Genefinder <--> AceView, Is AceView splitting or merging genes ?\n") ;
  
  freeOutf ("Some \"predicted_gene\"locusid have no NM or XM, they have an NC_ (44 cas/72 on human-ZH21)\n") ;
  freeOutf ("2005_08_30 we looked on CH22 at the case one pg touches 3 aceview: \n") ;
  freeOutf (" half the cases, we have partials (true splitting)\n") ;
  freeOutf (" half the cases, genefinder eats cloud genes which are inside introns\n") ;

  freeOutf ("This code does not work becuase of a bug in the query language\n") ; 
  goto done ;
  if (0)
    {
      freeOutf ("Type\tTotal\t0\t1\t2\t3\t4\t5\t6\t7\t8\t9\t>=10\n") ;
      /* count the main aceview genes with so many geneid */
      aks = ac_dbquery_keyset (db, "Find gene main && (!geneid || !genefinder)", h) ; 
      total = nn = ac_keyset_count (aks) ;
      ac_free (aks) ;
      aks = ac_dbquery_keyset (db, "Find gene main && geneid && genefinder", h) ; 
      nn1 = ac_keyset_count (aks) ;
      total += nn1 ;
      freeOutf ("AceView with # GeneID+Model\t%d\t%d", total, nn) ;
      for (ii = 1 ; ii < 10 ; ii++)
	{
	  nn = nn1 ;
	  aks1 = ac_ksquery_keyset (aks, messprintf ("COUNT ({>genefinder ; >geneid_pg } SETAND {>geneid}) > %d", ii), h) ;
	  nn1 = ac_keyset_count (aks1) ;
	  freeOutf ("\t%d", nn - nn1) ;
	  ac_free (aks) ;
	  aks = aks1 ;
	}
      freeOutf ("\t%d\n", nn1) ;
    }
  /* count the geneid with so many main aceview genes */
  aks = ac_dbquery_keyset (db, "Find geneid (!Predicted_gene || !gene)", h) ; 
  total = nn = ac_keyset_count (aks) ;
  ac_free (aks) ;
  aks = ac_dbquery_keyset (db, "Find geneid predicted_gene && gene", h) ; 
  nn1 = ac_keyset_count (aks) ;
  total += nn1 ;
  freeOutf ("AceView with # GeneID+Model\t%d\t%d", total, nn) ;
  for (ii = 1 ; ii < 10 ; ii++)
    {
      nn = nn1 ;
      aks1 = ac_ksquery_keyset (aks, messprintf ("COUNT ({>predicted_gene ; >geneid_pg } SETAND {>geneid}) > %d", ii), h) ;
      nn1 = ac_keyset_count (aks1) ;
      freeOutf ("\t%d", nn - nn1) ;
      ac_free (aks) ;
      aks = aks1 ;
    }
  freeOutf ("\t%d\n", nn1) ;
 
  n = ac_keyset_count (aks = ac_dbquery_keyset (db, "Find predicted_gene !gene && !nm_id && geneid_pg", h)) ;
  nn1 = ac_keyset_count (ac_ksquery_keyset (aks, "Find predicted_gene !nm_id && geneid_pg && pseudogene", h)) ;
  freeOutf ("In addition we have %d predicted_gene with geneid but no NM_ID (%d are pseudogenes)\n", n, nn1) ;
 done:
  ac_free (h) ;
  return 1 ;
} /* sLocusIdStats */

static int sGeneIdStats  (AC_DB db) /* added 2005_08_30 */
{
  AC_HANDLE h = handleCreate () ;
  AC_KEYSET aks1, aks2, aks3 ;
  int n1, n2, n3 ;
  int ngs, ngsid, ngu, nguid, nms, nmsid, nmu, nmuid, ncl ;
  ngs = ngsid = ngu = nguid = nms = nmsid = nmu = nmuid = ncl = 0 ;


  printf ("############################################################\n") ;
  printf ("## phase 54: sGeneIdStats  spliced genes with GeneId Refseq and support\tGene\tGeneID\n") ;

  aks1 = ac_dbquery_keyset (db, "Find gene nStandardIntrons && geneid", h) ;
  n1 = ac_keyset_count (aks1) ;
  aks2 = ac_ksquery_keyset (aks1, "COUNT {>transcribed_gene ; >read ; ! ref_seq } > 0", h) ;
  n2 = ac_keyset_count (aks2) ;
  aks3 = ac_ksquery_keyset (aks2, ">geneid", h) ;
  n3 = ac_keyset_count (aks3) ;
  freeOutf ("AceView spliced genes with GeneID and cDNA support\t%d\t%d\n", n2, n3) ;
  ngs += n2 ; ngsid += n3 ;

  ac_keyset_minus (aks1, aks2) ;
  n2 = ac_keyset_count (aks1) ;
  aks3 = ac_ksquery_keyset (aks1, ">geneid", h) ;
  n3 = ac_keyset_count (aks3) ;


  aks1 = ac_dbquery_keyset (db, "Find gene nStandardIntrons && !geneid", h) ;
  n1 = ac_keyset_count (aks1) ;
  freeOutf ("AceView spliced gene with cDNA support but no GeneID\t%d\t0\n", n1) ;
  ngs += n1 ; 
  freeOutf ("AceView spliced gene with GeneID but no cDNA support\t%d\t%d\n", n2, n3) ;
  nms += n2 ; nmsid += n3 ;

  aks1 = ac_dbquery_keyset (db, "Find gene !nStandardIntrons && geneid", h) ;
  n1 = ac_keyset_count (aks1) ;
  aks2 = ac_ksquery_keyset (aks1, "COUNT {>transcribed_gene ;  >read ; ! ref_seq} > 0", h) ;
  n2 = ac_keyset_count (aks2) ;
  aks3 = ac_ksquery_keyset (aks2, ">geneid", h) ;
  n3 = ac_keyset_count (aks3) ;
  freeOutf ("AceView unspliced gene with GeneID and cDNA support\t%d\t%d\n", n2, n3) ;
  ngu += n2 ; nguid += n3 ;

  aks1 = ac_dbquery_keyset (db, "Find gene !nStandardIntrons && main && !geneid", h) ;
  n1 = ac_keyset_count (aks1) ;
  aks2 = ac_ksquery_keyset (aks1, "COUNT {>transcribed_gene ;  >read ; ! ref_seq} > 0", h) ;
  n2 = ac_keyset_count (aks2) ;
  freeOutf ("AceView unspliced main gene with cDNA support but no GeneID\t%d\t0\n", n2) ;
  ngu += n2 ; nguid += 0 ;

  aks1 = ac_dbquery_keyset (db, "Find gene !nStandardIntrons && putative && !geneid", h) ;
  n1 = ac_keyset_count (aks1) ;
  aks2 = ac_ksquery_keyset (aks1, "COUNT {>transcribed_gene ;  >read ; ! ref_seq} > 0", h) ;
  n2 = ac_keyset_count (aks2) ;
  freeOutf ("AceView unspliced putative gene with cDNA support but no GeneID\t%d\t0\n", n2) ;
  ncl += n2 ;

  aks1 = ac_dbquery_keyset (db, "Find gene cloud", h) ;
  n1 = ac_keyset_count (aks1) ;
  freeOutf ("AceView cloud gene\t%d\t0\n", n1) ;
  ncl += n1 ;

  aks1 = ac_dbquery_keyset (db, "Find gene !nStandardIntrons && geneid", h) ;
  n1 = ac_keyset_count (aks1) ;
  aks2 = ac_ksquery_keyset (aks1, "COUNT {>transcribed_gene ;  >read ; ! ref_seq} > 0", h) ;
  ac_keyset_minus (aks1, aks2) ;
  n2 = ac_keyset_count (aks1) ;
  aks3 = ac_ksquery_keyset (aks1, ">geneid", h) ;
  n3 = ac_keyset_count (aks3) ;
  freeOutf ("AceView unspliced gene with GeneID but no cDNA support\t%d\t%d\n", n2, n3) ;
  nmu += n2 ; nmuid += n3 ;  

  aks1 = ac_dbquery_keyset (db, "Find predicted_gene !model_of_gene", h) ;
  n1 = ac_keyset_count (aks1) ;
  aks2 = ac_ksquery_keyset (aks1, "COUNT Source_Exons > 1", h) ;
  n2 = ac_keyset_count (aks2) ;
  aks3 = ac_ksquery_keyset (aks1, ">geneid_pg", h) ;
  n3 = ac_keyset_count (aks2) ;
  freeOutf ("Predicted spliced model with GeneID but no support in AceView\t%d\t%d\n", n2, n3) ;
  nms += n2 ; nmsid += n3 ;  

  ac_keyset_minus (aks1, aks2) ;
  n2 = ac_keyset_count (aks1) ;
  aks3 = ac_ksquery_keyset (aks1, ">geneid_pg", h) ;
  n3 = ac_keyset_count (aks3) ;
  freeOutf ("Predicted unspliced model with GeneID but no support in AceView\t%d\t%d\n", n2, n3) ;
  nmu += n2 ; nmuid += n3 ;  

  freeOutf ("\nGenes with cDNA support\t%d\t%d\n", ngs + ngu, ngsid + nguid) ;
  freeOutf ("Spliced\t%d\t%d\n", ngs, ngsid) ;
  freeOutf ("Not spliced\t%d\t%d\n", ngu, nguid) ;

  freeOutf ("\nCloud or putative AceView genes\t%d\n", ncl) ;

  freeOutf ("\nModels with no support\t%d\t%d\n", nms + nmu, nmsid + nmuid) ;
  freeOutf ("Spliced\t%d\t%d\n", nms, nmsid) ;
  freeOutf ("Not spliced\t%d\t%d\n", nmu, nmuid) ;

  ac_free (h) ;
  return 1 ;
} /* sLocusIdStats */

static void sCountCloudGenes (AC_DB db)  /* added feb 25 */
{
  AC_HANDLE h = handleCreate () ;

  AC_KEYSET mrnas, ks, ksg, ksp , genes, cloud ;
  AC_KEYSET ks1, ks11, ks110, ks111, ks12, ks120, ks121, ks2, ks21, ks210, ks211, ks22, ks220, ks221 ;  
  int nn1, nn11, nn110, nn111, nn12, nn120, nn121, nn2, nn21, nn210, nn211, nn22, nn220, nn221 ;  
  int nn3, ng, np, nMrna ;

  printf ("############################################################\n") ;
  printf ("## phase 50: sCountCloudGenes\n") ;

  /*
    date du download des data
    nb NM       total ali  from_gene > 1
    .  mrna      idem
    .  est       idem
    
    nb gene cloud !cloud
       tg cloud      classic introns    100<=prot<300   prot >=  300
                     non classic        100<=prot<300   prot >=  300
          !cloud     classic introns    100<=prot<300   prot >=  300
                     non classic        100<=prot<300   prot >=  300

    tg with classic intron
      -> # mrna  with classic introns ->    gap | !gap
                 no classic introns ->    gap | !gap

    among the 168000 mrna with introns how many have:
         one product     -> same num of products
         more than one   -> total nb of products
         cds || mrna covered_by  1, 2, > 2
  */
  genes = ac_dbquery_keyset (db , "find gene ! shedded_from && transcribed_gene", h) ;
  ng = ac_keyset_count (genes) ;
  cloud = ac_ksquery_keyset (genes, "cloud_gene", h) ; 
  ks1 = ac_ksquery_keyset (cloud, ">transcribed_gene", h) ;
  nn1 = ac_keyset_count (ks1) ;
  ks11 = ac_ksquery_keyset (ks1, ">transcribed_gene ; gt_ag || gc_ag", h) ;
  nn11 = ac_keyset_count (ks11) ;
  ks110 = ac_ksquery_keyset (ks11, " COUNT { >mrna ; !gap_length && longest_CDS >= 300} > 0", h) ;  
  ks111 = ac_ksquery_keyset (ks110, " COUNT { >mrna ; !gap_length && longest_CDS >= 900} > 0", h) ; 
  nn111 = ac_keyset_count (ks111) ;
  nn110 = ac_keyset_count (ks110) - nn111 ;
  ks12 = ac_copy_keyset (ks1, h) ;
  ac_keyset_minus (ks12, ks11) ;
  nn12 = ac_keyset_count (ks12) ;
  ks120 = ac_ksquery_keyset (ks12, " COUNT { >mrna ; !gap_length && longest_CDS >= 300} > 0", h) ;  
  ks121 = ac_ksquery_keyset (ks120, " COUNT { >mrna ; !gap_length && longest_CDS >= 900} > 0", h) ; 
  nn121 = ac_keyset_count (ks121) ;
  nn120 = ac_keyset_count (ks120) - nn121 ;

  ks2 = ac_copy_keyset (genes, h) ; 
  ac_keyset_minus (ks2, cloud) ;
  nn2 = ac_keyset_count (ks2) ;
  ks2 = ac_ksquery_keyset (ks2, ">transcribed_gene", h) ;
  ks21 = ac_ksquery_keyset (ks2, "gt_ag || gc_ag", h) ;
  ks210 = ac_ksquery_keyset (ks21, " COUNT { >mrna ; !gap_length && longest_CDS >= 300} > 0", h) ;  
  ks211 = ac_ksquery_keyset (ks210, " COUNT { >mrna ; !gap_length && longest_CDS >= 900} > 0", h) ; 
  nn211 = ac_keyset_count (ks211) ;
  nn210 = ac_keyset_count (ks210) - nn211 ;
  nn21 = ac_keyset_count (ks21) ;
  ks22 = ac_copy_keyset (ks2, h) ;
  ac_keyset_minus (ks22, ks21) ;
  nn22 = ac_keyset_count (ks22) ;
  ks220 = ac_ksquery_keyset (ks22, " COUNT { >mrna ; !gap_length && longest_ORF >= 300} > 0", h) ;  
  ks221 = ac_ksquery_keyset (ks220, " COUNT { >mrna ; !gap_length && longest_ORF >= 900} > 0", h) ; 
  nn221 = ac_keyset_count (ks221) ;
  nn220 = ac_keyset_count (ks220) - nn221 ;

  /* COUNT { >mrna ; !gap_length && longest_ORF >= 900} > 0 */

  freeOutf ("%d genes are defined in the database\n", ng) ;
  if (ng)
    {
      freeOutf ("\t%d\t%6.2f%% are cloud genes\n", nn1, (100.0 * nn1)/ng) ;  
      freeOutf ("\t\t%d\t%6.2f%% of the cloud genes have gt_ag or gc_ag introns\n", nn11, (100.0 * nn11)/nn1) ;
      freeOutf ("\t\t\t%d\t%6.2f%% have a CDS between 100 and 299 AA\n", nn110,  (100.0 * nn110)/nn11) ;
      freeOutf ("\t\t\t%d%6.2f%% have a CDS >= 300 AA\n", nn111,  (100.0 * nn111)/nn11) ;
      freeOutf ("\t\t%d\t%6.2f%% of the cloud genes do not\n", nn12, (100.0 * nn12)/nn1) ;
      freeOutf ("\t\t\t%d\t%6.2f%% have an ORF between 100 and 299 AA\n", nn120,  (100.0 * nn120)/nn12) ;
      freeOutf ("\t\t\t%d\t%6.2f%% have an ORF >= 300 AA\n", nn121,  (100.0 * nn121)/nn12) ;
      freeOutf ("\t%d\t%6.2f%% are not cloud genes\n", nn2, (100.0 * nn2)/ng) ;  
      freeOutf ("\t\t%d\t%6.2f%% of them have gt_ag or gc_ag introns\n", nn21, (100.0 * nn21)/nn2) ;
      freeOutf ("\t\t\t%d\t%6.2f%% have a CDS between 100 and 299 AA\n", nn210,  (100.0 * nn210)/nn21) ;
      freeOutf ("\t\t\t%d\t%6.2f%% have a CDS >= 300 AA\n", nn211,  (100.0 * nn211)/nn21) ;
      freeOutf ("\t\t%d\t%6.2f%% of them do not\n", nn22, (100.0 * nn22)/nn2) ;
      freeOutf ("\t\t\t%d\t%6.2f%% have an ORF between 100 and 299 AA\n", nn220,  (100.0 * nn220)/nn22) ;
      freeOutf ("\t\t\t%d\t%6.2f%% have an ORF >= 300 AA\n", nn221,  (100.0 * nn221)/nn22) ;
    }

  mrnas = ac_dbquery_keyset (db , "find mrna ; dna", h) ;
  nMrna = ac_keyset_count (mrnas) ;

  ks1 = ac_copy_keyset (mrnas, h) ;
  ks = ac_ksquery_keyset (ks1, "COUNT mrna_covered_by = 1", h) ;
  nn1 = ac_keyset_count (ks) ;
  ac_keyset_minus (ks1, ks) ;

  ks = ac_ksquery_keyset (ks1, "COUNT mrna_covered_by = 2", h) ;
  nn2 = ac_keyset_count (ks) ;
  ac_keyset_minus (ks1, ks) ;
  nn3  = ac_keyset_count (ks1) ;

  freeOutf ("%d mrna in the database have DNA\n", nMrna) ;
  if (nMrna)
    {
      freeOutf ("\t%d\t%6.2f%% are covered by a single clone\n", nn1, (100.0 * nn1)/nMrna) ;
      freeOutf ("\t%d\t%6.2f%% are covered by a 2 clones\n", nn2, (100.0 * nn2)/nMrna) ;
      freeOutf ("\t%d\t%6.2f%% are covered by a > 2 clones\n", nn3, (100.0 * nn3)/nMrna) ;
    }

  ks1 = ac_copy_keyset (mrnas, h) ;
  ks = ac_ksquery_keyset (ks1, "COUNT CDS_covered_by = 1", h) ;
  nn1 = ac_keyset_count (ks) ;
  ac_keyset_minus (ks1, ks) ;

  ks = ac_ksquery_keyset (ks1, "COUNT CDS_covered_by = 2", h) ;
  nn2 = ac_keyset_count (ks) ;
  ac_keyset_minus (ks1, ks) ;
  nn3  = ac_keyset_count (ks1) ;

  freeOutf ("%d mrna in the database have DNA\n", nMrna) ;
  if (nMrna)
    {
      freeOutf ("\t%d\t%6.2f%% have their CDS covered by a single clone\n", nn1, (100.0 * nn1)/nMrna) ;
      freeOutf ("\t%d\t%6.2f%% have their CDS covered by a 2 clones\n", nn2, (100.0 * nn2)/nMrna) ;
      freeOutf ("\t%d\t%6.2f%% have their CDS covered > 2 clones\n", nn3, (100.0 * nn3)/nMrna) ;
    }
 
  ksg = ac_dbquery_keyset (db, "find gene !cloud_gene && transcribed_gene", h) ;
  ks1 = ac_ksquery_keyset (ksg, "> transcribed_gene ; ! shedded_from ", h) ;
  ng = ac_keyset_count (ks1) ;
  ks = ac_ksquery_keyset (ks1, "gt_ag || gc_ag", h) ;
  nn1 = ac_keyset_count (ks) ;

  freeOutf ("%d non_cloud non_shed genes are defined\n", ng) ;
  if (ng)
    {
      freeOutf ("\t%d\t%6.2f%%\t of those have classic introns\n", nn1,  (100.0 * nn1)/ng) ;
      freeOutf ("\t%d\t%6.2f%%\t do not\n", ng -  nn1,  (100.0 * (ng-nn1))/ng) ;
    }

  
  ksp = ac_dbquery_keyset (db, "find mrna from_gene && dna && product", h) ;
  np = ac_keyset_count (ksp) ;
  ks1 = ac_ksquery_keyset (ksp, " COUNT product > 1", h) ;
  nn1 = ac_keyset_count (ks1) ;
  freeOutf ("%d mRNA have a product\n", np) ;
  if (np)
    {
      freeOutf ("\t%d\t%6.2f%%\t of those have more than 1 product\n", nn1,  (100.0 * nn1)/np) ;
    }
  
  ac_free (h) ;
}

/*************************************************************************/
/*************************************************************************/

static int cloudOrder (const void *a, const void *b)
{
  const CLOUD *ca = (const CLOUD *)a,  *cb = (const CLOUD *)b ;
  int n ;

  n = ca->map - cb->map ;  if (n) return n ;
  n = ca->a1 - cb->a1 ;  if (n) return n ;
  n = ca->a2 - cb->a2 ;  if (n) return n ;
  n = ca->isDown - cb->isDown ;  if (n) return n ;
  return ca->ii - cb->ii ;
}

static void sGeneCloud (AC_DB db, const char *dbName)
{
  FILE *f;
  char *cp ;
  int n, ii, jj, level, cumulUp = 0, cumulDown = 0 ;
  int nGene, nGeneUp, nGeneDown, nCoding, nIntron, nNeither ;
  CLOUD *up, *vp ;
  AC_HANDLE h = handleCreate () ; 
  DICT *dict = dictHandleCreate (100, h) ;
  Array cloud = arrayHandleCreate (100000, CLOUD, h) ;

  printf ("############################################################\n") ;
  printf ("## phase 34: sGeneCloud\n") ;


  /* create a huge table of all genes
     get the coordinates the strand the existence of introns
     and of coding potential

     we then compare the short non coding genes to the big good ones
  */
  f = filopen ("d14.geneCloud","txt","r") ;
  if (!f)
    messcrash ("cannot open %s/d14.geneCloud.txt", dbName) ;
  level = freesetfile (f, 0) ;
  ii = 0 ;
  nGene = nGeneUp = nGeneDown = nCoding = nIntron = nNeither = 0 ;
  while (freecard(level))
    {
      cp = freeword () ;
      if (!cp)
	continue ;
      up = arrayp (cloud, ii, CLOUD) ;
      up->ii = ii ;
      if ((cp = freeword ()))
	{
	  dictAdd (dict, cp, &n) ;
	  up->map = n ;
	}
      else
	continue ;
      if (freeint (&n))
	up->a1 = n ;
      else
	continue ;
      if (freeint (&n))
	up->a2 = n ;
      else
	continue ;
      if ((cp = freeword()))
	{
	  if (strstr (cp, "NULL"))
	    up->hasIntron = FALSE  ;
	  else
	    {
	      up->hasIntron = TRUE ;
	      nIntron++ ;
	    }
	}
      else
	continue ;

      if ((cp = freeword()))
	{
	  if (!strstr (cp, "NULL"))
	    up->isCoding = FALSE  ;
	  else
	    {
	      up->isCoding = TRUE ;
	      nCoding++ ;
	    }
	}
      else
	continue ;

      if (up->a1 < up->a2)
	{
	  nGeneDown++ ;
	  up->isDown = TRUE ;
	  if (up->hasIntron)
	    cumulDown += up->a2 - up->a1 ;
	}
      else
	{
	  int aa = up->a1 ;
	  up->a1 = up->a2 ;
	  up->a2 = aa ;
	  nGeneUp++ ;
	  up->isDown = FALSE ;
	  if (up->hasIntron)
	    cumulUp += up->a2 - up->a1 ;	  
	}
      if (!up->isCoding && !up->hasIntron)
	nNeither++ ;
      ii++ ; /* increment when garanteed complete */
    }
  freeclose (level) ;
  nGene = arrayMax (cloud) = ii ; /* may kill the last one if incomplete */
  
  printf ("Found\t%d genes,\t%d up,\t%d down\n"
	  , nGene, nGeneUp, nGeneDown) ;
  printf ("\t%d have classic introns,\t%d do not\n"
	  , nIntron, nGene - nIntron) ;
  printf ("\t%d are coding,\t%d are not\t%d are neither\n"
	  , nCoding, nGene - nCoding, nNeither) ;
  printf ("The genes with classic introns occupy\n") ;
  printf ("\t%d bp on the up strand\n", cumulUp) ;
  printf ("\t%d bp on the down strand\n", cumulDown) ;
  
  arraySort (cloud, cloudOrder) ;
  for (ii = 0, up = arrp (cloud, 0, CLOUD) ; ii < nGene ; ii++, up++)
    {
      if (!up->hasIntron) /* look for big gene */
	continue ;
      /* gobble small fish */
      for (jj = ii + 1, vp = arrp (cloud, jj, CLOUD) ; jj < nGene ; jj++, vp++)
	{
	  if (up->map != vp->map) break ;
	  if (vp->a1 > up->a2) break ;
	  if (vp->a2 > up->a2) continue ;
	  if (vp->isCoding || vp->hasIntron) continue ;
	  if (up->isDown == vp->isDown) vp->inSameStrand = TRUE ;
	  if (up->isDown != vp->isDown) vp->inOtherStrand = TRUE ;
	}
    }
  {
    /* report */
    int nInBoth = 0, nInSame = 0, nInOther = 0, nInNeither = 0 ;
    int tBoth = 0, tSame = 0, tOther = 0, tNeither = 0 ;
    
    for (ii = 0, up = arrp (cloud, 0, CLOUD) ; ii < nGene ; ii++, up++)
      { 
	if (up->isCoding || up->hasIntron) continue ;
	if (up->inSameStrand && up->inOtherStrand)
	  { nInBoth++ ; tBoth += up->a2 - up->a1 + 1 ; }
	else if (up->inSameStrand)
	  { nInSame++ ; tSame += up->a2 - up->a1 + 1 ; }
	else if (up->inOtherStrand)
	  { nInOther++ ; tOther += up->a2 - up->a1 + 1 ; }
	else 
	  { nInNeither++ ; tNeither += up->a2 - up->a1 + 1 ; }
      }
    printf ("The %d neither genes decompose in 4 exclusive classes\n", nNeither) ;
    printf ("\t%d are contained in gene with introns on each strand, they occupy\t%d bp\n", nInBoth, tBoth) ;
    printf ("\t%d are contained in gene with introns on same strand, they occupy\t%d bp\n", nInSame, tSame) ;
    printf ("\t%d are contained in gene with introns on other strand, they occupy\t%d bp\n", nInOther, tOther) ;
    printf ("\t%d are contained in gene with introns on neitherstrand, they occupy\t%d bp\n", nInNeither, tNeither) ;
    
  }
  
  messfree (h) ;
} /* sGeneCloud */

static void sGetMrnaIntrons (Array hits, AC_OBJ mrna, int m1)
{
  int ii, jj, p1, p2 ;
  HIT *up ;
  const char *ccp ;
  BOOL isCoding ;
  AC_OBJ product ;
  AC_TABLE splicing, products ;
  AC_HANDLE h = ac_new_handle () ;

  products = ac_tag_table (mrna, "Product", h) ;
  splicing = ac_tag_table (mrna, "Splicing", h) ;
  for (ii = 1 ; ii < splicing->rows - 1; ii++)
    {
      ccp = ac_table_printable (splicing, ii, 4, "") ;
      if (!strstr (ccp, "ntron"))
	continue ;
      ccp = ac_table_printable (splicing, ii, 5, "") ;
      if (strcmp (ccp, "gt_ag") && strcmp (ccp, "gt_ag"))
	continue ;
      /* we have a good classic intron */
      isCoding = FALSE ;
      up = arrayp (hits, arrayMax(hits), HIT) ;
      up->a1 = m1 + ac_table_int (splicing, ii, 0, 0) ;
      up->a2 = m1 + ac_table_int (splicing, ii, 1, 0) ;
      up->x1 = ac_table_int (splicing, ii - 1, 3, 0) ;
      up->x2 = ac_table_int (splicing, ii + 1, 2, 0) ;
      for (jj = 0 ; products && jj < products->rows ; jj++)
	{
	  product = ac_table_obj (products, jj, 0, h) ;
	  if (! ac_has_tag (product, "Best_product"))
	    continue ;
	  p1 = ac_table_int (products, jj, 1, 0) ;
	  p2 = ac_table_int (products, jj, 2, 0) ;
	  if (up->x1 >= p1 && up->x1 <= p2 &&
	      up->x2 >= p1 && up->x2 <= p2)
	    { isCoding = TRUE ; break ; }
	}
      up->x1 = isCoding ? 1 : 2 ; /* will sort coding before non coding */
      up->x2 = 3 ;
    }

  ac_free (h) ;
}

static BOOL sCountCodingIntrons (AC_OBJ tg, int *nCodingp, int *nNonCodingp, BOOL justVariantA)
{
  int ii, jj, m1 ;
  Array hits = arrayCreate (12, HIT) ;
  HIT *up, *vp ;
  BOOL isCoding ;
  AC_HANDLE h = handleCreate () ;
  AC_OBJ mrna ;
  AC_TABLE mrnas = ac_tag_table (tg, "mRNA", h) ;

  *nCodingp = *nNonCodingp = 0 ;
  /* accumulate all mrna introns */
  for (ii = 0 ; mrnas && ii < (justVariantA ? 1 : mrnas->rows)  ; ii++)
    {
      mrna = ac_table_obj (mrnas, ii, 0, h) ;
      if (ac_has_tag (mrna, "Gap")) 
	continue ;
      m1 = ac_table_int (mrnas, ii, 1, 0) ;
      if (m1 && mrna)
	sGetMrnaIntrons (hits, mrna, m1 - 1) ;
    }
  arraySort (hits, orderByA1) ; /* cdnaclone then est then a1 then x1 */
  arrayCompress (hits) ;
  for (ii = 0 ; ii < arrayMax(hits) ; ii++)
    {
      up = arrp (hits, ii, HIT) ;
      isCoding = up->x1 == 1 ? TRUE : FALSE ;
      if (isCoding)
	{
	  (*nCodingp)++ ;
	  /* a coding intron masks the same non coding in a different mrna */
	  for (jj = ii+1, vp = up+1 ; jj < arrayMax(hits) ; vp++, jj++)
	    {
	      if (vp->a1 != up->a1 || vp->a2 != up->a2)
		break ;
	      { up++ ; ii++ ; }
	    }
	}
      else
	(*nNonCodingp)++ ;
    }

  ac_free (h) ;
  return TRUE ;
}

static void sTgNonCodingIntrons (AC_DB db, BOOL doExportAceFile, BOOL justVariantA)
{
  int nTg = 0, n, ii,  nCoding = 0, nNonCoding = 0 ; 
  Array aa = arrayCreate (100, int) ;
  char *tgNames[300] ;
  AC_ITER iter ;
  AC_OBJ tg ;
  AC_HANDLE h = handleCreate () ;

  printf ("############################################################\n") ;
  if (! justVariantA)
    {
      printf ("## phase 31: sTgNonCodingIntrons\n") ;
      printf ("## computed on all introns of the gene, see is coding in at least one best product\n") ;
    }
  else
    {
      printf ("## phase 32: sTgNonCodingIntrons just variant a\n") ;
      printf ("## computed on all introns of variant a, see is coding in at least one best product\n") ;
    }

  memset (tgNames, 0, 100 * sizeof (char *)) ;
  iter = ac_query_iter (db, TRUE, "Find tg   mrna &&  ! shedded_from && Intron_boundaries", 0, h) ;
  while ((tg = ac_next_obj (iter)))
    { 
      sCountCodingIntrons (tg, &nCoding, &nNonCoding, justVariantA) ;
      if (nCoding + nNonCoding)
	{
	  nTg++ ;
	  if (nNonCoding)
	    {
	      array (aa, nNonCoding, int)++ ;	  
	      if (doExportAceFile)
		printf ("Transcribed_gene \"%s\"\nNb_non_coding_introns\t%d\n\n"
			, ac_name(tg), nNonCoding) ;
	    }
	  if (nNonCoding < 300)
	    {
	      if (!tgNames [nNonCoding])
		tgNames [nNonCoding] = strnew (ac_name(tg), h) ;
	    }
	}
      ac_free (tg) ;
    }
   if (! justVariantA)
     printf ("\nTable of non coding introns among %d non-shed genes with confirmed introns\n", nTg) ;
   else
     printf ("\nTable of non coding introns among %d .a mrna with introns\n", nTg) ;
  for (ii = 0 ; ii < arrayMax(aa) ; ii++)
    {
      n = arr (aa, ii, int) ;
      if (n)
	printf ("%6d genes with\t%3d non coding introns:\t%s\n"
		, n, ii, ii < 300 ?tgNames [ii] : "..." ) ;
    }
  printf ("//\n") ;
  arrayDestroy (aa) ;
  ac_free (h) ;
}
 
static void stMrnaMultiProteins (AC_OBJ mrna, int *n1p, int *n2p, int *n3p)
{
  int ii, jj, p1, p2, q1, q2, x1, x2, x3, upStop ;
  AC_TABLE tt ;
  AC_OBJ prod ;
  AC_HANDLE h = handleCreate () ;

  x1 = x2 = x3 = 0 ;
  tt = ac_tag_table (mrna, "Product", h) ;
  if (tt && tt->rows > 1)
    for (ii = 0 ; tt && ii < tt->rows - 1 ; ii++)
      {
	p1 = ac_table_int (tt, ii, 1, -1) ;
	p2 = ac_table_int (tt, ii, 2, -1) ;
	for (jj = ii + 1 ; jj < tt->rows ; jj++)
	  {
	    q1 = ac_table_int (tt, jj, 1, -1) ;
	    q2 = ac_table_int (tt, jj, 2, -1) ;
	    if (q1 > p1 + 2 * (p2 - p1) && q1 < p2 &&
		q2 > p2 + 2 * (p2 - q1))
	      { x3 = 1 ; break ; } /* frameshift */
	    else if (q1 > p1)
	      {
		prod = ac_table_obj (tt, jj, 0, 0) ;
		upStop = ac_tag_int (prod, "Up_stop", 1000) ; /* good if negative */
		ac_free (prod) ;
		q1 += upStop ; 
		if (!((q1 - p1)%3) && q1 > p2 && q1 < p2+7)
		  { x1 = 1 ; break ; } /* leaky stop */ 
		else if (q1 > p1 + 2 * (p2 - q1) && q1 < p2 &&
			 q2 > p2 + 2 * (p2 - q1))
		  { x3 = 1 ; break ; } /* frameshift */
		else
		  { x2 = 1 ; break ; } /* separated */
	      }
	  }
      }
  if (x1) (*n1p)++ ;
  if (x2) (*n2p)++ ;
  if (x3) (*n3p)++ ;
  
  ac_free (h) ;
} /* stMrnaMultiProteins */

static void sTgMultiProteins (AC_DB db, BOOL doShowAceFile) 
{
  int ntg_leaky, nmrna_leaky, ntg_sep, nmrna_sep, ntg_fs, nmrna_fs ;
  int ii, n1, n2, n3 ;
  char *tgNames[3] ;
  AC_TABLE tt ;
  AC_ITER iter ;
  AC_OBJ tg, mrna ;
  AC_HANDLE h = handleCreate () ;

  printf ("############################################################\n") ;
  printf ("## phase 33: sTgMultiProteins\n") ;

  tgNames[0] = tgNames[1] = tgNames[2] = 0 ;
  ntg_leaky = nmrna_leaky = ntg_sep = nmrna_sep = ntg_fs = nmrna_fs = 0 ;

  iter = ac_query_iter (db, TRUE, "Find tg  mrna &&  ! shedded_from &&  (gt_ag || gc_ag) && Product", 0, h) ;
  while ((tg = ac_next_obj (iter)))
    { 
      n1 = n2 = n3 = 0 ;
      tt = ac_tag_table (tg, "mrna", 0) ;
      for (ii = 0 ; tt && ii < tt->rows ; ii++)
	{
	  mrna = ac_table_obj (tt, ii, 0, h) ;
	  stMrnaMultiProteins (mrna, &n1, &n2, &n3) ;
	  ac_free (mrna) ;
	}
      if (n1)
	{
	  ntg_leaky++ ;
	  nmrna_leaky += n1 ;
	  if (doShowAceFile) 
	    printf ("leaky %s\n", ac_name(tg)) ;
	  if (!tgNames[0])
	    tgNames[0] = strnew (ac_name(tg), 0) ;
	}
      if (n2)
	{
	  ntg_sep++ ;
	  nmrna_sep += n2 ;
	  if (doShowAceFile) 
	    printf ("separated %s\n", ac_name(tg)) ;
	  if (!tgNames[1])
	    tgNames[1] = strnew (ac_name(tg), 0) ;
	}
      if (n3)
	{
	  ntg_fs++ ;
	  nmrna_fs += n3 ;
	  if (doShowAceFile) 
	    printf ("frame_shift %s\n", ac_name(tg)) ;
	  if (!tgNames[2])
	    tgNames[2] = strnew (ac_name(tg), 0) ;
	}
      ac_free (tt) ;
      ac_free (tg) ;
    }

  printf ("Found\t%d mrna from\t%d genes with leaky stop %s\n"
	  , nmrna_leaky, ntg_leaky, tgNames[0]) ;
  printf ("Found\t%d mrna from\t%d genes with well separated products\t%s\n"
	  , nmrna_sep, ntg_sep, tgNames[1]) ;
  printf ("Found\t%d mrna from\t%d genes with possible frame shift %s\n"
	  , nmrna_fs, ntg_fs, tgNames[2]) ;

  ac_free (h) ;
} /* stProdMultiProteins */

typedef struct { char title[60]; int count, nm, mrna, pfam, blast, psort, blastDate, pfamDate, psortDate ; } BPF ;
static void sTgXIntronsCount (AC_DB db, BOOL hasIntrons, BOOL isAntisens)
{
  int ii, nTg, iCds, iCds2, iTg  ;
  AC_KEYSET tgs, products, ks ;
  AC_ITER iter ;
  AC_OBJ tg = 0 ;
  BPF *bpf, *bpf2 ;
  BOOL hasBestProduct = FALSE ;
  char *cp ;
  AC_HANDLE h = handleCreate () ;
  Array aa = arrayHandleCreate (120, BPF, h) ;

  printf ("############################################################\n") ;
  if (hasIntrons && ! isAntisens)
    printf ("## phase 25: sTgXIntronsCount   CDS/Introns\n") ;
  else if (hasIntrons && isAntisens)
    printf ("## phase 26: sTgXIntronsCount   Antisens CDS/intron\n") ;
  else
    printf ("## phase 27: sTgXIntronsCount   CDS/No intron\n") ;

  printf ("## we count the CDS of the best product (Met to stop) or (first to stop if ! nh2-complete)\n") ; 
  if (hasIntrons && ! isAntisens)
    {
      tgs = ac_dbquery_keyset (db, "find tg  mrna &&  ! shedded_from && ( gt_ag || gc_ag) ", h) ;
      nTg = ac_keyset_count (tgs) ;
      printf ("%9d %s\n",nTg, "Genes non-shed with gt_ag || gc_ag introns") ;
    }
  else if (hasIntrons && isAntisens)
    {
      tgs = ac_dbquery_keyset (db, "find tg  mrna &&  ! shedded_from && ( gt_ag || gc_ag) && COUNT { >antisens_to ;  ( gt_ag || gc_ag)} > 0 ", h) ;
      nTg = ac_keyset_count (tgs) ;
      printf ("%9d %s\n",nTg, "Genes non-shed with gt_ag || gc_ag introns antisens to genes with gt_ag || gc_ag ") ;
    }
  else
    {
      tgs = ac_dbquery_keyset (db, "find tg mrna &&  ! shedded_from &&  !gt_ag && !gc_ag ", h) ;
      nTg = ac_keyset_count (tgs) ;
      printf ("%9d %s\n",nTg, "Genes non-shed without gt_ag || gc_ag introns") ;
    }

  bpf = arrayp (aa, 4, BPF) ; strcpy (bpf->title, ">= 300 AA" ) ;
  bpf = arrayp (aa, 3, BPF) ; strcpy (bpf->title, ">= 200 AA" ) ;
  bpf = arrayp (aa, 2, BPF) ; strcpy (bpf->title, ">= 100 AA" ) ;
  bpf = arrayp (aa, 1, BPF) ; strcpy (bpf->title, " < 100 AA" ) ;
  bpf = arrayp (aa, 0, BPF) ; strcpy (bpf->title, "any other") ;

  ks = ac_dbquery_keyset (db, "Find product best_product", 0) ;
  if (ac_keyset_count (ks))
    hasBestProduct = TRUE ;
  ac_free (ks) ;

  iter = ac_keyset_iter (tgs, TRUE, h) ;
  iTg = 0 ;
  while (ac_free (tg), iTg++, (tg = ac_next_obj (iter)))
    {
      iCds = 0 ;
      for (ii = 300 ; !iCds && ii >= 0 ; ii -= 100)
	{
	  cp =  messprintf (">mrna !gap_length ; >product %s  coding_length > %d", hasBestProduct ? "best_product &&" : "", 3*ii) ;
	  ks = ac_objquery_keyset (tg, cp,  0) ;
	  if (ac_keyset_count (ks)) iCds = 1 + ii/100 ;
	  ac_free (ks) ;
	}
      
      iCds2 = 0 ;
      for (ii = 300 ; !iCds2 && ii >= 0 ; ii -= 100)
	{
	  cp =  messprintf (">mrna !gap_length ; >product best_product && ( (! (First_Kozak && First_atg > 1 && NH2_complete) && coding_length > %d)  ||  ( (First_Kozak && First_atg > 1 && NH2_complete) && [coding_length - 3 * first_atg] > %d) ) ", 3*ii, 3*ii) ;

	  ks = ac_objquery_keyset (tg, cp,  0) ;
	  if (ac_keyset_count (ks)) iCds2 = 101 + ii/100 ;
	  ac_free (ks) ;
	}
      
      bpf = arrayp (aa, iCds, BPF) ;
      bpf->count++ ;

      bpf2 = arrayp (aa, iCds2, BPF) ;
      bpf2->count++ ;

      ks = ac_objquery_keyset (tg, ">read ref_seq", 0) ;
      if (ac_keyset_count (ks)) { bpf->nm++ ; bpf2->nm++ ; }
      ac_free (ks) ;

      ks = ac_objquery_keyset (tg, ">read ref_mrna && ! ref_seq", 0) ;
      if (ac_keyset_count (ks)) { bpf->mrna++ ; bpf2->mrna++ ; }
      ac_free (ks) ;

      products = ac_objquery_keyset (tg, ">mrna !gap_length ; >Product best_product", 0) ;

      ks = ac_ksquery_keyset (products, "PFam", 0) ;
      if (ac_keyset_count (ks)) bpf->pfam++ ;
      ac_free (ks) ;

      ks = ac_ksquery_keyset (products, ">Kantor ; PFam_date", 0) ;
      if (ac_keyset_count (ks)) bpf->pfamDate++ ;
      ac_free (ks) ;

      ks = ac_ksquery_keyset (products, "BlastP", 0) ;
      if (ac_keyset_count (ks)) bpf->blast++ ;
      ac_free (ks) ;

      ks = ac_ksquery_keyset (products, ">Kantor ; Blastp_date", 0) ;
      if (ac_keyset_count (ks)) bpf->blastDate++ ;
      ac_free (ks) ;

      ks = ac_ksquery_keyset (products, "Psort_title", 0) ;
      if (ac_keyset_count (ks)) bpf->psort++ ;
      ac_free (ks) ;


      ks = ac_ksquery_keyset (products, ">Kantor ; Psort_date", 0) ;
      if (ac_keyset_count (ks)) bpf->psortDate++ ;
      ac_free (ks) ;

      ac_free (products) ;
    }

 

  printf ("%15s\t%8s\t%15s\t%9s/%s\t%13s/%s\t%12s/%s\n","CDS  \\ "
	  ,"#Gene", "NM in gene", "Pfam", "Date", "Blastp", "Date", "Psort", "Date") ;
  
  for (ii = 4 ; ii >= 0 ; ii--)
    {
      bpf = arrp (aa, ii, BPF) ;
      printf ("%15s\t%8d\t%8d(%2.1f%%)\t%8d/%d(%2.0f%%)%8d/%d(%2.0f%%)%8d/%d(%2.0f%%)\n"
	      , bpf->title, bpf->count 
	      , bpf->nm, 100.0 * bpf->nm/(bpf->count? bpf->count: 1)
	      , bpf->pfam, bpf->pfamDate, 100.0 * bpf->pfam/(bpf->pfamDate ? bpf->pfamDate : 1)
	      , bpf->blast, bpf->blastDate, 100.0 * bpf->blast/(bpf->blastDate?bpf->blastDate: 1)
	      , bpf->psort, bpf->psortDate, 100.0 * bpf->psort/(bpf->psortDate?bpf->psortDate: 1)
	      ) ;
      
    }
 
  printf ("%15s\t%12s\t%12s\t%12s,  removing the Kozak to ATG segment\n in 5prime complete open mRNAs\n","CDS  \\ "
	  ,"#Gene", "NM in gene", "  mRNA in gene") ;
  
  for (ii = 104 ; ii >= 100 ; ii--)
    {
      bpf2 = arrp (aa, ii, BPF) ;
      bpf = arrp (aa, ii - 100, BPF) ;
      printf ("%15s\t%12d\t%12d\t%12d\n"
	      , bpf->title, bpf2->count 
	      , bpf2->nm, bpf2->mrna
	      ) ;
      
    }  

  printf ("%15s\t%12s\t%12s\t%12s,  removing the Kozak to ATG segment\n in 5prime complete open mRNAs\n","CDS  \\ "
	  ,"#Gene", "NM in gene", "  mRNA in gene") ;
  
  for (ii = 104 ; ii >= 100 ; ii--)
    {
      bpf2 = arrp (aa, ii, BPF) ;
      bpf = arrp (aa, ii - 100, BPF) ;
      printf ("%15s\t%12d\t%12d\t%12d\n"
	      , bpf->title, bpf2->count 
	      , bpf2->nm, bpf2->mrna
	      ) ;
      
    }  
  printf ("any other means gap or no best_product or no CDS_length in best_product\n" ) ;
  ac_free (h) ;
} /* sTgNoIntronsCount */

#ifdef JUNK
static void sMrnaIntronsCount (AC_DB db)
{
  int ii, nMrna ; 
  int in5, in3, inBoth, inNeither, myIn5 = 0, myIn3 = 0 ;
  AC_KEYSET mrnas ;
  AC_ITER iter ;
  AC_OBJ mrna ;
  AC_TABLE tt ;
  AC_HANDLE h = handleCreate () ;

  mrnas = ac_dbquery_keyset (db, "find mRNA", h) ;
  nMrna = ac_keyset_count (mrnas) ;
  printf ("%20s\t%9d\n","mRNAs", nMrna) ;

  mrnas = ac_dbquery_keyset (db, "find mRNA gt_ag || gc_ag", h) ;
  nMrna = ac_keyset_count (mrnas) ;
  printf ("%20s\t%9d\n","mRNAs with introns", nMrna) ;

  iter = ac_keyset_iter (mrnas, TRUE, h) ;
  in5 = in3 = inBoth = inNeither = 0 ;
  while ((mrna = ac_next_obj (iter)))
    {
      tt = ac_tag_table (mrna, "Splicing", h) ;
      for (ii = 0 ; tt && ii < tt->rows ; ii++)
	{
	}
      ac_free (tt) ;
      ac_free (mrna) ;
      if (myIn5) in5++ ;
      if (myIn3) in3++ ;
      if (myIn3 && myIn5) inBoth++ ;
      if (!myIn3 && !myIn5) inNeither++ ;
    }

  if (!nMrna) nMrna = 1 ;
  printf ("%25s\t%6d  %5.2f", "Introns in 5'UTR", in5, 100.0 * ((float)in5)/nMrna) ;
  printf ("%25s\t%6d  %5.2f", "Introns in 3'UTR", in3, 100.0 * ((float)in3)/nMrna) ;
  printf ("%25s\t%6d  %5.2f", "Introns in both", inBoth, 100.0 * ((float)inBoth)/nMrna) ;
  printf ("%25s\t%6d  %5.2f", "Introns in neither", inNeither, 100.0 * ((float)inNeither)/nMrna) ;

  ac_free (h) ;
} /* sMrnaIntronsCount */
#endif
static void sEstQuality (AC_DB db)
{
  int n1, nn, ll, ali, err, ii, jj, ll1, err1, ali1 ;
  AC_ITER iter ;
  AC_OBJ est, tg ;
  AC_TABLE tt ;
  char *myquery[] =
  {
    "find EST ref_seq && from_gene" ,
    "find EST  !(ref_seq) && (ref_mrna || is_mrna) && from_gene" ,
    "find EST !(ref_seq) && !(ref_mrna || is_mrna) && from_gene",
    "find EST  from_gene"
  } ;
  char *mytitle[] = { "NM", "mRNA", "EST", "EST" } ;
  float myAliRate, aliRate [] = { 99.9, 99.0, 95.0, 90.0, 80.0, 70.0, 60.0, 50.0, 0.0, -1 } ;
  float myErrRate, errRate [] = { 0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 100,  -1} ;
  float errRate1 [] = { 0, .5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 9.0, 10.0, 100, -1 } ;
  int nAliRate, nErrRate, nErrRate1, nAliRateMax ;
  int *res, *res1 ;
  AC_HANDLE h = handleCreate () ;
  
  printf ("############################################################\n") ;
  printf ("## phase 20: sEstQuality\n") ;

  for (nAliRate = 0 ; aliRate [nAliRate] > -1 ; nAliRate++) ; /* count aliRates */
  for (nErrRate = 0 ; errRate [nErrRate] > -1 ; nErrRate++) ; /* count errRates */
  for (nErrRate1 = 0 ; errRate1 [nErrRate1] > -1 ; nErrRate1++) ; /* count errRates */
  printf ("%12s\t%9s\t%9s\t%12s\t%12s\t%12s\n","Type","Number", "#ali","Length","Ali","Err") ;
  nAliRateMax = nAliRate ;
  /* #est #long #ali #err */
  for (jj = isWorm ? 3 : 0 ; jj < (isWorm ? 4 : 3) ; jj++)
    {
      res = (int *) messalloc ( nAliRate * nErrRate * sizeof(int)) ;
      res1 = (int *) messalloc ( nAliRate * nErrRate1 * sizeof(int)) ;
      n1 = nn = ali = ll = err = 0 ;
      iter = ac_query_iter (db, TRUE, myquery [jj], 0, h) ;
      while ((est = ac_next_obj (iter)))
	{
	  n1++ ;
	  tt = ac_tag_table (est, "From_gene", h) ;
	  for (ii = 0 ; tt && ii < tt->rows ; ii++)
	    { 
	      ll1 = ac_table_int (tt, ii, 1, 0) ; if (ll1 < 1) continue ;
	      nn++ ;
	      tg = ac_table_obj (tt, ii, 0, h) ;
	      ll1 = ac_table_int (tt, ii, 1, 0) ; ll += ll1 ;
	      ali1 = ac_table_int (tt, ii, 3, 0) ; ali += ali1 ;
	      err1 = ac_table_int (tt, ii, 4, 0) ; err += err1 ;
	      if (0) break ; /* single ali per est ! */
	      if (0) printf ("\"%s\"\t\"%s\"\t%d\t%d\t%d\n", ac_name(est), ac_name(tg), ll1, ali1, err1) ;
	      
	      if (ll1 <= 0) ll1 = 1 ; 
	      if (ali1 <= 0) ali1 = 1 ; 
	      if (ali1 > ll1) ali1 = ll1 ;
	      myAliRate = 100.0 * ((float)ali1)/ll1; myErrRate = 100.0 * ((float) err1)/ali1 ;
	      for (nAliRate = 0 ; aliRate [nAliRate] >  myAliRate ; nAliRate++) ; /* count aliRates */
	      for (nErrRate = 0 ; errRate [nErrRate] < myErrRate ; nErrRate++) ; /* count errRates */
	      for (nErrRate1 = 0 ; errRate1 [nErrRate1]  < myErrRate ; nErrRate1++) ; /* count errRates */
	      res [nAliRate + nAliRateMax * nErrRate]++ ;
	      res1 [nAliRate + nAliRateMax * nErrRate1]++ ;
	    }
	  ac_free (tt) ;
	  ac_free (est) ;
	  if (n1 >= 200000) break ;
	}
      if (!nn) nn = 1 ;
      if (nn > 1)
	{
	  printf ("%12s\t%9d\t%9d\t%12d\t%12d\t%12d\n"
		  , mytitle[jj], n1, nn, ll, ali, err) ;
	  if (!ali) ali = 1 ;
	  printf ("%10s\t%31d\t%12d\t%11.2f%%\n"
		  ,"average", ll/nn, ali/nn, ((float)err)/ali) ;
	}
      
      if (nn > 1)
	{
	   int z, zz = 0 ;
	  
	   printf ("Qualities\t%s\n%7s", mytitle[jj], "errMax") ;
	  for (nAliRate = 0 ; aliRate [nAliRate] > -1 ; nAliRate++)
	    printf ("%6.2f",aliRate [nAliRate]) ;
	  printf ("\n") ;
	  for (nErrRate = 0 ; errRate [nErrRate] > -1; nErrRate++)
	    {
	      printf ("\t%6.2f", errRate [nErrRate]) ;
	      for (nAliRate = 0 ; aliRate [nAliRate] > -1 ; nAliRate++)
		{
		  z = res [nAliRate + nAliRateMax * nErrRate] ; 
		  printf ("\t%5d", z) ;
		  zz += z ;
		}
	      printf ("\n") ;
	    }
	  printf ("\n total = %d\n", zz) ;
	}
      if (nn > 1)
	{
	  int z, zz = 0 ;
	  
	  printf ("Qualities\t%s\n%7s", mytitle[jj], "errMax") ;
	  for (nAliRate = 0 ; aliRate [nAliRate] > -1 ; nAliRate++)
	    printf ("\t%6.2f",aliRate [nAliRate]) ;
	  printf ("\n") ;
	  for (nErrRate1 = 0 ; errRate1 [nErrRate1] > -1; nErrRate1++)
	    {
	      if ( errRate1 [nErrRate1 + 1] == -1)
		printf ("\t%6s", "more") ;
	      else
		printf ("\t%6.2f", errRate1 [nErrRate1]) ;
	      for (nAliRate = 0 ; aliRate [nAliRate] > -1 ; nAliRate++)
		{
		  z = res1 [nAliRate + nAliRateMax * nErrRate1] ; 
		  printf ("\t%5d", z) ;
		  zz += z ;
		}
	      printf ("\n") ;
	    }
	  printf ("\n total = %d\n", zz) ;
	}
      if (nn > 1)
	{
	  float z, zz = 0 ;
	  
	  printf ("Qualities\t%s percent\n%6s", mytitle[jj], "errMax") ;
	  for (nAliRate = 0 ; aliRate [nAliRate] > -1 ; nAliRate++)
	    printf ("\t%6.2f",aliRate [nAliRate]) ;
	  printf ("\n") ;
	  for (nErrRate = 0 ; errRate [nErrRate] > -1; nErrRate++)
	    {
	      if ( errRate [nErrRate + 1] == -1)
		printf ("\t%6s", "more") ;
	      else
		printf ("\t%6.2f", errRate [nErrRate]) ;
	      for (nAliRate = 0 ; aliRate [nAliRate] > -1 ; nAliRate++)
		{
		  z = 100.0 * ((float)res [nAliRate + nAliRateMax * nErrRate])/nn ; 
		  printf ("\t%6.2f", z) ;
		  zz += z ;
		}
	      printf ("\n") ;
	    }
	  printf ("\n total = %6.2f\n", zz) ;
	}
      if (nn > 1)
	{
	  float z, zz = 0 ;

	  printf ("Qualities\t%s percent\n%6s", mytitle[jj], "errMax") ;
	  for (nAliRate = 0 ; aliRate [nAliRate] > -1 ; nAliRate++)
	    printf ("\t%6.2f",aliRate [nAliRate]) ;
	  printf ("\n") ;
	  for (nErrRate1 = 0 ; errRate1 [nErrRate1] > -1; nErrRate1++)
	    {
	      if ( errRate1 [nErrRate1 + 1] == -1)
		printf ("\t%6s", "more") ;
	      else
		printf ("\t%6.2f", errRate1 [nErrRate1]) ;
	      for (nAliRate = 0 ; aliRate [nAliRate] > -1 ; nAliRate++)
		{
		  z = 100.0 * ((float)res1 [nAliRate + nAliRateMax * nErrRate1])/nn ;
		  printf ("\t%6.2f", z) ;
		  zz += z ;
		}
	      printf ("\n") ;
	    }
	  printf ("\n total = %6.2f\n", zz) ;
	}
      messfree (res) ; 
      messfree (res1) ; 
    }
  
  ac_free (h) ;
} /* sEstQuality */

static void sTgCount (AC_DB db)
{
  int nG, nShed, nShed2Gene, nMrna ;
  AC_KEYSET ks = 0 ;
  AC_HANDLE h = handleCreate () ;

  printf ("############################################################\n") ;
  printf ("## phase 10: sTgCount\n") ;

  ks = ac_dbquery_keyset (db, "find gene transcribed_gene", h) ;
  nG = ac_keyset_count (ks) ;
  ac_free (ks) ;
  ks = ac_dbquery_keyset (db, "find tg  mrna &&  shedded_from && mrna ; > mrna", h) ;
  nShed = ac_keyset_count (ks) ;
  ac_free (ks) ;
  ks = ac_dbquery_keyset (db, "find tg  mrna &&  shedded_from && mrna ; > gene", h) ;
  nShed2Gene = ac_keyset_count (ks) ;
  ac_free (ks) ;
  ks = ac_dbquery_keyset (db, "find mRNA", h) ;
  nMrna = ac_keyset_count (ks) ;
  ac_free (ks) ;
 
  printf ("Total number of genes with clones, shed genes, mRNAs(including those from shed genes)\n") ;
  printf ("%9s\t%12s\t%18s\t%9s\n","Genes", "mRNAs", "Gene with shed", "Shed mRNA") ;
  printf ("%9d\t%12d\t%18d\t%9d\n", nG, nMrna, nShed2Gene, nShed) ;

  ac_free (h) ;
} /* sTgCount */

static void sGeneIntronsCount (AC_DB db, BOOL hasIntrons)
{
  int nTg, nPfam, nPfamNm, nPsort, nPsortNm, nNm ;
  AC_KEYSET tgs, genes, pfams, nms, psorts ;
  AC_HANDLE h = handleCreate () ;

  printf ("############################################################\n") ;
  if (hasIntrons)
    printf ("## phase 10.bis: sGeneIntronsCount   Introns/NM \n") ;
  else
    printf ("## phase 10.ter: sGeneIntronsCount   No intron/NM \n") ;

  if (hasIntrons)
    tgs = ac_dbquery_keyset (db, "find tg  mrna &&  ! shedded_from && (gt_ag || gc_ag)", h) ;
  else
    tgs = ac_dbquery_keyset (db, "find tg  mrna &&  ! shedded_from && ! (gt_ag || gc_ag)", h) ;

  genes = ac_ksquery_keyset (tgs, " > Gene", h) ;
  /* AC_KEYSET products = ac_ksquery_keyset (genes, " > Product", h) ; */
  nTg = ac_keyset_count (genes) ;
 
  nms = ac_ksquery_keyset (tgs, "COUNT {>read ; ref_seq} > 0 ; > gene", h) ;
  nNm = ac_keyset_count (nms) ;
 

  pfams = ac_ksquery_keyset (genes, "COUNT {>Product ; PFAM } > 0", h) ;
  nPfam = ac_keyset_count (pfams) ;
  nPfamNm = ac_keyset_and (pfams, nms) ;

  psorts = ac_ksquery_keyset (genes, "COUNT {>Product ; Psort_title } > 0", h) ;
  nPsort = ac_keyset_count (psorts) ;
  nPsortNm = ac_keyset_and (psorts, nms) ;

  printf ("%20s\t%s\t%9d\n","Genes non-shed ", hasIntrons ? "with introns" : "without intron", nTg) ;
  printf ("\t%19s\t%9d\n", "with NM", nNm) ;
  printf ("\t%19s\t%9d\t%9s\t%9d\t%9s\t%9d\n", "with PFAM", nPfam,"and NM", nPfamNm, "No NM", nPfam - nPfamNm ) ;
  printf ("\t%19s\t%9d\t%9s\t%9d\t%9s\t%9d\n", "with PSORT", nPsort, "and NM", nPsortNm, "No NM", nPsort - nPsortNm) ;

  ac_free (h) ;
} /* sGeneIntronsCount */
static void sGeneMrnaIntronsCount (AC_DB db, BOOL hasIntrons)
{
  int nTg, nPfam, nPfamNm, nPsort, nPsortNm, nNm ;
  AC_KEYSET tgs, genes, pfams, nms, psorts ;
  AC_HANDLE h = handleCreate () ;

  printf ("############################################################\n") ;
  if (hasIntrons)
    printf ("## phase 10.bis: sMrnaIntronsCount   Introns/NM \n") ;
  else
    printf ("## phase 10.ter: sMrnaIntronsCount   No intron/NM \n") ;

  if (hasIntrons)
    tgs = ac_dbquery_keyset (db, "find tg  mrna &&  ! shedded_from && (gt_ag || gc_ag)", h) ;
  else
    tgs = ac_dbquery_keyset (db, "find tg  mrna &&  ! shedded_from && ! (gt_ag || gc_ag)", h) ;

  genes = ac_ksquery_keyset (tgs, " > Gene", h) ;
  /* products = ac_ksquery_keyset (genes, " > Product", h) ; */
  nTg = ac_keyset_count (genes) ;
 
  nms = ac_ksquery_keyset (tgs, "COUNT {>read ; ref_seq} > 0 ; > gene", h) ;
  nNm = ac_keyset_count (nms) ;
 

  pfams = ac_ksquery_keyset (genes, "COUNT {>Product ; PFAM } > 0", h) ;
  nPfam = ac_keyset_count (pfams) ;
  nPfamNm = ac_keyset_and (pfams, nms) ;

  psorts = ac_ksquery_keyset (genes, "COUNT {>Product ; Psort_title } > 0", h) ;
  nPsort = ac_keyset_count (psorts) ;
  nPsortNm = ac_keyset_and (psorts, nms) ;

  printf ("%20s\t%s\t%9d\n","Genes non-shed ", hasIntrons ? "with introns" : "without intron", nTg) ;
  printf ("\t%19s\t%9d\n", "with NM", nNm) ;
  printf ("\t%19s\t%9d\t%9s\t%9d\t%9s\t%9d\n", "with PFAM", nPfam,"and NM", nPfamNm, "No NM", nPfam - nPfamNm ) ;
  printf ("\t%19s\t%9d\t%9s\t%9d\t%9s\t%9d\n", "with PSORT", nPsort, "and NM", nPsortNm, "No NM", nPsort - nPsortNm) ;

  {
    int nm, nmi, nmu, nig, nug ;
    AC_KEYSET mrnas, mi, mu, mig, mug ;

    mrnas = ac_dbquery_keyset (db, "follow mrna", h) ;
    nm = ac_keyset_count (mrnas) ;
    mi = ac_ksquery_keyset (mrnas, "gt_ag || gc_ag", h) ;
    nmi = ac_keyset_count (mi) ;
    mu = ac_ksquery_keyset (mrnas, "! (gt_ag || gc_ag)", h) ;
    nmu = ac_keyset_count (mu) ;
    mig = ac_ksquery_keyset (mi, "Gap", h) ;
    nig = ac_keyset_count (mig) ;
    mug = ac_ksquery_keyset (mu, "Gap", h) ;
    nug = ac_keyset_count (mug) ;
    printf ("\n%15s\t%12s\t%12s\t%12s\n", "These genes contain", "mRNA", "no-gap", "with-gap") ; 
    printf ("%15s\t%12d %12d %12d\n", "gt-ag/gc-ag", nmi, nmi - nmi, nig) ; 
    printf ("%15s\t%12d %12d %12d\n", "no gn_ag", nmu, nmu - nug, nug) ; 
    printf ("%15s\t%12d %12d %12d\n", "Total", nm, nm - nig - nug, nig + nug) ; 
  }
    
    ac_free (h) ;
} /* sGeneMrnaIntronsCount */

#ifdef JUNK
static int sAliCountAlignments (AC_KEYSET ks)
{
  AC_HANDLE h = handleCreate () ;
  AC_ITER iter ;
  AC_OBJ est = 0 ;
  AC_KEYSET ks2 ;
  int nAli = 0 ;

  iter = ac_keyset_iter (ks,TRUE,  h) ;
  while (ac_free (est), est = ac_iter_obj (iter))
    {
      ks2 = ac_objquery_keyset (est, ">From_gene", 0) ;
      nAli += ac_keyset_count (ks2) ;
      ac_free (ks2) ;
    }
  ac_free (h) ;

  return nAli ;
} /* sAliCountAlignments */
#endif

static void sAliCount (AC_DB db, BOOL isMulti)
{
  int nTotal=0, nNm=0, nMrna=0, nEst = 0, nTr = 0 ;
  int nTotalAli=0, nNmAli=0, nMrnaAli=0, nEstAli = 0 ;
  int nTotalInMrna = 0, nNmInMrna = 0, nMrnaInMrna = 0, nEstInMrna = 0 ;
  int nTotalInMrnaSpliced = 0, nNmInMrnaSpliced = 0, nMrnaInMrnaSpliced = 0, nEstInMrnaSpliced = 0 ;
  int nTotalFromGene = 0, nNmFromGene = 0, nMrnaFromGene = 0, nEstFromGene = 0 ;
  int nTotalFromGeneSpliced = 0, nNmFromGeneSpliced = 0, nMrnaFromGeneSpliced = 0, nEstFromGeneSpliced = 0 ;
  int nTotalFromGeneCoding = 0, nNmFromGeneCoding  = 0, nMrnaFromGeneCoding  = 0, nEstFromGeneCoding  = 0 ;
  int nTotalFromGeneSplicedNonCoding = 0, nNmFromGeneSplicedNonCoding = 0, nMrnaFromGeneSplicedNonCoding = 0, nEstFromGeneSplicedNonCoding = 0 ;
  int nTotalFromGeneSplicedMargeCoding = 0, nNmFromGeneSplicedMargeCoding = 0, nMrnaFromGeneSplicedMargeCoding = 0, nEstFromGeneSplicedMargeCoding = 0 ;
  AC_KEYSET ks = 0, ks1 = 0, ks2 = 0, ks3 = 0, ksAli = 0, ksMrna = 0, ksMrnaSpliced = 0, ksGenes = 0, ksGenesSpliced = 0, ksGenesCoding = 0, ksGenesSplicedNonCoding = 0, ksGenesSplicedMargeCoding = 0  ;
  AC_HANDLE h = handleCreate () ;


  printf ("\n#############################################################\n") ;
   if (!isMulti)
     printf ("## phase 3: sAliCount   Is Aceview aligning all reads ?\n") ;
   else
     printf ("## phase 4: sAliCount   Multi alignments ?\n") ;

  if (1)
    {
      AC_OBJ clo ;
      AC_ITER iter = ac_dbquery_iter (db, "find clone Nb_RefSeq && Nb_RefMrna", h) ;

      /* get from db if available */
      if (iter && (clo = ac_iter_obj (iter)))
	{
	  nNm = ac_tag_int (clo, "Nb_RefSeq", nNm) ;
	  nMrna = ac_tag_int (clo, "Nb_RefMrna", nMrna) ;
	  nEst = ac_tag_int (clo, "Nb_EST", nEst) ;
	  nTr = ac_tag_int (clo, "Nb_Trace", nTr) ;
	  ac_free (clo) ;
	}
      if (nTr > 1) nEst += nTr ;
      if (isMulti)
	ks = ac_dbquery_keyset (db, "find est ref_seq && COUNT from_gene > 1", h) ;
      else
	ks = ac_dbquery_keyset (db, "find sequence ref_seq && (from_gene || is_buried_under)", h) ;
      nNmAli = ac_keyset_count (ks) ;
      ksAli = ks ; 
      ks1 = ac_ksquery_keyset (ks, ">In_mRNA", h) ;
      nNmInMrna = ac_keyset_count (ks1) ;
      ks2 = ac_ksquery_keyset (ks1, "gt_ag || gc_ag", h) ; 
      nNmInMrnaSpliced = ac_keyset_count (ks2) ;   
      ksMrnaSpliced = ks2 ;
      ksMrna = ks1 ;
      ks1 = ac_ksquery_keyset (ks, ">From_gene", h) ;
      ksGenes = ks1 ;
      nNmFromGene = ac_keyset_count (ks1) ;
      ks2 = ac_ksquery_keyset (ks1, "gt_ag || gc_ag", h) ;
      nNmFromGeneSpliced = ac_keyset_count (ks2) ; 
      ksGenesSpliced = ks2 ;
      ks2 = ac_ksquery_keyset (ks1, ">gene", h) ;

      ks3 = ac_ksquery_keyset (ks2, "pastille_spliced_non_coding", h) ;
      nNmFromGeneSplicedNonCoding = ac_keyset_count (ks3) ; 
      ksGenesSplicedNonCoding =  ks3 ;

      ks3 = ac_ksquery_keyset (ks2, "Pastille_marginally_coding", h) ;
      nNmFromGeneSplicedMargeCoding = ac_keyset_count (ks3) ; 
      ksGenesSplicedMargeCoding = ks3 ;

      ks3 = ac_ksquery_keyset (ks2, "pastille_coding", h) ;
      nNmFromGeneCoding = ac_keyset_count (ks3) ; 
      ksGenesCoding = ks3 ;

      ac_free (ks2) ;
      
      if (isMulti)
	ks = ac_dbquery_keyset (db, "find est !ref_seq  && ref_mrna && COUNT from_gene > 1", h) ;
      else
	ks = ac_dbquery_keyset (db, "find est !ref_seq  && ref_mrna && (from_gene || is_buried_under)", h) ;
      nMrnaAli = ac_keyset_count (ks) ;
      ac_keyset_or (ksAli, ks) ;

      ks1 = ac_ksquery_keyset (ks, ">In_mRNA", h) ;
      nMrnaInMrna = ac_keyset_count (ks1) ;
      ks2 = ac_ksquery_keyset (ks1, "gt_ag || gc_ag", h) ; 
      nMrnaInMrnaSpliced = ac_keyset_count (ks2) ;   
      ac_keyset_or (ksMrnaSpliced, ks2) ;
      ac_free (ks2) ;
      ac_keyset_or (ksMrna, ks1) ;
      ac_free (ks1) ;
      ks1 = ac_ksquery_keyset (ks, ">From_gene", h) ;
      ac_keyset_or (ksGenes, ks1) ;
      nMrnaFromGene = ac_keyset_count (ks1) ;
      ks2 = ac_ksquery_keyset (ks1, "gt_ag || gc_ag", h) ;
      nMrnaFromGeneSpliced = ac_keyset_count (ks2) ; 
      ac_keyset_or (ksGenesSpliced, ks2) ;
      ac_free (ks2) ;
      ks2 = ac_ksquery_keyset (ks1, ">gene", h) ;

      ks3 = ac_ksquery_keyset (ks2, "pastille_spliced_non_coding", h) ;
      nMrnaFromGeneSplicedNonCoding = ac_keyset_count (ks3) ; 
      ac_keyset_or (ksGenesSplicedNonCoding, ks3) ;
      ac_free (ks3) ;

      ks3 = ac_ksquery_keyset (ks2, "Pastille_marginally_coding", h) ;
      nMrnaFromGeneSplicedMargeCoding = ac_keyset_count (ks3) ; 
      ac_keyset_or (ksGenesSplicedMargeCoding, ks3) ;
      ac_free (ks3) ;

      ks3 = ac_ksquery_keyset (ks2, "pastille_coding", h) ;
      nMrnaFromGeneCoding = ac_keyset_count (ks3) ; 
      ac_keyset_or (ksGenesCoding, ks3) ;
      ac_free (ks3) ;

      ac_free (ks2) ;
      ac_free (ks1) ;
      ac_free (ks) ;
      
      if (isMulti)
	ks = ac_dbquery_keyset (db, "find sequence !ref_seq  && !ref_mrna && COUNT from_gene > 1", h) ;
      else
	ks = ac_dbquery_keyset (db, "find sequence !ref_seq  && !ref_mrna && (is_buried_under || from_gene)", h) ;
      nEstAli = ac_keyset_count (ks) ;
      ac_keyset_or (ksAli, ks) ;
      ks1 = ac_ksquery_keyset (ks, ">In_mRNA", h) ;
      nEstInMrna = ac_keyset_count (ks1) ;
      ks2 = ac_ksquery_keyset (ks1, "gt_ag || gc_ag", h) ; 
      nEstInMrnaSpliced = ac_keyset_count (ks2) ;   
      ac_keyset_or (ksMrnaSpliced, ks2) ;
      ac_free (ks2) ;
      ac_keyset_or (ksMrna, ks1) ;
      ac_free (ks1) ;
      ks1 = ac_ksquery_keyset (ks, ">From_gene", h) ;
      ac_keyset_or (ksGenes, ks1) ;
      nEstFromGene = ac_keyset_count (ks1) ;
      ks2 = ac_ksquery_keyset (ks1, "gt_ag || gc_ag", h) ;
      nEstFromGeneSpliced = ac_keyset_count (ks2) ; 
      ac_keyset_or (ksGenesSpliced, ks2) ;
      ac_free (ks2) ;
      ks2 = ac_ksquery_keyset (ks1, ">gene", h) ;

      ks3 = ac_ksquery_keyset (ks2, "pastille_spliced_non_coding", h) ;
      nEstFromGeneSplicedNonCoding = ac_keyset_count (ks3) ; 
      ac_keyset_or (ksGenesSplicedNonCoding, ks3) ;
      ac_free (ks3) ;

      ks3 = ac_ksquery_keyset (ks2, "Pastille_marginally_coding", h) ;
      nEstFromGeneSplicedMargeCoding = ac_keyset_count (ks3) ; 
      ac_keyset_or (ksGenesSplicedMargeCoding, ks3) ;
      ac_free (ks3) ;

      ks3 = ac_ksquery_keyset (ks2, "pastille_coding", h) ;
      nEstFromGeneCoding = ac_keyset_count (ks3) ; 
      ac_keyset_or (ksGenesCoding, ks3) ;
      ac_free (ks3) ;

      ac_free (ks2) ;
      ac_free (ks1) ;
      ac_free (ks) ;
    }

  nTotal = nNm + nMrna + nEst ;
  nTotalAli = ac_keyset_count (ksAli) ;
  nTotalInMrna = ac_keyset_count (ksMrna) ;
  nTotalInMrnaSpliced = ac_keyset_count (ksMrnaSpliced) ;
  nTotalFromGene = ac_keyset_count (ksGenes) ;
  nTotalFromGeneSpliced = ac_keyset_count (ksGenesSpliced) ;
  nTotalFromGeneCoding = ac_keyset_count (ksGenesCoding) ;
  nTotalFromGeneSplicedMargeCoding = ac_keyset_count (ksGenesSplicedMargeCoding) ;
  nTotalFromGeneSplicedNonCoding = ac_keyset_count (ksGenesSplicedNonCoding) ;

  printf ("Alignment efficiency, RefSeq, mRNA the other reads in MrnaInfo, EST the rest\n") ;
  printf ("The first column is the total number of reads considered\n")  ;
  
  printf ("The second column counts the aligned or buried reads of the different categories\n") ;
  printf ("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s"
	  , "Type", "Total", isMulti ? "multi-aligned" : "aligned","%","In # mRNAs", "In # spliced mRNAs"
	  , "in # genes", "in # spliced genes", "in # coding genes", "in # spliced marginally coding genes", "in # spliced non coding genes") ;

  printf ("\n%s\t%d\t%d\t%5.2f%%\t%d\t%d\t%d\t%d\t%d\t%d\t%d"
	  ,"NM", nNm, nNmAli, (100.0 * nNmAli)/nNm
	  , nNmInMrna, nNmInMrnaSpliced
	  , nNmFromGene, nNmFromGeneSpliced, nNmFromGeneCoding, nNmFromGeneSplicedMargeCoding, nNmFromGeneSplicedNonCoding
	  ) ;

  printf ("\n%s\t%d\t%d\t%5.2f%%\t%d\t%d\t%d\t%d\t%d\t%d\t%d"
	  ,"mRNA", nMrna, nMrnaAli, (100.0 * nMrnaAli)/nMrna
	  , nMrnaInMrna, nMrnaInMrnaSpliced
	  , nMrnaFromGene, nMrnaFromGeneSpliced, nMrnaFromGeneCoding, nMrnaFromGeneSplicedMargeCoding, nMrnaFromGeneSplicedNonCoding
	  ) ;

  printf ("\n%s\t%d\t%d\t%5.2f%%\t%d\t%d\t%d\t%d\t%d\t%d\t%d"
	  ,"EST+Trace", nEst, nEstAli, (100.0 * nEstAli)/nEst
	  , nEstInMrna, nEstInMrnaSpliced
	  , nEstFromGene, nEstFromGeneSpliced, nEstFromGeneCoding, nEstFromGeneSplicedMargeCoding, nEstFromGeneSplicedNonCoding
	  ) ;

  printf ("\n%s\t%d\t%d\t%5.2f%%\t%d\t%d\t%d\t%d\t%d\t%d\t%d"
	  ,"Total", nTotal, nTotalAli, (100.0 * nTotalAli)/nTotal
	  , nTotalInMrna, nTotalInMrnaSpliced
	  , nTotalFromGene, nTotalFromGeneSpliced, nTotalFromGeneCoding, nTotalFromGeneSplicedMargeCoding, nTotalFromGeneSplicedNonCoding
	  ) ;
  printf ("\n") ;
  ac_free (h) ;
} /* sAliCount */

static void sBaliseCount (AC_DB db)
{
  AC_KEYSET ks, ks1, ks2 ;
  AC_HANDLE h = handleCreate () ;
  int nBalise,  nBaliseClones, nBaliseEsts, nNms ;

  printf ("\n#############################################################\n") ;
  printf ("## phase 5: sBaliseCount   Needed to evaluate the gene expression level in fichelogic.c:line 1953\n") ;

  /* balise genes */
  if (1)
    {
      ks = ac_dbquery_keyset (db, "find gene balise && transcribed_gene", h) ;
      nBalise = ac_keyset_count (ks) ;

      ks1 = ac_ksquery_keyset (ks, "> transcribed_gene  ; > cdna_clone", h) ;
      nBaliseClones  = ac_keyset_count (ks1) ;

      ks1 = ac_ksquery_keyset (ks, "> transcribed_gene  ; > read", h) ;
      nBaliseEsts  = ac_keyset_count (ks1) ;

      ks2 = ac_ksquery_keyset (ks1, "Ref_seq", h) ;
      nNms = ac_keyset_count (ks2) ;

      printf ("The %d balise genes contain (directly or in their shed variants)\n %d cdna_clones\t%d reads including\t%d NM/NR\n"
	      , nBalise,  nBaliseClones, nBaliseEsts, nNms) ;

    }
  ac_free (h) ;
} /* sBaliseCount */

static void sRefseqContribution (AC_DB db, BOOL isRefSeq)
{
  AC_KEYSET ks, ksTg, ks0, ks1;
  AC_HANDLE h = handleCreate () ;
  int nNm, nNmAli, nNm2Mrna, nNm2Tg, nNm2MrnaSpecific ;
  const char *type ;

  printf ("\n#############################################################\n") ;
  if (isRefSeq)
    {
      ks0 = ac_dbquery_keyset (db, "Find est ref_seq", h) ;
      printf ("## phase 6: RefSeq contribution of the NM/NR\n") ;
      type = "RefSeq (NM/NR)" ;
    }
  else
    {
      ks0 = ac_dbquery_keyset (db, "Find est ref_mrna && !ref_seq", h) ;
      printf ("## phase 7: mRNA contribution of the Genbank mRNAs\n") ;
      type = "mRNA" ;
    }
  /* nms */
  nNm = ac_keyset_count (ks0) ;
  
  ks = ac_ksquery_keyset (ks0,  "from_gene", h) ;
  nNmAli = ac_keyset_count (ks) ;
  ac_free (ks) ;
  
  ks = ac_ksquery_keyset (ks0,  ">in_mrna", h) ;
  nNm2Mrna = ac_keyset_count (ks) ;
  ks1 = ac_ksquery_keyset (ks, "COUNT cdna_clone = 1", h) ;
  nNm2MrnaSpecific = ac_keyset_count (ks1) ;
  ac_free (ks) ;
  
  ksTg = ac_ksquery_keyset (ks0,  ">from_gene", h) ;
  nNm2Tg = ac_keyset_count (ksTg) ;
  
  printf ("There are\t%d\t%s,\t%d aligned, in\t%d mrna in\t%d genes\n"
	  , nNm, type, nNmAli, nNm2Mrna, nNm2Tg) ;
  printf ("\t%d of these mRNAs are supported by a single clone\n"
	  , nNm2MrnaSpecific) ;
  
  /* count gt_ag gc_ag ct_ac supported by NMs */
  {
    AC_OBJ tg = 0, clo = 0, read = 0 ;
    AC_TABLE gIntrons ;
    AC_ITER iter = ac_keyset_iter (ksTg, TRUE, 0) ;
    int ir, x1, x2, y1, y2, dc ;
    int nGt_ag = 0, nGc_ag = 0, nAt_ac = 0, nCt_ac = 0, nOther = 0 ;
    const char *ccq ;
    
    while (ac_free (tg), (tg = ac_iter_obj (iter)))
      {
	gIntrons = ac_tag_table (tg, "Intron_boundaries", 0) ;
	x1 = x2 = y1 = y2 = 0 ;
	for (ir = 0; gIntrons && ir < gIntrons->rows ; ir++)
	  {
	    ccq = ac_table_printable (gIntrons, ir, 0, "") ;
	    if (! strcasecmp (ccq, "other"))
	      dc = 1 ;
	    else
	      dc = 0 ;
	    x1 = ac_table_int (gIntrons, ir, 2 + dc, 0) ;
	    x2 = ac_table_int (gIntrons, ir, 3 + dc, 0) ;
	    if (x1 == y1 && x2 == y2)
	      continue ;
	    clo = ac_table_obj (gIntrons, ir, 4 + dc, h) ; 
	    read = ac_tag_obj (clo, "Read", h) ;
	    if (
		(isRefSeq && ac_has_tag (read, "Ref_seq")) ||
		(!isRefSeq && ! ac_has_tag (read, "Ref_seq") && ac_has_tag (read, "Ref_mrna")) 
		)
	      {
		y1 = x1 ; y2 = x2 ;
		ccq = ac_table_printable (gIntrons, ir, 0, "") ;
		if (! strcasecmp (ccq, "gt_ag"))
		  nGt_ag++ ;
		if (! strcasecmp (ccq, "gc_ag"))
		  nGc_ag++ ;
		if (! strcasecmp (ccq, "at_ac"))
		  nAt_ac++ ;
		if (! strcasecmp (ccq, "ct_ac"))
		  nCt_ac++ ;
		if (! strcasecmp (ccq, "other"))
		  nOther++ ;
	      }
	    ac_free (clo) ;
	    ac_free (read) ;
	    }
	ac_free (gIntrons) ;
      }
    printf ("Introns supported by %s\n", type) ;
    printf ("%8d gt-ag\n", nGt_ag) ;
    printf ("%8d gc-ag\n", nGc_ag) ;
    printf ("%8d at-ac\n", nAt_ac) ;
    printf ("%8d ct-ac\n", nCt_ac) ;
    printf ("%8d other\n", nOther) ;
    printf ("%8d total\n", nGt_ag + nGc_ag + nAt_ac + nCt_ac + nOther) ;
  }
  ac_free (ksTg) ;

  ac_free (h) ;
} /* sBaliseCount */

static void sCountAllIntrons (AC_DB db)
{
  AC_KEYSET ks0, ksTg ;
  AC_HANDLE h = handleCreate () ;
  int nGene, nMrna, nEst ;

  printf ("\n#############################################################\n") ;
  printf ("## phase 8: Counting introns\n") ;
  ks0 = ac_dbquery_keyset (db, "Find gene transcribed_gene", h) ;
  nGene = ac_keyset_count (ks0) ;
  ks0 = ac_dbquery_keyset (db, "Find mrna cdna_clone", h) ;
  nMrna = ac_keyset_count (ks0) ;
  ks0 = ac_dbquery_keyset (db, "Find est from_gene", h) ;
  nEst = ac_keyset_count (ks0) ;

  printf ("There are\t%d\tsequences actually aligned and not buried, in\t%d\tmrna in %d genes\n"
	  , nEst, nMrna, nGene) ;
  ks0 = ac_dbquery_keyset (db, "Find transcribed_gene gt_ag || gc_ag", h) ;
  nGene = ac_keyset_count (ks0) ;
  ks0 = ac_dbquery_keyset (db, "Find mrna gt_ag || gc_ag", h) ;
  nMrna = ac_keyset_count (ks0) ;
  printf ("There are\t%d\ttranscipts with confirmed introns in\t%d\tgenes\n"
	  , nMrna, nGene) ;
  ksTg = ac_dbquery_keyset (db, "Find tg Intron_boundaries", h) ;
  
  /* count gt_ag gc_ag ct_ac */
  {
    AC_OBJ tg = 0 ;
    AC_TABLE gIntrons ;
    AC_ITER iter = ac_keyset_iter (ksTg, TRUE, 0) ;
    int ir, x1, x2, y1, y2, dc ;
    int nGt_ag = 0, nGc_ag = 0, nAt_ac = 0, nCt_ac = 0, nOther = 0 ;
    const char *ccq ;
    
    while (ac_free (tg), (tg = ac_iter_obj (iter)))
      {
	gIntrons = ac_tag_table (tg, "Intron_boundaries", 0) ;
	x1 = x2 = y1 = y2 = 0 ;
	for (ir = 0; gIntrons && ir < gIntrons->rows ; ir++)
	  {
	    ccq = ac_table_printable (gIntrons, ir, 0, "") ;
	    if (! strcasecmp (ccq, "other"))
	      dc = 1 ;
	    else
	      dc = 0 ;
	    x1 = ac_table_int (gIntrons, ir, 2 + dc, 0) ;
	    x2 = ac_table_int (gIntrons, ir, 3 + dc, 0) ;
	    if (x1 == y1 && x2 == y2)
	      continue ;
	    if (1)
	      {
		y1 = x1 ; y2 = x2 ;
		ccq = ac_table_printable (gIntrons, ir, 0, "") ;
		if (! strcasecmp (ccq, "gt_ag"))
		  nGt_ag++ ;
		if (! strcasecmp (ccq, "gc_ag"))
		  nGc_ag++ ;
		if (! strcasecmp (ccq, "at_ac"))
		  nAt_ac++ ;
		if (! strcasecmp (ccq, "ct_ac"))
		  nCt_ac++ ;
		if (! strcasecmp (ccq, "other"))
		  nOther++ ;
	      }
	    }
	ac_free (gIntrons) ;
      }
    printf ("Introns supported by EST mRNA or NM\n") ;
    printf ("%8d\tgt-ag\n", nGt_ag) ;
    printf ("%8d\tgc-ag\n", nGc_ag) ;
    printf ("%8d\tat-ac\n", nAt_ac) ;
    printf ("%8d\tct-ac\n", nCt_ac) ;
    printf ("%8d\tother\n", nOther) ;
    printf ("%8d\ttotal\n", nGt_ag + nGc_ag + nAt_ac + nCt_ac + nOther) ;
  }
  ac_free (ksTg) ;

  ac_free (h) ;
} /* sCountAllIntrons */

static void sAli (const char *ici, int phase)
{
  AC_DB db ;
  const char *s = "ok" ;

  db = ac_open_db (messprintf("%s",ici), &s);
  if (!db)
    messcrash ("Failed to open db %s, error %s", ici, s) ;
  printf ("\n#############################################################\n") ;
  printf ("## phase %d   %s   %s \n", phase, ici, timeShowNow()) ;

  switch (phase)
    {
    case 3: /* count ali */
      sAliCount (db, FALSE) ;
      break ;
    case 4: /* count multi ali */
      sAliCount (db, TRUE) ;
      break ;
    case 5: /* count balise */
      sBaliseCount (db) ;
      break ;
    case 6: /* count refseq contrib */
      sRefseqContribution (db, TRUE) ;
      break ;
    case 7: /* count refseq contrib */
      sRefseqContribution (db, FALSE) ;
      break ;
    case 8: /* count all introns */
      sCountAllIntrons (db) ;
      break ;
    case 10: /* count tg */
      sTgCount (db) ;
      sGeneIntronsCount (db, TRUE) ;
      sGeneIntronsCount (db, FALSE) ;
      sGeneMrnaIntronsCount (db, TRUE) ;
      sGeneMrnaIntronsCount (db, FALSE) ;
      break ;
    case 20:
      sEstQuality (db) ;
      break ;
    case 21: /* aligned bp and errors added feb 25, may be same as above ? */
      sAliStats (db) ; 
      break ;
    case 25: /* count mrna with introns sorted by coding length + pfam/psort */
      sTgXIntronsCount (db, TRUE, FALSE) ;
      break ;
    case 26: /* count mrna with introns sorted by coding length + pfam/psort */
      sTgXIntronsCount (db, TRUE, TRUE) ;
      break ;
    case 27: /* count mrna without introns sorted by coding length + pfam/psort */
      sTgXIntronsCount (db, FALSE, FALSE) ;
      break ;
    case 31: /* tg non coding introns histogram, exports an ace file  */
      sTgNonCodingIntrons (db, TRUE, FALSE) ;
      break ;
    case 32: /* idem just variant a */
      sTgNonCodingIntrons (db, TRUE, TRUE) ;
      break ;
    case 33: /* count multi proteines , aceExport the tag */
      sTgMultiProteins (db, TRUE) ;
      break ;
    case 34: /* analyse cloud of small genes */
      sGeneCloud (db, ici) ;
      break ; 
    case 40: /* multi count de la frequence d'alternatifs */
      sCountAlternatifs (db) ;
      break ;
    case 50: /* cloud genes */
      sCountCloudGenes (db) ;
      break ;
    case 52: /* alter stats genes, feb 25 */
      sAlterStats (db) ;
      break ;
    case 53: /* alter stats genes, feb 25 */
      sLocusIdStats (db) ;
      break ;
    case 54: /* gene stats genes, nov 2007 */
      sGeneIdStats (db) ;
      break ;
    case 99:
      sAliCount (db, FALSE) ;
      sAliCount (db, TRUE) ; 
      if (0) sAliStats (db) ;
      if (0) sEstQuality (db) ;
      sBaliseCount (db) ;
      sRefseqContribution (db, TRUE) ;
      sRefseqContribution (db, FALSE) ; 
      sTgCount (db) ;
      sGeneIntronsCount (db, TRUE) ;
      sGeneIntronsCount (db, FALSE) ;
      sTgXIntronsCount (db, TRUE, FALSE) ;
      sTgXIntronsCount (db, TRUE, TRUE) ;
      sTgXIntronsCount (db, FALSE, FALSE) ;
      sTgNonCodingIntrons (db, FALSE, FALSE) ; /* show stats do NOT export the ace file */
      sTgNonCodingIntrons (db, FALSE, TRUE) ; /* show stats do NOT export the ace file */
      sAlterStats (db) ;
      sCountAlternatifs (db) ;
      sTgMultiProteins (db, FALSE) ;
      sCountCloudGenes (db) ;
      sLocusIdStats (db) ;
      sGeneIdStats (db) ;
      if (0) 
	{
	  sGeneCloud (db, ici) ;
	  sAliStats (db) ;
	  sEstQuality (db) ;
	}
      break ;      
    default:
      messcrash ("sAli: bad phase %d\n", phase) ;
    }
      
  ac_db_close (db) ;
} /* sAli */

static void sGlobal (const char *ici, int phase)
{
  AC_DB db ;
  AC_KEYSET ks ;
  const char *s = "ok" ;
  int nTotal=0, nNm=0, nMrna=0, nEst = 0 ;
  int nTotalAli=0, nNmAli=0, nMrnaAli=0, nEstAli = 0 ;
  AC_HANDLE h = handleCreate () ;

  printf ("phase %d\n", phase) ;
  switch (phase)
    {
    case 1:
      db = ac_open_db (messprintf("%s/%s",ici,"MrnaInfo"), &s);
      if (!db)
	messcrash ("Failed to open db %s/MrnaInfo, error %s", ici, s) ;
      ks = ac_dbquery_keyset (db, "find est (Ref_Seq) && (is_mrna || ref_mrna)", h) ;
      nNm = ac_keyset_count (ks) ;
      ac_free (ks) ;
      ks = ac_dbquery_keyset (db, "find est (! Ref_Seq) && (is_mrna || ref_mrna)", h) ;
      nMrna = ac_keyset_count (ks) ;
      ac_free (ks) ;
      ac_db_close (db) ;
      break ;

    case 2: /* count EST */
      db = ac_open_db (messprintf("%s/%s",ici,"EstData"), &s);
      if (!db)
	messcrash ("Failed to open db %s/EstData, error %s", ici, s) ;
      ks = ac_dbquery_keyset (db, "find est !(ref_seq) && !(is_mrna || ref_mrna)", h) ;
      nEst = ac_keyset_count (ks) ;
      ac_free (ks) ;
      ac_db_close (db) ;
      break ;

    default:
      messcrash ("sGlobal: bad phase %d\n", phase) ;
    }
      
  nTotal = nNm + nMrna + nEst ;
  nTotalAli = nNmAli + nMrnaAli + nEstAli ;
  printf ("%8s\t%9d\t%9d\n","NM", nNm, nNmAli) ;
  printf ("%8s\t%9d\t%9d\n","mRNA", nMrna, nMrnaAli) ;
  printf ("%8s\t%9d\t%9d\n","EST", nEst, nEstAli) ;
  printf ("%8s\t%9d\t%9d\n","Total", nTotal, nTotalAli) ;
  ac_free (h) ;
} /* sGlobal */

static void usage (void)
{
  printf ("// Usage :blystats -db $ACEDB -p phase [-worm] \n") ;
  printf ("//                 =  1:global count NM mRNA\n") ;
  printf ("//                 =  2:global count EST\n") ;
  printf ("//                 =  3:count aligned NM, mRNA, EST\n") ;
  printf ("//                 =  4:count multi aligned NM, mRNA, EST\n") ;
  printf ("//                 =  5:count (!worm) balise genes and their clones\n") ;
  printf ("//                 =  6:Refseq Contribution\n") ;
  printf ("//                 =  7:Mrna Contribution\n") ;
  printf ("//                 =  8:Count all introns\n") ;
  printf ("//                 = 10:tg/mRNA: count the genes and their alternative forms\n") ;
  printf ("//                 = 20:table of quality of all alignments\n") ;
  printf ("//                 = 21:table of quality of all alignments, feb 25\n") ;
  printf ("//                 = 25: count mrna with introns sorted by coding length + pfam/psort\n") ;
  printf ("//                 = 26: count antisens gene with introns sorted by coding length + pfam/psort\n") ;
  printf ("//                 = 27: count mrna without intron sorted by coding length + pfam/psort\n") ;
  printf ("//                 = 31: ace file and histogram of genes witn nn non coding introns\n") ;
  printf ("//                 = 32: idem just variant a\n") ;
  printf ("//                 = 33:count multi proteins\n") ;
  printf ("//                 = 34:analyse cloud of small genes\n") ;
  printf ("//                 = 40:multi count de la frequence d'alternatifs\n") ;
  printf ("//                 = 50:cloud genes feb 25\n") ;
  printf ("//                 = 52:alternatives feb 25\n") ;
  printf ("//                 = 53:locusId<->gene stats aug 30 2005\n") ;
  printf ("//                 = 54:geneId<->gene stats nov 16 2007\n") ;
  printf ("//                 = 99:all\n") ;
  
  exit (1) ;
}


int main (int argc, const char **argv)
{
  int nn = 0 ;
  const char *ici, *phase ;

  freeinit () ;
  freeOutInit () ;

  ici = getArgV (&argc, argv, "-db") ;
  phase = getArgV (&argc, argv, "-p") ;
  isWorm = getArg (&argc, argv, "-worm") ;
  
  if (! ici)
    usage () ;
  if (phase) 
    sscanf (phase,"%d", &nn) ;
  switch (nn)
    {
    case 1: case 2: /* data gatherings */
      sGlobal (ici, nn) ;
      break ;
    case  3: /* count ali */
    case  4: /* count multi ali */
    case  5: /* count balise */
    case  6: /* RefSeq contrib */
    case  7: /* MnraContrib */
    case  8: /* Count introns */
    case 10: /* count tg mrna */
    case 20: /* count est quality */
    case 21: /* aligned bp and errors added feb 25, may be same as above ? */
    case 25: /* count mrna with introns sorted by coding length + pfam/psort */
    case 26: /* count antisens genes with introns sorted by coding length + pfam/psort */
    case 27: /* count mrna without introns sorted by coding length + pfam/psort */ 
    case 31: /* tg non coding introns histogram */
    case 32: /* idem, just variant a */
    case 33: /* count multi proteines */
    case 34: /* analyse cloud of small genes */
    case 40: /* multi count de la frequence d'alternatifs */
    case 52: /* alter stats genes, feb 25 */
    case 50: /* cloud genes */
    case 53: /* locusId <-> genes */
    case 54: /* locusId <-> genes */
    case 99: /* all */
      sAli (ici, nn) ;
      break ;
    default: 
      fprintf (stderr, "bad phase %d\n", nn) ;
      usage () ;
      break ;
    }

  sshowHits(0) ;/* for compiker happiness */
  printf ("// done\n") ;
  return 0 ;
}

