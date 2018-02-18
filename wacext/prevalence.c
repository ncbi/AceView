#include "wh/ac.h"
/*
 * Authors: Danielle and Jean Thierry-Mieg, NCBI, 
 * February 2014
 * prevalence.c, this code is part of the NCBI Magic RNA-seq pipeline
 * It his open access, with no warranty whatsoever
 * Aim
 *   In how many individual are these genes expressed for the union
 * Input: a concatenated set of tab delimited files
 * selected genes/mrnas ZZZZZ info about mrna<->gene ZZZZZ info about gene-types ZZZZZ index->expressed ZZZZZ read-count->touched
 */

/* #define ARRAY_CHECK  */

#include "regular.h"
#include <ac.h>
#include <aceio.h>
#include <dict.h>
#include <vtxt.h>
#include <math.h>

######################################

typdef prevalenceStruct {
  DICT *selectDict ;
  DICT *geneDict ;
  DICT *mrnaDict ; 
  Associator m2g ;
  const char *selectFileName, *mrnaMapFileName
} PRV ;

static void prevalenceGetSelection (PRV *prv)
{
  AC_HANDLE h = ac_new_handle ()
  ACEIN ai = aceInCreate (prv->selectFileName, FALSE, h) ;
  DICT *dict = prv->selectDict ;
  while (aceInCard (ai))
    {
      ccp = aceInWord () ; if (!ccp || *ccp == '#') continue ;
      dictAdd (dict, ccp, 0) ;
    }
  ac_free (h) ;
}

######################################

static void parseMrna_map_ln_gc_gene_geneid (PRV *prv, int wantGeneId)
{
  AC_HANDLE h = ac_new_handle ()
  ACEIN ai = aceInCreate (prv->mrnaMapFileName, FALSE, h) ;
  DICT *dict = prv->selectDict ;
  DICT *geneDict = prv->geneDict ;
  DICT *mrnaDict = prv->mrnaDict ;
  int gm, gene, mrna ;
  Associator m2g = prv->m2g ;
  char geneBuf[1024], mrnaBuf[1024] ;

  while (aceInCard (ai))
    {
      mrna = gene = 0 ;
      ccp = aceInWord () ; if (!ccp || *ccp == '#') continue ;
      strncpy (mrnaBuf, ccp, 1000) ; 


      aceInStep ('\t') ; aceInWord () ; /* map, ignore */
      aceInStep ('\t') ; aceInWord () ; /* mrna length, ignore */
      aceInStep ('\t') ; aceInWord () ; /* % GC,  ignore */
      
      aceInStep ('\t') ; aceInWord () ; if (!ccp || *ccp == '#') continue ;
      strncpy (geneBuf, ccp, 1000) ; 
      aceInStep ('\t') ; aceInWord () ;  /* geneId */

      /* register the gene mrna correspondence */
      dictAdd (mrnaDict, mrnaBuf, &mrna) ;
      dictAdd (geneDict, geneBuf, &gene) ; 
      assInsert (m2g, assVoid(m), assVoid(g)) ;

      /* do we have a GeneId as desired */
      switch (wantGeneId)
	{
	case 0: if (ccp) continue ; break ; /* please NO geneId */
	case 1: break ; /* doest not matter */
	case 2: if (!ccp) continue ; break ; /* please i NEED a geneId */
	}

      /* is this gene in the selection */
      if (! dictFind (dict, geneBuf, 0) && ! dictFind (dict, mrnaBuf, 0))
	continue ;

      keySet (geneok, gene) = 1 ;

      gm = gene 
      if (GM == _GENE || GM2 == _m2g)
	gm = gene ;
      else
	gm = mrna ;
      
      keySet (ggok, gm) = 1 ;
      keySet (gg, gm) = -1 ;
    }
  ac_free (h) ;
}

######################################




#ifdef JUNK
/^ZZZZZ/{
  zz++;
  next;
}
{if(zz<1) {
  selected_gene[$1] = 1 ;
  next ;
}}
######################################

/^#Run/{nnr++;if(nnr==1){for (i=1;i<=30 && i<=NF;i++){if($i=="_SumOfAllReadsInProject")imin=i+1;}for(i=imin;i<=NF;i++)if(substr($i,1,3)!="Rhs"){imax=i;next}imax=NF;}}
/^#/{next}

{if(zz<1) {
  selected_gene[$1] = 1 ;
  next ;
}}

{if(zz==1) {
  m = $1 ; g = $5 ; 
  m2g[m] = g ;  
  if (selected_gene[g] + selected_gene[m] < 1) next ;
  if((hasGid=="with" && $6=="")||(hasGid=="no" && $6)) next;
  if(length(g)>=1) 
    {  
      geneok[g] = 1 ;
      if (GM == "GENE" || GM2 == "m2g")
	{
	  ggok[g] = 1 ; gg[g] = -1 ; 
	}
      else
	{ ggok[m] = 1 ; gg[m] = -1 ; }
    }
  next ;
}}
{if(zz==2) {
  g=$1 ;
  if (geneok[g] == 1)
    {
      g2t[g]=$4 ; tt[$4]=1;
    }
  next;
}}


/AceView/{ if (index($2,"AceView")<1) next ;
  g=$1;
  if (GM2 == "m2g") 
    {
      m = g ;
      g = m2g[m] ;
      if (length(g)<1) printf("m=%s no Gene\n", m) ;
    }
  if (0+ggok[g] < 1) 
    {
      # printf("%s\t%s\trejected\n", g, m) ;
      next ;
    }

  n = 0 ;
  for (i = imin ; i <= NF && i < imax ; i++)
    {
      if ($i + 0 > 0)
	n++ ;
    }
  # printf("%s\t%s\taccepted\t%d\t%d\t%d\t%d\n", g, m, zz, n, imin, imax) ;

  if (n > nmax) nmax = n ;
  if (zz == 3 && n > 0+gg[g]) gg[g] = n ; 
  if (zz == 4 && n > 0 + touched[g]) touched[g] = n ; 

  next ;
}
END {
  for(g in gg)
    for(n=-1;n<=nmax;n++)
      {
	if (GM == "GENE" || GM2 == "m2g") 
	  t=g2t[g];
	else 
	  t = g2t[m2g[g]] ;
	if (length(t)<1) 
	  {
	    t="no_type" ; 
	    if(0)printf("ZZZ %s\n",g) ;
	  }

	if (n == -1 || (n == 0 && 0+touched[g]>0) || (n > 0 && 0+gg[g]>=n))
	  { 
	    tt[t]=1;t2n[t,n]++;ttn[n]++;
	  }
      }
  
  printf("A\tType\t%s\tTouched",title);
  for (n=1;n<=nmax;n++)printf("\t%d",n) ;
  printf("\t\tType\t%s\tTouched",title);
  for (n=1;n<=nmax;n++)printf("\t%.1f",100.0*n/(nmax+.000001)) ;
  for(t in tt)
    {
      if(0+t2n[t,0]<1)continue;
      printf("\nC\t%s",t);
      for(n=-1;n<=nmax;n++)
	printf("\t%d",t2n[t,n]);
      printf ("\t\t%s",t) ;
      for(n=-1;n<=nmax;n++)
	{
	  z=100*t2n[t,n]/(.0001+t2n[t,-1]) ;
	  printf("\t%.2f",z);
	}
    }
  printf("\nT\tAny type %s %s", GM2, title);
  for(n=-1;n<=nmax;n++)
    printf("\t%d",ttn[n]);
  printf ("\t\tAny type %s %s", GM2, title) ;
  for(n=-1;n<=nmax;n++)
    {
       z=100*ttn[n]/(.0001+ttn[-1]) ;
        printf("\t%.2f",z);
    }
  printf("\n");
  printf("a\nb\tType\tAnnotated\tTouched");
  for(n=1;n<=nmax;n++)
    printf("\t%d",n);
  for(t in tt)
    {
      printf("\np\t%% %s",t);
      for(n=-1;n<=nmax;n++)
	printf("\t%.1f",100*t2n[t,n]/(.0001+ttn[n]));
    }
  printf("\n");
}
#endif
