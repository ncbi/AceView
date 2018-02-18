/*
 Exportation des statistiques acembly pour la conference japon 2003
*/
/* 
   beautiful gawk code
   gawk '/^Variant/{aa[$2]++;}
        END{for (v in aa) printf("%s %d\n",v,aa[v]);}'
      variant.ace | wc


anatomy of a gene
tbly . << EOF
  query find gene structure
  show -a -f gene_structure.ace structure
  quit
EOF

#count the mrna
cat ZH* /d14.structure.*.ace | 
gawk '/^Gene /{ngene++;next;}/^$/{next;}{aa[$1]+=$2;}END{printf("using %d genes\n", ngene); for (v in aa) printf("%s %d\n",v,aa[v]);}'

#count by cumul
cat ZH* /d14.structure.*.ace | 
gawk '/^Gene /{ngene++;next;}/^$/{next;}{aa[$1]+=1;}END{printf("using %d genes\n", ngene); for (v in aa) printf("%s %d\n",v,aa[v]);}' gene_structure.ace
#count by cumul truque pour les promoteurs
gawk '/^Gene /{ngene++;next;}/^$/{next;}{if($2>1)aa[$1]+=1;}END{printf("using %d genes\n", ngene); for (v in aa) printf("%s %d\n",v,aa[v]);}' 
*/

#include "../wac/ac.h"
#include <errno.h>
#include "dict.h"
#include "freeout.h"
#include <math.h>

/*
static double mean (AC_ITER iter, char *tag)
{
  int nn = 0 ;
  double cumul = 0 ;
  AC_OBJ (obj) ;
  
  ac_iter_rewind (iter) ;

  while ((obj = ac_next_obj (iter)))
    {
      nn++ ;
      cumul += ac_tag_int (obj, tag, 0) ;
    }
  return nn ? cumul/nn : 0 ;
}
*/
/**********************************************************************/

typedef struct llStruct { int nP ; AC_DB db ; Array a ; AC_KEYSET ks ; AC_ITER iter ; AC_HANDLE h ; } LL ;

typedef struct Pstruct {
  BOOL man, woman, placebo, knee, hip, aceto ;
  int age, stratum, duration ;
  int pill_count ;
  float function, pain, height, weight, bmi, global, pills_per_day ;
  float auc_pain, auc_function, auc_global ;
  float pain2, global2, function2, pain4, global4, function4, pain6, global6, function6 ;
} P ;

static Array getAllData (LL *ll)
{
  int nn = 0, ir ;
  Array a = arrayHandleCreate (1000, P, ll->h) ;
  AC_OBJ (obj) ;
  AC_TABLE tbl ;
  P* pp ;
  
  ac_iter_rewind (ll->iter) ;

  while ((obj = ac_next_obj (ll->iter)))
    {
      pp = arrayp (a, nn++, P) ;
      pp->man = ac_has_tag (obj, "Man") ;
      pp->woman = ac_has_tag (obj, "Woman") ;
      pp->age = ac_tag_int (obj, "Age", -1) ;
      pp->stratum = ac_tag_int (obj, "Stratum", -1) ;
      pp->duration = ac_tag_int (obj, "Duration", -1) ;
      pp->placebo = ac_has_tag (obj, "Placebo") ;
      pp->knee  = ac_has_tag (obj, "Knee") ;
      pp->hip  = ac_has_tag (obj, "Hip") ;
      pp->global = ac_tag_int (obj, "Global", -1) ;
      pp->function = ac_tag_float (obj, "Function", -1) ;
      pp->pain = ac_tag_float (obj, "Pain", -1) ;
      pp->auc_global = ac_tag_float (obj, "auc_Global", -1) ;
      pp->auc_function = ac_tag_float (obj, "auc_Function", -1) ;
      pp->auc_pain = ac_tag_float (obj, "auc_Pain", -1) ;
      pp->height = ac_tag_float (obj, "Height", -1) ;
      pp->weight = ac_tag_float (obj, "Weight", -1) ;
      pp->bmi = ac_tag_float (obj, "BMI", -1) ;
      pp->aceto = ac_has_tag (obj, "Acetaminophen_user") ;
      pp->pill_count = ac_tag_int (obj, "Pill_count", -1) ;
      pp->pills_per_day = ac_tag_float (obj, "pills_per_day", -1) ;
      pp->pain2 = pp->pain4 = pp->pain6 = -9999 ;
      tbl = ac_tag_table (obj, "P_week", 0) ;
      for (ir = 0 ; tbl && ir < tbl->rows ; ir++)
	{
	  if (ac_table_int (tbl, ir, 0, -1) == 2)
	    pp->pain2 = ac_table_float (tbl, ir, 1, -9999) ;
	  if (ac_table_int (tbl, ir, 0, -1) == 4)
	    pp->pain4 = ac_table_float (tbl, ir, 1, -9999) ;
	  if (ac_table_int (tbl, ir, 0, -1) == 6)
	    pp->pain6 = ac_table_float (tbl, ir, 1, -9999) ;
	}
      ac_free (tbl) ;
      pp->function2 = pp->function4 = pp->function6 = -9999 ;
      tbl = ac_tag_table (obj, "F_week", 0) ;
      for (ir = 0 ; tbl && ir < tbl->rows ; ir++)
	{
	  if (ac_table_int (tbl, ir, 0, -1) == 2)
	    pp->function2 = ac_table_float (tbl, ir, 1, -9999) ;
	  if (ac_table_int (tbl, ir, 0, -1) == 4)
	    pp->function4 = ac_table_float (tbl, ir, 1, -9999) ;
	  if (ac_table_int (tbl, ir, 0, -1) == 6)
	    pp->function6 = ac_table_float (tbl, ir, 1, -9999) ;
	}
      ac_free (tbl) ;
      pp->global2 = pp->global4 = pp->global6 = -9999 ;
      tbl = ac_tag_table (obj, "G_week", 0) ;
      for (ir = 0 ; tbl && ir < tbl->rows ; ir++)
	{
	  if (ac_table_int (tbl, ir, 0, -1) == 2)
	    pp->global2 = ac_table_float (tbl, ir, 1, -9999) ;
	  if (ac_table_int (tbl, ir, 0, -1) == 4)
	    pp->global4 = ac_table_float (tbl, ir, 1, -9999) ;
	  if (ac_table_int (tbl, ir, 0, -1) == 6)
	    pp->global6 = ac_table_float (tbl, ir, 1, -9999) ;
	}
      ac_free (tbl) ;
      ac_free (obj) ;
    }
  return a ;
}

/**********************************************************************/
/**********************************************************************/
 
static void showAge (LL *ll) 
{
  int ip, ii, nn, nna ;
  double x, x2 ;
  double eps, m[2], s[2] ;
  P *pp ;

  for (ip = 1 ; ip >= 0 ; ip-- )
    {
      if (ip == 1)
	printf ("Placebo:\n") ;
      else
	printf ("Active:\n") ;
      x = x2 = 0 ;
      nn = nna = 0 ;
      for (ii = 0, pp = arrp (ll->a, ii, P) ; ii < arrayMax (ll->a) ; pp++, ii++)
	{
	  if ((ip && !pp->placebo) || (!ip && pp->placebo))
	      continue ;
	  nn++ ;
	  if (pp->age >= 0)
	    { x += pp->age ; x2 += pp->age * pp->age ; nna++ ; }
	}
      
      printf ("Found %d patients\n", nn) ;
      
      m[ip] = x/nna ; s[ip] = (x2 - x * x/nna)/(nna*(nna-1)) ; 
      printf ("\tWe recorded the age of %d,", nna) ;
      printf (" average %.2f years,", (float)(x/nna)) ;
      printf (" standard deviation %.2f years\n", (float)sqrt(x2 - (x*x)/nna)/(nna - 1)) ;
    }
  eps = (m[1] - m[0])/sqrt (s[0] + s[1]) ;
  if (eps < 0) eps = - eps ;
  if (eps < 1.96) printf ("The 2 observed means are compatible : eps = %.2f < 1.96\n", eps) ;
  else printf ("The 2 observed means are NOT compatible : eps = %.2f > 1.96\n", eps) ;

} /* showAge */

/**********************************************************************/
 
static void showSex (LL *ll) 
{
  int ip, ii, nn = 0, nMan = 0, nWoman = 0 ;
  int nm[2], nw[2] ;
  double chi2 ;
  P *pp ;

  for (ip = 1 ; ip >= 0 ; ip-- )
    {
      if (ip == 1)
	printf ("Placebo:\n") ;
      else
	printf ("Active:\n") ;

      nn = nMan = nWoman = 0 ;
      for (ii = 0, pp = arrp (ll->a, ii, P) ; ii < arrayMax (ll->a) ; pp++, ii++)
	{
	  if ((ip && !pp->placebo) || (!ip && pp->placebo))
	      continue ;
	  nn++ ;
	  if (pp->man) nMan++ ;
	  else if (pp->woman) nWoman++ ;
	}
      
      printf ("Found %d patients\n", nn) ;
      
      printf ("\t%d (%.2f%%) are men\n", nMan, (100.0*(float)nMan)/nn) ;
      printf ("\t%d (%.2f%%) are women\n", nWoman, (100.0*(float)(nWoman))/nn) ;
      printf ("\t%d (%.2f%%) are androgyne\n", nn - nMan - nWoman, (100.0*(float)(nn - nMan - nWoman))/nn) ;
      nm[ip] = nMan ;
      nw[ip] = nWoman ;
    }

  {
    double tm[2], tw[2], l1, l2, c1, c2, tot ;
    l1 = nm[0] + nm[1] ;
    l2 = nw[0] + nw[1] ;
    c1 = nm[0] + nw[0] ;
    c2 = nm[1] + nw[1] ;
    tot = l1 + l2 ;
    tm[0] = l1 * c1 / tot ; 
    tm[1] = l1 * c2 / tot ; 
    tw[0] = l2 * c1 / tot ; 
    tw[1] = l2 * c2 / tot ; 
    chi2 = (nm[0]-tm[0])*(nm[0]-tm[0])/tm[0] 
      + (nm[1]-tm[1])*(nm[1]-tm[1])/tm[1] 
      + (nw[0]-tw[0])*(nw[0]-tw[0])/tw[0] 
      + (nw[1]-tw[1])*(nw[1]-tw[1])/tw[1] ;

    printf ("chi2 = %.2f\n", (float)chi2) ;
  }
  } /* showSex */

/**********************************************************************/
 
static void showBmi (LL *ll) 
{
  int ip, ii, nn, nna ;
  double x, x2 ;
  double eps, m[2], s[2] ;
  P *pp ;

  for (ip = 1 ; ip >= 0 ; ip-- )
    {
      if (ip == 1)
	printf ("Placebo:\n") ;
      else
	printf ("Active:\n") ;
      x = x2 = 0 ;
      nn = nna = 0 ;
      for (ii = 0, pp = arrp (ll->a, ii, P) ; ii < arrayMax (ll->a) ; pp++, ii++)
	{
	  if ((ip && !pp->placebo) || (!ip && pp->placebo))
	      continue ;
	  nn++ ;
	  if (pp->bmi >= 0)
	    { x += pp->bmi ; x2 += pp->bmi * pp->bmi ; nna++ ; }
	}
      
      printf ("Found %d patients\n", nn) ;
      
      m[ip] = x/nna ; s[ip] = (x2 - x * x/nna)/(nna*(nna-1)) ; 
      printf ("\tWe recorded the bmi of %d,", nna) ;
      printf (" average bmi %.2f kg/m2,", (float)(x/nna)) ;
      printf (" standard deviation %.2f kg/m2\n", (float)sqrt(x2 - (x*x)/nna)/(nna - 1)) ;
    }
  eps = (m[1] - m[0])/sqrt (s[0] + s[1]) ;
  if (eps < 0) eps = - eps ;
  if (eps < 1.96) printf ("The 2 observed means are compatible : eps = %.2f < 1.96\n", eps) ;
  else printf ("The 2 observed means are NOT compatible : eps = %.2f > 1.96\n", eps) ;

} /* showBmi */

/**********************************************************************/
 
static void showDuration (LL *ll) 
{
  int ip, ii, nn, nna ;
  double x, x2 ;
  double eps, m[2], s[2] ;
  P *pp ;

  for (ip = 1 ; ip >= 0 ; ip-- )
    {
      if (ip == 1)
	printf ("Placebo:\n") ;
      else
	printf ("Active:\n") ;
      x = x2 = 0 ;
      nn = nna = 0 ;
      for (ii = 0, pp = arrp (ll->a, ii, P) ; ii < arrayMax (ll->a) ; pp++, ii++)
	{
	  if ((ip && !pp->placebo) || (!ip && pp->placebo))
	      continue ;
	  nn++ ;
	  if (pp->duration >= 0)
	    { x += pp->duration ; x2 += pp->duration * pp->duration ; nna++ ; }
	}
      
      printf ("Found %d patients\n", nn) ;
      
      m[ip] = x/nna ; s[ip] = (x2 - x * x/nna)/(nna*(nna-1)) ; 
      printf ("\tWe recorded the duration of %d,", nna) ;
      printf (" average duration %.2f years,", (float)(x/nna)) ;
      printf (" standard deviation %.2f years\n", (float)sqrt(x2 - (x*x)/nna)/(nna - 1)) ;
    }
  eps = (m[1] - m[0])/sqrt (s[0] + s[1]) ;
  if (eps < 0) eps = - eps ;
  if (eps < 1.96) printf ("The 2 observed means are compatible : eps = %.2f < 1.96\n", eps) ;
  else printf ("The 2 observed means are NOT compatible : eps = %.2f > 1.96\n", eps) ;

} /* showDuration */

/**********************************************************************/
 
static void showJoint (LL *ll) 
{
  int ip, ii, nn = 0, nKnee = 0, nHip = 0 ;
  P *pp ;

  for (ip = 1 ; ip >= 0 ; ip-- )
    {
      if (ip == 1)
	printf ("Placebo:\n") ;
      else
	printf ("Active:\n") ;

      nn = nKnee = nHip = 0 ;
      for (ii = 0, pp = arrp (ll->a, ii, P) ; ii < arrayMax (ll->a) ; pp++, ii++)
	{
	  if ((ip && !pp->placebo) || (!ip && pp->placebo))
	      continue ;
	  nn++ ;
	  if (pp->knee) nKnee++ ;
	  else if (pp->hip) nHip++ ;
	}
      
      printf ("Found %d patients\n", nn) ;
      
      printf ("\t%d (%.2f%%) are knee\n", nKnee, (100.0*(float)nKnee)/nn) ;
      printf ("\t%d (%.2f%%) are hip\n", nHip, (100.0*(float)(nHip))/nn) ;
      printf ("\t%d (%.2f%%) are unknown\n", nn - nKnee - nHip, (100.0*(float)(nn - nKnee - nHip))/nn) ;
    }

} /* showJoint */

/**********************************************************************/
 
static void showPain (LL *ll) 
{
  int ip, ii, nn, nna ;
  double x, x2 ;
  double eps, m[2], s[2] ;
  P *pp ;

  for (ip = 1 ; ip >= 0 ; ip-- )
    {
      if (ip == 1)
	printf ("Placebo:\n") ;
      else
	printf ("Active:\n") ;
      x = x2 = 0 ;
      nn = nna = 0 ;
      for (ii = 0, pp = arrp (ll->a, ii, P) ; ii < arrayMax (ll->a) ; pp++, ii++)
	{
	  if ((ip && !pp->placebo) || (!ip && pp->placebo))
	      continue ;
	  nn++ ;
	  if (pp->pain >= 0)
	    { x += pp->pain ; x2 += pp->pain * pp->pain ; nna++ ; }
	}
      
      printf ("Found %d patients\n", nn) ;
      
      m[ip] = x/nna ; s[ip] = (x2 - x * x/nna)/(nna*(nna-1)) ; 
      printf ("\tWe recorded the pain of %d,", nna) ;
      printf (" average pain %.2f VAS,", (float)(x/nna)) ;
      printf (" standard deviation %.2f VAS\n", (float)sqrt(x2 - (x*x)/nna)/(nna - 1)) ;
    }
  eps = (m[1] - m[0])/sqrt (s[0] + s[1]) ;
  if (eps < 0) eps = - eps ;
  if (eps < 1.96) printf ("The 2 observed means are compatible : eps = %.2f < 1.96\n", eps) ;
  else printf ("The 2 observed means are NOT compatible : eps = %.2f > 1.96\n", eps) ;

} /* showPain */

/**********************************************************************/
 
static void showFunction (LL *ll) 
{
  int ip, ii, nn, nna ;
  double x, x2 ;
  double eps, m[2], s[2] ;
  P *pp ;

  for (ip = 1 ; ip >= 0 ; ip-- )
    {
      if (ip == 1)
	printf ("Placebo:\n") ;
      else
	printf ("Active:\n") ;
      x = x2 = 0 ;
      nn = nna = 0 ;
      for (ii = 0, pp = arrp (ll->a, ii, P) ; ii < arrayMax (ll->a) ; pp++, ii++)
	{
	  if ((ip && !pp->placebo) || (!ip && pp->placebo))
	      continue ;
	  nn++ ;
	  if (pp->function >= 0)
	    { x += pp->function ; x2 += pp->function * pp->function ; nna++ ; }
	}
      
      printf ("Found %d patients\n", nn) ;
      
      m[ip] = x/nna ; s[ip] = (x2 - x * x/nna)/(nna*(nna-1)) ; 
      printf ("\tWe recorded the function of %d,", nna) ;
      printf (" average WOMAC function %.2f,", (float)(x/nna)) ;
      printf (" standard deviation %.2f\n", (float)sqrt(x2 - (x*x)/nna)/(nna - 1)) ;
    }
  eps = (m[1] - m[0])/sqrt (s[0] + s[1]) ;
  if (eps < 0) eps = - eps ;
  if (eps < 1.96) printf ("The 2 observed means are compatible : eps = %.2f < 1.96\n", eps) ;
  else printf ("The 2 observed means are NOT compatible : eps = %.2f > 1.96\n", eps) ;

} /* showFunction */

/**********************************************************************/
 
static void showGlobal (LL *ll) 
{
  int ip, ii, nn, nna ;
  double x, x2 ;
  double eps, m[2], s[2] ;
  P *pp ;

  for (ip = 1 ; ip >= 0 ; ip-- )
    {
      if (ip == 1)
	printf ("Placebo:\n") ;
      else
	printf ("Active:\n") ;
      x = x2 = 0 ;
      nn = nna = 0 ;
      for (ii = 0, pp = arrp (ll->a, ii, P) ; ii < arrayMax (ll->a) ; pp++, ii++)
	{
	  if ((ip && !pp->placebo) || (!ip && pp->placebo))
	      continue ;
	  nn++ ;
	  if (pp->global >= 0)
	    { x += pp->global ; x2 += pp->global * pp->global ; nna++ ; }
	}
      
      printf ("Found %d patients\n", nn) ;
      
      m[ip] = x/nna ; s[ip] = (x2 - x * x/nna)/(nna*(nna-1)) ; 
      printf ("\tWe recorded the global of %d,", nna) ;
      printf (" average patient global assessment at baseline global %.2f VAS,", (float)(x/nna)) ;
      printf (" standard deviation %.2f VAS\n", (float)sqrt(x2 - (x*x)/nna)/(nna - 1)) ;
    }
  eps = (m[1] - m[0])/sqrt (s[0] + s[1]) ;
  if (eps < 0) eps = - eps ;
  if (eps < 1.96) printf ("The 2 observed means are compatible : eps = %.2f < 1.96\n", eps) ;
  else printf ("The 2 observed means are NOT compatible : eps = %.2f > 1.96\n", eps) ;

} /* showGlobal */

/**********************************************************************/
 
static void getTable1 (LL *ll) 
{
  printf ("\n##################################################################\n") ;
  printf ("TABLE 1\n") ;

  printf ("\n############# Age\n") ;
  showAge (ll) ;
  printf ("\n############# Sex\n") ;
  showSex (ll) ;
  printf ("\n############# BMI\n") ;
  showBmi (ll) ;
  printf ("\n############# Duration\n") ;
  showDuration (ll) ;
  printf ("\n############# Joint\n") ;
  showJoint (ll) ;
  printf ("\n############# Pain\n") ;
  showPain (ll) ;
  printf ("\n############# Function\n") ;
  showFunction (ll) ;
  printf ("\n############# Global\n") ;
  showGlobal (ll) ;
} /* getTable1 */

/**********************************************************************/
/**********************************************************************/
 
static void showAcetoAge (LL *ll) 
{
  int ip, ii, nn, nna ;
  double x, x2 ;
  double eps, m[2], s[2] ;
  P *pp ;

  for (ip = 1 ; ip >= 0 ; ip-- )
    {
      if (ip == 1)
	printf ("Pills:\n") ;
      else
	printf ("No pill:\n") ;
      x = x2 = 0 ;
      nn = nna = 0 ;
      for (ii = 0, pp = arrp (ll->a, ii, P) ; ii < arrayMax (ll->a) ; pp++, ii++)
	{
	  if ((ip && !pp->pill_count) || (!ip && pp->pill_count))
	      continue ;
	  nn++ ;
	  if (pp->age >= 0)
	    { x += pp->age ; x2 += pp->age * pp->age ; nna++ ; }
	}
      
      printf ("Found %d patients\n", nn) ;
      
      m[ip] = x/nna ; s[ip] = (x2 - x * x/nna)/(nna*(nna-1)) ; 
      printf ("\tWe recorded the age of %d,", nna) ;
      printf (" average %.2f years,", (float)(x/nna)) ;
      printf (" standard deviation %.2f years\n", (float)sqrt(x2 - (x*x)/nna)/(nna - 1)) ;
    }
  eps = (m[1] - m[0])/sqrt (s[0] + s[1]) ;
  if (eps < 0) eps = - eps ;
  if (eps < 1.96) printf ("The 2 observed means are compatible : eps = %.2f < 1.96\n", eps) ;
  else printf ("The 2 observed means are NOT compatible : eps = %.2f > 1.96\n", eps) ;

} /* showAcetoAge */

/**********************************************************************/
 
static void showAcetoSex (LL *ll) 
{
  int ip, ii, nn = 0, nMan = 0, nWoman = 0 ;
  double chi2, l1, l2, c1, c2, ntot ;
  double obs[4], th[4] ;
  P *pp ;

  for (ip = 1 ; ip >= 0 ; ip-- )
    {
      if (ip == 1)
	printf ("Pills:\n") ;
      else
	printf ("No pill:\n") ;

      nn = nMan = nWoman = 0 ;
      for (ii = 0, pp = arrp (ll->a, ii, P) ; ii < arrayMax (ll->a) ; pp++, ii++)
	{
	  if ((ip && !pp->pill_count) || (!ip && pp->pill_count))
	      continue ;
	  nn++ ;
	  if (pp->man) nMan++ ;
	  else if (pp->woman) nWoman++ ;
	}
      
      printf ("Found %d patients\n", nn) ;
      obs[2*ip] = nMan ;       obs[2*ip + 1] = nWoman ; 
      printf ("\t%d (%.2f%%) are men\n", nMan, (100.0*(float)nMan)/nn) ;
      printf ("\t%d (%.2f%%) are women\n", nWoman, (100.0*(float)(nWoman))/nn) ;
      printf ("\t%d (%.2f%%) are androgyne\n", nn - nMan - nWoman, (100.0*(float)(nn - nMan - nWoman))/nn) ;
    }
  ntot = obs[0] + obs[1] + obs[2] + obs[3] ;
  l1 = obs[0] + obs[1] ; l2 = obs[2] + obs[3] ;
  c1 = obs[0] + obs[2] ; c2 = obs[1] + obs[3] ;
  th[0] = l1 * c1 / ntot ;
  th[1] = l1 * c2 / ntot ;
  th[2] = l2 * c1 / ntot ;
  th[3] = l2 * c2 / ntot ;
  for (ip = 0, chi2 = 0 ; ip < 4 ; ip++)
    if (th[ip] < 5)
      chi2 = -9999 ;
    else
      chi2 += (obs[ip] - th[ip]) *  (obs[ip] - th[ip])/ th[ip] ;
  
  if (chi2< 1.96) printf ("The chi2 (1ddl) shows compatibility: chi2 = %.2f < 1.96\n", (float) chi2) ;
  else printf ("The chi2 (1ddl) shows NON compatibility: chi2 = %.2f >= 1.96\n", (float) chi2) ;

} /* showAcetoSex */

/**********************************************************************/
 
static void showAcetoBmi (LL *ll) 
{
  int ip, ii, nn, nna ;
  double x, x2 ;
  double eps, m[2], s[2] ;
  P *pp ;

  for (ip = 1 ; ip >= 0 ; ip-- )
    {
      if (ip == 1)
	printf ("Pills:\n") ;
      else
	printf ("No pill:\n") ;
      x = x2 = 0 ;
      nn = nna = 0 ;
      for (ii = 0, pp = arrp (ll->a, ii, P) ; ii < arrayMax (ll->a) ; pp++, ii++)
	{
	  if ((ip && !pp->pill_count) || (!ip && pp->pill_count))
	      continue ;
	  nn++ ;
	  if (pp->bmi >= 0)
	    { x += pp->bmi ; x2 += pp->bmi * pp->bmi ; nna++ ; }
	}
      
      printf ("Found %d patients\n", nn) ;
      
      m[ip] = x/nna ; s[ip] = (x2 - x * x/nna)/(nna*(nna-1)) ; 
      printf ("\tWe recorded the bmi of %d,", nna) ;
      printf (" average bmi %.2f kg/m2,", (float)(x/nna)) ;
      printf (" standard deviation %.2f kg/m2\n", (float)sqrt(x2 - (x*x)/nna)/(nna - 1)) ;
    }
  eps = (m[1] - m[0])/sqrt (s[0] + s[1]) ;
  if (eps < 0) eps = - eps ;
  if (eps < 1.96) printf ("The 2 observed means are compatible : eps = %.2f < 1.96\n", eps) ;
  else printf ("The 2 observed means are NOT compatible : eps = %.2f > 1.96\n", eps) ;

} /* showAcetoBmi */

/**********************************************************************/
 
static void showAcetoDuration (LL *ll) 
{
  int ip, ii, nn, nna ;
  double x, x2 ;
  double eps, m[2], s[2] ;
  P *pp ;

  for (ip = 1 ; ip >= 0 ; ip-- )
    {
      if (ip == 1)
	printf ("Pills:\n") ;
      else
	printf ("No pill:\n") ;
      x = x2 = 0 ;
      nn = nna = 0 ;
      for (ii = 0, pp = arrp (ll->a, ii, P) ; ii < arrayMax (ll->a) ; pp++, ii++)
	{
	  if ((ip && !pp->pill_count) || (!ip && pp->pill_count))
	      continue ;
	  nn++ ;
	  if (pp->duration >= 0)
	    { x += pp->duration ; x2 += pp->duration * pp->duration ; nna++ ; }
	}
      
      printf ("Found %d patients\n", nn) ;
      
      m[ip] = x/nna ; s[ip] = (x2 - x * x/nna)/(nna*(nna-1)) ; 
      printf ("\tWe recorded the duration of %d,", nna) ;
      printf (" average duration %.2f years,", (float)(x/nna)) ;
      printf (" standard deviation %.2f years\n", (float)sqrt(x2 - (x*x)/nna)/(nna - 1)) ;
    }
  eps = (m[1] - m[0])/sqrt (s[0] + s[1]) ;
  if (eps < 0) eps = - eps ;
  if (eps < 1.96) printf ("The 2 observed means are compatible : eps = %.2f < 1.96\n", eps) ;
  else printf ("The 2 observed means are NOT compatible : eps = %.2f > 1.96\n", eps) ;

} /* showAcetoDuration */

/**********************************************************************/
 
static void showAcetoJoint (LL *ll) 
{
  int ip, ii, nn = 0, nKnee = 0, nHip = 0 ;
  double chi2, l1, l2, c1, c2, ntot ;
  double obs[4], th[4] ;
  P *pp ;

  for (ip = 1 ; ip >= 0 ; ip-- )
    {
      if (ip == 1)
	printf ("Pills:\n") ;
      else
	printf ("No pill:\n") ;

      nn = nKnee = nHip = 0 ;
      for (ii = 0, pp = arrp (ll->a, ii, P) ; ii < arrayMax (ll->a) ; pp++, ii++)
	{
	  if ((ip && !pp->pill_count) || (!ip && pp->pill_count))
	      continue ;
	  nn++ ;
	  if (pp->knee) nKnee++ ;
	  else if (pp->hip) nHip++ ;
	}
      
      printf ("Found %d patients\n", nn) ;
      obs[2*ip] = nKnee ;       obs[2*ip + 1] = nHip ; 
      printf ("\t%d (%.2f%%) are knee\n", nKnee, (100.0*(float)nKnee)/nn) ;
      printf ("\t%d (%.2f%%) are hip\n", nHip, (100.0*(float)(nHip))/nn) ;
      printf ("\t%d (%.2f%%) are unknown\n", nn - nKnee - nHip, (100.0*(float)(nn - nKnee - nHip))/nn) ;
    }
  ntot = obs[0] + obs[1] + obs[2] + obs[3] ;
  l1 = obs[0] + obs[1] ; l2 = obs[2] + obs[3] ;
  c1 = obs[0] + obs[2] ; c2 = obs[1] + obs[3] ;
  th[0] = l1 * c1 / ntot ;
  th[1] = l1 * c2 / ntot ;
  th[2] = l2 * c1 / ntot ;
  th[3] = l2 * c2 / ntot ;
  for (ip = 0, chi2 = 0 ; ip < 4 ; ip++)
    if (th[ip] < 5)
      chi2 = -9999 ;
    else
      chi2 += (obs[ip] - th[ip]) *  (obs[ip] - th[ip])/ th[ip] ;
  
  if (chi2< 1.96) printf ("The chi2 (1ddl) shows compatibility: chi2 = %.2f < 1.96\n", (float) chi2) ;
  else printf ("The chi2 (1ddl) shows NON compatibility: chi2 = %.2f >= 1.96\n", (float) chi2) ;

} /* showAcetoJoint */

/**********************************************************************/
 
static void showAcetoPain (LL *ll) 
{
  int ip, ii, nn, nna ;
  double x, x2 ;
  double eps, m[2], s[2] ;
  P *pp ;

  for (ip = 1 ; ip >= 0 ; ip-- )
    {
      if (ip == 1)
	printf ("Pills:\n") ;
      else
	printf ("No pill:\n") ;
      x = x2 = 0 ;
      nn = nna = 0 ;
      for (ii = 0, pp = arrp (ll->a, ii, P) ; ii < arrayMax (ll->a) ; pp++, ii++)
	{
	  if ((ip && !pp->pill_count) || (!ip && pp->pill_count))
	      continue ;
	  nn++ ;
	  if (pp->pain >= 0)
	    { x += pp->pain ; x2 += pp->pain * pp->pain ; nna++ ; }
	}
      
      printf ("Found %d patients\n", nn) ;
      
      m[ip] = x/nna ; s[ip] = (x2 - x * x/nna)/(nna*(nna-1)) ; 
      printf ("\tWe recorded the pain of %d,", nna) ;
      printf (" average pain %.2f VAS,", (float)(x/nna)) ;
      printf (" standard deviation %.2f VAS\n", (float)sqrt(x2 - (x*x)/nna)/(nna - 1)) ;
    }
  eps = (m[1] - m[0])/sqrt (s[0] + s[1]) ;
  if (eps < 0) eps = - eps ;
  if (eps < 1.96) printf ("The 2 observed means are compatible : eps = %.2f < 1.96\n", eps) ;
  else printf ("The 2 observed means are NOT compatible : eps = %.2f > 1.96\n", eps) ;

} /* showAcetoPain */

/**********************************************************************/
 
static void showAcetoFunction (LL *ll) 
{
  int ip, ii, nn, nna ;
  double x, x2 ;
  double eps, m[2], s[2] ;
  P *pp ;

  for (ip = 1 ; ip >= 0 ; ip-- )
    {
      if (ip == 1)
	printf ("Pills:\n") ;
      else
	printf ("No pill:\n") ;
      x = x2 = 0 ;
      nn = nna = 0 ;
      for (ii = 0, pp = arrp (ll->a, ii, P) ; ii < arrayMax (ll->a) ; pp++, ii++)
	{
	  if ((ip && !pp->pill_count) || (!ip && pp->pill_count))
	      continue ;
	  nn++ ;
	  if (pp->function >= 0)
	    { x += pp->function ; x2 += pp->function * pp->function ; nna++ ; }
	}
      
      printf ("Found %d patients\n", nn) ;
      
      m[ip] = x/nna ; s[ip] = (x2 - x * x/nna)/(nna*(nna-1)) ; 
      printf ("\tWe recorded the function of %d,", nna) ;
      printf (" average WOMAC function %.2f,", (float)(x/nna)) ;
      printf (" standard deviation %.2f\n", (float)sqrt(x2 - (x*x)/nna)/(nna - 1)) ;
    }
  eps = (m[1] - m[0])/sqrt (s[0] + s[1]) ;
  if (eps < 0) eps = - eps ;
  if (eps < 1.96) printf ("The 2 observed means are compatible : eps = %.2f < 1.96\n", eps) ;
  else printf ("The 2 observed means are NOT compatible : eps = %.2f > 1.96\n", eps) ;

} /* showAcetoFunction */

/**********************************************************************/
 
static void showAcetoGlobal (LL *ll) 
{
  int ip, ii, nn, nna ;
  double x, x2 ;
  double eps, m[2], s[2] ;
  P *pp ;

  for (ip = 1 ; ip >= 0 ; ip-- )
    {
      if (ip == 1)
	printf ("Pills:\n") ;
      else
	printf ("No pill:\n") ;
      x = x2 = 0 ;
      nn = nna = 0 ;
      for (ii = 0, pp = arrp (ll->a, ii, P) ; ii < arrayMax (ll->a) ; pp++, ii++)
	{
	  if ((ip && !pp->pill_count) || (!ip && pp->pill_count))
	      continue ;
	  nn++ ;
	  if (pp->global >= 0)
	    { x += pp->global ; x2 += pp->global * pp->global ; nna++ ; }
	}
      
      printf ("Found %d patients\n", nn) ;
      
      m[ip] = x/nna ; s[ip] = (x2 - x * x/nna)/(nna*(nna-1)) ; 
      printf ("\tWe recorded the global of %d,", nna) ;
      printf (" average patient global assessment at baseline global %.2f VAS,", (float)(x/nna)) ;
      printf (" standard deviation %.2f VAS\n", (float)sqrt(x2 - (x*x)/nna)/(nna - 1)) ;
    }
  eps = (m[1] - m[0])/sqrt (s[0] + s[1]) ;
  if (eps < 0) eps = - eps ;
  if (eps < 1.96) printf ("The 2 observed means are compatible : eps = %.2f < 1.96\n", eps) ;
  else printf ("The 2 observed means are NOT compatible : eps = %.2f > 1.96\n", eps) ;

} /* showAcetoGlobal */

/**********************************************************************/
 
static void getTable2 (LL *ll) 
{
  AC_KEYSET ks ;
  int nn ;

  printf ("\n##################################################################\n") ;
  printf ("TABLE 2\n") ;

  ks = ac_dbquery_keyset (ll->db, "Find Patient pills", 0) ;
  nn = ac_keyset_count (ks) ;
  ac_free (ks) ;

  printf ("We have the pill count for %d patients\n\n", nn) ;
} /* getTable2 */

/**********************************************************************/
 
static void getTable3 (LL *ll) 
{
  printf ("\n##################################################################\n") ;
  printf ("TABLE 3\n") ;

  printf ("sorry there is no data to complete table 3\n\n") ;
} /* getTable3 */

/**********************************************************************/

static void pill2relativechangeinpain (LL *ll)
{
  int ip, ii, nn = 0 ;
  double xx, yy, xy, x2, y2, x, y, r, t ;
  P *pp ;
  
  for (ip = 0 ; ip < 3 ; ip++)
    {
      switch (ip)
	{
	case 0:
	  printf ("Pill/pain Whole:\n") ;
	  break ;
	case 1:
	  printf ("Pill/pain Active:\n") ;
	  break ;
	case 2:
	  printf ("Pill/pain Placebo:\n") ;
	  break ;
	}

      nn = 0 ;
      xx = yy = xy = x2 = y2 = 0 ;
      for (ii = 0, pp = arrp (ll->a, ii, P) ; ii < arrayMax (ll->a) ; pp++, ii++)
	{
	  switch (ip)
	    {
	    case 0:
	      break ;
	    case 1:
	      if (pp->placebo)
		continue ;
	      break ;
	    case 2:
	      if (!pp->placebo)
		continue ;
	      break ;
	    }
	  x = pp->pill_count ;
	  y = pp->pain6 * 100.0/pp->pain ;
	  if (x < 0 || y < -999)
	    continue ;
	  xx += x ; yy += y ; 
	  x2 += x * x ; y2 += y * y ; xy += x * y ;
	  nn++ ;
	}
      
      r = ( xy - (xx * yy)/nn) / sqrt ((x2 - (xx * xx)/nn) * ( y2 - (yy * yy)/nn)) ;
      t = r * sqrt (nn - 2) / sqrt (1 - r*r) ;
      printf ("Found %d patients: rho = %.2f t= %.2f\n", nn, (float) r, (float) t) ;
      if (0) printf ("xx=%lf yy=%lf xy=%lf x2 = %lf y2 = %lf\n"
	      , (float)xx
	      , (float)yy
	      , (float)xy
	      , (float)x2
	      , (float)y2
	      ) ;
    }

} /* pill2relativechangeinpain */

/**********************************************************************/

static void pill2relativechangeinglobal (LL *ll)
{
  int ip, ii, nn = 0 ;
  double xx, yy, xy, x2, y2, x, y, r, t ;
  P *pp ;
  
  for (ip = 0 ; ip < 3 ; ip++)
    {
      switch (ip)
	{
	case 0:
	  printf ("Pill/global Whole:\n") ;
	  break ;
	case 1:
	  printf ("Pill/global Active:\n") ;
	  break ;
	case 2:
	  printf ("Pill/global Placebo:\n") ;
	  break ;
	}

      nn = 0 ;
      xx = yy = xy = x2 = y2 = 0 ;
      for (ii = 0, pp = arrp (ll->a, ii, P) ; ii < arrayMax (ll->a) ; pp++, ii++)
	{
	  switch (ip)
	    {
	    case 0:
	      break ;
	    case 1:
	      if (pp->placebo)
		continue ;
	      break ;
	    case 2:
	      if (!pp->placebo)
		continue ;
	      break ;
	    }
	  x = pp->pill_count ;
	  y = pp->global6 * 100.0/pp->global ;
	  if (x < 0 || y < -999)
	    continue ;
	  xx += x ; yy += y ; 
	  x2 += x * x ; y2 += y * y ; xy += x * y ;
	  nn++ ;
	}
      
      r = ( xy - (xx * yy)/nn) / sqrt ((x2 - (xx * xx)/nn) * ( y2 - (yy * yy)/nn)) ;
      t = r * sqrt (nn - 2) / sqrt (1 - r*r) ;
      printf ("Found %d patients: rho = %.2f  t= %.2f\n", nn, (float) r, (float) t) ;
      if (0) printf ("xx=%lf yy=%lf xy=%lf x2 = %lf y2 = %lf\n"
	      , (float)xx
	      , (float)yy
	      , (float)xy
	      , (float)x2
	      , (float)y2
	      ) ;
    }

} /* pill2relativechangeinglobal */

/**********************************************************************/

static void pill2relativechangeinfunction (LL *ll)
{
  int ip, ii, nn = 0 ;
  double xx, yy, xy, x2, y2, x, y, r, t ;
  P *pp ;
  
  for (ip = 0 ; ip < 3 ; ip++)
    {
      switch (ip)
	{
	case 0:
	  printf ("Pill/function Whole:\n") ;
	  break ;
	case 1:
	  printf ("Pill/function Active:\n") ;
	  break ;
	case 2:
	  printf ("Pill/function Placebo:\n") ;
	  break ;
	}

      nn = 0 ;
      xx = yy = xy = x2 = y2 = 0 ;
      for (ii = 0, pp = arrp (ll->a, ii, P) ; ii < arrayMax (ll->a) ; pp++, ii++)
	{
	  switch (ip)
	    {
	    case 0:
	      break ;
	    case 1:
	      if (pp->placebo)
		continue ;
	      break ;
	    case 2:
	      if (!pp->placebo)
		continue ;
	      break ;
	    }
	  x = pp->pill_count ;
	  y = pp->function6 * 100.0/pp->function ;
	  if (x < 0 || y < -999)
	    continue ;
	  xx += x ; yy += y ; 
	  x2 += x * x ; y2 += y * y ; xy += x * y ;
	  nn++ ;
	}
      
      r = ( xy - (xx * yy)/nn) / sqrt ((x2 - (xx * xx)/nn) * ( y2 - (yy * yy)/nn)) ;
      t = r * sqrt (nn - 2) / sqrt (1 - r*r) ;
      printf ("Found %d patients: rho = %.2f t= %.2f\n", nn, (float) r, (float) t) ;
      if (0) printf ("xx=%lf yy=%lf xy=%lf x2 = %lf y2 = %lf\n"
	      , (float)xx
	      , (float)yy
	      , (float)xy
	      , (float)x2
	      , (float)y2
	      ) ;
    }

} /* pill2relativechangeinfunction */

/**********************************************************************/
 
static void getTable4 (LL *ll) 
{
  printf ("\n##################################################################\n") ;
  printf ("TABLE 4\n") ;
  printf ("Schwartz 4e edition page 212\n") ;
  pill2relativechangeinpain (ll) ;
  pill2relativechangeinglobal (ll) ;
  pill2relativechangeinfunction (ll) ;
} /* getTable4 */

/**********************************************************************/
 
static void getTable5 (LL *ll) 
{
  printf ("\n##################################################################\n") ;
  printf ("TABLE 5\n") ;

  /*
  trmt2pain (ll) ;
  trmt2global (ll) ;
  trmt2function (ll) ;
  */
} /* getTable5 */

/**********************************************************************/
 
static void showTable6 (LL *ll, int type, int nPill) 
{
  int ip, ii, nn, ntot ;
  double ss, chi2, m[2], l1, l2, c1, c2 ;
  double obs[4], th[4] ;
  P *pp ;

  for (ip = 1 ; ip >= 0 ; ip-- )
    {
      if (ip == 1)
	printf ("Placebo:\n") ;
      else
	printf ("Active:\n") ;

      nn = ntot = 0 ;
      for (ii = 0, pp = arrp (ll->a, ii, P) ; ii < arrayMax (ll->a) ; pp++, ii++)
	{
	  if ((ip && !pp->placebo) || (!ip && pp->placebo))
	      continue ;
	  if (pp->pain6 < -999)
	    continue ;
	  ntot++ ;
	  if (pp->pills_per_day < 0 || pp->pills_per_day > nPill/1000.0)
	    continue ;
	  if (100.0 * pp->pain6/pp->pain > -30.0)
	    continue ;
	  nn++ ;
	}
      
      printf ("Found %d patients\n", nn) ;
      
      m[ip] = ((double)nn)/ntot ;
      printf ("\tChange in pain < -30 in  %.3f (%d/%d) of the patients,", m[ip], nn, ntot) ;
      obs[2*ip] = nn ; obs[2*ip + 1] = ntot - nn ;
     }
  printf ("Treatment effect %.2f %%\n", (float) 100.0 * (m[0] - m[1])) ;
  ntot = obs[0] + obs[1] + obs[2] + obs[3] ;
  l1 = obs[0] + obs[1] ; l2 = obs[2] + obs[3] ;
  c1 = obs[0] + obs[2] ; c2 = obs[1] + obs[3] ;
  th[0] = l1 * c1 / ntot ;
  th[1] = l1 * c2 / ntot ;
  th[2] = l2 * c1 / ntot ;
  th[3] = l2 * c2 / ntot ;
  for (ip = 0, chi2 = 0 ; ip < 4 ; ip++)
    if (th[ip] < 5)
      chi2 = -9999 ;
    else
      chi2 += (obs[ip] - th[ip]) *  (obs[ip] - th[ip])/ th[ip] ;
  
  if (chi2< 1.96) printf ("The chi2 (1ddl) shows compatibility: chi2 = %.2f < 1.96\n", (float) chi2) ;
  else printf ("The chi2 (1ddl) shows NON compatibility: chi2 = %.2f >= 1.96\n", (float) chi2) ;
  ss = 13.0 * (m[1]*(1.0 - m[1]) + m[0]*(1.0 - m[0]))/((m[1] - m[0])*(m[1] - m[0])) ;
  printf ("Sample size N %d \n(N = C(beta) (p1.q1 + p2.q2)/(p2 - p1)^2,  with C(beta)=13, q = 1-p)\n", (int) ss) ;
  printf ("######################\n") ;
} /* showTable6 */

/**********************************************************************/
 
static void getTable6 (LL *ll) 
{
  int nn, nPill ;
  Array aa = arrayCreate (100, int) ;
  P *pp ;
  int ii, n, n1, n2, n3 ;
  printf ("\n##################################################################\n") ;
  printf ("TABLE 6\n") ;
  
  nn = nPill = 0 ;
  for (ii = 0, pp = arrp (ll->a, ii, P) ; ii < arrayMax (ll->a) ; pp++, ii++)
    {
      if (pp->pills_per_day >= 0)
	{
	  nn++ ; nPill += 1000 * pp->pills_per_day ;
	  array (aa, 1000 * pp->pills_per_day, int)++ ;
	}
    }
  
  n = 0 ; 
  n1 = n2 = n3 = -1 ;
  for (ii = 0 ; ii < arrayMax (aa) ; ii++)
    {
      n += array (aa, ii, int) ;
      if (n1 < 0 && n > nn/4) n1 = ii ;
      if (n2 < 0 && n > nn/2) n2 = ii ;
      if (n3 < 0 && n > nn * .75) n3 = ii ;
    }

  for (ii = 0 ; ii < 6 ; ii++)
    {
      switch (ii)
	{
	case 0:
	  printf ("%% patients with change in relative pain\n") ;
	  showTable6 (ll, ii, 999999) ;
	  break ;
	case 1:
	  printf ("idem && aceto intake < %.3f pills per day (first quartile)\n", n1/1000.0) ;
	  showTable6 (ll, ii, n1) ;
	  break ;
	case 2:
	  printf ("idem && aceto intake <  %.3f pills per day (median)\n", n2/1000.0) ;
	  showTable6 (ll, ii, n2) ;
	  break ;
	case 3:
	  printf ("idem && aceto intake <  %.3f pills per day (third quartile)\n", n3/1000.0) ;
	  showTable6 (ll, ii, n3) ;
	  break ;
	case 4:
	  printf ("%% patients with change in pain and aceto < 2 days/week\n") ;
	  printf ("je sais pas\n") ;
	  break ;
	case 5:
	  printf ("%% patients with change in pain and aceto < 4 days/week\n") ;
	  printf ("je sais pas\n") ;
	  break ;
	}
    }
} /* getTable6 */

/**********************************************************************/
 
static void getTable7 (LL *ll) 
{
  Array aa = 0 ;
  P *pp ;
  int jj, ii, n, n1, n2, n3, nn, nnP, nnA, nn0, nnP0, nnA0 ;
  
  printf ("\n##################################################################\n") ;
  printf ("TABLE 7\n") ;
  printf ("We only count patient where the pils_per_day is known\n") ;

  for (jj = 0 ; jj < 3 ; jj++)
    {
      /* fabricate the histo of pills per day */
      nn = nnA = nnP = 0 ; aa = arrayReCreate (aa, 100, int) ;
      for (ii = 0, pp = arrp (ll->a, ii, P) ; ii < arrayMax (ll->a) ; pp++, ii++)
	{
	  if (pp->pills_per_day < 0)
	    continue ;
	  
	  switch (jj)
	    {
	    case 0: if (pp->placebo) continue ; break ;  /* active */
	    case 1: if (!pp->placebo) continue ; break ; /* placebo */
	    case 2: break ; /* whole */
	    }

	  if (1)
	    {
	      nn++ ;
	      array (aa, 1000.0 * pp->pills_per_day, int)++ ;
	    }
	}
      
      /* locate the quartiles */
      n = 0 ; 
      n1 = n2 = n3 = -1 ;
      for (ii = 0 ; ii < arrayMax (aa) ; ii++)
	{
	  n += array (aa, ii, int) ;
	  if (n1 < 0 && n > nn/4) n1 = ii ;
	  if (n2 < 0 && n > nn/2) n2 = ii ;
	  if (n3 < 0 && n > nn * .75) n3 = ii ;
	}
      /* export the limits of the quartiles */
	  
      switch (jj)
	{
	case 0: /* active */
	  printf ("\nFor the %d active patients the quartiles are %.3f ;  %.3f ;  %.3f pills per day\n"
		  , nn
		  , (float) (n1/1000.0)
		  , (float) (n2/1000.0)
		  , (float) (n3/1000.0)
		  ) ;
	  break ;
	case 1: /* placebo */
	  printf ("\nFor the %d placebo patients the quartiles are %.3f ;  %.3f ;  %.3f pills per day\n"
		  , nn
		  , (float) (n1/1000.0)
		  , (float) (n2/1000.0)
		  , (float) (n3/1000.0)
		  ) ;
	  break ;
	case 2:  /* whole */
	  printf ("\nFor the %d whole patients the quartiles are %.3f ;  %.3f ;  %.3f pills per day\n"
		  , nn
		  , (float) (n1/1000.0)
		  , (float) (n2/1000.0)
		  , (float) (n3/1000.0)
		  ) ;
	  break ;
	}      
    }
  /* note that we came out with the limits for whole */
  /* count patients in different quartiles */
  printf ("%30s%8s%8s%8s\n", "Type", "Whole", "Active", "Placebo") ;
  nn0 = nnA0 = nnP0 = 0 ;
  for (ii = 0, pp = arrp (ll->a, ii, P) ; ii < arrayMax (ll->a) ; pp++, ii++)
    {
      if (pp->pills_per_day < 0)
	continue ;
      nn0++ ;
      if (pp->placebo)
	nnP0++ ; 
      else
	nnA0++ ;
    }
  printf ("%30s%8d%8d%8d\n","# patient", nn0, nnA0, nnP0) ;

  nn = nnA = nnP = 0 ;
  for (ii = 0, pp = arrp (ll->a, ii, P) ; ii < arrayMax (ll->a) ; pp++, ii++)
    {
      if (pp->pills_per_day <0)
	continue ;
      if (1000.0 * pp->pills_per_day > n1)
	continue ;
      nn++ ;
      if (pp->placebo)
	nnP++ ; 
      else
	nnA++ ;
    }
  printf ("%30s (%.3f pills_per_day) %.1f%% (%d/%d)\t%.1f%% (%d/%d)\t%.1f%% (%d/%d)\n","1st quartile intake"
	  , (float) (n1/1000.0)
	  , (float)(100.0 * nn)/nn0, nn, nn0
	  , (float)(100.0 * nnA)/nnA0, nnA, nnA0
	  , (float)(100.0 * nnP)/nnP0, nnP, nnP0
	  ) ;
  nn = nnA = nnP = 0 ;
  for (ii = 0, pp = arrp (ll->a, ii, P) ; ii < arrayMax (ll->a) ; pp++, ii++)
    {
      if (pp->pills_per_day <0)
	continue ;
      if (1000.0 * pp->pills_per_day > n2)
	continue ;
      nn++ ;
      if (pp->placebo)
	nnP++ ; 
      else
	nnA++ ;
    }
  printf ("%30s (%.3f pills_per_day) %.1f%% (%d/%d)\t%.1f%% (%d/%d)\t%.1f%% (%d/%d)\n","2nd quartile intake"
	  , (float) (n2/1000.0)
	  , (float)(100.0 * nn)/nn0, nn, nn0
	  , (float)(100.0 * nnA)/nnA0, nnA, nnA0
	  , (float)(100.0 * nnP)/nnP0, nnP, nnP0
	  ) ;
  nn = nnA = nnP = 0 ;
  for (ii = 0, pp = arrp (ll->a, ii, P) ; ii < arrayMax (ll->a) ; pp++, ii++)
    {
      if (pp->pills_per_day <0)
	continue ;
      if (1000.0 * pp->pills_per_day > n3)
	continue ;
      nn++ ;
      if (pp->placebo)
	nnP++ ; 
      else
	nnA++ ;
    }
  printf ("%30s (%.3f pills_per_day) %.1f%% (%d/%d)\t%.1f%% (%d/%d)\t%.1f%% (%d/%d)\n","3rd quartile intake"
	  , (float) (n3/1000.0)
	  , (float)(100.0 * nn)/nn0, nn, nn0
	  , (float)(100.0 * nnA)/nnA0, nnA, nnA0
	  , (float)(100.0 * nnP)/nnP0, nnP, nnP0
	  ) ;

  printf ("\nmean number of days of intake, we do not have those data\n") ;
  printf ("which test should we use to prove normality of the distribs ? \n") ;
  arrayDestroy (aa) ;
} /* getTable7 */

/**********************************************************************/
 
static void getTable8 (LL *ll) 
{
  printf ("\n##################################################################\n") ;
  printf ("TABLE 8\n") ;

  showAcetoAge (ll) ;
  showAcetoSex (ll) ;
  showAcetoBmi (ll) ;
  showAcetoDuration (ll) ;
  showAcetoJoint (ll) ;
  showAcetoPain (ll) ;
  showAcetoFunction (ll) ;
  showAcetoGlobal (ll) ;
} /* getTable8 */

/**********************************************************************/

static void pill2age (LL *ll)
{
  int ip, ii, nn = 0 ;
  double xx, yy, xy, x2, y2, x, y, r, t ;
  P *pp ;
  
  for (ip = 0 ; ip < 3 ; ip++)
    {
      switch (ip)
	{
	case 0:
	  printf ("Pill/age Whole:\n") ;
	  break ;
	case 1:
	  printf ("Pill/age Active:\n") ;
	  break ;
	case 2:
	  printf ("Pill/age Placebo:\n") ;
	  break ;
	}

      nn = 0 ;
      xx = yy = xy = x2 = y2 = 0 ;
      for (ii = 0, pp = arrp (ll->a, ii, P) ; ii < arrayMax (ll->a) ; pp++, ii++)
	{
	  switch (ip)
	    {
	    case 0:
	      break ;
	    case 1:
	      if (pp->placebo)
		continue ;
	      break ;
	    case 2:
	      if (!pp->placebo)
		continue ;
	      break ;
	    }
	  x = pp->pills_per_day ;
	  y = pp->age ;
	  if (x < 0 || y < 0)
	    continue ;
	  xx += x ; yy += y ; 
	  x2 += x * x ; y2 += y * y ; xy += x * y ;
	  nn++ ;
	}
      
      r = ( xy - (xx * yy)/nn) / sqrt ((x2 - (xx * xx)/nn) * ( y2 - (yy * yy)/nn)) ;
      t = r * sqrt (nn - 2) / sqrt (1 - r*r) ;
      printf ("Found %d patients: rho = %.2f t = %.2f\n", nn, (float) r, (float) t) ;
      if (0) printf ("xx=%lf yy=%lf xy=%lf x2 = %lf y2 = %lf\n"
	      , (float)xx/nn
	      , (float)yy/nn
	      , (float)xy/nn
	      , (float)x2/nn
	      , (float)y2/nn
	      ) ;
    }

} /* pill2age */

/**********************************************************************/

static void pill2bmi (LL *ll)
{
  int ii, nn = 0 ;
  double xx, yy, xy, x2, y2, x, y, r, t ;
  P *pp ;
  
  printf ("Pill/bmi Whole:\n") ;

  nn = 0 ;
  xx = yy = xy = x2 = y2 = 0 ;
  for (ii = 0, pp = arrp (ll->a, ii, P) ; ii < arrayMax (ll->a) ; pp++, ii++)
    {
      x = pp->pills_per_day ;
      y = pp->bmi ;
      if (x < 0 || y < 0)
	continue ;
      xx += x ; yy += y ; 
      x2 += x * x ; y2 += y * y ; xy += x * y ;
      nn++ ;
    }
  
  r = ( xy - (xx * yy)/nn) / sqrt ((x2 - (xx * xx)/nn) * ( y2 - (yy * yy)/nn)) ;
  t = r * sqrt (nn - 2) / sqrt (1 - r*r) ;
  printf ("Found %d patients: rho = %.2f t= %.2f\n", nn, (float) r, (float) t) ;
  if (0) printf ("xx=%lf yy=%lf xy=%lf x2 = %lf y2 = %lf\n"
		 , (float)xx
		 , (float)yy
		 , (float)xy
		 , (float)x2
		 , (float)y2
		 ) ;
} /* pill2bmi */

/**********************************************************************/

static void pill2baselinepain (LL *ll)
{
  int ii, nn = 0 ;
  double xx, yy, xy, x2, y2, x, y, r, t ;
  P *pp ;
  
  printf ("Pill/baselinepain Whole:\n") ;

  nn = 0 ;
  xx = yy = xy = x2 = y2 = 0 ;
  for (ii = 0, pp = arrp (ll->a, ii, P) ; ii < arrayMax (ll->a) ; pp++, ii++)
    {
      x = pp->pills_per_day ;
      y = pp->pain ;
      if (x < 0 || y < 0)
	continue ;
      xx += x ; yy += y ; 
      x2 += x * x ; y2 += y * y ; xy += x * y ;
      nn++ ;
    }
  
  r = ( xy - (xx * yy)/nn) / sqrt ((x2 - (xx * xx)/nn) * ( y2 - (yy * yy)/nn)) ;
  t = r * sqrt (nn - 2) / sqrt (1 - r*r) ;
  printf ("Found %d patients: rho = %.2f t= %.2f\n", nn, (float) r, (float) t) ;
  if (0) printf ("xx=%lf yy=%lf xy=%lf x2 = %lf y2 = %lf\n"
		 , (float)xx
		 , (float)yy
		 , (float)xy
		 , (float)x2
		 , (float)y2
		 ) ;
} /* pill2baselinepain */

/**********************************************************************/

static void pill2baselinefunction (LL *ll)
{
  int ii, nn = 0 ;
  double xx, yy, xy, x2, y2, x, y, r, t ;
  P *pp ;
  
  printf ("Pill/baselinefunction Whole:\n") ;

  nn = 0 ;
  xx = yy = xy = x2 = y2 = 0 ;
  for (ii = 0, pp = arrp (ll->a, ii, P) ; ii < arrayMax (ll->a) ; pp++, ii++)
    {
      x = pp->pills_per_day ;
      y = pp->function ;
      if (x < 0 || y < 0)
	continue ;
      xx += x ; yy += y ; 
      x2 += x * x ; y2 += y * y ; xy += x * y ;
      nn++ ;
    }
  
  r = ( xy - (xx * yy)/nn) / sqrt ((x2 - (xx * xx)/nn) * ( y2 - (yy * yy)/nn)) ;
  t = r * sqrt (nn - 2) / sqrt (1 - r*r) ;
  printf ("Found %d patients: rho = %.2f t= %.2f\n", nn, (float) r, (float) t) ;
  if (0) printf ("xx=%lf yy=%lf xy=%lf x2 = %lf y2 = %lf\n"
		 , (float)xx
		 , (float)yy
		 , (float)xy
		 , (float)x2
		 , (float)y2
	      ) ;
} /* pill2baselinefunction */

/**********************************************************************/

static void pill2baselineglobal (LL *ll)
{
  int ii, nn = 0 ;
  double xx, yy, xy, x2, y2, x, y, r, t ;
  P *pp ;
  
  printf ("Pill/baselineglobal Whole:\n") ;

  nn = 0 ;
  xx = yy = xy = x2 = y2 = 0 ;
  for (ii = 0, pp = arrp (ll->a, ii, P) ; ii < arrayMax (ll->a) ; pp++, ii++)
    {
      x = pp->pills_per_day ;
      y = pp->global ;
      if (x < 0 || y < 0)
	continue ;
      xx += x ; yy += y ; 
      x2 += x * x ; y2 += y * y ; xy += x * y ;
      nn++ ;
    }
  
  r = ( xy - (xx * yy)/nn) / sqrt ((x2 - (xx * xx)/nn) * ( y2 - (yy * yy)/nn)) ;
  t = r * sqrt (nn - 2) / sqrt (1 - r*r) ;
  printf ("Found %d patients: rho = %.2f t= %.2f\n", nn, (float) r, (float) t) ;
  if (0) printf ("xx=%lf yy=%lf xy=%lf x2 = %lf y2 = %lf\n"
		 , (float)xx
		 , (float)yy
		 , (float)xy
		 , (float)x2
		 , (float)y2
	      ) ;
} /* pill2baselineglobal */

/**********************************************************************/

static void pill2aucpain (LL *ll)
{
  int ii, nn = 0 ;
  double xx, yy, xy, x2, y2, x, y, r, t ;
  P *pp ;
  
  printf ("Pill/AUC pain Whole:\n") ;
  nn = 0 ;
  xx = yy = xy = x2 = y2 = 0 ;
  for (ii = 0, pp = arrp (ll->a, ii, P) ; ii < arrayMax (ll->a) ; pp++, ii++)
    {
      x = pp->pills_per_day ;
      y = pp->auc_pain ;
      if (x < 0 || y < 0)
	continue ;
      xx += x ; yy += y ; 
      x2 += x * x ; y2 += y * y ; xy += x * y ;
      nn++ ;
    }
  
  r = ( xy - (xx * yy)/nn) / sqrt ((x2 - (xx * xx)/nn) * ( y2 - (yy * yy)/nn)) ;
  t = r * sqrt (nn - 2) / sqrt (1 - r*r) ;
  printf ("Found %d patients: rho = %.2f t= %.2f\n", nn, (float) r, (float) t) ;
  if (0) printf ("xx=%lf yy=%lf xy=%lf x2 = %lf y2 = %lf\n"
		 , (float)xx
		 , (float)yy
		 , (float)xy
		 , (float)x2
		 , (float)y2
		 ) ;
} /* pill2aucpain */

/**********************************************************************/

static void pill2aucfunction (LL *ll)
{
  int ii, nn = 0 ;
  double xx, yy, xy, x2, y2, x, y, r, t ;
  P *pp ;
  
  printf ("Pill/AUC function Whole:\n") ;
  nn = 0 ;
  xx = yy = xy = x2 = y2 = 0 ;
  for (ii = 0, pp = arrp (ll->a, ii, P) ; ii < arrayMax (ll->a) ; pp++, ii++)
    {
      x = pp->pills_per_day ;
      y = pp->auc_function ;
      if (x < 0 || y < 0)
	continue ;
      xx += x ; yy += y ; 
      x2 += x * x ; y2 += y * y ; xy += x * y ;
      nn++ ;
    }
  
  r = ( xy - (xx * yy)/nn) / sqrt ((x2 - (xx * xx)/nn) * ( y2 - (yy * yy)/nn)) ;
  t = r * sqrt (nn - 2) / sqrt (1 - r*r) ;
  printf ("Found %d patients: rho = %.2f t= %.2f\n", nn, (float) r, (float) t) ;
  if (0) printf ("xx=%lf yy=%lf xy=%lf x2 = %lf y2 = %lf\n"
		 , (float)xx
		 , (float)yy
		 , (float)xy
		 , (float)x2
		 , (float)y2
		 ) ;
} /* pill2aucfunction */

/**********************************************************************/

static void pill2aucglobal (LL *ll)
{
  int ii, nn = 0 ;
  double xx, yy, xy, x2, y2, x, y, r, t ;
  P *pp ;
  
  printf ("Pill/AUC global Whole:\n") ;
  nn = 0 ;
  xx = yy = xy = x2 = y2 = 0 ;
  for (ii = 0, pp = arrp (ll->a, ii, P) ; ii < arrayMax (ll->a) ; pp++, ii++)
    {
      x = pp->pills_per_day ;
      y = pp->auc_global ;
      if (x < 0 || y < 0)
	continue ;
      xx += x ; yy += y ; 
      x2 += x * x ; y2 += y * y ; xy += x * y ;
      nn++ ;
    }
  
  r = ( xy - (xx * yy)/nn) / sqrt ((x2 - (xx * xx)/nn) * ( y2 - (yy * yy)/nn)) ;
  t = r * sqrt (nn - 2) / sqrt (1 - r*r) ;
  printf ("Found %d patients: rho = %.2f t= %.2f\n", nn, (float) r, (float) t) ;
  if (0) printf ("xx=%lf yy=%lf xy=%lf x2 = %lf y2 = %lf\n"
		 , (float)xx
		 , (float)yy
		 , (float)xy
		 , (float)x2
		 , (float)y2
		 ) ;
} /* pill2aucglobal */

/**********************************************************************/
 
static void pill2sex (LL *ll) 
{
  int ii, nn = 0, nMan = 0, nWoman = 0 ;
  double eps, sm, sw, xm, xw, xm2, xw2, zm, zw ;
  P *pp ;

  nn = nMan = nWoman = 0 ;
  xm = xm2 = xw = xw2 = 0 ;
  for (ii = 0, pp = arrp (ll->a, ii, P) ; ii < arrayMax (ll->a) ; pp++, ii++)
    {
      nn++ ;
      if (pp->man) { xm += pp->pills_per_day ; xm2 += pp->pills_per_day * pp->pills_per_day ; nMan++ ; }
      else if (pp->woman) { xw += pp->pills_per_day ; xw2 += pp->pills_per_day * pp->pills_per_day ; nWoman++ ; }
    }
      
  printf ("Found %d patients\n", nn) ;
  sm = sqrt (xm2 - (xm*xm)/nMan)/(nMan - 1) ;
  zm = (xm2 - (xm*xm)/nMan)/(nMan * (nMan - 1)) ;
  printf ("\t%d (%.2f%%) are men m=%.2f sd = %.2f\n"
	  , nMan, (100.0*(float)nMan)/nn, xm/nMan
	  , (float)sm
	  ) ;
  sw = sqrt (xw2 - (xw*xw)/nWoman)/(nWoman - 1) ;
  zw = (xw2 - (xw*xw)/nWoman)/(nWoman * (nWoman - 1)) ;
  printf ("\t%d (%.2f%%) are women m=%.2f sd = %.2f\n"
	  , nWoman, (100.0*(float)(nWoman))/nn, xw/nWoman
	  , (float)sw
	  ) ;
  printf ("\t%d are androgyne\n", nn - nMan - nWoman) ;

  eps = (xm/nMan - xw/nWoman)/sqrt (zm + zw) ;
  if (eps < 0) eps = - eps ;
  if (eps < 1.96) printf ("The 2 observed means are compatible : eps = %.2f < 1.96\n", (float)eps) ;
  else printf ("The 2 observed means are NOT compatible : eps = %.2f > 1.96\n", (float)eps) ;
  if (0) printf ("xm=%f xm2 = %f, xw=%f nMan = %d nWoman = %d sm = %f sw= %f\n"
	  , xm, xm2, xw, nMan, nWoman, sm, sw) ;
} /* pill2sex */

/**********************************************************************/
 
static void getTable8b (LL *ll) 
{
  printf ("\n##################################################################\n") ;
  printf ("TABLE 8b\n") ;

  pill2age (ll) ;
  pill2bmi (ll) ;
  pill2baselinepain (ll) ;
  pill2baselinefunction (ll) ;
  pill2baselineglobal (ll) ;
  pill2aucpain (ll) ;
  pill2aucfunction (ll) ;
  pill2aucglobal (ll) ;
  pill2sex (ll) ;
} /* getTable8b */

/**********************************************************************/
 
static void getTable9 (LL *ll) 
{
  printf ("\n##################################################################\n") ;
  printf ("TABLE 9\n") ;
  printf ("Je ne comprends pas bien  pas ce que tu voudrais calculer\n") ;
} /* getTable9 */

/**********************************************************************/
 
static void showTable10 (LL *ll, int type, int nPill) 
{
  int ip, ii, nn, ntot ;
  double ss, chi2, m[2], l1, l2, c1, c2 ;
  double obs[4], th[4] ;
  P *pp ;

  for (ip = 1 ; ip >= 0 ; ip-- )
    {
      if (ip == 1)
	printf ("Placebo:\n") ;
      else
	printf ("Active:\n") ;

      nn = ntot = 0 ;
      for (ii = 0, pp = arrp (ll->a, ii, P) ; ii < arrayMax (ll->a) ; pp++, ii++)
	{
	  if ((ip && !pp->placebo) || (!ip && pp->placebo))
	      continue ;
	  if (pp->pain6 < -999)
	    continue ;
	  ntot++ ;
	  if (pp->pill_count < 0 || pp->pill_count > nPill)
	    continue ;
	
	  /* NO test on pain change 
	     if (100.0 * pp->pain6/pp->pain > -30.0)
	     continue ;
	  */
	  nn++ ;
	}
      
      printf ("ZZFound %d patients\n", nn) ;
      
      m[ip] = ((double)nn)/ntot ;
      printf ("\t%.2f (%d/%d) of the patients take this amount of pills,", m[ip], nn, ntot) ;
      obs[2*ip] = nn ; obs[2*ip + 1] = ntot - nn ;
     }
  printf ("Treatment effect %.2f\n", (float) (m[1] - m[0])) ;
  ntot = obs[0] + obs[1] + obs[2] + obs[3] ;
  l1 = obs[0] + obs[1] ; l2 = obs[2] + obs[3] ;
  c1 = obs[0] + obs[2] ; c2 = obs[1] + obs[3] ;
  th[0] = l1 * c1 / ntot ;
  th[1] = l1 * c2 / ntot ;
  th[2] = l2 * c1 / ntot ;
  th[3] = l2 * c2 / ntot ;
  for (ip = 0, chi2 = 0 ; ip < 4 ; ip++)
    if (th[ip] < 5)
      chi2 = -9999 ;
    else
      chi2 += (obs[ip] - th[ip]) *  (obs[ip] - th[ip])/ th[ip] ;
  
  if (chi2< 1.96) printf ("The chi2 (1ddl) shows compatibility: chi2 = %.2f < 1.96\n", (float) chi2) ;
  else printf ("The chi2 (1ddl) shows NON compatibility: chi2 = %.2f >= 1.96\n", (float) chi2) ;
  ss = 13.0 * (m[1]*(1.0 - m[1]) + m[0]*(1.0 - m[0]))/((m[1] - m[0])*(m[1] - m[0])) ;
  printf ("Sample size N %d \n(N = C(beta) (p1.q1 + p2.q2)/(p2 - p1)^2,  with C(beta)=13, q = 1-p)\n", (int) ss) ;

  printf ("######################\n") ;
} /* showTable10 */

/**********************************************************************/
 
static void getTable10 (LL *ll) 
{
  int nn, nPill ;
  Array aa = arrayCreate (100, int) ;
  P *pp ;
  int ii, n, n1, n2, n3 ;

  printf ("\n##################################################################\n") ;
  printf ("TABLE 10\n") ;
  printf ("discriminant capacity of concomitant therapy in itself\n") ;
if (0)   printf ("je DUPLIQUE ici la table 6, puisque je vois pas la difference\n") ;

  nn = nPill = 0 ;
  for (ii = 0, pp = arrp (ll->a, ii, P) ; ii < arrayMax (ll->a) ; pp++, ii++)
    {
      if (pp->pill_count >= 0)
	{
	  nn++ ; nPill += pp->pill_count ;
	  array (aa, pp->pill_count, int)++ ;
	}
    }
  
  n = 0 ; 
  n1 = n2 = n3 = -1 ;
  for (ii = 0 ; ii < arrayMax (aa) ; ii++)
    {
      n += array (aa, ii, int) ;
      if (n1 < 0 && n > nn/4) n1 = ii ;
      if (n2 < 0 && n > nn/2) n2 = ii ;
      if (n3 < 0 && n > nn * .75) n3 = ii ;
    }

  for (ii = 0 ; ii < 4 ; ii++)
    {
      switch (ii)
	{
	case 0:
	  printf ("%% patients taking any # of pills (expect 100%%\n") ;
	  showTable10 (ll, ii, 999999) ;
	  break ;
	case 1:
	  printf ("aceto intake < %d (first quartile)\n", n1) ;
	  showTable10 (ll, ii, n1) ;
	  break ;
	case 2:
	  printf ("aceto intake < %d (median)\n", n2) ;
	  showTable10 (ll, ii, n2) ;
	  break ;
	case 3:
	  printf ("aceto intake < %d (third quartile)\n", n3) ;
	  showTable10 (ll, ii, n3) ;
	  break ;
	case 4:
	  printf ("%% patients with change in pain and aceto < 2 days/week\n") ;
	  printf ("je sais pas\n") ;
	  break ;
	case 5:
	  printf ("%% patients with change pain and aceto < 4 days/week\n") ;
	  printf ("je sais pas\n") ;
	  break ;
	}
    }
} /* getTable10 */

/**********************************************************************/

static void sGlobal (LL *ll, int phase)
{
  printf ("phase %d\n", phase) ;

  ll->ks = ac_dbquery_keyset (ll->db, "Find patient ", ll->h) ;
  ll->iter = ac_keyset_iter (ll->ks, TRUE, ll->h) ;
  ll->nP = ac_keyset_count (ll->ks) ;
  ll->a = getAllData (ll) ;

  printf ("Found %d patients\n", ll->nP) ;

  if (1) getTable1 (ll) ;
  if (1) getTable2 (ll) ;
  if (1) getTable3 (ll) ;
  if (1) getTable4 (ll) ;
  if (1) getTable5 (ll) ;
  if (1) getTable6 (ll) ;
  if (1) getTable7 (ll) ;
  if (1) getTable8 (ll) ;
  if (1) getTable8b (ll) ;
  if (1) getTable9 (ll) ;
  if (1) getTable10 (ll) ;
} /* sGlobal */

/**********************************************************************/
/**********************************************************************/

static void exportAll (LL *ll)
{
  AC_OBJ obj ;
  AC_TABLE tb ;
  int ir ;

  printf ("AN\tSex\tAge\tHeight\tWeight\tBMI") ;
  printf ("Trmt\t") ;
  printf ("\tJoint\tDuration\tStratum") ;
  printf ("\tPill_count\tPill_duration\tPills_per_day") ;
  printf ("\tGlobal\tFunction\tPain") ;
  printf ("\tAUC_Global\tAUC_Function\tAUC_Pain") ;
  printf ("\tGlobal_2\tGlobal_4\tGlobal_6") ;
  printf ("\tFunction_2\tFunction_4\tFunction_6") ;
  printf ("\tPain_2\tPain_4\tPain_6") ;
  printf ("\tCountry\n") ;

  ll->ks = ac_dbquery_keyset (ll->db, "Find patient ", ll->h) ;
  ll->iter = ac_keyset_iter (ll->ks, TRUE, ll->h) ;

  while ((obj = ac_next_obj (ll->iter)))
    {
      printf ("%s", ac_name(obj) + 3) ;
      printf ("\t%c", ac_has_tag (obj, "Man") ? 'M' : 'F') ;
      printf ("\t%d", ac_tag_int (obj, "Age", -1)) ;
      printf ("\t%.2f", ac_tag_float (obj, "Height", -1)) ;
      printf ("\t%.2f", ac_tag_float (obj, "Weight", -1)) ;
      printf ("\t%.2f", ac_tag_float (obj, "BMI", -1)) ;
      printf ("\t%d", ac_has_tag(obj, "Placebo") ? 1 : ac_tag_int (obj, "Active", -1)) ;
      printf ("\t%s", ac_has_tag(obj, "Knee") ? "Knee" : (ac_has_tag (obj, "Hip") ? "Hip" : "0")) ;
      printf ("\t%d", ac_tag_int (obj, "Duration", -1)) ;
      printf ("\t%s", ac_has_tag (obj, "NSAID_user") ? "NSAID" : (ac_has_tag (obj, "Acetaminophen_user") ? "Acetaminophen_user": "0")) ;
      printf ("\t%d\t%d\t%.2f", ac_tag_int (obj, "Pill_count", 0), ac_tag_int (obj, "Pill_duration", 0), ac_tag_float (obj, "Pills_per_day", 0)) ;
      printf ("\t%d\t%.2f\t%.2f", ac_tag_int (obj, "Global", 0), ac_tag_float (obj, "Function", 0), ac_tag_float (obj, "Pain", 0)) ;
      printf ("\t%.2f\t%.2f\t%.2f", ac_tag_float (obj, "AUC_Global", 0), ac_tag_float (obj, "AUC_Function", 0), ac_tag_float (obj, "AUC_Pain", 0)) ;

      tb = ac_tag_table (obj, "G_week", 0) ;
      for (ir = 0 ; ir < 3 ; ir++)
	{
	  int week = tb ? ac_table_int (tb, ir, 0, -1) : -1 ;
	  float x = tb ? ac_table_float (tb, ir, 1, -999) : -999;
	  if (week != 2 * (ir + 1))
	    x = -9999 ;
	  printf ("\t%.2f", x) ;
	}
      ac_free (tb) ;

      tb = ac_tag_table (obj, "F_week", 0) ;
      for (ir = 0 ; ir < 3 ; ir++)
	{
	  int week = tb ? ac_table_int (tb, ir, 0, -1) : -1 ;
	  float x = tb ? ac_table_float (tb, ir, 1, -999) : -999;
	  if (week != 2 * (ir + 1))
	    x = -9999 ;
	  printf ("\t%.2f", x) ;
	}
      ac_free (tb) ;

      tb = ac_tag_table (obj,"P_week", 0) ;
      for (ir = 0 ; ir < 3 ; ir++)
	{
	  int week = tb ? ac_table_int (tb, ir, 0, -1) : -1 ;
	  float x = tb ? ac_table_float (tb, ir, 1, -999) : -999;
	  if (week != 2 * (ir + 1))
	    x = -9999 ;
	  printf ("\t%.2f", x) ;
	}
      ac_free (tb) ;

      printf ("\t%s\n", ac_has_tag (obj, "US") ? "US" : (ac_has_tag (obj, "Multinational") ? "Multinational": "0")) ;
    }
} /* exportAll */

/**********************************************************************/
/**********************************************************************/

static void exportAUC (LL *ll)
{
  AC_OBJ obj ;
  AC_TABLE tbl ;
  float x, x0 ;
  int ir, week ;

  ll->ks = ac_dbquery_keyset (ll->db, "Find patient ", ll->h) ;
  ll->iter = ac_keyset_iter (ll->ks, TRUE, ll->h) ;

  while ((obj = ac_next_obj (ll->iter)))
    {
      printf ("Patient %s\n", ac_name(obj) ) ;
      x0 = ac_tag_float (obj, "Pain", -1) ;
      tbl = ac_tag_table (obj, "P_week", 0) ;
      for (x = 0, ir = 0 ; tbl && ir < tbl->rows ; ir++)
	{
	  week = ac_table_int (tbl, ir, 0, -1) ;
	  if (week == 2 || week == 4 || week == 6)
	    x += x0 + ac_table_float (tbl, ir, 1, -9999) ;
	}
      ac_free (tbl) ;
      if (x > 0)
	printf ("AUC_pain %f\n", x) ;


      x0 = ac_tag_float (obj, "Function", -1) ;
      tbl = ac_tag_table (obj, "F_week", 0) ;
      for (x = 0, ir = 0 ; tbl && ir < tbl->rows ; ir++)
	{
	  week = ac_table_int (tbl, ir, 0, -1) ;
	  if (week == 2 || week == 4 || week == 6)
	    x += x0 + ac_table_float (tbl, ir, 1, -9999) ;
	}
      ac_free (tbl) ;
      if (x > 0)
	printf ("AUC_function %f\n", x) ;


      x0 = ac_tag_int (obj, "Global", -1) ;
      tbl = ac_tag_table (obj, "G_week", 0) ;
      for (x = 0, ir = 0 ; tbl && ir < tbl->rows ; ir++)
	{
	  week = ac_table_int (tbl, ir, 0, -1) ;
	  if (week == 2 || week == 4 || week == 6)
	    x += x0 + ac_table_float (tbl, ir, 1, -9999) ;
	}
      ac_free (tbl) ;
      if (x > 0)
	printf ("AUC_global %f\n", x) ;

      printf ("\n") ;
      ac_free (obj) ;
    }
} /* exportAll */

/**********************************************************************/
/**********************************************************************/

static void helpUsage (void)
{
  printf ("// Usage :laure $ACEDB [-p phase \n") ;
  printf ("//                 =  1:global \n") ;
}

/**********************************************************************/

int main (int argc, const char **argv)
{
  int nn = 0 ;
  const char *dbName,  *phase ;
  const char *error = "" ;
  LL ll0 ;
  LL *ll = &ll0 ;

  freeinit () ;
  freeOutInit () ;

  phase = getArgV (&argc, argv, "-p") ;
  nn = phase ? atoi (phase) : 1 ;
  dbName = argc>=2 ? argv[1] : "." ;
  
  ll->h = ac_new_handle () ;
  ll->db = ac_open_db (dbName, &error) ;
  if (!ll->db)
    messcrash ("cannot open acedb database %s, %s sorry", dbName, error) ;

  switch (nn)
    {
    case 1: /* data gatherings */
      sGlobal (ll, nn) ;
      break ;
    case 2:
      exportAll (ll) ;
      break ;
    case 3:
      exportAUC (ll) ;
      break ;
    default: 
      helpUsage () ;
      messcrash ("bad phase %d\n", nn) ;
      break ;
    }

  ac_free (ll->h) ;

  printf ("// done\n") ;
  return 0 ;
}

/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
