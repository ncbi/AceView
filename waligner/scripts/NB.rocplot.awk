/^#/ {next ;}
{ c = kk ; # column with relevant phenotype
  if ($c == "###") 
    { w = 1 ; f++; nn++;a[nn]=1; alpha[nn] = $ii ; } 
  if ($c == "---") 
    { w = 0; u++ ; s += f; u2f[u] = f;  nn++;a[nn]=-1; alpha[nn] = $ii ;}  # ii is the column with index
  # if (nn > 90 && nn < 96) { print "MCC",nn,a[nn], $1,$2,$7,$20;}
}
END {
  if (w == 1)
    {  # register the last value
      u++ ; u2f[u] = f ; 
    }
   beta = s ;
   if (u*f > 0) beta = s/(u*f) ; 
   printf("s=%d uf=%d u=%d f=%d AUC=%.2f\t",s,u*f,u,f,100*beta) ;
   printf("\nFalse positive counts\tFalse positive counts") ; 
   for (i = 0 ; i <= u ; i++)
     {
       printf ("\t%d", i) ;
     }
   printf ("\nTrue positive counts\t%s",title) ;
   for (i = 0 ; i <= u; i++)
     {
       printf("\t%d", u2f[i]);
     }
   printf("\nFalse positive rate\tFalse positive rate") ; 
   for (x = 0 ; x <= 100; x++)
     printf ("\t%.2f", x/100) ;
   printf("\nTrue positive rate\t%s",title) ;
   u2f[u+1] = u2f[u] ;
   for (x = 0 ; x <= 100 ; x++)
     {
       i = int(u*x/100) + 1 ;
       if (x == 100)
	 y = 1 ;
       else
	 {
	   y1 = u2f[i] ;
	   dy = u2f[i+1] - y1 ;
	   y = y1 ;
	   dx = u * x / 100 - i + 1 ;
	   if (dx>0)
	     y = y + dx * dy / 100;
	   if (i < -2)
	     printf("ERROR dx=%f x=%d u=%d ux=%d i=%d dy=%.2f y1=%.2f y=%.2f\n",dx,x,u,u*x,i,dy,y1,y);
	   if (u2f[u] > 0) y = y / u2f[u];
	 }
       printf ("\t%.3f",y);
     }
   printf("\n");

# look for the threshold with optimal MCC Matthews correlation coefficient
   bestjj = 0 ; bestmcc = 0 ;
   for (jj = 1 ; jj < nn ; jj++)
     {
# TP all positive below jj
# FP all positive above jj
# TN all positive above jj
# FN all positive below jj
       tp = 0 ; fp = 0 ; tn = 0 ; fn = 0 ;
       for (i = 0 ; i < jj ; i++)
	 {
	   if (a[i] == 1) tp++ ;
	   if (a[i] == -1) fp++ ;
	 }
       for (i = jj ; i <= nn ; i++)
	 {
	   if (a[i] == 1) fn++ ;
	   if (a[i] == -1) tn++ ;
	 }
       ttp[jj] = tp ;
       ttn[jj] = tn ;
       tfp[jj] = fp ;
       tfn[jj] = fn ;
       sn[jj] = 0 ; if (tp>0)   sn[jj] = tp/(tp+fn) ; # sensitivity or TPR true positive rate or Recall or Hit-rate
       sp[jj] = 0 ; if (tn>0)   sp[jj] = tn/(tn+fp) ; # specifity

# Youden's J statistics = intercept of ROC curve with anti diagonal  = sn + sp -1  # not recommended
# Gini coefficint = 2 * AUC - 1
       ppv[jj] = 0 ;  if (tp + fp > 0) tp/(tp + fp) ;   # positive predictive value or Precision
       npv[jj] = 0 ;  if (tn + fn > 0) tn/(tn + fn) ;   # negative predictive value
       fdr[jj] = 0 ;  if (fp + tp > 0) fp/(fp + tp) ;   # FDR false discovery rate
       fpr[jj] = 0 ;  if (fp + tn > 0) fp (fp + tn) ;   # false positive rate or fall-out

       accu[jj] = 0 ; if (tn + tp > 0) accu[jj] = (tn + tp)/(tn+tp+fn+fp) ;
       mcc[jj] = 0 ; if ((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn) > 0) mcc[jj] = (tp*tn - fp*fn)/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)) ;
       if(mcc[jj]<0)mcc[jj] = -mcc[jj] ;
       if (mcc[jj] > bestmcc) { bestjj = jj ; bestmcc = mcc[jj] ; }
     }
   jj = bestjj ;
   printf ("MCC sort on %d %s clf=%s Bestj=%d alpha=%.3f AUC=%.2f TP=%d TN=%d FP=%d FN=%d SENS=%.2f SP=%.2f ACCURACY=%.2f MCC=%.2f\n", kk, title, clf,jj, (alpha[jj]+ alpha[jj-1])/2,100*beta,ttp[jj],ttn[jj], tfp[jj], tfn[jj],100*sn[jj], 100*sp[jj], 100*accu[jj], 100*mcc[jj]);
   if (0)
     {
       printf ("MCC\tThreshold\tAlpha\tROC\tTP\tFP\tFN\tFP\tSensibility\tSpecificity\tAccuracy\tMCC\tSP+SN\tSP*SN\n") ;
       for (jj=1;jj<nn;jj++)
	 printf("MCC\t%d\t%.3f\t%d\t%d\t%d\t%d\t%d\t%.2f\t%.2f\t%.2f\t%.6f\t%.6f\t%.2f\n", jj, alpha[jj],a[jj],ttp[jj],tfp[jj],tfn[jj],ttn[jj],sn[jj], sp[jj], accu[jj],mcc[jj], sn[jj]+sp[jj],sn[jj]*sp[jj]);
     }
}

