
/^ZZZZZ/ { zz++ ; next ; }

{ if (zz + 0 < 1) {r2OS[$1] = $2; next ; }}


{ if (zz == 1)
  { # register High-MYCN status of all runs
    z = substr ($1,6) ;
    gsub(/\"/,"",z) ; 
    HM[z] = 1 ;
    if (0) printf( "z=%s\n",z);
    next ;
  }
}

/^#Sample/ { 
  for (i = 3 ; i <= NF ; i++)
    {
      if (justLowMYCN == 0 || HM[r[i]] < 1) # register the sample name of each run
        {
	  x = 1 ;    # disregard sex misclassified samples
	  if (substr($i,1,1) == "X" || substr($i,1,1) == "Y") x = 2 ;
	  s[i] = substr($i,x) ;
	}
    }
next;
}

/^#Run\t/ { # register the run name of each column
  for (i = 3 ; i <= NF ; i++)
    {
      r[i] = $i ;
    } 
  next;
}

/^#Title\t/ { # register the run name of each column, and the parity of the title
  for (i = 3 ; i <= NF ; i++)
    {
      t[i] = $i ;
      z = substr(t[i], length(t[i])) ; z1 = z ;
# print "XXXX",$i,z
      if (z == "a" || z == "b")
        z = substr(t[i], length(t[i])-1,1) ;
# print "XXXXY",z
      pp[i] = (z+0) % 2 ; if (pp[i] == 0) pp[i] = 2 ;
	if (z1 == "a1" || z1 == "b") pp[i] = -1 ;
    } 
  next;
}
{iter = 0 ;}
/^Ratio:0/ {iter=0;}
/^Ratio:1/ {iter=1;}
/^Ratio:[0-9] Favorable_NB_178/        {for(i=3;i<=NF;i++)ff0[i,iter]=-$i;nf=NF;next;}
/^Ratio:[0-9] Favorable_TrainingSet_87/   {for(i=3;i<=NF;i++)fft[i,iter]=-$i;nf=NF;next;}
/^Ratio:[0-9] EventFree_TrainingSet/    {for(i=3;i<=NF;i++)ffef[i,iter]=-$i;nf=NF;next;}
/^Ratio:[0-9] OS_Died_TrainingSet/    {for(i=3;i<=NF;i++)ffos[i,iter]=$i;nf=NF;next;}
/^Ratio:[0-9] HR_Died_TrainingSet_43/  {for(i=3;i<=NF;i++)ffhr[i,iter]=$i;nf=NF;next;}
/^Ratio:[0-9] MYCN_high/        {for(i=3;i<=NF;i++)ff1[i,iter]=$i;nf=NF;next;}
/^Ratio:[0-9] MYCN2_high/       {for(i=3;i<=NF;i++)ff2[i,iter]=$i;nf=NF;next;}
/^Ratio:[0-9] MYCN3_high/       {for(i=3;i<=NF;i++)ff3[i,iter]=$i;nf=NF;next;}
/^Ratio:[0-9] SG1_/             {for(i=3;i<=NF;i++)sg1[i,iter]=$i;nf=NF;next;}
/^Ratio:[0-9] SG2_/             {for(i=3;i<=NF;i++)sg2[i,iter]=$i;nf=NF;next;}
/^Ratio:[0-9] SG3_/             {for(i=3;i<=NF;i++)sg3[i,iter]=$i;nf=NF;next;}
/^Ratio:[0-9] HR_EF1/           {for(i=3;i<=NF;i++)hr_ef12[i,iter]-=$i;nf=NF;next;}
/^Ratio:[0-9] HR_EF2/           {for(i=3;i<=NF;i++)hr_ef12[i,iter]-=$i;nf=NF;next;}
/^Ratio:[0-9] HR_EF/            {for(i=3;i<=NF;i++)hr_ef[i,iter]-=$i;nf=NF;next;}
/^Ratio:[0-9] OS1_/             {for(i=3;i<=NF;i++)os12[i,iter]+=$i;nf=NF;for(i=3;i<=NF;i++)if(r2OS[r[i]]=="OSx2")os1or2[i,iter]+=$i;next;}
/^Ratio:[0-9] OS2_/             {for(i=3;i<=NF;i++)os12[i,iter]+=$i;nf=NF;for(i=3;i<=NF;i++)if(r2OS[r[i]]=="OS1")os1or2[i,iter]+=$i;next;}
/^Ratio:[0-9] HR1_30_die\/HR1_30_surv/ {for(i=3;i<=NF;i++)hr_os123[i,iter]+=$i;nf=NF;next;}
/^Ratio:[0-9] HR1_30_die\/HR2_30_surv/ {for(i=3;i<=NF;i++)hr_os123[i,iter]+=$i;nf=NF;next;}
/^Ratio:[0-9] HR1_30_die\/HR3_30_surv/ {for(i=3;i<=NF;i++)hr_os123[i,iter]+=$i;nf=NF;next;}
/^Ratio:[0-9] HR1_30_surv\/HR2_30_die/ {for(i=3;i<=NF;i++)hr_os123[i,iter]-=$i;nf=NF;next;}
/^Ratio:[0-9] HR1_30_surv\/HR3_30_die/ {for(i=3;i<=NF;i++)hr_os123[i,iter]-=$i;nf=NF;next;}
/^Ratio:[0-9] HR2_30_die\/HR2_30_surv/ {for(i=3;i<=NF;i++)hr_os123[i,iter]+=$i;nf=NF;next;}
/^Ratio:[0-9] HR2_30_die\/HR3_30_surv/ {for(i=3;i<=NF;i++)hr_os123[i,iter]+=$i;nf=NF;next;}
/^Ratio:[0-9] HR2_30_surv\/HR3_30_die/ {for(i=3;i<=NF;i++)hr_os123[i,iter]-=$i;nf=NF;next;}
/^Ratio:[0-9] Gender_AnnotAs/               {for(i=3;i<=NF;i++)female[i,iter]-=$i;nf=NF;next;}

END {
if (parity < 1) printf ("# The ROC curves are evaluated only against all samples (test set + training set)\n") ;
if (parity == 1) printf ("# The ROC curves are evaluated only against the odd samples (training set)\n") ;
if (parity == 2) printf ("# The ROC curves are evaluated only against the even samples (test set)\n") ;
  printf ("# Run\tSample\tTitle") ;

     printf ("\tFemina\tTraining f/u\tTraining Event free\tTraining Overall survival\tMYCN april\tMYCN april iterated\tMYCN april 2 iter\tHR EF\tHR EF stratified\tHigh Risk\tHR OS stratified") ;

     printf ("\tFemina iterated\tTraining f/u iterated\tTraining Event free iterated\tTraining Overall survival iterated\tMYCN april iterated\tMYCN april iterated iterated\tMYCN april 2 iter iterated\tHR EF iterated\tHR EF stratified iterated\tHigh Risk iterated\tHR OS stratified  iterated") ;

  printf ("\tf/u\tEvent free\tOverall Survival\tMale#/Female-\n") ;

  for (i = 3 ; i <= nf ; i++)
    {
      if (substr (r[i], 1, 3) != "Rhs") 
        continue ;
      if (parity > 0 && parity != pp[i]) 
        continue ;
      if (pp[i] == -1) continue ;
      uu1 = "" ;      uu2 = "" ;      uu3 = "" ;
      if (1)   # u/f outcome
        { 
           if (substr (s[i],1,1) == "f") uu1 = "###" ;
           if (substr (s[i],1,1) == "u") uu1 = "---" ;
        }
      if (2)   # event free outcome
        { 
           x = 0 ;
           if (substr (s[i], 1 , 1) == "f" || substr(s[i],1,1) == "u") x = 1 ;
           if (substr (s[i], 1+x, 1) == "0") uu2 = "###" ;
           if (substr (s[i], 1+x, 1) == "1") uu2 = "---" ;
        }
      if (3)   # overall survival outcome
        { 
           x = 0 ;
           if (substr (s[i], 1 , 1) == "f" || substr(s[i],1,1) == "u") x = 1 ;
           if (substr (s[i], 2+x, 1) == "0") uu3 = "###" ;
           if (substr (s[i], 2+x, 1) == "1") uu3 = "---" ;
        }
      if (4)   # sex outcome
        { 
           x = 0 ;
           if (index(s[i],"_M_") > 0) uu4 = "###" ;
           if (index(s[i],"_F_") > 0) uu4 = "---" ;
        }
      osinternal = "" ;
      if (r2OS[r[i]]=="OS1") osinternal = "OSINTERNAL" ;
      printf ("%s\t%s\t%s", r[i],t[i],s[i]) ;

      for (iter = 0 ; iter < 2 ; iter++)
        printf ("\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f",
              female[i,iter],fft[i,iter], ffef[i,iter], ffos[i,iter], ff1[i,iter], ff2[i,iter], ff3[i,iter], hr_ef[i,iter], hr_ef12[i,iter], ffhr[i,iter], hr_os123[i,iter]) ;

      printf ("\t%s\t%s\t%s\t%s\t%s\n",
               uu1, uu2, uu3, uu4, osinternal);}
}
