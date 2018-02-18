/^ZZZZZ/{pp++; next; }
/Title_line/ {
   printf("Method") ;
   printf("\tSex rained on gender\tSex iter\tSex V\tSex iter V") ;
   printf("\tu/f trained on u/f\tu/f on u/f iter\tu/f trained on u/f V\tu/f on u/f iter V") ;
   printf("\tu/f trained on OS\tu/f on OS iter\tu/f trained on OS V\tu/f on OS iter V") ;

   printf("\tEF trained on EF\tEF on EF iter\tEF trained on EF V\tEF on EF iter V") ;
   printf("\tEF trained on OS\tEF on OS iter\tEF trained on OS V\tEF on OS iter V") ;

   printf("\tOS trained on OS\tOS on OS iter\tOS trained on OS V\tOS on OS iter V") ;

   printf("\tOS among HR\tOS among HR iter") ;
   printf("\tEF among HR\tEF among HR iter") ;
   printf("\n") ;
   pp = 99 ;
   next ;
 } 
{ pp = pp + 0 ; auc = "AUC=" ; mcc = "MCC=" ;}
{ iter = 0 ;}
/iterated/ { iter = 1 ;}
/Outcome Sex\tClassifier: Femina/                    { i=index($4,auc); z[1,pp,iter]= 0 + substr($4,i+4); next;}
/Outcome u\/f\tClassifier: Training f\/u/            { i=index($4,auc); z[2,pp,iter]= 0 + substr($4,i+4); next;}
/Outcome u\/f\tClassifier: Training Overall survival/{ i=index($4,auc); z[3,pp,iter]= 0 + substr($4,i+4); next;}
/Outcome EF\tClassifier: Training Event free/        { i=index($4,auc); z[4,pp,iter]= 0 + substr($4,i+4); next;}
/Outcome EF\tClassifier: Training Overall survival/  { i=index($4,auc); z[5,pp,iter]= 0 + substr($4,i+4); next;}


/Outcome OS\tClassifier: Training Overall survival/  { i=index($4,auc); z[6,pp,iter]= 0 + substr($4,i+4); next;}


/Outcome OS\tClassifier: High Risk/                  { i=index($0,auc); if (pp < 1)z[7,0,iter]= 0 + substr($0,i+4); next;}
/Outcome EF\tClassifier: HR EF/                      { i=index($0,auc); if (pp < 1)z[7,0,iter]= 0 + substr($0,i+4); next;}



/sex.*clf=Femina/                                             { i=index($1,mcc);zmcc[1,pp,iter]= substr($1,i+4); next;}
/favorable versus unfavorable.*clf=Training f\/u/             { i=index($1,mcc);zmcc[2,pp,iter]= substr($1,i+4); next;}
/favorable versus unfavorable.*clf=Training Overall survival/ { i=index($1,mcc);zmcc[3,pp,iter]= substr($1,i+4); next;}
/event free survival clf=Training Event free/                           { i=index($1,mcc);zmcc[4,pp,iter]= substr($1,i+4); next;}
/event free survival clf=Training Overall survival/                     { i=index($1,mcc);zmcc[5,pp,iter]= substr($1,i+4); next;}
/overal survival clf=Training Overall survival/                         { i=index($1,mcc);zmcc[6,pp,iter]= substr($1,i+4); next;}
/event free survival clf=High Risk/                                     { i=index($1,mcc);zmcc[7,pp,iter]= substr($1,i+4); next;}
/overal survival clf=High Risk /                                        { i=index($1,mcc);zmcc[8,pp,iter]= substr($1,i+4); next;}

END {
  if (pp < 99)
    {
      printf ("%s", m) ;
      for (i =1 ; i <= 7 ; i++)
	for (pp = 0 ; pp < 2 ; pp++)
	  for (iter = 0 ; iter < 2 ; iter++)
	    if (mm == "MCC") printf("\t%.2f", zmcc[i,pp,iter]); 
            else printf("\t%.2f", z[i,pp,iter]); 
      printf ("\n") ;
    }
}



