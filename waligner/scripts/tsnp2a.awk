
/^#/{ next ; }
{ 
   if (substr($3,1,4) == "Sub:")
   {
       x = $2 ; x1 = x - 1 ; x2 = x + 1 ; dx = 0 ;
       c = $9 ; m = $10 ; w = $11 ; 
       mp = $12 ; wp = $13 ; mm = $14 ; wm = $15 ;
       v = $1 ":" $2 ":" $3 ;


       split ($3,aa,":") ; tag = aa[2] "2" aa[3] ;
   }
   else
       next ;
   printf ("%s\t%s\t", v, run) ;
   printf ("tttiiiiiiittt") ;
   printf ("\t%d\t%d\t%d",x1,x2,dx) ;
   printf ("\t%d\t%d\t%d\t%d\t%d\t%d\t%d",c,m,w,mp,wp,mm,wm) ;
   printf ("\tBRS\t%s\t%d\t%s\n", tag, dx, insert) ;

}




