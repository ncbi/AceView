
/^#/{ next ; }
{ 
   type = $3 ;
   split (type, aa, ":") ;
   split (aa[1], bb,  "/") ;

   split (type, aa, ":") ;
   ref = aa[2] ;
   var = aa[3] ;

   if (substr(type,1,4) == "Sub:")
   {
       x = $2 ; x1 = x - 1 ; x2 = x + 1 ; dx = 0 ;
       v = $1 ":" $2 ":" $3 ;
       tag = ref "2" var ;
   }
   else if (substr(type,1,4) == "Del_")
   {
       da = 0 + substr (bb[1],5) ;
       x = $2 ; x1 = x - 1 ; x2 = x + da ; dx = -da ;

       ref = substr(ref,1,1+da) ;
       var = substr(ref,1,1) ;
       motif = substr(ref,2) ;

       v = $1 ":" x1 ":Del:" ref ":" var ;
       if (da == 1)
	   tag = "Del" motif ;
       else
	   tag = "Multi_deletion\t" da "\t"  motif ;
   }
   else if (substr(type,1,4) == "Ins_")
   {
       da = 0 + substr (bb[1],5) ;
       x = $2 ; x1 = x - 1 ; x2 = x + da ; dx = da ;

       ref = substr(ref,1,1) ;
       var = substr(var,1,1+da) ;
       motif = substr(var,2) ;
       
       v = $1 ":" x1 ":Ins:" ref ":" var ;
       if (da == 1)
	   tag = "Ins" motif ;
       else
	   tag = "Multi_insertion\t" dx "\t"  motif ;
   }
   else
       next ;

   printf ("%s\t%s\t", v, run) ;
   printf ("tttiiiiiiittt") ;
   printf ("\t%d\t%d\t%d",x1,x2,dx) ;


   c = $9 ; m = $10 ; w = $11 ; 
   mp = $12 ; wp = $13 ; mm = $14 ; wm = $15 ;
   printf ("\t%d\t%d\t%d\t%d\t%d\t%d\t%d",c,m,w,mp,wp,mm,wm) ;
   printf ("\tBRS\t%s\n", tag) ;

}




