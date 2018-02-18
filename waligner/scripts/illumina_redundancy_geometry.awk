{ 
  for (i=2 ; i<=NF ; i++)
    {
      split ($i, aa, ":") ; c[i]=aa[5] ; x[i]=aa[6] ; y[i]=aa[7] ; if  (0)print "cxy", i, c[i], x[i], y[i],  $1, aa[2] ; 
    }
  for  (i=2 ; i<NF ; i++)
   {
     for (j=i+1 ; j<=NF ; j++)
       {
         if (c[i]==c[j])
          {
            dx=x[i]-x[j] ; dy=y[i]-y[j] ; 
            if (dx<0)dx=-dx ; 
            if (dy<0)dy=-dy ; 
            if (dxmax<dx)dxmax=dx ; 
            if (dymax<dy)dymax=dy ; 
            r=int (.5+sqrt (dx*dx+dy*dy)) ; 
            if (rmax<r)rmax=r ; 
            if (r>30)r=30 ; if (dx>30)dx=30  ; if (dy>30)dy=30 ; 
            nnn++ ; nn[dx,dy]++ ; 
            nx[dx]++ ; ny[dy]++ ; 
            rrt++ ; rr[r]++ ; 

            if (0) print i, j, y[i], y[j], dx, dy ; 
          }
       }
   }
}
END {
  printf ("Lane %s Number of pairs of identical reads at given radial distance,  saturated at 30\n", lane) ; 
  for (i=1 ; i<=30 ; i++)
    printf ("\t%d", i) ; 
  printf ("\n%s", lane) ; 
  for (i=1 ; i<=30 ; i++)
    printf ("\t%d", rr[i]) ; 

  printf ("\t%d\t%g", rrt-rr[30], 100*rrt/ (rr[30]+rrt)) ;
  printf ("\n\nLane %s Number of pairs of identical reads at given  (x,y) distance, saturated at 30, in reality dxmax=%d  dymax=%d\n", lane, dxmax, dymax) ; 
  for (i=1 ; i<=30 ; i++)
    printf ("\t%d", i) ; printf ("\tTotal") ; 
  for (i=1 ; i<=30 ; i++)
    {
      printf ("\n%d", i) ; 
      for (j=1 ; j<=30 ; j++)
        printf ("\t%d", nn[i,j]) ; 
      printf ("\t%d", nx[i]) ; 
    }
  printf ("\nTotal") ; 
  for (j=1 ; j<=30 ; j++)
    printf ("\t%d", ny[j]) ; 
  printf ("\n")  ; 
}
