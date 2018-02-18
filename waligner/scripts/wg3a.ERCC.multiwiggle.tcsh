#!bin/tcsh -f

set ici=`pwd`

if (! -d tmp/WIGGLE) mkdir tmp/WIGGLE 
if (! -d tmp/WIGGLE/SpikeIn) mkdir tmp/WIGGLE/SpikeIn
foreach ERCC (ERCC1 ERCC2)
  set out=tmp/WIGGLE/SpikeIn/multi_wiggle.$ERCC.txt
  echo "#" > $out.1
  foreach run (`cat  tmp/WIGGLE/SpikeIn/$ERCC.runs.list | gawk '/^Run/{gsub(/\"/,"",$0);print $2;}'`)

    if (-e  tmp/WIGGLERUN/$run/R.chrom.SpikeIn.u.f.BF.gz) then
      echo "RUN $run" >> $out.1
      foreach fr (f r)
        echo "STRAND $fr" >> $out.1
        gunzip -c  tmp/WIGGLERUN/$run/R.chrom.SpikeIn.u.$fr.BF.gz  >> $out.1
      end
    endif
  end

  date > $out
  echo "Multi wiggle for the $ERCC SpikeIn in project $MAGIC" >> $out 

  cat $out.1 | gawk '/^RUN/{nr++;run[nr]=$2;next;}/^STRAND/{strand=$2;next;}/^#/{next;}/^track/{next}/^fixedStep/{split($2,aa,"=");target=aa[2];split($3,aa,"=");x0=aa[2];split($4,aa,"=");step=aa[2];x0 -= step; tt=t2i[target];if(tt<1){nt++;i2t[nt]=target;t2i[target]=nt;tt=nt;};next;}{x0+=step;yy[nr,tt,strand,x0]=$1;if(xmax[tt]<x0);xmax[tt]=x0;nnn[tt]+=$1;}END{for (tt=1;tt<=nt;tt++)nnn[tt]/=(1+xmax[tt]);printf("Target\tPosition");for (r=1;r<=nr;r++)printf("\tf:%s\tr:%s",run[r],run[r]);for (t1=1;t1<=nt;t1++){n=0;for (t2=1;t2<=nt;t2++){if(nnn[t2]>n){n=nnn[t2];tt=t2;}}nnn[tt]=0;if(0)printf("t1=%d tt=%d n=%d max=%d %s\n",t1,tt,n,xmax[tt],i2t[tt]);for (x=0;x<=xmax[tt];x+=10){if(x==0)kk=5;else kk=1;for(k=0;k<kk;k++){printf("\n%s\t%d",i2t[tt],x);for (r=1;r<=nr;r++)printf("\t%d\t%d",yy[r,tt,"f",x],yy[r,tt,"r",x]);}}}printf("\n");}' >> $out

 echo $out
 \cp $out RESULTS/$MAGIC.multi_wiggle.$ERCC.txt

  set out1=tmp/WIGGLE/SpikeIn/multi_wiggle.$ERCC.normed.txt
  date > $out1
  echo "Multi wiggle for the $ERCC SpikeIn in project $MAGIC" >> $out1 

  cat $out.1 | gawk '/^RUN/{nr++;run[nr]=$2;next;}/^STRAND/{strand=$2;next;}/^#/{next;}/^track/{next}/^fixedStep/{split($2,aa,"=");target=aa[2];split($3,aa,"=");x0=aa[2];split($4,aa,"=");step=aa[2];x0 -= step; tt=t2i[target];if(tt<1){nt++;i2t[nt]=target;t2i[target]=nt;tt=nt;};next;}{x0+=step;yy[nr,tt,strand,x0]=$1;if(xmax[tt]<x0);xmax[tt]=x0;nnn[tt]+=$1;nnnr[nr,tt]+=$1;}END{for (tt=1;tt<=nt;tt++)nnn[tt]/=(1+xmax[tt]);printf("Target\tPosition");for (r=1;r<=nr;r++)printf("\tf:%s\tr:%s",run[r],run[r]);for (t1=1;t1<=nt;t1++){n=0;for (t2=1;t2<=nt;t2++){if(nnn[t2]>n){n=nnn[t2];tt=t2;}}nnn[tt]=0;if(0)printf("t1=%d tt=%d n=%d max=%d %s\n",t1,tt,n,xmax[tt],i2t[tt]);for (x=0;x<=xmax[tt];x+=10){if(x==0)kk=5;else kk=1;for(k=0;k<kk;k++){printf("\n%s\t%d",i2t[tt],x);for (r=1;r<=nr;r++)printf("\t%d\t%d",yy[r,tt,"f",x],1000*yy[r,tt,"r",x]/(1+nnnr[r,tt]));}}}printf("\n");}' >> $out1

 echo $out1
 \cp $out1 RESULTS/$MAGIC.multi_wiggle.$ERCC.txt

end

