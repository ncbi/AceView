#!bin/tcsh -f

set run=$1
set lane=$2

if ($2 == "") goto c2


if (! -e tmp/COUNT/$lane.err.hack,done) then

  mv tmp/COUNT/$lane.errorProfile.gz tmp/COUNT/$lane.errorProfile.gz_old
  bin/bestali -i tmp/COUNT/$lane.hits.gz -errorProfile -strategy RNA_seq -o  tmp/COUNT/$lane.new -gzo
  mv  tmp/COUNT/$lane.new.errorProfile.gz  tmp/COUNT/$lane.errorProfile.gz
  \rm  tmp/COUNT/$lane.new.*
  touch tmp/COUNT/$lane.err.hack,done

endif
exit 0

c2:

echo hello $run
set MAGIC_COUNT_DIR=COUNT

if (! -e tmp/COUNT/$run/c2hack.err.ace3) then

  set toto=tmp/COUNT/$run/c2hack.err
  echo "Ali $run" > $toto.ace 
  echo "-D Error_position" >> $toto.ace 
  echo "-D Ambiguous_position" >> $toto.ace 
  echo "-D First_base_aligned" >> $toto.ace 
  echo "-D Last_base_aligned" >> $toto.ace 
  echo "-D Error_profile" >> $toto.ace  
  echo "-D Error_spike" >> $toto.ace  



echo toto >  $toto.pos
echo toto >  $toto.amb
echo toto > $toto.FB
echo toto > $toto.LB
echo toto > $toto.err_spike

      gunzip -c  tmp/$MAGIC_COUNT_DIR/$run/*.errorProfile.gz   | grep RPOS >> $toto.err_spike
      gunzip -c tmp/$MAGIC_COUNT_DIR/$run/*.errorProfile.gz |  gawk -F '\t' '/^POS/{prefix=$2;x=$3;prefixes[prefix]=1;xx[x]=1;p[$2,x]+=$4;}END{for(prefix in prefixes)for(x in xx)if(0+p[prefix,x]>0 && x+0<=1000)printf("Error_position  %s %d %d \n", prefix, x, p[prefix,x]) ;}'  | sort -k 1,1 -k 2,2 -k 3,3n >> $toto.pos
      gunzip -c tmp/$MAGIC_COUNT_DIR/$run/*.errorProfile.gz |  gawk -F '\t' '/^FB/{prefix=$2;x=$3;prefixes[prefix]=1;xx[x]=1;p[$2,x]+=$4;}END{for(prefix in prefixes)for(x in xx)if(0+p[prefix,x]>0 && x+0<=1000)printf("FB %s %d %d \n", prefix, x, p[prefix,x]) ;}'  | sort -k 1,1 -k 2,2 -k 3,3n >> $toto.FB
      gunzip -c tmp/$MAGIC_COUNT_DIR/$run/*.errorProfile.gz |  gawk -F '\t' '/^LB/{prefix=$2;x=$3;prefixes[prefix]=1;xx[x]=1;p[$2,x]+=$4;}END{for(prefix in prefixes)for(x in xx)if(0+p[prefix,x]>0 && x+0<=1000)printf("LB %s %d %d \n", prefix, x, p[prefix,x]) ;}'  | sort -k 1,1 -k 2,2 -k 3,3n >> $toto.LB
      
gunzip -c tmp/$MAGIC_COUNT_DIR/$run/*.errorProfile.gz |  gawk -F '\t' '/^NNPOS/{prefix=$2;x=$3;prefixes[prefix]=1;xx[x]=1;p[$2,x]+=$4;}END{for(prefix in prefixes)for(x in xx)if(0+p[prefix,x]>0 && x+0<=1000)printf("Ambiguous_position %s  %s %d \n", prefix, x, p[prefix,x]) ;}' >>  $toto.amb

gunzip -c tmp/$MAGIC_COUNT_DIR/$run/*.errorProfile.gz | gawk -F '\t' '/^ERRALI/{if(NF==4){s=$2;nkb[s]+=$3;}if(NF==3){s="f1";nkb[s]+=$2;}next;}/^ERRS/{s=$2;ss[s]=1;t=$3;if(t>tmax)tmax=t;tt[t]=$4;ns[s,t]+=$5;nns[s]+=$5;next;}/^ERR/{s="f1";ss[s]=1;t=$2;if(t>tmax)tmax=t;tt[t]=$3;ns[s,t]+=$4;nns[s]+=$4;next;}END{tmax++;tt[tmax]="Any";for(s in ss){ns[s,tmax]=nns[s];for(t=1;t<=tmax;t++)if(nkb[s]>0)if(t==tmax || ns[s,t]>0)printf("Error_profile %s %s %d %d %.1f\n", s, tt[t],ns[s,t],nkb[s]/1000,1000* ns[s,t]/(nkb[s]));}}' >> $toto.ace


cat $toto.err_spike | gawk '/RPOS/{ z = $2 " " $3 " \"" $6 "\"" ; nt[z] += $5 ; ne[z] += $7 ;}END {for (z in nt)if(ne[z]>0 && 3*ne[z] > nt[z])printf ("Error_spike %s %d %d\n", z, ne[z], nt[z]);  }'  >>   $toto.ace

cat $toto.pos | sort  -k 1,1 -k 2,2 -k 3,3n | gawk '/^Error_position/{nf[$2]++;nn[$3,$2]=$4;if($3>mx)mx=$3;}END{for (x=1;x<=mx;x++){printf("Error_position %d", x); if(nf["f"]>0)printf(" %d",nn[x,"f"]); if(nf["f1"]>0)printf(" %d",nn[x,"f1"]); if(nf["f2"]>0)printf(" %d",nn[x,"f2"]); if(nf["f3"]>0)printf(" %d",nn[x,"f3"]);printf("\n");}}' >>   $toto.ace

cat $toto.FB | sort  -k 1,1 -k 2,2 -k 3,3n | gawk '/^FB/{nf[$2]++;nx[$3]++;nn[$3,$2]=$4;if($3>mx)mx=$3;}END{for (x=1;x<=mx;x++){if (nx[x]>0){printf("First_base_aligned %d", x); if(nf["f"]>0)printf(" %d",nn[x,"f"]); if(nf["f1"]>0)printf(" %d",nn[x,"f1"]); if(nf["f2"]>0)printf(" %d",nn[x,"f2"]); if(nf["f3"]>0)printf(" %d",nn[x,"f3"]);printf("\n");}}}' >>   $toto.ace

cat $toto.LB | sort  -k 1,1 -k 2,2 -k 3,3n | gawk '/^LB/{nf[$2]++;nx[$3]++;nn[$3,$2]=$4;if($3>mx)mx=$3;}END{for (x=1;x<=mx;x++){if (nx[x]>0){printf("Last_base_aligned %d", x); if(nf["f"]>0)printf(" %d",nn[x,"f"]); if(nf["f1"]>0)printf(" %d",nn[x,"f1"]); if(nf["f2"]>0)printf(" %d",nn[x,"f2"]); if(nf["f3"]>0)printf(" %d",nn[x,"f3"]);printf("\n");}}}' >>   $toto.ace

cat $toto.amb | gawk '/^Ambiguous_position/{nf[$2]++;nx[$3]++;nn[$2,$3]=$4;if($3>mx)mx=$3;}END{for (x=1;x<=mx;x++){if (nx[x]>0){printf("Ambiguous_position %d", x); if(nf["f"]>0)printf(" %d",nn["f",x]); if(nf["f1"]>0)printf(" %d",nn["f1",x]); if(nf["f2"]>0)printf(" %d",nn["f2",x]); if(nf["f3"]>0)printf(" %d",nn["f3",x]);printf("\n");}}}' >>   $toto.ace

echo >> $toto.ace

  touch tmp/COUNT/$lane.err.hack.done

endif

exit 0
