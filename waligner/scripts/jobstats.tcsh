#!bin/tcsh -f

set run=$1
set lane=$2
set inside=$3

# sync targets2 between MAGIC and jobstats.tcsh

setenv targets2 "any genome seqc av RefSeq UCSC mito rnaGene rrna SpikeIn DNASpikeIn gdecoy EBI smallRNA"
setenv targets2 "any $RNAtargets $DNAtargets"

set toto=tmp/COUNT/$lane.jobstats
set totopreace=tmp/COUNT/$lane.jobstats.preace
if (-e $toto) \rm $toto
if (-e $totopreace) \rm $totopreace

    printf "jobstats c1: $lane\n"
      printf "$lane" > $toto

        foreach target ($targets2)
          if ($target == any) continue

          set ff=tmp/PHITS_$target/$lane.err
          if (-e $ff) then
            gawk 'BEGIN{z=-1}/^Command exited with non-zero status/{z=$6 ;}/^user/{if(z==-1)z=0;}/^Status/{n++;z=$4;}END{printf("\t%d",z);}' $ff >> $toto
          else
            printf "\t-1" >> $toto
          endif

        end

          set ff=tmp/COUNT/$lane.err
          if (-e $ff) then
             gawk 'BEGIN{z=-1}/^Command exited with non-zero status/{z=$6 ;}/^user/{if(z==-1)z=0;}/Status/{n++;z=$4;}END{printf("\t%d",z);}' $ff >> $toto
          else
            printf "\t-1" >> $toto
          endif

          set ff=tmp/COUNT/$lane.a123.err
          if (-e  $ff) then      
             gawk 'BEGIN{z=-1}/^Command exited with non-zero status/{z=$6 ;}/^user/{if(z==-1)z=0;}/Status/{n++;z=$4;}END{printf("\t%d",z);}' $ff >> $toto
          else
            if ($inside == 1) then
              printf "\t0" >> $toto
            else 
              printf "\t-1" >> $toto
            endif
          endif

        printf "\t\t$lane" >> $toto

        set total=0
        foreach target ($targets2)
          if ($target == any) continue

          set ff=tmp/PHITS_$target/$lane.err
          if (-e $ff) then
            set cpu=`cat $ff | gawk '/^user/{n++;s+=$2;}/^sys/{n++;s+=$2;}END{if(n==0)s=-10;printf("%d",s);}' $ff`
            printf "\t$cpu" >> $toto
            if ($cpu != -10) @ total = $total + $cpu
            if ($cpu != -10) printf "CPU\t$run\t$lane\t$target\t$cpu\n" >> $totopreace
          else
            printf "\t-10" >> $toto
          endif

        end

          set ff=tmp/COUNT/$lane.err
          if (-e $ff) then
            set cpu=`cat  $ff | gawk '/^user/{n++;s+=$2;}/^sys/{n++;s+=$2;}END{if(n==0)s=-10;printf("%d",s);}' $ff`
            printf "\t$cpu" >> $toto
            if ($cpu != -10) @ total = $total + $cpu
            if ($cpu != -10) printf "CPU\t$run\t$lane\tSelectBest\t$cpu\n" >> $totopreace
          else
            printf "\t-10" >> $toto
          endif

	printf "\t$total" >> $toto

        printf "\t\t$lane" >> $toto
        set sk=`gawk '/Sequence_kept/{printf("%s",$2);}' Fastc/$lane.count`
        set sr=`gawk '/Sequence_rejected/{printf("%s",$2);}' Fastc/$lane.count`
        printf "\t$sk\t$sr" >> $toto
        echo toto | gawk '{if(sk==0)sk=1;printf("\t%.2f",100*sr/(sr+sk));}' sr=$sr sk=$sk >> $toto
        printf "\t$lane" >> $toto
        foreach target ($targets2)

          set ff=tmp/COUNT/$lane.count
          if (-e $ff) then
            gawk '/HITS/{if($3=="any"){if (index($2,target)>0 && index($2,target)<5){n=1;printf("\t%s",$4);}}}END{if(n==0)printf("\t0");}' target=$target $ff >> $toto
          else
            printf "\t0" >> $toto
          endif

        end
        printf "\t\t$lane" >> $toto
        foreach target ($targets2)

          set ff=tmp/COUNT/$lane.count
          if (-e $ff) then
            gawk '/HITS/{if($3=="any"){if(sk==0)sk=1;if (index($2,target)>0 && index($2,target)<5){n=1;printf("\t%.2f",100*$4/sk);}}}END{if(n==0)printf("\t0");}' sk=$sk target=$target $ff >> $toto
          else
            printf "\t0" >> $toto
          endif

        end

        printf "\t\t$lane" >> $toto
        set tk=`gawk '/Tags_kept/{printf("%s",$2);}' Fastc/$lane.count`
        set tr=`gawk '/Tags_rejected/{printf("%s",$2);}' Fastc/$lane.count`
        printf "\t$tk\t$tr" >> $toto
        echo toto | gawk '{if(tk==0)tk=1;printf("\t%.2f",100*tr/(tr+tk));}' tr=$tr tk=$tk >> $toto
        printf "\t\t$lane" >> $toto
        foreach target ($targets2 )

          set ff=tmp/COUNT/$lane.count
          if (-e $ff) then
            gawk '/HITS/{if($3=="any"){if (index($2,target)>0 && index($2,target)<5){n=1;printf("\t%s",$8);}}}END{if(n==0)printf("\t0");}' target=$target $ff >> $toto
          else
            printf "\t0" >> $toto
          endif

        end
        printf "\t\t$lane" >> $toto
        foreach target ($targets2 )

          set ff=tmp/COUNT/$lane.count
          if (-e $ff) then
            gawk '/HITS/{if($3=="any"){if(tk==0)tk=1;if (index($2,target)>0 && index($2,target)<5){n=1;printf("\t%.2f",100*$8/tk);}}}END{if(n==0)printf("\t0");}' tk=$tk target=$target $ff >> $toto
          else
            printf "\t0" >> $toto
          endif

        end

        printf "\t\t$lane" >> $toto
        set tk=`gawk '/Bases_tags_kept/{printf("%s",$2);}' Fastc/$lane.count`
        set tr=`gawk '/Bases_tags_rejected/{printf("%s",$2);}' Fastc/$lane.count`
        printf "\t$tk\t$tr" >> $toto
        echo toto | gawk '{if(tk==0)tk=1;printf("\t%.2f",100*tr/(tr+tk));}' tr=$tr tk=$tk >> $toto
        printf "\t\t$lane" >> $toto
        foreach target ($targets2 )

          set ff=tmp/COUNT/$lane.count
          if (-e $ff) then
            gawk '/HITS/{if($3=="any"){if (index($2,target)>0 && index($2,target)<5){n=1;printf("\t%s",$9);}}}END{if(n==0)printf("\t0");}' target=$target $ff >> $toto
          else
            printf "\t0" >> $toto
          endif

        end
        printf "\t\t$lane" >> $toto
        foreach target ($targets2 )

          set ff=tmp/COUNT/$lane.count
          if (-e $ff) then
            gawk '/HITS/{if($3=="any"){if(tk==0)tk=1;if (index($2,target)>0 && index($2,target)<5){n=1;printf("\t%.2f",100*$9/tk);}}}END{if(n==0)printf("\t0");}' tk=$tk target=$target $ff >> $toto
          else
            printf "\t0" >> $toto
          endif

        end

        printf "\t\t$lane" >> $toto
        foreach target ($targets2)
          if ($target == any) continue

          set ff=tmp/PHITS_$target/$lane.err
          if (-e $ff) then
            set mem=`gawk -F '\n' '/max memory/{i = index($1,"max memory");s=substr($1,i+11);split(s,aa," ");mem=aa[1];}END{printf("\t%d",mem);}' $ff`
            printf "\t$mem" >> $toto
            printf "Memory\t$run\t$lane\t$target\t$mem\n" >> $totopreace
          else
            printf "\t-1" >> $toto
          endif

        end
          set ff=tmp/COUNT/$lane.err
          if (-e $ff) then
            set mem=`gawk -F '\n' '/max memory/{i = index($1,"max memory");s=substr($1,i+11);split(s,aa," ");mem=aa[1];}END{printf("\t%d",mem);}' $ff`
            printf "\t$mem" >> $toto
            printf "Memory\t$run\t$lane\tSelectBest\t$mem\n" >> $totopreace
          else
            printf "\t-1" >> $toto
          endif

        echo  >> $toto
