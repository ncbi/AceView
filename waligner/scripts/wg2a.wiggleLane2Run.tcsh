#!bin/tcsh -f

set run=$1
set uu=$2

set out_step="-out_step 10"
if ($?wiggle_step) then
   set out_step="-out_step $wiggle_step"
endif

foreach chrom (mito SpikeIn $chromSetAll)
  echo "... wg2a $run $chrom"
  if (! -d tmp/WIGGLERUN/$run/$chrom) mkdir tmp/WIGGLERUN/$run/$chrom

      foreach fr (f r ELF ELR ERF ERR)

          # 2013_03_25  the genes cases is not programmed in WIGGLELANE, the -ventilate option does not support it
          #   tmp/WIGGLERUN/$run/$chrom/R.genes.$uu.$fr.BF.gz

          if (! -e  tmp/WIGGLERUN/$run/$chrom/R.chrom.$uu.$fr.BF.gz) then
            set n=`ls -ls tmp/WIGGLELANE/$run/*/K.*mapped.$chrom.$uu.$fr.BV.gz | gawk '{if($6==20)next;n++;}END{print 0+n;}'`
            if ($n > 0) then
              echo "Constructing tmp/WIGGLERUN/$run/$chrom/R.chrom.$uu.$fr.BF.gz"
              echo     "tmp/WIGGLELANE/$run/*/K.*mapped.$chrom.$uu.$fr.BV.gz | bin/wiggle -I BV -O BF $out_step -gzo -o tmp/WIGGLERUN/$run/$chrom/R.chrom.$uu.$fr"
              gunzip -c tmp/WIGGLELANE/$run/*/K.*mapped.$chrom.$uu.$fr.BV.gz | bin/wiggle -I BV -O BF $out_step -gzo -o tmp/WIGGLERUN/$run/$chrom/R.chrom.$uu.$fr
	    endif
          endif
      end
      if ((-e  tmp/WIGGLERUN/$run/$chrom/R.chrom.$uu.f.BF.gz || -e  tmp/WIGGLERUN/$run/$chrom/R.chrom.$uu.r.BF.gz) && ! -e  tmp/WIGGLERUN/$run/$chrom/R.chrom.frns.$uu.BF.gz) then
        gunzip -c tmp/WIGGLERUN/$run/$chrom/R.chrom.$uu.[fr].BF.gz | bin/wiggle -I BF -O BF $out_step -gzo -o  tmp/WIGGLERUN/$run/$chrom/R.chrom.frns.$uu -cumul
     endif
      if (( -e  tmp/WIGGLERUN/$run/$chrom/R.chrom.$uu.ELF.BF.gz || -e  tmp/WIGGLERUN/$run/$chrom/R.chrom.$uu.ELR.BF.gz) && ! -e  tmp/WIGGLERUN/$run/$chrom/R.chrom.EL.$uu.BF.gz) then
        gunzip -c tmp/WIGGLERUN/$run/$chrom/R.chrom.$uu.EL[FR].BF.gz | bin/wiggle -I BF -O BF $out_step -gzo -o  tmp/WIGGLERUN/$run/$chrom/R.chrom.EL.$uu -cumul
     endif
      if (( -e  tmp/WIGGLERUN/$run/$chrom/R.chrom.$uu.ERF.BF.gz || -e  tmp/WIGGLERUN/$run/$chrom/R.chrom.$uu.ERR.BF.gz) && ! -e  tmp/WIGGLERUN/$run/$chrom/R.chrom.ER.$uu.BF.gz) then
        gunzip -c tmp/WIGGLERUN/$run/$chrom/R.chrom.$uu.ER[FR].BF.gz | bin/wiggle -I BF -O BF $out_step -gzo -o  tmp/WIGGLERUN/$run/$chrom/R.chrom.ER.$uu -cumul
     endif

  # 2015_06_12  Nettoyage des fuites de strand, est ce necessaire ?


end


foreach chrom ($chromSetAll)
    foreach cover (10)
        gunzip -c  tmp/WIGGLERUN/$run/$chrom/R.chrom.frns.$uu.BF.gz | bin/wiggle -I BF -O BV  -gauss 20 -minCover $cover -peaks -o tmp/WIGGLERUN/$run/$chrom/coverome.$cover.$uu
    end

    if ($uu == u) then
      cat  tmp/WIGGLEGROUP/$run/$chrom/coverome.$minCoveron.u.peaks | gawk '/^#/{next;}{nn++;chrom=$1;a1=$2;a2=$3;if(a1>a2){a0=a1;a1=a2;a2=a0;}printf("%s\t%09d\t%09d\tCOVERON\tCoveron.%s.%s_%d\t%d\t%s\t%s\t%s\n",chrom,a1,a2,run,chrom,nn,$4,$5,$6,$7); }' run=$run >> toto.sp1
      cat toto.sp1 | sort -k 1,1 -k2,2n -k3,3n |  gawk -F '\t' 'BEGIN{printf("#Target\ta1\ta2\tCoveron\tGene\tLength bp\tmax cover\taverage cover\taligned bp\n");}/^#/{next;}{chrom=$1;a1=$2;a2=$3;type=$4;if(type=="GENE"){g++;gnam[g]=$5;gch[g]=chrom;g1[g]=a1;g2[g]=a2;gz=0;if(a2>amax){amax=a2;top=g;}next;}if(type=="COVERON"){gz=0;for(i=g;i>=top;i--){if(gch[i]==chrom && a1<g2[i] && a2>g1[i]){gz=gnam[i];break;}}if(gz==0){ngz++;gz="New_" chrom "_"  ngz;}printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$5,gz,$6,$7,$8,$9);}}'  > tmp/WIGGLEGROUP/$run/$chrom/coverome2gene.$minCoveron.u.peaks

    \rm toto.sp[01]
    endif
end


touch tmp/WIGGLERUN/$run/wg2a.$uu.done
\rm -rf tmp/WIGGLELANE/$run/*/*.$uu.*.BV.gz




       
