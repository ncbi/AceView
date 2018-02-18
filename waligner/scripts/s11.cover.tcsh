#!bin/tcsh -f

set zone=$1
set MAGIC2=$2

set cover=10
set detect=count

echo "s11 start project:$MAGIC"
set ok=1

    if (! -e tmp/SNPH/$zone/$MAGIC.snp && ! -e tmp/SNPH/$zone/$MAGIC.snp.sorted) then

      echo "s11 hello 1A chrom=$zone"
      set ok=0
      if (-e tmp/SNPH/$zone/$MAGIC.fList) \rm  tmp/SNPH/$zone/$MAGIC.fList
      touch  tmp/SNPH/$zone/$MAGIC.fList
      foreach run (`cat MetaDB/$MAGIC/RunsList  MetaDB/$MAGIC/GroupSnpList `)
   
        if (-e tmp/SNP/$run/$MAGIC2.$zone.count.u.snp.gz) then
          echo  tmp/SNP/$run/$MAGIC2.$zone.count.u.snp.gz >>   tmp/SNPH/$zone/$MAGIC.fList
        endif

      end
      if (-e tmp/SNPH/$zone/$MAGIC.snp) \rm tmp/SNPH/$zone/$MAGIC.snp

    endif

    if (! -e tmp/SNPH/$zone/$MAGIC.snp && ! -e tmp/SNPH/$zone/$MAGIC.snp.sorted) then
      foreach ff (`cat tmp/SNPH/$zone/$MAGIC.fList`)
        gunzip -c $ff | gawk -F '\t' '{if($9 >= m)print}' m=$cover | sort -k 1,1 -k 2,2n -k 3,3 >>  tmp/SNPH/$zone/$MAGIC.snp
      end
    endif

    wc tmp/SNPH/$zone/$MAGIC.fList
    ls -ls  tmp/SNPH/$zone/$MAGIC.snp
    if (-e tmp/SNPH/$zone/$MAGIC.snp && ! -e tmp/SNPH/$zone/$MAGIC.snp.sorted) then

      # nic 2016_09_30  why remove snp.parseChrom.done ? this causes snp.parse_remap.tcsh to reparse huge static acefiles!
      if (0 && -e tmp/SNP_DB/$zone/snp.parseChrom.done) \rm  tmp/SNP_DB/$zone/snp.parseChrom.done
      if (-e tmp/SNP_DB/$zone/snp.parseSnp.done) \rm  tmp/SNP_DB/$zone/snp.parseSnp.done
      if (-e tmp/SNP_DB/$zone/snp.remap.done) \rm  tmp/SNP_DB/$zone/snp.remap.done
      if (-e tmp/SNP_DB/$zone/snp.translate.done) \rm  tmp/SNP_DB/$zone/snp.translate.done

      echo "sorting  tmp/SNPH/$zone/$MAGIC.snp"
      cat  tmp/SNPH/$zone/$MAGIC.snp | sort -k 1,1 -k 2,2n -k 3,3  | gawk -F '\t' '/^#/{nt[$0]++;if(nt[$0]<2)print;next;}{if(NF > nfMax)nfMax=NF;z= $1 "\t" $2 "\t" $3 ; if(z != oldz){ if(index(bestP bestS,"NN")==0){for(n=1;n<=0+nn;n++){printf("%s\t%s\t%s",oldz,bestP,bestS);for(i=6;i<=nfMax;i++)printf("\t%s",dd[n,i]);printf("\n");}}nn=0;bestP=$4;bestS=$5;oldz=z;pHasN=split(bestP,aa,"N")-1;sHasN=split(bestS,aa,"N")-1;}nn++;for(i=6;i<=nfMax;i++)dd[nn,i]=$i;if(pHasN>0){pHasN2=split($4,aa,"N")-1;if(pHasN2<pHasN){pHasN=pHasN2;bestP=$4;}}if(sHasN>0){sHasN2=split($5,aa,"N")-1;if(sHasN2<sHasN){sHasN=sHasN2;bestS=$5;}}}' > tmp/SNPH/$zone/$MAGIC.snp.sorted
      cat tmp/SNPH/$zone/$MAGIC.snp.sorted | gawk -F '\t' '{n++;}/^#/{next;}{if (NF != 23 || length($0)>1000)print n,NF; else n23++;}END{print 23,n23;}' >  tmp/SNPH/$zone/$MAGIC.snp.sorted.NF
      \rm tmp/SNPH/$zone/$MAGIC.snp
    endif 

     if (-e tmp/SNPH/$zone/$MAGIC.snp.sorted && ! -e tmp/SNPH/$zone/$MAGIC.snp.sorted.homozygous) then
       cat  tmp/SNPH/$zone/$MAGIC.snp.sorted |  gawk -F '\t' '{r=$6;c=$9;m=$10;if (c >= 10 && (100*m > 95*c || 100*m < 5 * c))print }' >  tmp/SNPH/$zone/$MAGIC.snp.sorted.homozygous 
     endif

     # 2016_08   disable differential_snps , we now prefer to use the characteristic snps prepard in s21 after filtering
     if (0 && -e tmp/SNPH/$zone/$MAGIC.snp.sorted && -e tmp/GENERUNS/$MAGIC.snp.SNP.info.ace && ! -e tmp/SNPH/$zone/$MAGIC.differential_snps) then
      echo " bin/geneindex -deepSNP tmp/SNPH/$zone/$MAGIC.snp.sorted -u -runList MetaDB/$MAGIC/RunListSorted -runAce tmp/GENERUNS/$MAGIC.snp.SNP.info.ace -pA -snpEval >  tmp/SNPH/$zone/$MAGIC.differential_snps"
             bin/geneindex -deepSNP tmp/SNPH/$zone/$MAGIC.snp.sorted -u -runList MetaDB/$MAGIC/RunListSorted -runAce tmp/GENERUNS/$MAGIC.snp.SNP.info.ace -pA -snpEval >  tmp/SNPH/$zone/$MAGIC.differential_snps
    endif

    if ( -e  tmp/SNPH/$zone/$MAGIC.differential_snps && ! -e tmp/SNPH/$zone/$MAGIC.snp.sorted.differential) then
	cat tmp/SNPH/$zone/$MAGIC.differential_snps ZZZZZ  tmp/SNPH/$zone/$MAGIC.snp.sorted | gawk -F '\t'  '/ZZZZZ/{zz++;next;}/^#/{next;}{if(zz<1){split($1,aa,":");ok[aa[1] ":" aa[2] ":" aa[3]]=1;next;}}{t=$3;gsub(">","2",t);z=$1 ":" $2 ":" t;if(ok[z]==1)print;}' | sort -u  >  tmp/SNPH/$zone/$MAGIC.snp.sorted.differential
    endif

exit 0

echo aaa
cat run2title.txt | gawk '{gsub("\"","",$0);r=$1;t=$2;printf("gunzip -c tmp/Unaligned/%s/f2.*.fastc.gz | gawk \'/^ > /{s=$1 ; next}{if ( length ( \$1 ) == 202 ) {print s ; print ; }}\' | dna2dna -I fastc -O fasta -o ~/ftp-private/NB_unali/%s.unaligned\n",r,t);}' > _m
echo ccc
cat run2title.txt | gawk '{gsub("\"","",$0);r=$1;t=$2;printf("gunzip -c tmp/Unaligned/%s/f2.*.fastc.gz | gawk \'/^ > /{s=$1 ; next}{if ( length ( \$1 ) == 202 ) {print s ; print ; }}\'\n",r,t);}' > _m
echo bbb
cat run2title.txt | gawk '{gsub("\"","",$0);r=$1;t=$2;printf("gunzip -c tmp/Unaligned/%s/f2.*.fastc.gz\t%s\t%s\n",r,r,t);}' > _m
cat _m | gawk -F '\t' '{printf("%s|gawk -f titi.awk |  dna2dna -I fastc -O fasta -o ~/ftp-private/NB_unali/%s.unaligned\n",$1,$3);}' > _m2

## count the number of SNPs seen per gene


 
exit 0

################
################

if (0) then

foreach bb (`ls -d B[0-9]*`)
  setenv MAGIC $bb
  echo $bb
  pushd $bb
    MAGIC s11
  popd
end

foreach ff (`ls tmp/SNPH/*/*.snp.sorted `)
  cat $ff | gawk -F '\t' '{n++;}/^#/{next;}{if (NF != 23 || length($0)>1000)print n,NF; else n23++;}END{print 23,n23;}' >  $ff.NF
end

foreach ff (`ls B[0-9]*/tmp/SNPH_july2014/*/*.snp.sorted `)
  cat $ff | gawk -F '\t' '{n++;}/^#/{next;}{if (NF != 23 || length($0)>1000)print n,NF; else n23++;}END{print 23,n23;}' >  $ff.NF
end

foreach bb (`ls -d B[0-9]`)
  date ; echo $bb
  tar cf ALL_SNPH.2014_04_07/$bb.snp.cover.10.sorted.tar $bb/tmp/SNPH/*/*.snp.sorted &
end

sleep 10000
foreach bb (`ls -d B[0-9]*`)
  gzip ALL_SNPH.2014_04_07/$bb.snp.cover.10.sorted.tar &
end

sleep 7200
rsync -avh --whole-file  ALL_SNPH.2014_04_07/*.snp.cover.10.sorted.tar.gz ~/ftp-private/SNPH

