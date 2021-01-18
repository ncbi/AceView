#!bin/tcsh -f

set run=$1
set ff0=$2
set PAIRED=0
if ($3 == 1) then
  set PAIRED=1
endif

echo aa $ff0
if (! -e $ff0) then
  echo "file not found $ff0"
  goto done
endif

set ff=Fastc/$run/f2
if (-e $ff.1.fastc.gz) goto done

set ff=Fastc/$run/tmp/f2
if (-e $ff.1.fasta.gz) goto done
if (-e $ff.fasta.gz) goto done

if (! -d Fastc) then
  echo "Missing directory Fastc"
  goto done
endif
if (! -d Fastc/$run) then
  echo "Missing directory Fastc/$run"
  goto done
endif
if (! -d Fastc/$run/tmp) then
  mkdir Fastc/$run/tmp
endif
if (! -d Fastc/$run/tmp) then
  echo "Missing directory Fastc/$run/tmp"
  goto done
endif

setenv LANG C

echo -n "... samtools: bam -> raw "
date
samtools view $ff0 -f 16 -F 256  | bin/dna2dna -I raw10 -complement -O raw -o $ff.r -gzo -keepNameBam 
samtools view $ff0 -F 272  | bin/dna2dna -I raw10 -O raw -o $ff.f -gzo  -keepNameBam 


if ($PAIRED == 1) then

  echo -n "... alphabetical sort by sequence identifer: raw -> raw.sorted "
  date
  gunzip -c $ff.f.raw.gz | grep '/1' | gawk '{printf("%s\t%s\n",$2,$1);}' | sort > $ff.f.1.raw.sorted 
  gunzip -c $ff.f.raw.gz | grep '/2' | gawk '{printf("%s\t%s\n",$2,$1);}' | sort > $ff.f.2.raw.sorted 
  gunzip -c $ff.r.raw.gz | grep '/1' | gawk '{printf("%s\t%s\n",$2,$1);}' | sort > $ff.r.1.raw.sorted 
  gunzip -c $ff.r.raw.gz | grep '/2' | gawk '{printf("%s\t%s\n",$2,$1);}' | sort > $ff.r.2.raw.sorted 

  echo -n "... merge the .f and the .r reads [f  r].raw.sorted -> fr.raw.sorwted "
  date
  cat  $ff.f.1.raw.sorted $ff.r.1.raw.sorted | sort >  $ff.fr.1.raw.sorted
  cat  $ff.f.2.raw.sorted $ff.r.2.raw.sorted | sort >  $ff.fr.2.raw.sorted

  echo -n "... verif against the eventual orphans: raw.sorted -> verified.raw "
  date
  # it may happen that some reads are absent on one side
  echo ZZZZZ > ZZZZZ
  cat  $ff.fr.1.raw.sorted  $ff.fr.2.raw.sorted ZZZZZ $ff.fr.1.raw.sorted | gawk '/^ZZZZZ/{zz++;next;}{z=substr($1,1,length($1)-1);if(zz<1){nn[z]++;next;}if(nn[z]==2)print;}' >  $ff.fr.verified.1.raw
  cat  $ff.fr.1.raw.sorted  $ff.fr.2.raw.sorted ZZZZZ $ff.fr.2.raw.sorted | gawk '/^ZZZZZ/{zz++;next;}{z=substr($1,1,length($1)-1);if(zz<1){nn[z]++;next;}if(nn[z]==2)print;}' >  $ff.fr.verified.2.raw

  echo -n "... raw -> fasta.gz " 
  date
  cat  $ff.fr.verified.1.raw | gawk '{printf(">%s\n%s\n",$1,$2);}' | gzip > $ff.1.fasta.gz 
  cat  $ff.fr.verified.2.raw | gawk '{printf(">%s\n%s\n",$1,$2);}' | gzip > $ff.2.fasta.gz 
else
  gunzip -c $ff.[fr].raw.gz | sort >  $ff.sorted.raw
  # cat $ff.sorted.raw | gawk '{printf(">%s\n%s\n",$2,$1);}' | gzip  > $ff.fasta.gz
  cat $ff.sorted.raw ZZZZZ | gawk '{s=$1;if(s == old){n++;next;}if(n+0>0){ss += 10 + length(old); nf = 1+int(ss/(1000000*Mb));kk++;printf (">n.%d#%d\n%s\n",kk,n,old) > fil"."nf".fastc";}n=1;old=$1;}' fil=$ff Mb=200
  foreach g (`ls $ff.*.fastc | sed -e 's/\.fastc//'`)
    dna2dna -i $g.fastc -I fastc -count > $g.count
    gzip $g.fastc 
  end
  mv $ff.*.fastc.gz $ff.*.count Fastc/$run/
  \rm -rf  Fastc/$run/tmp
endif


# \rm $ff.*raw*

done:
ls -ls  Fastc/$run
echo -n "... done "
date


