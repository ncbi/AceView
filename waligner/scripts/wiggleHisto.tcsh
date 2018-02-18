#!bin/tcsh -ef

set chrom=$1
set tissue=$2
set uu=$3

echo "wiggleHisto.tcsh $chrom $tissue"
if ($tissue == any || $tissue == BodyMapSEQC) then

  if (-e tmp/WIGGLE/$chrom/$tissue.$uu.E.f.bf.gz ) then 
    gunzip -c  tmp/WIGGLE/$chrom/$tissue.$uu.E.[fr].bf.gz |  gawk -F '\t' '/[a-z_]/{next}{x=$1 ; nn[x]++;if(x>xmax)xmax=x;}END{for(x in nn)printf("%s\t%d\n",x,nn[x]);}'  | sort -k 1n > tmp/WIGGLE/$chrom/$tissue.$uu.histo.fr.txt
  endif
  if (-e  tmp/WIGGLE/$chrom/$tissue.$uu.E.ns.bf.gz) then
    gunzip -c  tmp/WIGGLE/$chrom/$tissue.$uu.E.ns.bf.gz | gawk -F '\t' '/[a-z_]/{next}{x=$1 ; nn[x]++;if(x>xmax)xmax=x;}END{for(x in nn)printf("%s\t%d\n",x,nn[x]);}' | sort -k 1n > tmp/WIGGLE/$chrom/$tissue.$uu.histo.ns.txt
  endif

else

  gunzip -c  tmp/WIGGLE/$chrom/$tissue.any_manip.$uu.[fr].bf.gz |  gawk -F '\t' '/[a-z_]/{next}{x=$1 ; nn[x]++;if(x>xmax)xmax=x;}END{for(x in nn)printf("%s\t%d\n",x,nn[x]);}'  | sort -k 1n > tmp/WIGGLE/$chrom/any.$tissue.$uu.histo.fr.txt
  gunzip -c  tmp/WIGGLE/$chrom/$tissue.any_manip.$uu.ns.bf.gz | gawk -F '\t' '/[a-z_]/{next}{x=$1 ; nn[x]++;if(x>xmax)xmax=x;}END{for(x in nn)printf("%s\t%d\n",x,nn[x]);}' | sort -k 1n > tmp/WIGGLE/$chrom/any.$tissue.$uu.histo.ns.txt

endif



