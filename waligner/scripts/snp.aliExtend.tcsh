#!bin/tcsh -ef

set run=$1

if (! -e  tmp/SNP/snp.u.complete_list) then
  cat tmp/SNP/*.u.snp.list > tmp/SNP/snp.u.complete_list
endif

if (-e tmp/COUNT2/$run/g.1.hits.gz && ! -e tmp/SNP/$run/aliExtend.txt) then
  gunzip -c tmp/COUNT2/$run/g.1.hits.gz | head -100 | bin/snp -aliExtend -snp_list tmp/SNP/snp.u.complete_list -o tmp/SNP/$run/aliExtend
endif
