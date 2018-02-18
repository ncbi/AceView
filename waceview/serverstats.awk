/^2014/ { mm=0+substr($1,6,2); dd =0+ substr($1,9,2) ; if (mm == mmm) { tothits++; nnhits[dd]++ ;} }
/^2014.*REMOTE_ADDR/ { mm=0+substr($1,6,2); dd =0+ substr($1,9,2) ; if (mm == mmm) {uuu=dd $9; client[$9]++;dayclient[uuu]++; if(client[$9]==1)totClt++;if(dayclient[uuu]==1)clt[dd]++;} }
/^2014.*c=product.*a=fiche/ { mm=0+substr($1,6,2); dd =0+ substr($1,9,2) ; if (mm == mmm) { totfmrna++; nnfmrna[dd]++ ;}next; }
/^2014.*open/ { mm=0+substr($1,6,2); dd =0+ substr($1,9,2) ; if (mm == mmm) { totopen++; nnopen[dd]++ ;} next;}

/^2014.*c=mrna.*a=fiche/ { mm=0+substr($1,6,2); dd =0+ substr($1,9,2) ; if (mm == mmm) { totfmrna++; nnfmrna[dd]++ ;}next; }
/^2014.*a=fiche/ { mm=0+substr($1,6,2); dd =0+ substr($1,9,2) ; if (mm == mmm) { totgene++; nngene[dd]++ ;} }
/^2014.*a=fgene/ { mm=0+substr($1,6,2); dd =0+ substr($1,9,2) ; if (mm == mmm) { totgene++; nngene[dd]++ ;} }
/^2014.*a=fmol/ { mm=0+substr($1,6,2); dd =0+ substr($1,9,2) ; if (mm == mmm) { totfmol++; nnfmol[dd]++ ;} }
/^2014.*a=fexp/ { mm=0+substr($1,6,2); dd =0+ substr($1,9,2) ; if (mm == mmm) { totfexp++; nnfexp[dd]++ ;} }
/^2014.*a=ffunc/ { mm=0+substr($1,6,2); dd =0+ substr($1,9,2) ; if (mm == mmm) { totffunc++; nnffunc[dd]++ ;} }
/^2014.*a=vmrna/ { mm=0+substr($1,6,2); dd =0+ substr($1,9,2) ; if (mm == mmm) { totvmrna++; nnvmrna[dd]++ ;} }
/^2014.*&s&/ { mm=0+substr($1,6,2); dd =0+ substr($1,9,2) ; if (mm == mmm) { totflash++; nnflash[dd]++ ;} }
/^2014.*&G&/ { mm=0+substr($1,6,2); dd =0+ substr($1,9,2) ; if (mm == mmm) { totgif++; nngif[dd]++ ;} }
/^2014.*a=vgene/ { mm=0+substr($1,6,2); dd =0+ substr($1,9,2) ; if (mm == mmm) { totvgene++;nnvgene[dd]++ ;} }
/^2014.*a=vmrna/ { mm=0+substr($1,6,2); dd =0+ substr($1,9,2) ; if (mm == mmm) { totvmrna++;nnvmrna[dd]++ ;} }
/^2014.*a=hmrna/ { mm=0+substr($1,6,2); dd =0+ substr($1,9,2) ; if (mm == mmm) { tothmrna++;nnhmrna[dd]++ ;} }
/^2014.*a=locatorbig/ { mm=0+substr($1,6,2); dd =0+ substr($1,9,2) ; if (mm == mmm) { tothloc++;nnhloc[dd]++ ;} }
/^2014.*a=locator&/ { mm=0+substr($1,6,2); dd =0+ substr($1,9,2) ; if (mm == mmm) { tothlocbig++;nnhlocbig[dd]++ ;} }
/^2014.*a=clone/ { mm=0+substr($1,6,2); dd =0+ substr($1,9,2) ; if (mm == mmm) { totclo++; nnclo[dd]++ ;} }
/^2014.*a=fasta/ { mm=0+substr($1,6,2); dd =0+ substr($1,9,2) ; if (mm == mmm) { totfasta++; nnfasta[dd]++ ;} }
/^2014.*a=DNA/ { mm=0+substr($1,6,2); dd =0+ substr($1,9,2) ; if (mm == mmm) { totDNA++; nnDNA[dd]++ ;} }

/^2014.*c=PFam/ { mm=0+substr($1,6,2); dd =0+ substr($1,9,2) ; if (mm == mmm) { totpfam++; nnpfam[dd]++ ;} }
/2014.*blastn/ {if ($2==month) { dd = $3 ; totblast++;nnblast[dd]++;}}
/^2014.*org=/ { mm=0+substr($1,6,2); dd =0+ substr($1,9,2) ; if (mm == mmm) { totfromloc++; nnfromloc[dd]++ ;} }
/^2014.*ctx=ctx/ { mm=0+substr($1,6,2); dd =0+ substr($1,9,2) ; if (mm == mmm) { totctx++; nnctx[dd]++ ;} }
/^2014.*webpfamg/ { mm=0+substr($1,6,2); dd =0+ substr($1,9,2) ; if (mm == mmm) { totpfX++; nnpfX[dd]++ ;} }

END {
for (nam in client) print "XXXXX", client[nam], nam ;
  printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
         "day","client","hits","gene","hrna","loc","big","vgen","mrna","vrna","mol","exp","func", "DNA","fasta","pfam","open","ctx", "blst","pfX","flash","gif",":ll","clo") ;
			

    printf("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
      "cumul", totClt, tothits,totgene, tothmrna,tothloc,tothlocbig,totvgene,totfmrna,totvmrna,totfmol,totfexp,totffunc,totDNA, totfasta,totpfam,totopen,totctx,totblast,totpfX,totflash,totgif,totfromloc,totclo) ; 

    for (ii=0; ii <= 31; ii++) 
      { if (nnhits[ii]>0) printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", 
    ii,clt[ii],nnhits[ii],nngene[ii],nnhmrna[ii],nnhloc[ii],nnhlocbig[ii],nnvgene[ii],nnfmrna[ii],nnvmrna[ii],nnfmol[ii],nnfexp[ii],nnffunc[ii],nnDNA[ii], nnfasta[ii],nnpfam[ii],nnopen[ii],nnctx[ii],nnblast[ii],nnpfX[ii],nnflash[ii],nngif[ii],nnfromloc[ii],nnclo[ii]) ;
      } 
}
