
# rs1102 | human | 9606 | snp | genotype=NO | submitterlink=YES | updated 2004-10-04 13:37
# ss1127 | WIAF | WIAF-1571 | orient=+ | ss_pick=YES
# SNP | alleles='C/T' | het=? | se(het)=?
# VAL | validated=NO | min_prob=? | max_prob=? | notwithdrawn
# CTG | assembly=Celera | chr=18 | chr-pos=30414911 | NW_927095.1 | ctg-start=15083316 | ctg-end=15083316 | loctype=2 | orient=-
# CTG | assembly=GRCh37 | chr=18 | chr-pos=33607030 | NT_010966.14 | ctg-start=15096132 | ctg-end=15096132 | loctype=2 | orient=-
# LOC | RPRD1A | locus_id=55197 | fxn-class=missense | allele=C | frame=1 | residue=P | aa_position=208 | mrna_acc=NM_018170.3 | prot_acc=NP_060640.2
# LOC | RPRD1A | locus_id=55197 | fxn-class=reference | allele=T | frame=1 | residue=S | aa_position=208 | mrna_acc=NM_018170.3 | prot_acc=NP_060640.2
# CTG | assembly=HuRef | chr=18 | chr-pos=30466019 | NW_001838467.2 | ctg-start=10934783 | ctg-end=10934783 | loctype=2 | orient=+
# SEQ | 142385371 | source-db=remap | seq-pos=628 | orient=+

/^rs/ { if(ok>=3)printf("%s\n", old) ; old="" ; ok = 1 ; hetOk=0; rs=$1 ; printf("### rs=%s", rs) ; next ; } 
/^SNP/ { i=index($2,"alleles=");s=substr($2,i+9);gsub(/'/,"",s);allele=s;  gsub(/ /,"",allele);printf("### allele=%s\n", s) ;ok++ ;}
/^SNP/ { i=index($3,"het=");s=substr($3,i+5);i=index(s, " ");het=substr(s,1,i-1);if(het != "?")hetOk=het ;printf("### het=%s\n", het) ;}
/^CTG/ {if($2 == " assembly=GRCh37 "){ i=index($3,"chr=");chrom=substr($3,i+4);gsub(/ /,"",chrom);i=index($4,"chr-pos=");pos=substr($4,i+8);gsub(/ /,"",pos);printf("### pos=%s\n", pos) ;vnam = "chrom" chrom ":" pos allele ; ok++; old = "Variant \"" vnam "\"\n" ; if(hetOk > 0) old = old "Frequency " hetOk "\n" ;}}
END { if (ok >= 0) printf("%s\n", old) ; }
