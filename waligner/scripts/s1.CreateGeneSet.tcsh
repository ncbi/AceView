#!bin/tcsh -f

cat <<EOF > gene.set.list
TP53
BRCA1P1
CSMD3
NF1
CDK12
FAT3
GABRA6
BRCA2
RB1
ATM
ATR
C11orf30
FANCD2
RAD51C
CHEK2
FOXM1
CDC25B
AFMIDandBIRC5
PLK1
AURKB
CCNB1
CDKN2A
CCNE1
CCND1
CCND2
JAG1
NUDT14andJAG2
NOTCH3
MAML1
MAML2
MAML3
PTEN
PIK3CA
AKT1
AKT2
NF1
KRAS
BRAF
EOF

cat <<EOF > mrna.set.list
RPS11andSNORD35B.b	RPS11.b
FTL.d	FTL.d
FTL.a	FTL.a
RPL3_.a	RPL3.a
RPL15.e	RPL15.e
RPS3andSNORD15B.c	RPS3.c clipped at 839 as the rest of the UTR is not well supported, and not in refseq
RPL37andSNORD72.g	RPL37.g
RPL4_.a	RPL4.a
RPS14.b	RPS14.b 
RPS14.f	RPS14.f
ACTG1.b	ACTG1.b
EIF1.d	EIF1.d	
H3F3B.f	H3F3B.f
FTH1.b	FTH1.b
RPL6.a	RPL6.a
RPL5andSNORD21andSNORA66.a	RPL5.a
RPL13andSNORD68.b	RPL13.b
RPL13andSNORD68.a	RPL13.a
RPS7.e	RPS7.e
CALM2andC2orf61.i	variant I, subtract first 47 bases as supported only by NM
EEF2.a	variant a
HSP90AA1.c	variant c
ALDOA.f	variant f
ATP5BandSNORD59AandSNORD59B.a	variant a
DDX5.b	variant b (NM)
CALR.a	variant a
ARL6IP1andRPS15A.h	variant h
ATP5A1.c	variant c
MALAT1.a	FJ209305
ATP1A1.b	variant b
FLNA.a	variant a
FLNA.b	variant b
COX6C.d	variant d (NM)
ATP5J.j	variant j reduced 680 to end 1309
CTSD.b	variant b
CTSB.h	variant h, but long UTR is not supported: clip at base 1986 
ATP5I.a	variant a
DAZAP2.c	variant c, clipped at 2005
CCDC72.b	variant b (NM)
AHNAK.a	variant a
NEAT1.a	variant a
EOF


cat gene.set.list ZZZZZ TARGET/GENES/av.gene.ace | gawk '/^ZZZZZ/{zz=1;g=0;next;}{if(zz<1){gg[$1]=1;next;}}/^Gene/{gsub(/\"/,"",$0);if(gg[$2]==1)g=$2;next;}/^IntMap/{if(g!=0){gsub(/\"/,"",$0);ch=$2;a1=$3;a2=$4;da=a2-a1;if(da<0)da=-da;da++;printf("%s\t%d\t%d\t%s\t%d\t%d\n",g,1,da,ch,a1,a2);g=0;}}' > tmp/SNP/gene.set.remap
cat gene.set.list ZZZZZ TARGET/GENES/av.gene.ace | gawk '/^ZZZZZ/{zz=1;g=0;next;}{if(zz<1){gg[$1]=1;next;}}/^Gene/{gsub(/\"/,"",$0);if(gg[$2]==1)g=$2;next;}/^IntMap/{if(g!=0){gsub(/\"/,"",$0);ch=$2;a1=$3;a2=$4;da=a2-a1;if(da<0){da=-da;a0=a1;a1=a2;a2=a0;}da++;printf("%s\t%d\t%d\t999\n",ch,a1,a2);g=0;}}' > tmp/SNP_ZONE/zoneg.999.txt

bin/dna2dna -i TARGET/Targets/$species.genome.fasta.gz -select tmp/SNP_ZONE/zoneg.999.txt -I fasta -O fasta -gzo -o tmp/SNP_ZONE/zoneg.999



