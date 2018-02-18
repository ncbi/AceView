#!bin/tcsh -f

set run=$1
set zone=$2

    gunzip -c tmp/SNP_BRS/$run/$zone.BRS.u.gz  | gawk -F '\t' '/^# Target/{zz=1;next;}{if(zz<1)next;}{if($1!="~"){t=$1;}x=$2;if (x!=oldx){if(c>=20)printf("%s\t%d\t%d\t%d\t%d\n",oldt,oldx,c,m,int((100*m/c)));m=0;c=0;}oldx=x;oldt=t;if(($4 == "+" && $3 == "a>g") || ($4 == "-" && $3 == "t>c")){m=$6/100;c=$8/100;next;}if(c==0 && (($4== "+" && $3 == "a") || ($4== "-" && $3 == "t"))){c=$6/100;m=0;}}' | gzip >  tmp/SNP_BRS/$run/$zone.a2g.u.gz

# cut -f 5 | tags | sort -k 1n


