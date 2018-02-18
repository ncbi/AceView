#!bin/tcsh -f

# 2015_10_02
# Start from SraRunInfo.csv
# Using the SRA interface, select relevant runs
# use menu save-as file -> format RunInfo
#    the data should be set in excell to the format 12-Mar-15 (for 2015)
#    open the file in Excell on the PC, save as txt tab delimited, 
#    this transforms the format of the dates (which in csv contain comas)
#                  file -> format XML

# please edit the configuration
set phase = $1

if ($species == coli) then
  set date=2015_10_03
  set ff=/home/mieg/ACEVIEWHELP/E.coli/1734E.coliK12_Oct7_2015_SRApublic_SraRunInfo.txt
endif

if ($species == Dmelanogaster) then
  set date=2016_04_19
  set ff=/home/mieg/ACEVIEWHELP/Droso/DATA/Droso_allRNA.2016_04_18.SraRunInfo.txt
  set date=2016_10_12 
  set ff=/home/mieg/ACEVIEWHELP/Droso/DATA/SraRunInfo_projectsMissingFromSRXdb_oct12_2016_andInFBfile.txt
  set date=2016_10_14
  set ff=/home/mieg/ACEVIEWHELP/Droso/DATA/2016Oct15_Drosophila_melanogaster_Organism_NOT_biomol_dna_Properties_13863_SraRunInfo.txt
  set ff=/home/mieg/ACEVIEWHELP/Droso_DATA/DATA/DrosophilaMelano_2017Mar20_4640RNA_since2016Apr10_PreviousWasApr19_SraRunInfo.txt
endif
if ($species == hs) then
  set date=2016_10_14
  set ff='/home/mieg/ACEVIEWHELP/Human_DATA/*SraRunInfo.txt'
endif
if ($species == rn) then
  set date=2016_11_04
  set ff='/home/mieg/ACEVIEWHELP/Rat_DATA/2016Nov05_6365_Rattus_norvegicus_Organism_NOT_biomol_dna_PropertiesSraRunInfo.txt'
endif
if ($species == mm) then
  set date=2016_11_04
  set ff='/home/mieg/ACEVIEWHELP/Mouse_DATA/2016Nov5_131768_mus_musculus_Organism_NOT_biomol_dnaPropertiesSraRunInfo.txt'
endif
if ($species == dog) then
  set date=2016_11_04
  set ff='/home/mieg/ACEVIEWHELP/Dog_DATA/Dog_3strandedRunsForTestTarget.txt'
endif
if ($species == Campylobacter || $species == jejuni) then
  set date=2017_03_07
  set ff='/home/mieg/ACEVIEWHELP/Viruses_microbes_DATA/Campylobacter_Skesa_OXA_SraRunInfo.txt'
  set ff='/home/mieg/ACEVIEWHELP/Viruses_microbes_DATA/CampylobacterProject_SRP063302_1952entries_SraRunInfo.txt'
  set ff='/home/mieg/ACEVIEWHELP/Viruses_microbes_DATA/Campylobacter_jejuni_genome_Truseq_SraRunInfo.txt'
endif

echo "species=$species date=$date"
set today=`echo $date | sed -e 's/_/-/g'`
set dd=DATA/$date
setenv ici `pwd` 
set ff1=$dd/runInfo.tsv


# Create a work directory
if (! -d DATA) mkdir DATA
if (! -d DATA/$date) mkdir DATA/$date
if (! -d $dd) mkdir $dd

# copy the RunInfo file in Unix format and transform to tab delimited
cat $ff | sed -e 's/\r//' | grep -v "Synthetic-Long-Read" > $dd/runInfo.tsv
set ff1=$dd/runInfo.tsv


if ($phase == SRR) goto phaseSRR
if ($phase == SRP) goto phaseSRP
if ($phase == GEO) goto phaseGEO
if ($phase == Sample) goto phaseSample
if ($phase == SRX) goto phaseSRX
if ($phase == Files) goto phaseFiles
if ($phase == Titles) goto phaseTitles
if ($phase == srr2run) goto phase_srr2run
if ($phase == srr2srx) goto phase_srr2srx
if ($phase == Papers) goto phasePapers

if ($phase == sraDownload) goto phase_sraDownload
if ($phase == sraDownloadTest) goto phase_sraDownload


echo "usage: SRX_import.tcsh SRR SRP GEO Sample SRX Files Papers  Titles  srr2run |  sraDownload sraDownLoadTest | srr2srx"
goto phaseLoop

############
phaseSRR:
# parse the SRR number, use it as intermediary Run name
cat $ff1 | gawk -F '\t' -f scripts/SRX_import.1.awk today=$today | grep -v Submission_dateOK > $dd/run_info.ace

if (! -d SRX_DB) then
  ln -s ~/ace/waligner/metaData
  ln -s ~/ace/waligner/scripts
  ln -s ~/ace/bin.ICC_centos7 bin
  mkdir SRX_DB
  pushd  SRX_DB
    ln -s ../metaData/wspec.SRX wspec
    mkdir database
    echo y | ../bin/tacembly .
  popd
endif

bin/tacembly SRX_DB <<EOF
  read-models
  parse $dd/run_info.ace
  edit project $MAGIC
  save
  quit
EOF
goto phaseLoop

####### Parse the SRP objects i.e. title and abstracts for the runs
phaseSRP:
if (! -e $dd/SRP) mkdir $dd/SRP
bin/tacembly SRX_DB <<EOF
  query find project $MAGIC ; >SRR ; > SRP 
  spush
  query find srp srr && ! title
  sor
  spop
  // query ! title
  list -a -f $dd/srp.list
  quit
EOF

if (-e  $dd/SRP/_wget) \rm  $dd/SRP/_wget
foreach srp (`cat $dd/srp.list | gawk '/^SRP/{gsub(/\"/,"",$2);if(n[$2]<1){n[$2]=1;print $2;}}'`)
   if ( -e  $dd/SRP/$srp.html) continue 
   echo "wget -O "$dd"/SRP/"$srp".html  "'"http://trace.ncbi.nlm.nih.gov/Traces/sra/?study='$srp'"' >> $dd/SRP/_wget
end

if (-e  $dd/SRP/_wget) then
  wc  $dd/SRP/_wget
  source  $dd/SRP/_wget
endif

if (-e  $dd/srp.ace) \rm  $dd/srp.ace
foreach srp (`cat $dd/srp.list | gawk '/^SRP/{gsub(/\"/,"",$2);if(n[$2]<1){n[$2]=1;print $2;}}'`)
   if (! -e  $dd/SRP/$srp.html) continue 
   cat  $dd/SRP/$srp.html | gawk -f scripts/SRX_import.2.awk srp=$srp | sed -e 's/\\\"//g' -e s'/\\//g'  >> $dd/srp.ace
end

bin/tacembly SRX_DB <<EOF
  read-models
  parse $dd/srp.ace
  save
EOF
goto phaseLoop

####### Parse the GEO to find the geo->contibutors and affiliations
phaseGEO:
if (! -e $dd/GEO) mkdir $dd/GEO

bin/tacembly SRX_DB <<EOF
  query find project $MAGIC ; >SRR ; >SRP ; >GEO
  list -a -f $dd/geo.list
EOF

if (-e $dd/GEO/_wget) \rm $dd/GEO/_wget
foreach geo (`cat $dd/geo.list | gawk '/^GEO/{gsub(/\"/,"",$2);if(n[$2]<1){n[$2]=1;print $2;}}'`)
   if ( -e  $dd/GEO/$geo.html) continue 
   echo "wget -O "$dd"/GEO/"$geo".html  "'"http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='$geo'"' >>  $dd/GEO/_wget
end

wc $dd/GEO/_wget
source  $dd/GEO/_wget

if (-e  $dd/geo.ace) \rm  $dd/geo.ace
foreach geo (`cat $dd/geo.list | gawk '/^GEO/{gsub(/\"/,"",$2);if(n[$2]<1){n[$2]=1;print $2;}}'`)
   if (! -e  $dd/GEO/$geo.html) continue 
   cat  $dd/GEO/$geo.html | gawk -f scripts/SRX_import.3.awk geo=$geo  >> $dd/geo.ace
end

bin/tacembly SRX_DB <<EOF
  read-models
  parse $dd/geo.ace
  save 
EOF

goto phaseLoop

####### Parse the biosamples to get the attributes (organism, tissues etc)
phaseSample:
if (! -e $dd/BIOSAMPLE) mkdir $dd/BIOSAMPLE

bin/tacembly SRX_DB <<EOF
  query find biosample ! Species
  query find project $MAGIC ; > srr ; >biosample
  list -a -f $dd/biosample.list
EOF


if (-e $dd/BIOSAMPLE/_wget) \rm $dd/BIOSAMPLE/_wget
foreach bs (`cat $dd/biosample.list | gawk '/^Biosample/{gsub(/\"/,"",$2);if(n[$2]<1){n[$2]=1;print $2;}}'`)
   if ( -e  $dd/BIOSAMPLE/$bs.html) continue 
   echo "wget -O "$dd"/BIOSAMPLE/"$bs".html  "'"http://www.ncbi.nlm.nih.gov/biosample/?term='$bs'"' >>  $dd/BIOSAMPLE/_wget
end

wc $dd/BIOSAMPLE/_wget
source  $dd/BIOSAMPLE/_wget

if (-e  $dd/biosample.ace) \rm  $dd/biosample.ace
foreach biosample (`cat $dd/biosample.list | gawk '/^Biosample/{gsub(/\"/,"",$2);if(n[$2]<1){n[$2]=1;print $2;}}'`)
   if (! -e  $dd/BIOSAMPLE/$biosample.html) continue 
   cat  $dd/BIOSAMPLE/$biosample.html | gawk -f scripts/SRX_import.4.awk biosample=$biosample  >> $dd/biosample.ace
end

# cat toto.html |  gawk -f scripts/SRX_import.4.awk biosample=SAMN04584997
if (0) then 
  set dd=DATA/2016_04_19
  cat $dd/biosample.ace | gawk '/^Biosample /{print $2}' | head -500000 > _t1
  \rm _t5 
  foreach b ( `cat _t1` )
    cat $dd/BIOSAMPLE/$b.html | gawk -f scripts/SRX_import.4.awk biosample=$b >> _t5
  end
endif

bin/tacembly SRX_DB <<EOF
  read-models
  parse $dd/biosample.ace
  save 
EOF

goto phaseLoop

####### Parse the SRX get details on the preparation of the library
phaseSRX:
if (! -e $dd/SRX) mkdir $dd/SRX

bin/tacembly SRX_DB <<EOF
  query find SRX
  query find project $MAGIC ; > srr ; >srx
 list -a -f $dd/SRX.list
EOF

if (-e $dd/SRX/_wget) \rm $dd/SRX/_wget

if (1) then
  foreach srx (`cat $dd/SRX.list | gawk '/^SRX/{gsub(/\"/,"",$2);if(n[$2]<1){n[$2]=1;print $2;}}'`)
     if ( -e  $dd/SRX/$srx.html) continue 
     echo "wget -O "$dd"/SRX/"$srx".html  "'"http://www.ncbi.nlm.nih.gov/sra/' $srx'[accn]"' >>  $dd/SRX/_wget
  end
else
  cat  $dd/SRX.list | gawk '/^SRX/{gsub(/\"/,"",$2);if(n[$2]<1){n[$2]=1; srx=$2;gsub(/ /,"",srx);nn++ ;if(nn%200==1)printf("mkdir %s/SRX/%d\n",dd,int(nn/200)); printf ("wget -O %s/SRX/%d/%s.html  \"http://www.ncbi.nlm.nih.gov/sra/%s[accn]\"\n",dd,int(nn/200),srx,srx);}}' dd=$dd  > $dd/SRX/_wget
endif

wc  $dd/SRX/_wget
source  $dd/SRX/_wget

\rm  $dd/srx.ace
foreach ff (`ls  $dd/SRX/*.html $dd/SRX/*/*.html`)
  set srx=`echo $ff | gawk '{n=split ($1,aa,"/");s=aa[n];gsub(".html","",s);print s;}'` 
  cat $ff  | gawk -f scripts/SRX_import.5.awk srx=$srx  >> $dd/srx.ace
end

bin/tacembly SRX_DB <<EOF
  read-models
  parse $dd/srx.ace
  save 
  // query find SRX
  // list -a -f $dd/SRX.list
EOF

goto phaseLoop

######################
## add the file description

phaseFiles:
tbly SRX_DB <<EOF
  query find SRR !File && !Solid && Paired_end && project == $MAGIC
  list -a -f $dd/r2f.p
  query find SRR !File && !Solid && !Paired_end && project == $MAGIC
  list -a -f $dd/r2f.s
  query find SRR !File && Solid && Paired_end && project == $MAGIC
  list -a -f $dd/r2f.p.s
  query find SRR !File && Solid && !Paired_end && project == $MAGIC
  list -a -f $dd/r2f.s.s
EOF

cat $dd/r2f.p  | gawk  '/^SRR/{gsub(/\"/,"",$0);printf("SRR %s\n-D File\nFile fasta/1 %s/SRA/%s_1.fasta.gz\nFile fasta/2 %s/SRA/%s_2.fasta.gz\n\n",$2,dd,$2,dd,$2);}' dd=$dd > $dd/r2f.ace
cat $dd/r2f.s  | gawk  '/^SRR/{gsub(/\"/,"",$0);printf("SRR %s\n-D File\nFile fasta %s/SRA/%s.fasta.gz\n\n",$2,dd,$2,dd,$2);}' dd=$dd >> $dd/r2f.ace
cat $dd/r2f.p.s  | gawk  '/^SRR/{gsub(/\"/,"",$0);printf("SRR %s\n-D File\nFile fasta/1 %s/SRA/%s_1.csfasta.gz\nFile fasta/2 %s/SRA/%s_2.csfasta.gz\n\n",$2,dd,$2,dd,$2);}' dd=$dd >> $dd/r2f.ace
cat $dd/r2f.s.s  | gawk  '/^SRR/{gsub(/\"/,"",$0);printf("SRR %s\n-D File\nFile fasta %s/SRA/%s.csfasta.gz\n\n",$2,dd,$2,dd,$2);}' dd=$dd >> $dd/r2f.ace

echo "pparse  $dd/r2f.ace" | tbly SRX_DB -no_prompt

goto phaseLoop

#############
## add the super titles
phaseTitles:

bin/sra_metadata -db SRX_DB -a -dbEdit
bin/sra_metadata -db SRX_DB -s -dbEdit

tbly SRX_DB <<EOF
  date
  bql -a -o r2b2t.txt select r,b,m from r in class srr, b in r->biosample, m in b->magic_sample2 where m
  date
EOF
tbly SRX_DB <<EOF
  date
  bql -a -o b2t.txt select b,t from b in class biosample, t in b->title
  bql -a -o b2p.txt select b,t,t2 from b in class biosample, t in b->biosample_attribute, t2 in t[1]
  quit
EOF

cat r2b2t.txt | gawk -F '\t' '{printf("SRR %s\nTitle %s\n\n", $1,$3);}' > r2b2t.ace
tbly SRX_DB <<EOF
  pparse r2b2t.ace
  save
  quit
EOF

#############
## count the SRX spots cumulating the spots of the SRR sublibraries

if (! -e ZZZZZ) echo ZZZZZ > ZZZZZ

tace SRX_DB <<EOF
  select -a -o $dd/SRX_count_spots.txt srx,r,s1,s2,s3,s4,s5,s6,s7,s8,s9 from srx in ?srr, r in srx->sublibraries where r,s in r#spots where s,s1 in s[1],s2 in s[2],s3 in s[3],s4 in s[4],s5 in s[5],s6 in s[6],s7 in s[7],s8 in s[8],s9 in s[9]
EOF
cat  $dd/SRX_count_spots.txt ZZZZZ | gawk -F '\t' '{gsub(/\"/,"",$0);}{if($1 != srx) {if(ns>0){printf("SRX %s\nSpots %d ",srx,ns); if(nb>0) { printf(" bases_in_SRA %d ", nb); if(nA >0){printf(" Average_length %d ",nA); if (nSize > 0){printf(" Insert_size %s ", nSize); if( nMates > 0) printf(" spots_with_mate %d ", nMates);}}} printf("\n\n");} ns = 0 ; nb = 0 ; nA = 0 ; nSize = 0 ; nMates = 0 ;} srx =$1;if ($3 > 0) { ns += $3;nb += $5; nMates += $11; nA = (nA * (ns - $3) + $7 * $3)/(ns); nSize = (nSize * (ns - $3) + $9 * $3)/(ns);}}' >  $dd/SRX_count_spots.ace

tbly SRX_DB <<EOF
  parse  $dd/SRX_count_spots.ace
  save
  quit
EOF

goto phaseLoop

#############
## get the papers
phasePapers:

if (! -d $dd/PAPERS) then
  mkdir $dd/PAPERS
  pushd  $dd/PAPERS
  cvs checkout IMPORT_DATA/BIBLIO/PmImport
  mv IMPORT_DATA/BIBLIO .
  popd
endif

tbly SRX_DB <<EOF
  query find paper IS pm* AND NOT citation
  list -a -f $dd/PAPERS/a5.pmnocit.list
EOF
wc  $dd/PAPERS/a5.pmnocit.list
pushd  $dd/PAPERS
  perl  BIBLIO/PmImport/medlineGet.pl < a5.pmnocit.list >! a5.pmnocit.gb
  perl  BIBLIO/PmImport/medline2ace.pl < a5.pmnocit.gb >! a5.pmnocit.preace
popd

cat $dd/PAPERS/a5.pmnocit.preace | grep -v EC_symbol > $dd/PAPERS/a5.pmnocit.ace
tbly SRX_DB << EOF
    pparse   $dd/PAPERS/a5.pmnocit.ace
    save
    quit
EOF

#######  create the sublibraries

tbly SRX_DB << EOF
    query find srx COUNT { >SRR ; File}  > 1
    bql -a -o $dd/multisrx.txt select srx, srr from srx in @, srr in srx->srr where srr#File
    save
    quit
EOF
cat $dd/multisrx.txt | gawk -F '\t' '/^#/{next;} {printf("SRR %s\nSublibraries %s\n\n", $1,$2);}' > $dd/multisrx.ace
tbly SRX_DB << EOF
    pparse  $dd/multisrx.ace
    save
    quit
EOF
tbly SRX_DB << EOF
    query find SRR sublibrary_of
    show -a -f $dd/sublib.ace
    save
    quit
EOF

cat $dd/multisrx.txt ZZZZZ $dd/sublib.ace | gawk '/^ZZZZZ/{zz++;next;}{if(zz<1){srx[$2] = $1;next;}}/^$/{if(ok==1)print;ok=0;next;}/^SRR_download/{next;}/^SRR/{gsub(/\"/,"",$2);ok=0;z = srx[$2];if(length(z)<3)next;ok=1;printf("\nSRR %s\n",z);next;}/^File/{next;}/^Sublibrary_of/{next;}{if(ok==1) print}' >  $dd/sublib.ace2
tbly SRX_DB << EOF
    // pparse  $dd/sublib.ace2
    query find srr biosample
    bql -a -o $dd/title.txt  select srr,t from srr in @, b in srr->biosample , t in b->magic_sample2 
    save
    quit
EOF
cat  $dd/title.txt | gawk -F '\t' '/^#/{next;}{printf ("SRR %s\nTitle %s\n\n", $1,$2);}' > $dd/titles.ace
tbly SRX_DB << EOF
    pparse  $dd/titles.ace
    save
    quit
EOF

goto phaseLoop

#######  done

bly SRX_DB &

#############
## download the actual fasta files

phase_sraDownload:

if (! -d  $dd/SRA ) mkdir $dd/SRA
# check which we already know   // key toto.a2.list

tbly SRX_DB <<EOF
  select  -o $dd/SRA/srrP.txt1  r from p in ?project where p == $MAGIC, r in p->SRR where r#file
  bql -o $dd/SRA/srrP.txt2 select r,p,d from r in @, p in r#paired_end, d in r->SRR_download
EOF
wc $dd/SRA/srrP.txt[12] 
 cat  $dd/SRA/srrP.txt2 | gawk -F '\t' '{srr=$1;pair=$2;file=$3;if( index(file,"http")==1)printf ("%s___%s___%s\n",srr,file,pair);}' > $dd/SRA/srr.todo
wc $dd/SRA/srr.todo

if (-e toto88) \rm toto88
touch toto88
foreach ss (`cat  $dd/SRA/srr.todo`)
  set srr=`echo $ss | gawk -F '___' '{print $1}'`
  if (! -d Fastc/$srr &&  ! -e $dd/SRA/$srr && ! -e $dd/SRA/$srr.fasta.gz &&  ! -e $dd/SRA/$srr'_1'.fasta.gz && ! -e $dd/SRA/$srr.csfasta.gz &&  ! -e $dd/SRA/$srr'_1'.csfasta.gz) then
    echo $ss >> toto88
   if ($phase == sraDownload) scripts/submit $dd/SRA/$srr "scripts/SRX_import_fasta.tcsh $dd/SRA $ss" 
  endif
end
wc toto88

goto phaseLoop
ls DATA/2016_01_08/SRA/*fasta* | gawk '{split($1,aa,"/");split(aa[4],bb,"_");split(bb[1],cc,".");printf("Run %s\n", cc[1]);}' | sort -u > sraReady.list
ls DATA/2016_01_08/SRA/*fasta* | gawk '{split($1,aa,"/");split(aa[4],bb,"_");split(bb[1],cc,".");printf("SRR %s\n", cc[1]);}' | sort -u > sraReady.srr.list
######################
## clean up
######################
## Transfer to RunMasterDB

phase_srr2run:

tbly SRX_DB <<EOF
  query find geo reference && srp
  bql -a -o $dd/geo2srp2ref.txt  select geo,srp,ref from geo in @, srp in geo->srp, ref in geo->reference
  quit
EOF
 cat $dd/geo2srp2ref.txt | gawk '/^"/{printf ("SRP %s\nReference %s\n\n",$2,$3);}' > $dd/geo2srp2ref.ace
 echo "pparse $dd/geo2srp2ref.ace " | tbly  SRX_DB -no_prompt

tbly SRX_DB <<EOF
  query find srr  project == $MAGIC
  // key sraReady.srr.list
  kstore kk
  show -a -f $dd/srr.dump.ace
  query follow srp
  show -a -f $dd/srp.preace
  kget kk
  query follow srx
  show -a -f $dd/srx.preace
  kget kk
  query Nucleic_Acid_Extraction
  show -a -f $dd/srr.nuc.ace Nucleic_Acid_Extraction
  kget kk
  query sraNucleic_Acid_Extraction
  show -a -f $dd/srr.sranuc.ace sraNucleic_Acid_Extraction
  find paper
  show -a -f $dd/srr.paper.preace
  follow abstract
  show -a -f $dd/srr.paper_abstract.ace
  find biosample
  show -a -f $dd/srr.biosample.preace
  kget kk
  query  Sample_builder
  show -a -f $dd/srr.sample_builder.preace Sample_builder
  kget kk
  follow compare
   show -a -f $dd/compare.preace
  quit
EOF

echo ZZZZZ > ZZZZZ
# cat $dd/srr.dump.ace   ZZZZZ $dd/srr.dump.ace | gawk '/^ZZZZZ/{zz++;next;}{if(zz<1){if($1 == "SRR"){ srr=$2;srr2run[srr];}if($1=="Run")srr2run[srr]=srr;next;}}/^Run/{next}/^SRR /{printf("Run %s\nSRR %s\n",srr2run[$2],$2);next;}{print}' > $dd/srr2run.preace
cat $dd/srr.dump.ace > $dd/srr2run.preace

cat <<EOF > $dd/srr2run.awk
/^Adult/{next;}
/^Annotation_problem/{next;}
/^Biosample/{next;}
/^Body_site/{next;}
/^Cap_CAGE_RACE/{print;next;}
/^Cell_line/{next;}
/^Center_name/{next;}
/^Control/{next;}
/^Date_received/{print;next;}
/^Digestive/{next;}
/^ERROR/{next;}
/^Embryo/{next;}
/^Female/{next;}
/^File/{print;next;}
/^Forward/{next;}
/^Gene_selection/{next;}
/^Germline_and_development/{next;}
/^Group/{print ; next;}
/^Head/{next;}
/^Helicos/{print;next;}
/^Illumina/{print;next;}
/^Ion_Torrent/{print;next;}
/^L1/{next;}
/^L2/{next;}
/^L3/{next;}
/^L4/{next;}
/^Any_larva/{next;}
/^Larva/{next;}
/^Library_name/{next;}
/^Magic_author2/{gsub("Magic_author2" ,"Author", \$0);print; next;}
/^Male/{next;}
/^Microbiome/{next;}
/^Mixed_sex/{next;}
/^Nascent_RNA/{print;next;}
/^Nerve/{next;}
/^Oxford_nanopore/{print;next;}
/^Paired_end/{print;next;}
/^Project/{print;next;}
/^Pupa/{next;}
/^RIP_CLIP/{print;next;}
/^RNA/{print;next;}
/^Reference/{print;next;}
/^Roche_454/{print;next;}
/^Run/{next;printf("\n");print;next;}
/^SNP/{print;next;}
/^SOLiD/{print;next;}
/^SRP/{next;}
/^SRR_download/{print; next;}
/^SRR/{printf ("\n\nRun %s\n-D Nucleic_Acid_Extraction\n-D Origin\n-D group\n-D Title\n-D SRR\n-D Sample\n-D SRP\n",\$2); x=substr(\$2,2,3); if(x=="DRR" || x=="SRR" || x=="ERR")printf("SRR %s\n",\$2); next;}
/^SRX/{next;}
/^Sample_name/{next;}
/^Small_RNA/{print;next;}
/^SOLiD/{print;next;}
/^Sorting_title/{print;next;}
/^Species/{print;next;}
/^Spots/{print;next;}
/^Stranded/{next;print "Strand";next;}
/^Sublibraries/{print;next;}
/^Sublibrary_of/{print;next;}
/^Submission_date/{print;next;}
/^Title/{print;next;}
/^Treatment/{next;}
/^Total_RNA/{print;next;}
/^Total/{next;}
/^Union_of/{print;next;}
/^Unspecified_RNA/{print;next;}
/^Warning/{next;}
/^Whole_genome/{print;next;}
/^Whole_organism/{next;}
/^nonStranded/{print;next;}
/^polyA/{print;next;}
/^sraUnspecified_RNA/{next;}
/^sraGenesraGene_selection/{next;}
/^sraGenesraGenesraPolyA/{next;}
/^sraRIP_CLIP/{next;}
/^sraSmall_RNA/{next;}
/^sraUnspecified_RNA/{next;}

END {printf("\n");}
EOF

cat  $dd/srr2run.preace | gawk -f $dd/srr2run.awk >  $dd/srr2run.ace

cat <<EOF > $dd/srr2srr.awk
/^Adult/{next;}
/^Annotation_problem/{next;}
/^Biosample/{next;}
/^Body_site/{next;}
/^Cap_CAGE_RACE/{next;}
/^Cell_line/{next;}
/^Center_name/{next;}
/^Control/{next;}
/^Date_received/{next;}
/^Digestive/{next;}
/^ERROR/{next;}
/^Embryo/{next;}
/^Female/{next;}
/^File/{next;}
/^Forward/{next;}
/^Gene_selection/{next;}
/^Germline_and_development/{next;}
/^Group/{next;}
/^Head/{next;}
/^Helicos/{next;}
/^Illumina/{next;}
/^Ion_Torrent/{next;}
/^L1/{next;}
/^L2/{next;}
/^L3/{next;}
/^L4/{next;}
/^Any_larva/{next;}
/^Larva/{next;}
/^Library_name/{next;}
/^Magic_author2/{next;}
/^Male/{next;}
/^Mixed_sex/{next;}
/^Microbiome/{next;}
/^Nascent_RNA/{next;}
/^Nerve/{next;}
/^Oxford_nanopore/{next;}
/^Paired_end/{next;}
/^Project/{next;}
/^Pupa/{next;}
/^RIP_CLIP/{next;}
/^RNA/{next;}
/^Reference/{next;}
/^Roche_454/{next;}
/^Run/{next;}
/^SNP/{next;}
/^SOLiD/{next;}
/^SRP/{next;}
/^SRR_download/{next;}
/^SRR/{printf("\n");print;next;}
/^SRX/{next;}
/^Sample_name/{next;}
/^Small_RNA/{next;}
/^SOLiD/{next;}
/^Sorting_title/{next;}
/^Species/{print;next;}
/^Spots/{next;}
/^Stranded/{print;next;}
/^Sublibraries/{next;}
/^Sublibrary_of/{next;}
/^Submission_date/{next;}
/^Title/{next;}
/^Treatment/{next;}
/^Total_RNA/{next;}
/^Total/{next;}
/^Union_of/{next;}
/^Unspecified_RNA/{next;}
/^Warning/{next;}
/^Whole_genome/{next;}
/^Whole_organism/{next;}
/^nonStranded/{next;}
/^polyA/{next;}
/^sraUnspecified_RNA/{next;}
/^sraGenesraGene_selection/{next;}
/^sraGenesraGenesraPolyA/{next;}
/^sraRIP_CLIP/{next;}
/^sraSmall_RNA/{next;}
/^sraUnspecified_RNA/{next;}

END {printf("\n");}
EOF
cat  $dd/srr2run.preace | gawk -f $dd/srr2srr.awk >  $dd/srr2srr.ace

cat <<EOF > $dd/srp2run.awk
/^Abstract/{print;next;}
/^Description/{print;next;}
/^GEO/{next;}
/^Identifier/{print;next;}
/^Reference/{print;next;}
/^SRP/{printf("\n");print;print "-D T" ; next;}
/^SRR/{print;printf ("Run %s\n",\$2);next;}
/^Species/{print;next;}
/^Title/{print;next;}
END {printf("\n");}
EOF

cat  $dd/srp.preace | gawk -f $dd/srp2run.awk >  $dd/srp.ace

cat <<EOF > $dd/srx2run.awk
/^SRX/{printf("\n");print; print "-D T" ; next;}
/^Construction_protocol/{print;next;}
/^Design/{print;next;}
/^SRR/{print;next;}
/^Submitted_by/{print;next;}
/^Spots/{print;next;}
/^Title/{print;next;}
END {printf("\n");}
EOF

cat  $dd/srx.preace | gawk -f $dd/srx2run.awk >  $dd/srx.ace
cat  $dd/compare.preace |  sed -e 's/^SRR[ \t]/Runs /' > $dd/compare.ace
tbly MetaDB -no_prompt <<EOF
  read-models
  query find project $MAGIC ; > run
  edit -D group
  query find project $MAGIC ; > compare
  kill
  query find project $MAGIC 
  kill
  pparse  $dd/srr2run.ace
  pparse  $dd/srr2srr.ace
  pparse  $dd/srp.ace
  pparse  $dd/srx.ace
  pparse  $dd/srr.nuc.ace
  pparse  $dd/srr.sranuc.ace
  pparse  $dd/compare.ace
  save
  quit
EOF


cat <<EOF > $dd/srx2paper.awk
/^GEO/{next;}
{print}
END {printf("\n");}
EOF

cat  $dd/srr.paper.preace | gawk -f $dd/srx2paper.awk >  $dd/srr.paper.ace
echo "pparse  $dd/srr.paper.ace" | tbly MetaDB -no_prompt
echo "pparse  $dd/srr.paper_abstract.ace" | tbly MetaDB -no_prompt

cat <<EOF > $dd/srx2paper.awk
/^GEO/{next;}
{print}
END {printf("\n");}
EOF

cat $dd/srr.biosample.preace ZZZZZ  $dd/srr.sample_builder.preace | gawk '{line++;}/^ZZZZZ/{zz++;inside=0;next;}/^$/{inside=0;print;next;}/^Biosample_/{next;}/^Sample_name/{next;}/^Biosample/{bio=$2;if(length(bio)>0){printf("Sample %s\n-D T \n",$2);inside=1;}next;}/^Magic_sample2/{gsub("Magic_sample2","Title", $0);if(inside==1)print;next;}/^SRR/{if(zz<1){print; srr2bio[$2]=bio;}next;}{if(zz==1 && inside==1)print;next;}' > $dd/srr.sample.ace

echo "pparse  $dd/srr.sample.ace" | tbly MetaDB -no_prompt


time tbly MetaDB <<EOF
  query find project $MAGIC ; >run
  find run
  bql -a -o $dd/r2s2t.txt  select r,srr,s,t from r in @ , srr in r->srr, s in srr->sample, t in s->title 
  query find project $MAGIC ; >run
  find run
  bql -a -o $dd/r2s2t2.txt  select r,srr,s,t from r in @ , sub in r->sublibraries, srr in sub->srr, s in srr->sample, t in s->title 

  query find project $MAGIC ; >run
  find run
  bql -a -o $dd/r2s2t.txt1  select r,srr,s,t from r in @ , srr in r->srr, s in srr->sample, t in s->title 
  query find project $MAGIC ; >run
  find run
  bql -a -o $dd/r2s2t2.txt1  select r,srr,s,t from r in @ , sub in r->sublibraries, srr in sub->srr, s in srr->sample, t in s->title 

EOF

cat $dd/r2s2t.txt  $dd/r2s2t2.txt | gawk -F '\t' '{gsub(/SRR:/,"",$2);gsub(/Sample:/,"",$3);gsub(/Text:/,"",$4);printf("Run %s\nSample %s\nTitle %s\n\n", $1, $3, $4);}' | grep -v NULL > $dd/r2s2t.ace
echo "pparse  $dd/r2s2t.ace" | tbly MetaDB -no_prompt

tbly MetaDB <<EOF
  query find run ; sublibrary_of && Group
  bql -a -o rSubGr.txt  select r,sub,g from r in @ , sub in r->sublibrary_of, g in r->group
  undo
  bql -a -o rSubGr.txt1  select r,sub,g from r in @ , sub in r->sublibrary_of, g in r->group
EOF
cat  rSubGr.txt | gawk -F '\t' '{printf ( "Run %s\n-D Union_of %s\nUnion_of %s\n\n",$3,$1,$2 ) ; }' >  rSubGr.fix.ace

echo "pparse  rSubGr.fix.ace" | tbly MetaDB -no_prompt
\rm  rSubGr.fix.ace  rSubGr.txt

goto phaseLoop

#######################################################
## Inside SRX_db, transfer the manual annoations from SRR to SRX
phase_srr2srx:
tbly SRX_DB <<EOF
  query find project $MAGIC ; >srr ;
  query find SRR
  query file && sublibrary_of && ! union_of
  select -a -o srr2srx.txt  srr,srx from srr in @, srx in srr->sublibrary_of
  show -a -f $MAGIC.srr.ace
EOF

echo ZZZZZ > ZZZZZ
cat  srr2srx.txt  ZZZZZ  $MAGIC.srr.ace | gawk -f $dd/srr2srx.awk >  $dd/$MAGIC.srr2srx.ace
tbly SRX_DB <<EOF
  // pparse $dd/$MAGIC.srr2srx.ace
  save
  quit
EOF


goto phaseLoop

#######################################################


##########
# danielle score for droso runs
cat RESULTS/scoring.txt | sed -e 's/\r//' | gawk -F '\t' '{run=$2;k=0;x=$4;if(x>10000)k+=5;else if (x>8000)k+4;else if (x>5000)k+=2;else k++;x=$5;if(x<5)k+=3;else if (x<10)k+=2;else if (x<20)k+=.5;else if (x>50)k-=1;x=$6;if(x>4500)k+=4;else if (x>3000)k+=4;else if (x >2000)k+=1.5;else if (x>1000)k+=.5; x=$7;z="";if(x>99)z="++";else if (x >94)z="+";else if (x<1)z="--";else if (x<6)z="-";printf("Run %s\t%f\t%s\n",run,k,z);}' >  RESULTS/scoring.2.txt

cat RESULTS/scoring.2.txt | gawk '/^Rdm/{if($3=="++" && $2 >=11)print "Run " $1;}' > toto.Aplus.list
cat RESULTS/scoring.2.txt | gawk '/^Rdm/{if($3=="--" && $2 >=11)print "Run " $1;}' > toto.Aminus.list
cat RESULTS/scoring.2.txt | gawk -F '\t' '/^#/{next;}{printf("Run %s\nScore %s \"%s\"\n\n",$1,$2,$3);}' >  RESULTS/scoring.2.ace

#############


phaseLoop:
 date
