#!bin/tcsh -f

set phase=$1
set chrom=$2

if ($phase == k1) goto phasek1
if ($phase == k2) goto phasek2
if ($phase == k9) goto phasek9

if ($species == hs && $chrom != chr4) goto phaseLoop

############################################################
## Phase k1 : recover the available kantor

phasek1:

echo -n "Start of phase k1 chrom=$chrom "
date

# construct the list of existing kantor files, and parse them
  echo 'read-models' >  Kantor/$chrom/_r
 
  if (-d Kantor/$chrom/tmp) then
    pushd Kantor/$chrom/tmp/$species_kantor.data/$chrom
      pwd
      ls
      foreach dd (`ls`)
        cat $dd/ace.* > ../../../ace.$dd 
        # \rm -rf $dd
      end
    popd 
  endif

  set ok=0

  ls -ls Kantor/$chrom/ace.*
  foreach ff (`ls  Kantor/$chrom/ace.*`)
    # a hack to clean a current PSort bug 
    cat $ff | gawk '/^Psort Localization \"\"/{next;}{print}' > $ff.bis
    mv  $ff.bis $ff

    echo "pparse $ff" >>  Kantor/$chrom/_r
    set ok=1
  end
  echo save >> Kantor/$chrom/_r
  echo quit >> Kantor/$chrom/_r
echo aaa
  if ($ok == 1 && -e tmp/XH$chrom/database/log.wrm && ! -e tmp/XH$chrom/database/lock.wrm) then

    if (-e tmp/XH$chrom/k9.megaRun.done) \rm tmp/XH$chrom/k9.megaRun.done

    tacembly tmp/XH$chrom <  Kantor/$chrom/_r

    tacembly tmp/XH$chrom << EOF
      query find pfam IS Pox_polyA_pol
      kill
      query find Kantor psr ; > product ; >mrna 
      acembly
        cdna_kantor // kantorizes active mrna set
        quit
      // dna mmMrna.$chrom.dna
      save
      find model
      list -a -f tmp/XH$chrom/k1.kantor_parse.done
      quit
EOF
    if (! -d Kantor/$chrom/done) mkdir Kantor/$chrom/done
    mv Kantor/$chrom/ace.* Kantor/$chrom/done
  endif
touch tmp/XH$chrom/k1.kantor_parse.done
goto phaseLoop

############################################################
############################################################
## Phase k9 : megaRun this chromosome

#to change the kantor because of the model change

phasek9:

echo -n "Start of phase k9 megaRun chrom=$chrom "
date

  # setenv chrom is used by megaRun
  setenv chrom $chrom
  setenv ici `pwd`
  pushd tmp/XH$chrom

  setenv port 1235901

    if ($chrom == "1") setenv port 1235901
    if ($chrom == "2") setenv port 1235002
    if ($chrom == "2L") setenv port 2235002
    if ($chrom == "2R") setenv port 3235002
    if ($chrom == "3") setenv port 1235003
    if ($chrom == "3L") setenv port 2235003
    if ($chrom == "3R") setenv port 3235003
    if ($chrom == "4") setenv port 1235004
    if ($chrom == "5") setenv port 1235005
    if ($chrom == "6") setenv port 1235006
    if ($chrom == "7") setenv port 1235007
    if ($chrom == "8") setenv port 1235008
    if ($chrom == "9") setenv port 1235009
    if ($chrom == "10") setenv port 1235010
    if ($chrom == "11") setenv port 1235011
    if ($chrom == "12") setenv port 1235012
    if ($chrom == "13") setenv port 1235013
    if ($chrom == "14") setenv port 1235014
    if ($chrom == "15") setenv port 1235015
    if ($chrom == "16") setenv port 1235016
    if ($chrom == "17") setenv port 1235017
    if ($chrom == "18") setenv port 1235018
    if ($chrom == "19") setenv port 1235019
    if ($chrom == "20") setenv port 1235020
    if ($chrom == "21") setenv port 1235021
    if ($chrom == "22") setenv port 1235022
    if ($chrom == "X") setenv port 1235023
    if ($chrom == "Y") setenv port 1235024
    if ($chrom == "Un") setenv port 1235025
    if ($chrom == "MT") setenv port 1235026
    if ($chrom == "Pltd") setenv port 1235027

    if ($chrom == "chr1") setenv port 1236901
    if ($chrom == "chr2") setenv port 1236002
    if ($chrom == "chr3") setenv port 1236003
    if ($chrom == "chr4") setenv port 1236004
    if ($chrom == "chr5") setenv port 1236005
    if ($chrom == "chr6") setenv port 1236006
    if ($chrom == "chr7") setenv port 1236007
    if ($chrom == "chr8") setenv port 1236008
    if ($chrom == "chr9") setenv port 1236009
    if ($chrom == "chr10") setenv port 1236010
    if ($chrom == "chr11") setenv port 1236011
    if ($chrom == "chr12") setenv port 1236012

    set iiStart=0

  iciStart: 
    echo "warning: at Try 0 we expect the error message connect: Connection refused"

    echo -n "Try $iiStart "  
    $ici/bin/taceclient a:localhost:$port  <<EOF >! k9.test
      date
      quit
EOF
    @ iiStart = $iiStart + 1
    set startOk = `grep -c 'Active Objects' k9.test `
    echo " ok=$startOk"  
    if ($startOk > 0) goto iciStartDone 
    if ($iiStart > 1 && $iiStart < 60) then
      sleep 10
      goto iciStart
    endif


    echo "launching  taceserver XH$chrom $port"
    #    cat /etc/hosts | grep $HOST | gawk '{printf("%s w\n", $1)}' >>  wspec/acetcp_access.wrm

    $ici/bin/taceserver . $port 3600:3600:0 &
    echo 'sleep 10'
    sleep 10
    goto iciStart

  iciStartDone:

    popd
    touch  tmp/XH$chrom/k9.megaRun.start
    pushd Kantor/$chrom
      if (! -d tmp) mkdir tmp
      if (! -d tmp/$species_kantor.data) mkdir tmp/$species_kantor.data
      setenv megaRun ~/MEGA3/scripts/megaRun
      echo -n 'starting MEGA3/scripts/megaRun'
      date
      $megaRun a:localhost:$port  psort $species_kantor -nowrite 
      $megaRun a:localhost:$port acekog $species_kantor  -nowrite
      $megaRun a:localhost:$port   pfam $species_kantor -nowrite 
      $megaRun a:localhost:$port blastp $species_kantor  -nowrite
      # $megaRun a:localhost:$port oligo $species_kantor  -nowrite 
      # $megaRun a:localhost:$port acekog_n $species_kantor  -nowrite
      date
    popd

    $ici/bin/taceclient a:localhost:$port <<EOF > /dev/null
      shutdown now
      quit
EOF
    if (-e tmp/XH$chrom/k1.kantor_parse.done) \rm tmp/XH$chrom/k1.kantor_parse.done
    touch  tmp/XH$chrom/k9.megaRun.done
    \rm tmp/XH$chrom/k9.megaRun.start
  endif
  popd
end

echo -n 'End of phase k9'
date

goto phaseLoop

############################################################
phaseLoop:
 echo "phase $phase done"
 exit 0

