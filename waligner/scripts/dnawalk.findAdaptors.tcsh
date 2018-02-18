#!bin/tcsh -f

set run=$1
set read=$2

set r2=""
if ($read == 2) set r2='-read2'

echo " gunzip -c Fastc/$run/f*.fastc.gz  | bin/dnawalk -prefix -walk -wordLength 11 -vary 30 -run $run $r2 -o tmp/Adaptors/$run/$run.$read"
gunzip -c Fastc/$run/f*.fastc.gz  | bin/dnawalk -prefix -walk -wordLength 11 -vary 30 -run $run $r2 -o tmp/Adaptors/$run/$run.$read
