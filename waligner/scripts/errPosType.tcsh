#!bin/tcsh -ef

set target=$1
set type="$2"
set out="$3"
  
set minPlatform=$4
set minZ=$5

 set type2=`echo $type | sed -e 's/a/A/' -e 's/t/T/' -e 's/g/G/' -e 's/c/C/' -e 's/p2/insert/'  -e 's/m2/delete/'`

 date >! $out
 echo "$type2 modification seen N times by at least $minPlatform platform(s)" >>  $out
 gunzip -c tmp/ERR_$target/*/*/*errPosTypeCover.txt.gz | gawk -F '\t' -f bin/errPosMap.awk type=$type type2=$type2 minZ=$minZ minProtocol=$minPlatform >> $out
