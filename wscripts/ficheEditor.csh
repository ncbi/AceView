#!/bin/csh -f

echo "java -classpath /net/vesta/export/home/simonyan/yk/wscripts/fiche ficheEditor $2 $3 $4 $5 "

setenv MACHINE `uname`
if ($MACHINE == "OSF1") then
  setenv MACHINE "alpha"
  setenv java /usr/opt/java122/bin/java
else
  setenv MACHINE `mach`
  if ($MACHINE == "sparc") then
    setenv java /usr/java1.2/bin/java
  else
    setenv java java
  endif
endif

$java -classpath  ~/yk/wscripts/fiche ficheEditor $2 $3 $4 $5


# Finally, mv $2.b to $2.out and signal acedb
#mv $2.b $2.out 

if ($1 != "0") kill -USR1 $1

