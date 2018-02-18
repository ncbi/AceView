#!/bin/csh -f
# parse models.wrm and generate
# tags.tmp and classes.tmp from it
# Otto Ritter 3.11.93
#
#  $Id: models2tags.csh,v 1.1.1.1 2002/07/19 20:23:38 sienkiew Exp $

\rm classes.tmp
\rm tags.tmp

cat  $1 | nawk '{gsub(/\/\/.*$/,""); print}'  | \
grep -v '^$' | \
 tr -cs '\_\-A-Za-z0-9?' '\012'   | \
egrep -vi 'Text|Int|Float|XREF|ANY|FREE|UNIQUE|REPEAT'  |\
sort -u  | \
nawk 'BEGIN {while ((getline <"classes.wrm") >0) \
               if ($1=="#define"){ class[$2] = 1; ClassOffset =$3+1}\
             close("classes.wrm"); \
             while ((getline <"tags.wrm") >0) \
               if ($1=="#define"){ tag[$2] = 1; TagOffset = $3+1}\
             close("tags.wrm"); \
             while ((getline <"sysclasses.wrm") >0) \
               if ($1=="#define"){ class[$2] = 1;}\
             close("sysclasses.wrm"); \
             while ((getline <"systags.wrm") >0) \
               if ($1=="#define"){ tag[$2] = 1;}\
             close("systags.wrm"); \
} \
    {cl = "_V" $0; x = gsub(/\?/,"",cl); } \
    {tg = "_" $0;  x = gsub(/\?/,"",tg); } \
    /^\?/ { if (cl in class) ;else print ("#define",cl,ClassOffset++)>>"classes.tmp"}; \
    { if (tg  in tag) ;else print ("#define",tg,TagOffset++)>>"tags.tmp" }'






