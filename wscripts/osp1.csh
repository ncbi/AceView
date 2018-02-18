#!/bin/csh -f
# @(#)osp1.csh	1.3 6/25/97
# param 1 should be the sequence name; 2,3,4,the maxScore, min max lengths 
#                                      5,6 TmMin TmMax the oligos,
#                                      7 sequence || mRNA, used in the returned ace file
#                                      8 the fasta file
# ./osp1.pl -m $2 -l $3,$4 -t $5 -T $6 -s $8 | gawk -f ./osp1.awk zzseq=$1 seqMrna=$7 
# exit 0


# debug version of same
echo osp1.csh running >! /tmp/xxxaa

set dd="."
if (-d /home/mieg/ace/wscripts) set dd=/home/mieg/ace
if ($?ACEDB) set dd=$ACEDB
echo "$dd/wscripts/osp1.pl -m $2 -l $3,$4 -t $5 -T $6 -s $8"  >> /tmp/xxxaa 
#

$dd/wscripts/osp1.pl -m $2 -l $3,$4 -t $5 -T $6 -s $8   >! /tmp/xxx2
gawk -f $dd/wscripts/osp1.awk zzseq=$1 seqMrna=$7 /tmp/xxx2 >!  /tmp/xxx3
cat /tmp/xxx3

#
