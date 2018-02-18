#!/bin/csh -f
# @(#)osp1.csh	1.3 6/25/97
# this extra layer is for the nbenefit of cdscriptpaipe
# do not remove
./osp2.pl $*
#sleep 1
exit 0


# debug version of same
echo osp2.csh running >! /tmp/xxx2
echo $* >> /tmp/xxx21
#
./osp2.pl $*   >! /tmp/xxx31
cat /tmp/xxx31
sleep 1
#
