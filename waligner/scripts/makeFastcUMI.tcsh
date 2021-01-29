#!/bin/tcsh -f

set run=$1
set xx=$2

if (! -d Fastc/$run) goto done

if (! -d Fastc/$run/tmp) mkdir  Fastc/$run/tmp

# extract the bar code
foreach lane (`cat Fastc/$run/List | head -1`)
  zcat Fastc/$lane/