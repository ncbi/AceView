#!bin/tcsh -f
# do not use -ef in this scripts, some files are not always constructed, but we attempt to cp them

set action=$1
set target=$2
set nwait=0


# export can be done in parallel
# notice that the file name is relative to -db database dir
if ($action == export) then

  set tt=$target
  if ($target == seqc) set tt="per SEQC gene" 
  if ($target == RefSeq) set tt="per RefSeq gene" 
  if ($target == UCSC) set tt="per UCSC gene" 
  if ($target == av) set tt="per AceView gene" 
  if ($target == ERCC) set tt="per ERCC" 
  
  #\rm RESULTS/$MAGIC.geneTable.*

  if (! -d RESULTS/Expression)  mkdir  RESULTS/Expression

  \rm GeneIndexDB/*RefSeq.noGeneId*txt

  foreach uu (u nu )

  if ($target == ERCC && $uu == nu) continue 

  if ($uu == u) then
    \cp GeneIndexDB/$MAGIC.geneTable.$uu.ERCC.ERCC1.Index.txt    RESULTS/Expression/$MAGIC.ERCC1.index.$uu.txt
    \cp GeneIndexDB/$MAGIC.geneTable.$uu.ERCC.ERCC2.Index.txt    RESULTS/Expression/$MAGIC.ERCC2.index.$uu.txt
    \cp GeneIndexDB/$MAGIC.geneTable.$uu.ERCC.noGeneId.Index.txt    RESULTS/Expression/$MAGIC.ERCC.index.$uu.txt
    \cp GeneIndexDB/$MAGIC.geneTable.$uu.ERCC.noGeneId.Reads.txt    RESULTS/Expression/$MAGIC.ERCC.Reads.$uu.txt
    \cp GeneIndexDB/$MAGIC.geneTable.$uu.ERCC.noGeneId.kb.txt   RESULTS/Expression/$MAGIC.ERCC.kb.$uu.txt
  endif
 
endif


