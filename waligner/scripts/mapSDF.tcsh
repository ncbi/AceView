#!bin/tcsh -ef

set intronList=$1
set outFile=$2

set stranded=""
if ($3 == stranded) set stranded="-stranded"

set CDS=EBI_seqc

# genome gene_and_new_exon new_and_known_exon known_exon RefSeq

    set zcR="RefSeq_CDS:tmp/METADATA/RefSeq.ns.cds.sponge"
    set zcA="AceView_CDS:tmp/METADATA/av.ns.cds.sponge"
    set zcD="Discovery_CDS:tmp/METADATA/seqc.ns.cds.sponge"
    set zcE="EBI_CDS:tmp/METADATA/EBI.ns.cds.sponge"

    set ztR="RefSeq_transcripts:tmp/METADATA/RefSeq.ns.mrna.sponge"
    set ztA="AceView_transcripts:tmp/METADATA/av.ns.mrna.sponge"
    set ztD="Discovery_transcripts:tmp/METADATA/seqc.ns.mrna.sponge"
    set ztE="EBI_transcripts:tmp/METADATA/EBI.ns.mrna.sponge"

    set zgR="RefSeq_genes:tmp/METADATA/RefSeq.ns.gene.sponge"
    set zgA="AceView_genes:tmp/METADATA/av.ns.gene.sponge"
    set zgD="Discovery_genes:tmp/METADATA/seqc.ns.gene.sponge"
    set zgE="EBI_genes:tmp/METADATA/EBI.ns.gene.sponge"

    set zgenome="Genome:tmp/METADATA/genome.ns.sponge"

    set mask="$ztR,$ztA,$ztD,$zgR,$zgA,$zgD,$zgenome"
    
    set                       mask="$zcR,$zcA,$zcE,$zcD,$ztR,$ztA,$ztE,$ztD,$zgR,$zgA,$zgE,$zgD,$zgenome"
    if ($CDS == seqc_EBI) set mask="$zcR,$zcA,$zcD,$zcE,$ztR,$ztA,$ztD,$ztE,$zgR,$zgA,$zgD,$zgE,$zgenome"
    

    if (! -e $outFile.sdfPosition.txt) then
      echo "bin/geneelements -spongeFile $mask -mapSDF $intronList $stranded -o $outFile"
            bin/geneelements -spongeFile $mask -mapSDF $intronList $stranded -o $outFile
    endif
