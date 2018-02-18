#!bin/tcsh -ef

set intronList=$1
set outFile=$2

set stranded=""
if ($3 == stranded) set stranded="-stranded"

set CDS=EBI_seqc

# genome gene_and_new_exon new_and_known_exon known_exon RefSeq

set fr=ns
    set zcR="RefSeq_CDS:tmp/METADATA/RefSeq_cds.$fr.sponge"
    set zcA="AceView_CDS:tmp/METADATA/av_cds.$fr.sponge"
    set zcD="Discovery_CDS:tmp/METADATA/seqc_cds.$fr.sponge"
    set zcE="EBI_CDS:tmp/METADATA/EBI_cds.$fr.sponge"

    set ztR="RefSeq_transcripts:tmp/METADATA/RefSeq.$fr.sponge"
    set ztA="AceView_transcripts:tmp/METADATA/av.$fr.sponge"
    set ztD="Discovery_transcripts:tmp/METADATA/seqc.$fr.sponge"
    set ztE="EBI_transcripts:tmp/METADATA/EBI.$fr.sponge"

    set zgR="RefSeq_genes:tmp/METADATA/RefSeq.$fr.gene.sponge"
    set zgA="AceView_genes:tmp/METADATA/av.$fr.gene.sponge"
    set zgD="Discovery_genes:tmp/METADATA/seqc.$fr.gene.sponge"
    set zgE="EBI_genes:tmp/METADATA/EBI.$fr.gene.sponge"

    set zgenome="Genome:tmp/METADATA/genome.ns.sponge"

    set mask="$ztR,$ztA,$ztD,$zgR,$zgA,$zgD,$zgenome"
    
    set                       mask="$zcR,$zcA,$zcE,$zcD,$ztR,$ztA,$ztE,$ztD,$zgR,$zgA,$zgE,$zgD,$zgenome"
    if ($CDS == seqc_EBI) set mask="$zcR,$zcA,$zcD,$zcE,$ztR,$ztA,$ztD,$ztE,$zgR,$zgA,$zgD,$zgE,$zgenome"
    

    if (! -e $outFile.sdfPosition.txt) then
      echo "bin/geneelements -spongeFile $mask -mapSDF $intronList $stranded -o $outFile"
            bin/geneelements -spongeFile $mask -mapSDF $intronList $stranded -o $outFile
    endif
