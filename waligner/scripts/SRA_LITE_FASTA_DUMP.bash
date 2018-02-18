#!/bin/bash
SRA_TOOLKIT=/panfs/traces01/trace_software/toolkit/centos64
PATH="$SRA_TOOLKIT/bin:$PATH"
LD_LIBRARY_PATH="$SRA_TOOLKIT/lib:$LD_LIBRARY_PATH"
export PATH LD_LIBRARY_PATH
/panfs/traces01/trace_software/toolkit/centos64/bin/fastq-dump  $*
# end
