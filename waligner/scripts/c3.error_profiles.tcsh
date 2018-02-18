#!bin/tcsh -f

set title=$1
set toto_ace=$2
set toto_out=$3
set toto=$toto_out.mismatch.txt
echo $toto

    date > $toto
    echo  "Distribution of mismatches in best unique alignment in sample $title" >> $toto
    echo "cat $toto_ace | gawk  -f scripts/c3.error_profiles.awk title=$title"
    cat $toto_ace | gawk  -f scripts/c3.error_profiles.awk title="$title" >> $toto


