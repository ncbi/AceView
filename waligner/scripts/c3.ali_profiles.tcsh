#!bin/tcsh -f

set title=$1
set toto_ace=$2
set toto_out=$3
    
set toto=$toto_out.ali_quality.txt

# the awk script only uses: Cumulated_mismatches

    echo $toto

    date > $toto
    echo "Per cent mismatch $title" >> $toto
    echo -n "Percent mismatch" >> $toto
    cat  $toto_ace | gawk  -f scripts/c3.ali_profiles.awk BP=0 tag=Per_cent_mismatch >> $toto
    echo >> $toto
    date >> $toto
    echo "Count mismatch $title" >> $toto
    echo -n "Count mismatch" >> $toto
    cat  $toto_ace | gawk  -f scripts/c3.ali_profiles.awk BP=0 tag=Count_mismatch >> $toto
    echo >> $toto
    date >> $toto
    echo "Count ambiguous $title" >> $toto
    echo -n "Count ambiguous" >> $toto
    cat  $toto_ace | gawk  -f scripts/c3.ali_profiles.awk BP=0 tag=Count_ambiguous >> $toto
    echo >> $toto
    date >> $toto
    echo "Percent aligned $title" >> $toto
    echo -n "Per cent aligned" >> $toto
    cat  $toto_ace | gawk  -f scripts/c3.ali_profiles.awk BP=0 tag=Per_cent_aligned >> $toto
    echo >> $toto
    date >> $toto
    echo "Aligned length $title" >> $toto
    echo -n "Aligned length" >> $toto
    cat  $toto_ace | gawk  -f scripts/c3.ali_profiles.awk BP=1 tag=Aligned_length >> $toto
    echo >> $toto
    date >> $toto
    echo "Clipped length $title" >> $toto
    echo -n "Clipped length" >> $toto
    cat  $toto_ace | gawk -f scripts/c3.ali_profiles.awk BP=1 tag=Clipped_length >> $toto

    echo >> $toto
    date >> $toto






