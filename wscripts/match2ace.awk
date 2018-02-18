/^SCORE/ { score = $3 ; gsub (/,/,"",score) ; mm = $6 ;}
/^Match/ {
 x1 = $4 ; x2 = $6 ;
 y1 = $10 ; y2 = $12 ;
 seq1 = $8 ; seq2 = $14 ;

 printf("Sequence \"%s\"\n",seq1) ;
 printf("Homol \"%s\" ICA_Match  %d %d %d %d %d\n\n",
	seq2, score, x1, x2, y1, y2) ;
 printf("Sequence \"%s\"\n",seq2) ;
 printf("Homol \"%s\" ICA_Match  %d %d %d %d %d\n\n",
	seq1, score, y1, y2, x1, x2) ;
}
