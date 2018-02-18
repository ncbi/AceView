/^Sequence/ {sequence=$2; print ; next ; }
/^Accession/ {printf("!"); print ; next ;}
/^Database/ {print ; next ;}
/^From_database/ {printf("!"); print ; next ;}
/^Lab_host/ {printf("!"); print ; next ;}
/^Medline_acc/ {printf("!"); print ; next ;}
/^Nucleotide_id/ {printf("!"); print ; next ;}
/^Organ/ {printf("!"); print ; next ;}
/^Paper/ { print ; next ;}
/^Reference/ {  print ; next ;}
/^Site_/ {printf("!"); print ; next ;}
/^Tissue/ {printf("!Tissue"); print ; next ;}
/^gene/ {printf("!"); print ; next ;}
/^mat_peptide/ {printf("!"); print ; next ;}
/^misc_/ {printf("!"); print ; next ;}
/^polyA_signal/ {printf("!"); print ; next ;}
/^polyA_site/ {printf("PolyA_after_base %s\n", $2); next ;}
/^prim/ {printf("!"); print ; next ;}
/^protein_bind/ {printf("!"); print ; next ;}
/^repeat_/ {printf("!"); print ; next ;}
/^sig_peptide/ {printf("!"); print ; next ;}
/^terminator/ {printf("!"); print ; next ;}
/^transit_peptide/ {printf("!"); print ; next ;}
/^variation/ {printf("!"); print ; next ;}
/^MGC_Remark/ {printf("!"); print ; next ;}
{ print ;}
