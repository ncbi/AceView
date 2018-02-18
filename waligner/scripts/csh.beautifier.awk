/^#/{line++;printf("%4d %2d %2d\t",line,n1,n2);print;next;}
/ then/{if ($NF == "then") n2++;}
/else if/{if ($NF == "then") n2--;}
/endif/{if ($1 == "endif") n2--;}
/foreach/{n1++;}
/ while/{n1++;}
/end/{if ($1 == "end")n1--;}

{line++;printf("%4d %2d %2d\t",line,n1,n2);print}
