/^#/{next;}
{if ($8 == "U_introns"){intron = substr($11,4);n[intron]=1;if($10 ==1)n_u[intron]+=$3;else if ($10>0)n_nu[intron]+=$3/$10;}}END{for(intron in n) printf("%s\t%s\t%d\t%.1f\n", intron, run, n_u[intron],n_nu[intron]);
}
