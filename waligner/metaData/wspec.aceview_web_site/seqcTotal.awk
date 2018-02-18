/^Element/ {for (k in n1){if(n1[k]>0)printf("%s %d %d\n",k, n1[k],n2[k]);n1[k]=0;n2[k]=0;}printf("\n");print;next;}



/^ILM_A/ { k="X_pA_A"; n1[k]+= $2; n2[k] += $3; k="X_A"; n1[k]+= $2; n2[k] += $3;  next;}
/^ILM_B/ { k="X_pA_B"; n1[k]+= $2; n2[k] += $3; k="X_B"; n1[k]+= $2; n2[k] += $3;  next;}
/^ILMnoSt_A/ {  k="X_ns_A"; n1[k]+= $2; n2[k] += $3;  next;}
/^ILMnoSt_B/ {  k="X_ns_B"; n1[k]+= $2; n2[k] += $3;  next;}
/^ILM_totA/ { k="X_tot_A"; n1[k]+= $2; n2[k] += $3; k="X_A"; n1[k]+= $2; n2[k] += $3;  next;}
/^ILM_totB/ { k="X_tot_B"; n1[k]+= $2; n2[k] += $3; k="X_B"; n1[k]+= $2; n2[k] += $3;  next;}
/^ILMnoSt_totA/ {  k="X_tot_ns_A"; n1[k]+= $2; n2[k] += $3;  next;}
/^ILMnoSt_totB/ {  k="X_tot_ns_B"; n1[k]+= $2; n2[k] += $3;  next;}

/^LIF_A/ { k="X_pA_A"; n1[k]+= $2; n2[k] += $3; k="X_A"; n1[k]+= $2; n2[k] += $3;  next;}
/^Anti_ILM_totB/ { k="Anti_X_tot_B"; n1[k]+= $2; n2[k] += $3; k="Anti_X_B"; n1[k]+= $2; n2[k] += $3;  next;}
/^Anti_ILMnoSt_totA/ {  k="Anti_X_tot_ns_A"; n1[k]+= $2; n2[k] += $3;  next;}
/^Anti_ILMnoSt_totB/ {  k="Anti_X_tot_ns_B"; n1[k]+= $2; n2[k] += $3;  next;}

/^Anti_LIF_A/ { k="Anti_X_pA_A"; n1[k]+= $2; n2[k] += $3; k="Anti_X_A"; n1[k]+= $2; n2[k] += $3;  next;}
/^Anti_LIF_B/ { k="Anti_X_pA_B"; n1[k]+= $2; n2[k] += $3; k="Anti_X_B"; n1[k]+= $2; n2[k] += $3;  next;}
/^Anti_LIF_totB/ { k="Anti_X_tot_B"; n1[k]+= $2; n2[k] += $3; k="Anti_X_B"; n1[k]+= $2; n2[k] += $3;  next;}
/^Anti_HELdg_A/ { k="Anti_X_pA_A"; n1[k]+= $2; n2[k] += $3; k="Anti_X_A"; n1[k]+= $2; n2[k] += $3;  next;}
/^Anti_HELdg_B/ { k="Anti_X_pA_B"; n1[k]+= $2; n2[k] += $3; k="Anti_X_B"; n1[k]+= $2; n2[k] += $3;  next;}
/^Anti_HEL_A/ { k="Anti_X_pA_A"; n1[k]+= $2; n2[k] += $3; k="Anti_X_A"; n1[k]+= $2; n2[k] += $3;  next;}
/^Anti_HEL_B/ { k="Anti_X_pA_B"; n1[k]+= $2; n2[k] += $3; k="Anti_X_B"; n1[k]+= $2; n2[k] += $3;  next;}
/^Anti_R454_A/ { k="Anti_X_ns_A"; n1[k]+= $2; n2[k] += $3;  next;}
/^Anti_R454_B/ { k="Anti_X_ns_B"; n1[k]+= $2; n2[k] += $3;  next;}
/^Anti_S_cap_B/ { k="Anti_X_B"; n1[k]+= $2; n2[k] += $3;  next;}


