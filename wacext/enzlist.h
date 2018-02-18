/* enzlist.h */

typedef struct enz_struct {
  char *name;
  char *seq;
  char *dna ;
  int length;
  int overhang;
}ENZ;

Array enzInit (void)  ;
