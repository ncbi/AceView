






static int rank ;

#define MX(_m,_i,_j) _m[rank*(_i)+(_j)]


/* A(m/n) Cartan matrix */
static int *cartanA(int m, int n) 
{
  int i, j ;
  int *mm = malloc(rank * rank * sizeof(int)) ;

  memset (mm, 0, rank * rank * sizeof(int)) ;
  for (i = 0 ; i < m ; i++)
    { 
      MX(mm,i,i) = 2 ;
      if (i > 0) MX(mm, i-1, i) = -1 ;
      MX(mm,i,i+1) = -1 ;
    }
  for (i = m+1 ; i < m+n ; i++)
    { 
      MX(mm,i,i) = 2 ;
      if (i > 0) MX(mm, i-1, i) = -1 ;
      if (i < m+n-1)
	MX(mm,i,i+1) = -1 ;
    }
  MX(mm, i-1, m-1) = 1 ;
  MX(mm, i+1, m-1) = -1 ;
  
  return mm ;
}

int main ()
{
  int *cartan ;
  int m = 2 ;
  int n= 1 ;

  rank = m + n ;
  cartan = cartanA (m,n) ;
  matrixShow (cartan) ;
}
