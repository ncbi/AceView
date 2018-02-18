function show(aa) { for (i=5;i<=na && i < 217 ; i+=5){x=0;for(j=0;j<5;j++)x+=aa[i+j]; printf("\t%.2f",x) ;}}
function aamax(aa) { m=0;for (i=25;i<=na ; i++)if(aa[i] > m) m=aa[i]; return m;}
function aasum(aa) { s=0;for (i=5;i<=na ; i++) s += aa[i]; return s;}
function tame(aa,x) { z = 0 ;  for (i = 5 ; i <= 5+x ; i++){z += aa[i]; aa[i] = 0 ;} aa[x] = z ; }

function  dZZ1(_dx,_dy,_ratio_bound) {
  _z = _ratio_bound ; 
  if (1)
    { # methode additive, count area only if above 2 fold */
      _z = 1 ;
      if (_dy > 0)
	{ 
          _m = 2.5 ;
	  _z = _dx/_dy ;
	  if (_z < _m)
	    _z = 2.2 * _z - (2.2 * _m - 1) ;
	  if (_z < 0)
	    _z = 0 ;
	  if (_z > 1)
	    _z = 1 ; 
	}
    }
  else
    {    /* methode multiplicative */
      if (_dy > 0) 
	_z = _dx/_dy ;
      if (_z < 3)
	_z = 2 * _z - 3 ; 
      if (_z > 0) 
	{
	  if (_z > _ratio_bound) 
	_z = _ratio_bound ;
	} 
      if (_z < 0)
	_z = 0 ;
    }

  _dzu = (_dx - _dy) * _z ; 
  return _dzu ;
}
function  dZZ(dx,dy,aa,bb,ratio_bound,ii1,ii2) {
  ok = 1 ; test=0 ;
  dz = dZZ1(dx,dy,ratio_bound);
  if (test>1) printf ("###dZZ \t%s\tii1=%d\tii2=%d\tdx=%.2f\tdy=%.2f\tdz=%.1f\tz=%.1f\n",gene,ii1,ii2,dx,dy,dz,dz/(dx-dy)) ;
  if (0) return dz ;
  while (ok == 1 && ii1 < ii2)
    {
      ok = 0 ;
      dx1 = dx - aa[ii1] ;
      dy1 = dy - bb[ii1] ;
      dz1 = dZZ1(dx1,dy1,ratio_bound); 
      if (test>1) print gene,".........OK#",ok,ii1,dz,dz1 ;
      if (dz1 > dz)
	{
	  ok = 1 ; dz = dz1 ; dx1 = dx ; dy1 = dy ; ii1 = ii1+1 ; if (test>1) print gene,"OK=",ok,ii1 ;
	}
      if (test>1) printf ("###dZZ1 \t%s\tok=%d\tii1=%d\tii2=%d\tdx=%.2f\tdy=%.2f\tdz1=%.1f\tz=%.1f\n",gene,ok,ii1,ii2,dx1,dy1,dz1,dz1/(dx1-dy1)) ;
      dx1 = dx - aa[ii2] ;
      dy1 = dy - bb[ii2] ;
      dz1 = dZZ1(dx1,dy1,ratio_bound);
      if (test>1) print gene,".........OK#",ok,ii2,dz,dz1 ;
      if (dz1 > dz)
	{
	  ok = 1 ; dz = dz1 ; dx1 = dx ; dy1 = dy ; ii2=ii2-1 ; if (test>1) print gene,"OK=",ok,ii2 ;
	}
      if (test>1) printf ("###dZZ2\t%s\tok=%d\tii1=%d\tii2=%d\tdx=%.2f\tdy=%.2f\tdz1=%.1f\tz=%.1f\n",gene,ok,ii1,ii2,dx1,dy1,dz1,dz1/(dx1-dy1)) ;
    }
  if (test>1) if(x+y>0)printf ("###dZZ3\t%s\tii1=%d\tii2=%d\tdx=%.2f\tdy=%.2f\tdz=%.1f\tz=%.1f\n",gene,ii1,ii2,dx,dy,dz,dz/(dx-dy)) ;
  return dz ;
}

{ gene = $1 ; if (gene != old) { na = 0 ; nb = 0 ; nnnnn = 0 ;} old = gene ;  nnnnn++ ; }
{ for(i=1;i<NF;i++) {if(nnnnn == 1) { na = NF ; a[i] = $i ; } if(nnnnn==2) { nb = NF ; b[i] = $i ; } if(nnnnn==3) { nc = NF ; c[i] = $i ;}}}
{
  if (nnnnn < 3) next ;
   if (na != nb || nb < 1) printf("\tERROR na=%d nb=%d", na, nb) ;
  z = 0 ; # score

  # method 1: overlap
  s1 = 0 ; s2 = 0 ;
  for (i = 5 ; i <= na ; i++)
    { 
      x = b[i] ; y = c[i] ; 
      # print gene, "nn=", nn, i, x, "y=",y ;
      if (x > y) { s1 += x - y ; s2 += x ; }
      else { s1 += y - x ; s2 += y ; }
    } 
  idx[1] = 1 - s1/s2 ;

  # method 2

  #  x = 50 ;  tame(a,x) ; tame(b,x) ; tame(c,x) ;

  mb = aamax(b) ; mc = aamax(c) ;
  s1 = 0 ; s2 = 0 ; i1 = 0 ;
  for (i = 5 ; i <= na ; i++)
    { 
      x = b[i] ; y = c[i] ; 
      if (x + y > 0.5 || x + y > (mb+mc)/10.0) { if (i1 == 0) i1 = i ; i2 = i ;} /* [i1,i2] extent */
      if (x > y) { s1 += x - y ; s2 += x ; }
      else { s1 += y - x ; s2 += y ; }
    } 
  idx[2] = s1 + i2 - i1 ;

  # method 5

  #  x = 50 ;  tame(a,x) ; tame(b,x) ; tame(c,x) ;

  sb = aasum(b) ; sc = aasum(c) ;
  mb = aamax(b) * 100/sb ; mc = aamax(c) * 100/sc ;
  # sb = 100 ; sc = 100 ;
  s1 = 0 ; s2 = 0 ; i1 = 0 ; nx = 0 ; old = 0 ;
  for (i = 5 ; i <= na ; i++)
    { 
      x = b[i] ; y = c[i] ; 
      x = x * 100/sb ; y = y * 100/sc ;
      mmx = mb  ; mmy = mc ; 
      if (x > y && x > mmx/10) { if(old == -1) nx++ ; old = 1 ; }
      if (y > x && y > mmy/10) { if(old == 1) nx++ ; old = -1 ; }
      if (x > y) { mm = mmx ; mmx = mmy ; mmy = mm ; z = x ; x = y ; y = z }
      # if (((i == 5 && y > 10 * x) || (i > 5 && y > 2*x)) && y > 0) { s1 += ( y - x) * (1 + y/(x+3*y)) ; }
      if (i < 60)
	{
	  if (y > 10 * x && y > 0) { s1 += 3 * ( y - x) * (1 + y/(x+y)) ; }
	}
      else
	{
	  if (y > x && y > 0) { s1 += ( y - x) * (1 + y/(x+y)) ; }
	  if (i > 60 && i < 80 && y > 2*x && y > mmy/10 && x < mmx/1) { s2 += ( 1 ) * (1 + y/(x+y)) ; }
	  if (i >= 80 && i < 120 && y > 2*x && y > mmy/10 && x < mmx/1) { s2 += ( 2 ) * (1 + y/(x+y)) ; }
	  if (i >= 120  && y > 2*x && y > mmy/10 && x < mmx/1) { s2 += ( 2.0 ) * (1 + y/(x+y)) ; }
	}
      # print i,y,x,s1,s2 ;
    } 
  # s1 = s1 * 2 ;
  s2 = 2 * s2 ;
  idx[5] = s1 + s2 - 30 * (nx - 1) ;
  idx[3] = s1 ;
  idx[4] = s2 ;
  if (1)printf ("%s %s sb=%d sc=%d diff=%.1f step=%.1f cross=%d score=%.1f mx=%.1f my=%.1f\n", gene, NA, sb, sc, s1,s2,nx,idx[5],mb,mc);

# method 6 : from-right to left diff times ratio times step , down to 9 

  threshold = 7 ;
  ratio_bound  = 4 ; # back to 6 # mb bound from 6 top  3  added if x/y<1.5 do not count, all this favors the large guys
  test=0 ;
  # initialize the loop
  sx = 0 ; sy = 0 ; sk = 0 ; z = 0 ; 
  iMin = 5 + 10 * threshold  ;   # col 5 -> point zero, then we step by 1/10, hence we want 5 + 10 * desired_threshold
  nx = 0 ; old = 0 ;
  dx = 0 ; dy = 0 ; dz = 0 ;
  dnix = 0 ; dniy = 0 ; dnx = 0 ; dny = 0 ;
  x = 0 ; y = 0 ; xmi = 0 ; ymi = 0 ;
  ii1 = 0 ; ii2 = 0 ;

  sb = aasum(b) ; sc = aasum(c) ;
  mx = aamax(b) ; my = aamax(c) ;
  
  mmx = mb  ; mmy = mc ; 
  # mmx = mx * 100/sb ; mmy = my * 100/sc  ;
  for (i = 5 ; i <= na+1 ; i++)
    { 
      if (i <= iMin)
	{
	  x += b[i] ; y += c[i] ; 
	  if (i < iMin) continue ;
	}
      if (i > iMin)
	{ x = b[i] ; y = c[i] ; }

      # normalize :  x = x * 100/sb ; y = y * 100/sc ;

      if (i == iMin)
	{
if (test>0) if(x+y>0)printf ("###old=0\t%s\ti=%d\tx=%.3f\ty=%.3f\tdx=%.3f\tdy=%.3f\tz=%.2f\tsx=%.2f\tsy=%.2f\tnx=%d\n",gene,i,x,y,dx,dy,z,sx,sy,nx);
	  if (x > y) old = 1 ;
	  else old = -1 ;
	}
      # for the evaluation of the first and last zone
      if (i == na + 1)   
	{
	  if (old == -1) x = mmx + y ;
	  else y = mmy + x ;
	}

      # compute the sum weigthed by the bounded ratio
      if (y > x && y > mmy/10)
	{ 
	  if (old == 1) 
	    { 
	      dz = dZZ(dx,dy,b,c,ratio_bound,ii1,ii2);
	      if (dz > 3)
		{
		  nx++ ;
		  sx += dz ; 
		  xmi = dnix / dnx ;
		}
	      if (test>1) if(x+y>0)printf ("###old=1\t%s\ti=%d\tx=%.3f\ty=%.3f\tdx=%.3f\tdy=%.3f\tz=%.2f\tsx=%.2f\tsy=%.2f\tnx=%d\n",gene,i,x,y,dx,dy,z,sx,sy,nx);
	      dx = 0 ; dy = 0 ; dnix = 0 ; dniy = 0 ; dnx = 0 ; dny = 0 ; ii1 = i ; ii2 = i ;
	    }
	  old = -1 ;
	}

      if (x > y && x > mmx/10)
	{ 
	  if (old == -1) 
	    { 
	      dz = dZZ(dy,dx,c,b,ratio_bound,ii1,ii2);
	      if (dz > 3)
		{
		  nx++ ;
		  ymi = dniy / dny ;
		  sy += dz ; 
		}
	      if (test>1) if(x+y>0)printf ("###old=-1\t%s\ti=%d\tx=%.3f\ty=%.3f\tdx=%.3f\tdy=%.3f\tz=%.2f\tsx=%.2f\tsy=%.2f\tnx=%d\n",gene,i,x,y,dx,dy,z,sx,sy,nx);
	      dx = 0 ; dy = 0 ; dnix = 0 ; dniy = 0 ; dnx = 0 ; dny = 0 ; ii1 = i ; ii2 = i ;
	    }
	  old = 1 ;
	}
      if (test>1) if(x+y>0)printf ("###old=-1\t%s\ti=%d\tx=%.3f\ty=%.3f\tdx=%.3f\tdy=%.3f\tz=%.2f\tsx=%.2f\tsy=%.2f\tnx=%d\tmx=%.2f\tmy=%.2f\n",gene,i,x,y,dx,dy,z,sx,sy,nx,mmx,mmy);

      # if (x > 1.5 * y || y > 1.5 * x)  stupid way to define a zone
      dx += x ; dy += y ; dnix += i * x ; dniy += i * y ; dnx += x ; dny += y ; if (ii1 == 0) ii1 = i ; ii2 = i ;
      z = 0 ;
    }

  dz = 0 ;
  if (nx == 2 && ymi < xmi + 2 && ymi > xmi - 2) dz += 2 ; 
  if (nx == 3) dz += 3 ;
  dz = 10 * dz * (sb + sc)/100 ;
  if (nx > 3) dz = sx + sy ;

  if (test>0) if(x+y>0)printf ("###COUNT\t%s\tnx=%d dz=%.2f sx=%.2f sy=%.2f sb=%.1f sc=%.1f\n",gene,nx,dz,sx,sy,sb,sc) ;

  idx[5] = sx + sy - dz ;
  if (uf == "f") { idx[3] = sx ;  idx[4] = sy ;}
  else { idx[3] = sy ;  idx[4] = sx ;}


#  print gene, i1, i2, i2 -i1, idx[4] ;

  if (0)
    {
# method 3 sumA
# method 4 sumF
# method 5 sumU
      x = 0 ; y = 0 ; z = 0 ;
      for (i = 5 ; i <= na ; i++)
	{ 
	  x += a[i] ; y += b[i] ; z += c[i] ;
	}
      idx[3] = x ;
      idx[4] = y ;
      idx[5] = z ;
    }
     
#############################################
#############################################
# method 7 : try o define the zones

  threshold = 7 ;
  ratio_bound  = 4 ; # back to 6 # mb bound from 6 top  3  added if x/y<1.5 do not count, all this favors the large guys
  jjShift=3 ;

  test=1 ;
  if (MM > 0) test = 0 ;
  # initialize the loop
  iMin = 5 + 10 * threshold  ;   # col 5 -> point zero, then we step by 1/10, hence we want 5 + 10 * desired_threshold


  # compute the total surface and the max of the continuous zone (above index 2) 
  sx = aasum(b) ; sy = aasum(c) ;
  mx = aamax(b) ; my = aamax(c) ;

  printf ("=== %s mx=%.1f my=%.1f sx =%.1f sy=%.1f\n", gene, mx, my, sx, sy) ;
  for (jjj = 0 ; jjj < 2 ; jjj++)
    {
      sens[jjj] = 0 ;
      scores[jjj] = 0 ;
      scoreA[jjj] = 0 ;
      scoreB[jjj] = 0 ;
      boni[jjj] = 0 ;
      iiix1[jjj] = 0 ; iiix2[jjj] = 0 ; 
      iiiy1[jjj] = 0 ; iiiy2[jjj] = 0 ; 
    }

   for (jjj = 0 ; jjj < 2 ; jjj++)
    {
      peakx = 0; peaky = 0; dp = 0; bonus = 0 ;
      
      aBest1 = 0 ; aBest2 = 0 ; bBest1 = 0 ; bBest2 = 0 ;
      a1 = 0 ; a2 = 0 ; b1 = 0 ; b2 = 0 ;
      dx = 0 ; dy = 0 ;
      sx1 = sy1 = 0 ;
      ijjxa = 0 ; ijjxb = 0 ; ijjya = 0 ; ijjyb = 0 ;
      ixa = 0 ; ixn = 1 ; iya = 0 ; iyn = 1 ; 
      jxa = 0 ; jxn = 1 ; jya = 0 ; jyn = 1 ;
      ix1 = 0 ; ix2 = 0 ; iy1 = 0 ; iy2 = 0 ;
      jx1 = 0 ; jx2 = 0 ; jy1 = 0 ; jy2 = 0 ;
      
# look for the best zone from 0 to i
# collect at the same time the position of the dominant peaks
      for (i = 5 - jjShift ; i <= na + jjShift - 1 ; i = i+1)
	{ 
	  if (0)
	    { x = b[i] ; y = b[i] ; }
	  else
	    {
	      ijjx = i + jjShift * (1 - jjj) ;
	      ijjy = i + jjShift * jjj ;
	      
	      if (ijjx < 5 || ijjx > na) x = 0 ;
	      else x = b[ijjx] ;

	      if (ijjy < 5 || ijjy > na) y = 0 ;
	      else y = c[ijjy] ;
	    }
	  if (0)
	    {
	      /* normalize : */
	      x = x * 100.0 / sx ; y = y * 100.0 / sy ; 
	    }
	  dx += x ; dy += y ;
	  
	  if (i < iMin) continue ;
	  if (i == iMin) {  x = dx ; y = dy ; }
	  
	  # average of significant region 
	  if (x > 1.5 * y && x > mx/3 ) { ixa += ijjx * dx ; ixn += dx ; sx1 += x ; if(! ix1) ix1 = ijjx ; ix2 = ijjx ; }
	  if (y > 1.5 * x && y > my/3) 	{ iya += ijjy * dy ; iyn += dy ; sy1 += y ; if(! iy1) iy1 = ijjy ; iy2 = ijjy ; }
	  
	  # true average
	  if (x > 0 && x > mx/3 ) 	{ jxa += ijjx * dx ; jxn += dx ; if(! ijjxa) ijjxa = ijjx ; ijjxb = ijjx ; }
	  if (y > 0 && y > my/3) 	{ jya += ijjy * dy ; jyn += dy ; if(! ijjya) ijjya = ijjy ; ijjyb = ijjy ; }
	  
	  dza = dZZ1(dx,dy,ratio_bound) ;
	  if (dza > aBest1) { aBest1 = dza ; a2 = ijjx ; }
	  
	  dzb = dZZ1(dy,dx,ratio_bound) ;
	  if (dzb > bBest1) { bBest1 = dzb ; b2 = ijjy ; }

	  if (test > 1 && x+y>0)
	    printf ("+++ %s i=%d  x=%.2f y=%.2f dx=%.1f dy=%.1f dza=%.1f dzb=%.1f aBest1=%.1f aBest2=%.1f bBest1=%.1f bBest2=%.1f\n", gene, i, x,y,dx,dy,dza,dzb,aBest1, aBest2, bBest1,bBest2) ;
	}
      

# look for the best zone from na  to i
      dx = 0 ; dy = 0 ;
      for (i = na + jjShift - 1 ; i >= 5 - jjShift ; i = i - 1)
	{ 
	  if (0)
	    { x = b[i] ; y = b[i] ; }
	  else
	    {
	      ijjx = i + jjShift * (1 - jjj) ;
	      ijjy = i + jjShift * jjj ;
	      
	      if (ijjx < 5 || ijjx > na) x = 0 ;
	      else x = b[ijjx] ;

	      if (ijjy < 5 || ijjy > na) y = 0 ;
	      else y = c[ijjy] ;
	    }

	  dx += x ; dy += y ;
	  if (i < iMin) continue ;
	  
	  if (x > 1.5 * y && x > mx/3 ) { if(! jx2) jx2 = ijjx ; jx1 = ijjx ; }
	  if (y > 1.5 * x && y > my/3) 	{ if(! jy2) jy2 = ijjy ; jy1 = ijjy ; }
	  
	  dza = dZZ1(dx,dy,ratio_bound) ;
	  if (dza > aBest2) { aBest2 = dza ; a2 = ijjx ; }
	  
	  dzb = dZZ1(dy,dx,ratio_bound) ;
	  if (dzb > bBest2) { bBest2 = dzb ; b2 = ijjy ; }
	  
	  if (test > 1 && x+y>0)
	    printf ("--- %s i=%d  x=%.2f y=%.2f dx=%.1f dy=%.1f dza=%.1f dzb=%.1f aBest1=%.1f aBest2=%.1f bBest1=%.1f bBest2=%.1f\n", gene, i, x,y,dx,dy,dza,dzb,aBest1, aBest2, bBest1,bBest2) ;
	}
      
      dx = 0 ; dy = 0 ; 
      if (aBest1 > aBest2)  { dx = ix2 - ix1 ;} # D8  measure the zone width and add a bonus */
      if (aBest1 < aBest2)  { dx = jx2 - jx1 ;}
      if (bBest1 > bBest2)  { dy = iy2 - iy1 ;}
      if (bBest1 < bBest2)  { dy = jy2 - jy1 ;}
      

      # D8  if the interesting zone < 1 on each side, disfavor, else favor, 10points per index step */

      score = 0 ;
      if (aBest1 >= aBest2 && bBest2 >= bBest1) 
	{
	  sens[jjj] = 1 ;
	  scores[jjj] = aBest1 + bBest2 ;
	  scoreA[jjj] = aBest1 ; 	  
	  scoreB[jjj] = bBest2 ; 
	  iiix1[jjj] = ix1 ; iiix2[jjj] = ix2 ; 
	  iiiy1[jjj] = jy1 ; iiiy2[jjj] = jy2 ; 
	}
      if (aBest1 <= aBest2 && bBest2 <= bBest1)
	{
	  sens[jjj] = 2 ;
	  scores[jjj] = aBest2 + bBest1 + bonus ;
	  scoreA[jjj] = aBest2 ; 	  
	  scoreB[jjj] = bBest1 ; 
	  iiix1[jjj] = jx1 ; iiix2[jjj] = jx2 ; 
	  iiiy1[jjj] = iy1 ; iiiy2[jjj] = iy2 ; 
	}

      bonus = 0 ;

      if (1)
	{
	  peakx = ixa/ixn ;  peaky = iya/iyn ; dp = peaky - peakx ; 
	  if (dp < 0) dp = -dp ;
	  
	  if (0)  # danielle
	    {
	      if (100 * scoreA[jjj] > 15 * sx && 100 * scoreB[jjj] > 10 * sy)
		bonus += .5 * (dp - 20)/10 * (sx + sy)/85 ;
	      if (100 * scoreB[jjj] > 15 * sy && 100 * scoreA[jjj] > 10 * sx)
		bonus += .5 * ((dp - 20)/10) * (sx + sy)/85 ;
	    }

	  if (0 && 10 * sx1 > sx && 10 * sy1 > sy)
	    {
# peak distance in significant differential positions
	      peakx = ixa/ixn ;  peaky = iya/iyn ; dp = peaky - peakx ; if (dp < 0) dp = -dp ;
	      if (dp > 40) dp = 40 ;
	      bonus += 1 * 0.5 * 0.5 * dp * (sx + sy)/100 ;
	      
# peak distance in significant not necesseraly differential positions
	      peakx = jxa/jxn ;  peaky = jya/jyn ; dp = peaky - peakx ; if (dp < 0) dp = -dp ;
	      if (dp > 40) dp = 40 ;
	      dp -= 10 ;  
              bonus +=  7 * (dp/40)  * (sx + sy)/100 ; 
	      # printf("# %s Bonus= %.1f delta=%.1f = %.1f | %.1f\n", gene,  bonus, dp, peakx, peaky) ;
	    }

	  if (bonus < 0) bonus = 0 ;
	  if (0 && bonus >  0) printf ("##dp\t%s\t%d\n", gene, dp) ;
	}


      boni[jjj] = bonus ;
      scores[jjj] += bonus ;

      if (test > 0)
        printf ("==+ %s\tscore=%.1f\tsens %d\t\tA %.1f::%d:%d//%.1f::%d:%d//%.1f::%d:%d\tB %.1f::%d:%d//%.1f::%d:%d//%.1f::%d:%d\n", gene,scores[jjj], jjj, scoreA[jjj], ijjxa,ijjxb,aBest1,ix1,ix2,aBest2,jx1,jx2,scoreB[jjj],ijjya,ijjyb,bBest1,iy1,iy2,bBest2,jy1,jy2) ;
    }
   if (scores[0] < scores[1]) # take the least favorable score */
     jjj = 0 ;
   else
     jjj = 1 ;

   if (sens[0] == sens[1])
     {
       score = scores[jjj] ;
       bonus = boni[jjj] ;
       aBest1 = scoreA[jjj] ;
       bBest1 = scoreB[jjj] ;
       ix1 = iiix1[jjj] ; 
       ix2 = iiix2[jjj] ;   
       iy1 = iiiy1[jjj] ;      
       iy2 = iiiy2[jjj] ;
     }
   else
     {
       jjj = -1 ; 
       score = 0 ;
       bonus = 0 ;
       aBest1 = 0 ;
       bBest1 = 0 ;
       ix1 = 0 ;
       ix2 = 0 ;
       iy1 = 0 ;
       iy2 = 0 ;
     }

   if (test > 0)
     printf ("=== %s\tscore=%.1f\tsens %d\tbonus=%.1f\t\tA %.1f::%d:%d\tB %.1f::%d:%d\n", gene,score, jjj, bonus, aBest1, ix1, ix2, bBest1, iy1, iy2) ;
   
   if (uf == "f") jjj = 0 ;
   else jjj = 1 ;

   idx[3+jjj] = aBest1 "::" ix1 ":" ix2 ;
   idx[4-jjj] = bBest1 "::" iy1 ":" iy2 ;

   idx[5] = score ;

#############################################

   if (1)
     {
       printf ("MMMMMM\t%s\t%.2f\n", gene, idx[5]) ;
     }
   if (1)
     {
       if (firstPass < 1)
	 {
	   firstPass = 1 ; 
	   printf ("000000000\tGene\tType\tGene\tOverlap fraction\tExtent\tS1\tS2\tscore\tbonus\tdelta\tMix\tGene", ratio_bound) ; for(i=0;i<=na && i<217;i+=5) printf ("\t%d", i) ;
	 }
       a[4] = a[3]+a[4] ;
       printf ("\n%09d\t%s %s\tA",100000000-100*idx[5],gene,NA) ; printf ("\t%s %d/%d/%d", gene,idx[3],idx[4],bonus) ; for(i=1;i<=5;i++) printf ("\t%.2f", idx[i]) ; printf ("\t%.2f", idx[5] - idx[3] - idx[4]) ;for(i=2;i<=3;i++) printf ("\t%.2f", a[i]) ; printf ("\t%s %d/%d", gene,idx[3],idx[4]) ;show(a) ;
       printf ("\n%09d\t%s %s\tF",100000000-100*idx[5],gene,NA) ; printf ("\t%s %d/%d/%d", gene,idx[3],idx[4],bonus) ; for(i=1;i<=5;i++) printf ("\t%.2f", idx[i]) ; printf ("\t%.2f", idx[5] - idx[3] - idx[4]) ;for(i=2;i<=3;i++) printf ("\t%.2f", a[i]) ;  printf ("\t%s %d/%d", gene,idx[3],idx[4]) ;if (uf == "f") show(b) ; else show(c) ;
       printf ("\n%09d\t%s %s\tU",100000000-100*idx[5],gene,NA) ; printf ("\t%s %d/%d/%d", gene,idx[3],idx[4],bonus) ; for(i=1;i<=5;i++) printf ("\t%.2f", idx[i]) ; printf ("\t%.2f", idx[5] - idx[3] - idx[4]) ; for(i=2;i<=3;i++) printf ("\t%.2f", a[i]) ; printf ("\t%s %d/%d", gene,idx[3],idx[4]) ;if (uf == "f") show(c) ; else show(b) ;
       printf ("\n") ;
     }
}

      
