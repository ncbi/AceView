/*  File: gmapconvert.c
 *  Author: Simon Kelley (srk@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1991
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: assorted utility functions for extracting data from
 * Map objects  
 * Exported functions:
 *     str2tag
 *     gMapNeighbours
 *     gMapPositive
 *     gMapNegative
 *     gMapMultiPt
 *     gMap2Pt
 *     gMapGetMultiData
 *     gMapGet2PtData
 *     getPos
 *     setTestPos
 *     gMapGetMapObject
 *     gMapSymbols
 *     gMapKeyOnTag
 *     
 * HISTORY:
 * Last edited: Dec 21 11:32 1998 (fw)
 *-------------------------------------------------------------------
 */

/* $Id: gmapconvert.c,v 1.3 2007/03/21 21:02:52 mieg Exp $ */

#include "gmap.h"
#include "pick.h"
#include "freeout.h"


KEY gMapCacheKey(KEY map, char *suffix)
{ char buf[100];
  KEY ret;
  OBJ obj = bsCreate(map);
  BOOL noCache = FALSE;

  if (obj)
    { noCache = bsFindTag(obj, str2tag("No_cache"));
      bsDestroy(obj);
    }
  if (noCache)
    return 0;
  strcpy(buf, name(map));
  strcat(buf, ".");
  strcat(buf, suffix);
  if (lexword2key(buf, &ret, _VgMap))
    return ret;
  else
    return 0;
}

void gMapCacheStore(KEY map, char *suffix, Array a, char *format)
{ char buf[100];
  KEY key;
  strcpy(buf, name(map));
  strcat(buf, ".");
  strcat(buf, suffix);
  lexaddkey(buf, &key, _VgMap);
  arrayStore(key, a, format);
}

void gMapCacheKill(KEY map)
{ char buf[100];
  KEYSET kset;
  int i;

  strcpy(buf, "FIND gMap ");
  strcat(buf, name(map));
  strcat(buf, ".*");
  kset = query(0, buf);
  for (i = 0 ; i < keySetMax(kset); i++)
    arrayKill(keySet(kset, i));
}


BOOL gMapNeighbours(MAPCONTROL map, KEYSET *keyset, KEY key)
{ KEYSET a;
  int i, j;

  *keyset = keySetReCreate(*keyset);
 
  a = bsKeySet(key);
  if (!a) return FALSE;

  for(i=0, j=0 ; i<keySetMax(a) ; i++)
    if (class(keySet(a,i)))
      keySet(*keyset,j++) = keySet(a,i) ;
  keySetDestroy(a) ;

  return keySetMax(*keyset) != 0;
  

}

static BOOL gMapContainment(int level, MAPCONTROL map, KEYSET *ks, 
			    KEY key, KEY senseTag)
/* senseTag is _Positive or _Negative */
{ KEY mapk = map->key;
  OBJ obj, comp;
  KEY k, mo1, mo2;
  int typeOf1=0, typeOf2=0, i;
  Array loci = arrayCreate(50, BSunit);

	/* get info in object */
  if ((obj = bsCreate(key)))
    { if (bsGetKey(obj, str2tag("More_data"), &k))
	do 
	  { if (level>100)
	      messerror("Have chased 100 levels of More_data and got to %s, "
			"is there a loop?", name(k));
	    else
	      gMapContainment(level+1, map, ks, k, senseTag);
	  } 
	while (bsGetKey(obj, _bsDown, &k));
  
    }
      
  if (bsFindTag(obj, senseTag) && bsFlatten(obj, 2, loci))
    for (i = 1; i<arrayMax(loci); i+= 2)
      { k = arr(loci, i, BSunit).k;
	if (k)
	  keySet(*ks, keySetMax(*ks)) = k;
      }

	/* get info in Pos_neg_data */
  if (bsGetKey(obj, str2tag("Pos_neg_data"), &k))
    do 
      {
	comp = bsCreate(k);
	if (!bsFindTag(comp, senseTag))
	  { bsDestroy(comp);
	    continue;
	  }
	
	if (bsFindTag(comp, str2tag("Item_1")) &&
	    bsGetKeyTags(comp, _bsRight, &mo1) &&
	    bsGetKey(comp, _bsRight, &mo1))
	  typeOf1 = gMapGetMapObject(mo1, mapk, 0, 20, 0, 0, 0, 0);
	if (bsFindTag(comp, str2tag("Item_2")) &&
	    bsGetKeyTags(comp, _bsRight, &mo2) &&
	    bsGetKey(comp, _bsRight, &mo2))
	  typeOf2 = gMapGetMapObject(mo2, mapk, 0, 20, 0, 0, 0, 0);
	bsDestroy(comp);
	
	if (((typeOf1&FLAG_ANY_INTERVAL) && (typeOf2&FLAG_ANY_LOCUS)) ||
	    ((typeOf2&FLAG_ANY_INTERVAL) && (typeOf1&FLAG_ANY_LOCUS)))
	  {
	    if (mo1 == key)
	      keySet(*ks, keySetMax(*ks)) = mo2;
	    
	    if (mo2 == key)
	      keySet(*ks, keySetMax(*ks)) = mo1;
	  }
      }
    while (bsGetKey(obj, _bsDown, &k));
  
  
  keySetSort(*ks);
  bsDestroy(obj);
  arrayDestroy(loci);
  return keySetMax(*ks) != 0;
}

BOOL gMapPositive(MAPCONTROL map, KEYSET *ks, KEY key)
{ *ks = keySetReCreate(*ks);
  return gMapContainment(0, map, ks, key, str2tag("Positive"));
}

BOOL gMapNegative(MAPCONTROL map, KEYSET *ks, KEY key)
{ *ks = keySetReCreate(*ks);
  return gMapContainment(0, map, ks, key, str2tag("Negative"));
}

BOOL gMapMultiPt(MAPCONTROL map, KEYSET *ks, KEY key)
{
  OBJ obj ;
  KEY k ;
  
  *ks = keySetReCreate(*ks) ;
  if ((obj = bsCreate (key)))
    { if (bsGetKey (obj, _Multi_point, &k)) do
	keySet (*ks, keySetMax(*ks)) = k ;
      while (bsGetKey (obj, _bsDown, &k)) ;
      bsDestroy (obj) ;
      keySetSort(*ks);
    }

  return keySetMax(*ks) != 0;
}

BOOL gMap2Pt(MAPCONTROL map, KEYSET *ks, KEY key)
{
  OBJ obj ;
  KEY k ;
  
  *ks = keySetReCreate(*ks) ;
  if ((obj = bsCreate (key)))
    { if (bsGetKey (obj, _2_point, &k)) do
	keySet (*ks, keySetMax(*ks)) = k ;
      while (bsGetKey (obj, _bsDown, &k)) ;
      bsDestroy (obj) ;
      keySetSort(*ks);
    }

  return keySetMax(*ks) != 0;
}

static void addMultiData (OBJ obj, Array counts, Array loci)
{
  int i, count ;
  KEY locus, tag ;
  static BSMARK mark ;

  mark = bsMark (obj, mark) ;

  for (i = 0 ; ; ++i)
    { bsPushObj(obj);
      if (!bsGetKeyTags(obj, _bsRight, &tag) ||
	  !bsGetKey (obj, _bsRight, &locus))
	{ messerror ("No terminating locus in multi results") ;
	  goto done ;
	}
      if (i >= arrayMax(loci))
	array (loci, i, KEY) = locus ;
      else if (locus != array(loci,i,KEY))
	{ messerror ("Mismatching locus %s in multi results", 
		     name(locus)) ;
	  goto done ;
	}
      if (!bsGetData (obj, _bsRight, _Int, &count))
	goto done ;
      if (i >= arrayMax(counts))
	array(counts,i,int) = count ;
      else
	array(counts,i,int) += count ;
    }

 done:
  bsGoto (obj, mark) ;
}


BOOL gMapGetMultiPtData(MAPCONTROL map, 
			KEY multi, 
			AC_HANDLE handle, 
			MULTIPTDATA *ret)
{ Array counts, loci;
  float y, min, max;
  OBJ obj;
  MULTIPTDATA data;
  int i;
  
  if (!(obj = bsCreate (multi)))
    return FALSE;
  
  counts = arrayHandleCreate(20, int, handle);
  loci = arrayHandleCreate(20, KEY, handle);
  
  if (bsFindTag (obj, _Combined))
    addMultiData (obj, counts, loci) ;
  else 
    { if (bsFindTag (obj, _A_non_B))
	addMultiData (obj, counts, loci) ;
      if (bsFindTag (obj, _B_non_A))
	addMultiData (obj, counts, loci) ;
    }

  bsDestroy (obj);

  if (arrayMax(loci) <=2 && iskey(multi) > 1)  /* i.e. not just a name */
    {
      arrayDestroy(counts);
      arrayDestroy(loci);
      return FALSE;
    }
  
  min = 1000000 ; max = -1000000 ;
  for (i = 0 ; i < arrayMax(loci) ; ++i)
    if (getPos (map, arr(loci,i,KEY), &y))
      { if (y < min) min = y ;
	if (y > max) max = y ;
      }

  if (min > max)
    { arrayDestroy(counts);
      arrayDestroy(loci);
      return FALSE;
    }

  data = (MULTIPTDATA)handleAlloc(0, handle, sizeof(struct MultiPtData));
  data->counts = counts;
  data->loci = loci;
  data->min = min;
  data->max = max;
  *ret = data;
  return TRUE;
}


 
BOOL gMapGet2PtData(MAPCONTROL map, 
		    KEY key, 
		    AC_HANDLE handle, 
		    TWOPTDATA *ret) 
{ OBJ obj;
  KEY type;
  TWOPTDATA data;
  BOOL barf;
  int *p;
  int i;
  

  obj = bsCreate(key);
  if (!obj)
    return FALSE; 

  data = (TWOPTDATA)handleAlloc(0, handle, sizeof(struct TwoPtData));
  
  
  if (!bsGetKey(obj, _Locus_1, &data->loc1) || 
      !bsGetKey(obj, _Locus_2, &data->loc2) ||
      !getPos(map, data->loc1, &data->y1) ||
      !getPos(map, data->loc2, &data->y2))
    { bsDestroy(obj);
      messfree(data);
      return FALSE;
    }
  

  if (!bsGetKeyTags(obj, _Calculation, &type))
    { data->type = 0;
      if (!bsGetData(obj, _Distance, _Float, &data->distance))
	data->distance = 0;
      if (!bsGetData(obj, _Error, _Float, &data->error))
	data->error = 0;
    }
  else
    { p = &data->n1;
      for (i=0; bsGetData(obj, _bsRight, _Int, p) && i<4; p++, i++);
      data->count = i;
      data->type = type;
      
      barf = FALSE;

      if ((type == _Full ||		/* WT X Y XY */
	   type == _Backcross ||	/* WT X Y XY */
	   type == _Sex_full) &&        /* WT X Y XY */
	  i != 4)
	barf = TRUE;

      if ((type == _Recs_all ||		/* X Y ALL */
	   type == _Tetrad) &&  	/* PD NPD TT */
	  i != 3)
	barf = TRUE;

      if ((type == _One_recombinant ||	/* WT X */
	   type == _Selected ||		/* X XY */
	   type == _One_all ||		/* X ALL */
	   type == _One_let ||		/* X ALL */
	   type == _Tested ||		/* H X */
	   type == _Selected_trans ||	/* X XY */
	   type == _Back_one ||		/* WT X */
	   type == _Sex_one ||		/* WT X */
	   type == _Sex_cis ||		/* X ALL */
	   type == _Dom_one ||		/* WT nonWT */
	   type == _Dom_selected ||	/* WT X */
	   type == _Dom_semi ||		/* XD ALL */
	   type == _Dom_let ||		/* WT ALL */
	   type == _Direct ||		/* R T */
	   type == _Complex_mixed ||	/* X ALL */
	   type == _Centromere_segregation) && /* 1st 2nd */
	  i != 2)
	barf = TRUE;

      if (barf)
	{ messerror("Wrong number of data (%d) for 2 point calculation %s",
		    i, name(type));
	  messfree(data);
	  bsDestroy(obj);
	  return FALSE;
	}

    }  
 
  bsDestroy(obj);
  *ret = data;
  return TRUE;
  
}

/* This is a utility function which returns the map position of any
   locus-like object on the map. As side effects, it is responsible
   for setting the values of map->min, map->max etc and look->orderedLoci,
   which is a keyset of all the loci which have a Positive_clone entry.
   It may be called with key zero to get its side-effects */

BOOL getPos (MAPCONTROL map, KEY key, float *y)
{ OBJ Chrom, Locus;
  KEY chrom, locus, clone;
  KEY posCache, orderedCache;
  int i, j, flags;
  float xmin = -1, xmax = 1; /* minimum sensible range */
  float x, dx, width;
  BOOL flipped;
  GeneticMap look = (GeneticMap)map->look;
  Array cacheArray = 0, orderedArray = 0, loci = 0;
  union 
    { float f;
      void *vs;
     } u;
  struct cacheEntry
    { KEY key;
      float x;
    }; 
  Associator a = look->posAss;
  

  if (!a)
    { look->posAss = a = assHandleCreate(map->handle);
      look->orderedLoci = arrayHandleCreate(200, KEY, map->handle);
      if ((posCache = gMapCacheKey(map->key, "pos")) && 
	  (cacheArray = arrayGet(posCache, struct cacheEntry, "kf")) &&
	  (orderedCache = gMapCacheKey(map->key, "ordered")) &&
	  (orderedArray = arrayGet(orderedCache, KEY, "k")))
	{ for (i= 0; i <arrayMax(orderedArray); i++)
	    keySet(look->orderedLoci, i) = arr(orderedArray, i, KEY);
	  for (i = 2; i <arrayMax(cacheArray); i++)
	    { u.f = arr(cacheArray, i, struct cacheEntry).x;
	      assInsert(a, assVoid(arr(cacheArray, i,struct cacheEntry).key), 
			u.vs);
	    }
	  xmin = arr(cacheArray, 0, struct cacheEntry).x;
	  xmax = arr(cacheArray, 1, struct cacheEntry).x;
	}
      else
	{ 
	  loci = gMapMapObjects(map->key);
	  cacheArray = arrayCreate(200, struct cacheEntry);
	  j = 2; /* zero and one for min and max */
	  for (Locus = 0, i = 0; i <arrayMax(loci); bsDestroy (Locus) , i++)
	    { locus = arr(loci, i, KEYPAIR).obj;
	      chrom = arr(loci, i, KEYPAIR).map;
	      flags = gMapGetMapObject(locus, map->key, chrom,
				       20, &x, &dx, &Locus, 0);
	      if (!flags)
		continue;
	      if (flags & FLAG_ANY_LOCUS)
		{ /* put info into cacheArray */
		  array(cacheArray, j, struct cacheEntry).x = x;
		  array(cacheArray, j, struct cacheEntry).key = locus;
		  j++;
		  /* This is very bizarre, we need need polymorphic assocs ASAP */
		  u.f = x;
		  assInsert(a, assVoid(locus), u.vs);
		  if (x < xmin)
		    xmin = x ;
		  if (x > xmax)
		    xmax = x ;
		}
	      else /* Interval */
		{ if (x-dx < xmin)
		    xmin = x-dx;
		  if (x+dx > xmax)
		    xmax = x+dx;
		}
	      
	      if (bsGetKey (Locus,_Positive_clone,&clone))
		keySetInsert(look->orderedLoci, locus);
	    }
	  
	  array(cacheArray, 0, struct cacheEntry).x = xmin;
	  array(cacheArray, 1, struct cacheEntry).x = xmax;
	  orderedArray = arrayCopy(look->orderedLoci);

	  gMapCacheStore(map->key, "pos", cacheArray, "kf");
	  gMapCacheStore(map->key, "ordered", orderedArray, "k");
	}
      
      arrayDestroy(loci);
      arrayDestroy(cacheArray);
      arrayDestroy(orderedArray);
      /* If we can't find display tags in this map, 
	 we look in any it inherits from */
      map->min = xmin;
      map->max = xmax;
      chrom = map->key;
      while ((Chrom = bsCreate(chrom)))
	{ if (bsGetData(Chrom, _Extent, _Float, &map->min))
	    { bsGetData(Chrom, _bsRight, _Float, &map->max);
	      break;
	    }
	  if (bsGetKey(Chrom, _From_map, &chrom))
	    { bsDestroy(Chrom);
	      continue;
	    }
	  break;
	}
      bsDestroy(Chrom);
      
      map->centre = (map->max + map->min)/2.0; /* default values */
      width = (map->max - map->min)/3.0;
      flipped = FALSE;
      chrom = map->key;
      while ((Chrom = bsCreate(chrom)))
	{ if (bsFindTag(Chrom, _Flipped))
	    flipped = !flipped;
	  if (bsGetData(Chrom, _Centre, _Float, &map->centre))
	    { bsGetData(Chrom, _bsRight, _Float, &width);
	      break;
	    }
	  if (bsGetKey(Chrom, _From_map, &chrom))
	    { bsDestroy(Chrom);
	      continue;
	    }
	  break;
	}
      bsDestroy(Chrom);
      map->mag = (map->control->graphHeight-5) / width;
      if (flipped)
	map->mag = - map->mag;
      
    }
  
  if (key && assFind(look->posAss, assVoid(key), &u.vs))
    { if (y) 
	*y = u.f;
      return TRUE;
    }
  
  return FALSE;
}


/* override the position value for a Locus temporarily */

void setTestPos (MAPCONTROL map, KEY key, float pos)
{
  GeneticMap look = (GeneticMap)map->look;
  union { float f;
	  void *vs;
	} u;
  
  assRemove(look->posAss, assVoid(key));

  u.f = pos;
  assInsert(look->posAss, assVoid(key), u.vs);

}

/* This is a fundamental routine. It takes the key of 
   some object which has a Map tag followed by ?Map #Map_position
   and the key of the map on which the postition is required.
   It returns zero if something's missing,
   or FLAG_ANY_INTERVAL if it's an interval or FLAG_ANY_LOCUS if it's
   a Locus. 
   If the object doesn't have a position on the given map, but the map
   has a From_map tag, it searches recursively in the inherited from maps.
   For an interval the centre goes into *xret  and half the length
   into *dxret. For a Locus the position goes into *xret and the 
   error into *dxret.
   dxret and/or xret may be zero without segfaulting.
   The map key is needed to resolve the positions the objects which
   appear on more than one map. If objp is non-zero, the bsObject 
   for the object is not closed, but returned (a bsGoto(0) is performed on it).
   ind is the limit on how far one should recurse to resolve maps_with
   before giving up. 
   If count is non-null, multiple positions may be returned. In this
   case *count on calling gives the maximum number of positions to
   be returned (must be >=1) and xret and dxret must be suitably sized arrays.
   On return, *count if the actual number of positions found.

   Note that the object need not be mapped directly on the map given as a 
   parameter, it may be mapped on another map which the parameter map inherits
   from, or another map which it includes. See gMapMapObjects for details.

   
*/

int gMapGetMapObject
(KEY key, KEY map, KEY hint, int ind, 
 float *xret, float *dxret, OBJ *objp, int *count)
{
  OBJ item;	
  int flags;
  float offset;
  float x1, x2;
  int max, ourCount;
  BSMARK mark = 0;
  KEY _Multi_Ends = str2tag("Multi_Ends");
  float mapOffset, mapFactor; 

  if (count)
    max = *count;
  else
    max = 1;

  if (ind == 0)
    { messerror("Chased too many mapable objects trying to "
		"resolve \"maps with\". Is there a cycle involving %s?", 
		name(key));
      if (count)
	*count = 0;
      return 0;
    }

  if (!key || !(item = bsCreate(key)))
    { if (count)
	*count = 0;
      return 0;
    }

/* now find out what's what in the mapping dept. If we have a hint, and it's
   the same as the map, that's easy. A hint != map, means the locus is on an 
   inhertits-from map, or a sub-map. we find which, and determine the 
   transformation. If there is no hint, we need to first find the actual map 
   that the locus is on, using the same priority rules as gMapMapObjects
   (IMPORTANT), then proceed as above. */
 
  if (hint)
    {
      if (map != hint)
	/* find the map which includes the hint, or is the hint */
	{
	  OBJ Map=0;
	  while (hint != map &&
		 (Map = bsCreate(map)) &&
		 !bsFindKey(Map, str2tag("Includes"), hint))
	    /* cannot find mapping for this map, it may inherit from another */
	    { if (!bsGetKey(Map, str2tag("From_map"), &map))
	      { bsDestroy(Map);
	      bsDestroy(item);
	      if (count)
		*count = 0;
	      return 0;
	      }
	    bsDestroy(Map);
	    }
	  bsDestroy(Map);
	}
      
      if (!bsFindKey(item, str2tag("Map"), hint))
	{ bsDestroy(item);
	  if (count)
	    *count = 0;
	  return 0;
	}
    }
  else /* !hint */
    while  (!hint) 
      { if (bsFindKey(item, str2tag("Map"), map)) 
	  /* it's on the map we want, easy */
	  hint = map;
      else
	/* it's not, see if it's on a sub-map */
	{ OBJ Map = bsCreate(map);
	  KEY submap;
	  if (bsGetKey(Map, str2tag("Includes"), &submap)) 
	    do 
	      { if (bsFindKey(item, str2tag("Map"), submap))
		  hint = submap;
	      }
	  while (!hint && bsGetKey(Map, _bsDown, &submap));
	  bsDestroy(Map);
	}

      if (!hint) /* not found it yet, try inherits from maps. */
	{ OBJ Map = bsCreate(map);
	  if (!Map || !bsGetKey(Map, str2tag("From_map"), &map))
	    /* No more, failed to find it. */
	    { bsDestroy(Map);
	      bsDestroy(item);
	      if (count)
		*count = 0;
	      return 0;
	    }
	  bsDestroy(Map);
	}
      }
  
  /* Ok , we now have a sequence map, hint, item, if hint != map, 
     we can calculate the scaling factor and offset */
  
  if (hint == map)
    { mapFactor = 1; 
      mapOffset = 0;
    }
  else
    { float leftEnd, rightEnd, extentLeft, extentRight;
      OBJ Submap = bsCreate(hint);
      if (!bsFindKey(Submap, str2tag("Map"), map) || 
	  !bsPushObj(Submap) ||
	  !bsGetData(Submap, _Left, _Float, &leftEnd) ||
	  !bsGetData(Submap, _Right, _Float, &rightEnd))
	messerror("Failed to find map position for submap %s in map %s",
		  name(hint), name(map));
      if (leftEnd == rightEnd)
	messerror("Submap %s has zero length on map %s - cannot scale.",
		  name(hint), name(map));
      bsGoto(Submap, 0);
      if (!bsGetData(Submap, str2tag("Extent"), _Float, &extentLeft) ||
	  !bsGetData(Submap, _bsRight, _Float, &extentRight))
	messerror("Map %s is a submap of map %s, but does not have an extent",
		  name(hint), name(map));
      
      bsDestroy(Submap);
      
      if (leftEnd > rightEnd)
	{ float tmp = rightEnd;
	  rightEnd = leftEnd;
	  leftEnd = tmp;
	}

      mapOffset = leftEnd;
      mapFactor = (extentRight - extentLeft) / (leftEnd - rightEnd);
    }
      


  if (!bsPushObj(item)) /* to #Map_position */
    { bsDestroy(item);
      if (count)
	*count = 0;
      return 0;
    }

  if (bsGetData(item, _Multi_Position, _Float, &x1))
    { ourCount = 0;
      do 
	{ mark = bsMark(item, mark);
	  bsPushObj(item); /* To #Map_Error */
	  if (!bsGetData(item, _Error, _Float, &x2))
	    x2 = 0.0; /* default if no error */
	  
	  if (xret) 
	    *(xret++) = mapOffset + mapFactor*x1;
	  if (dxret)
	    *(dxret++) = mapFactor*x2;
	  
	  ourCount++;
	  bsGoto(item, mark);
	} 
      while (ourCount<max && bsGetData(item, _bsDown, _Float, &x1));
      
      bsMarkFree(mark);
      
      if (count) 
	*count = ourCount;
      
      if (objp)
	{ bsGoto(item, 0);
	  *objp = item;
	}
      else
	bsDestroy(item);
      
      return FLAG_ANY_LOCUS;
    }
  
  if (bsGetData(item, _Multi_Ends, _Float, &x1))
    { ourCount = 0;
      do 
	{ mark = bsMark(item, mark);
	  if (!bsGetData(item, _bsRight, _Float, &x2))
	    continue; /* lines with no right */
	  if (x1 < x2)
	    { float dx = (x2 - x1)/2.0;
	      if (xret)
		*(xret++) = mapOffset + mapFactor * (x1 + dx);
	      if (dxret)
		*(dxret++) = mapFactor*dx;
	    }
	  else
	    { float dx = (x1 - x2)/2.0;
	      if (xret)
		*(xret++) = mapOffset + mapFactor * (x2 + dx);
	      if (dxret)
		*(dxret++) =  mapFactor*dx; 
	    }
	  
	  ourCount++;
	  bsGoto(item, mark);
	} 
      while (ourCount<max && bsGetData(item, _bsDown, _Float, &x1));
      
      bsMarkFree(mark);
      
      if (count) 
	*count = ourCount;
      
      if (objp)
	{ bsGoto(item, 0);
	  *objp = item;
	}
      else
	bsDestroy(item);
      
      return FLAG_ANY_INTERVAL;
    }
  
  if (bsGetData(item, _Position, _Float, &x1)) /* Locus-like */
    { bsPushObj(item); /* to #Map_Error */
      if (!bsGetData(item, _Error, _Float, &x2))
	x2 = 0.0; /* default if no error */
      if (objp)
	{ bsGoto(item, 0);
	  *objp = item;
	}
      else
	bsDestroy(item);
      if (xret) 
	*xret = mapOffset + mapFactor*x1;
      if (dxret)
	*dxret = mapFactor*x2;
      if (count)
	*count = 1;
      return FLAG_ANY_LOCUS;
    }
   
  if (bsGetData(item , _Left, _Float, &x1) && /* Interval-like */
      bsGetData(item, _Right, _Float, &x2))
    { if (objp)
	{ bsGoto(item, 0);
	  *objp = item;
	}
      else
	bsDestroy(item);
      
      if (x1 < x2)
	{ float dx = (x2 - x1)/2.0;
	  if (xret)
	    *xret = mapOffset + mapFactor*(x1 + dx);
	  if (dxret)
	    *dxret = mapFactor*dx;
	}
      else
	{ float dx = (x1 - x2)/2.0;
	  if (xret)
	    *xret = mapOffset + mapFactor * (x2 + dx);
	  if (dxret)
	    *dxret = mapFactor * dx;
	}

      if (count)
	*count = 1;
      return FLAG_ANY_INTERVAL;
    }


  if (bsFindTag(item, str2tag("With")) && 
      bsGetKeyTags(item, _bsRight, &key) &&
      bsGetKey(item, _bsRight, &key))
    { int i;
      flags = gMapGetMapObject(key, map, 0, ind-1, xret, dxret, 0, count);
      if (flags && 
	  bsPushObj(item) && bsFindTag(item, str2tag("Relative")) && 
	  bsPushObj(item))
	{ if (bsGetData(item, _Position, _Float, &offset))
	    { if (xret) 
		for (i= 0; i < (count ? *count : 1); i++)
		  xret[i] += offset;
	      flags = FLAG_ANY_LOCUS ;
	    }
	  else if (bsGetData(item, _Left, _Float, &x1) &&
		   bsGetData(item, _Right, _Float, &x2))
	    { float dx ;
	      if (x2 > x1) 
		{ int tmp = x1 ; x1 = x2 ; x2 = tmp ; }
	      dx = (x1 + x2) / 2.0 ;
	      if (xret) 
		for (i= 0; i < (count ? *count : 1); i++)
		  xret[i] += dx ;
	      dx = (x2 - x1) / 2.0 ;
	      if (dxret)
		for (i= 0; i < (count ? *count : 1); i++)
		  dxret[i] = dx ;
	      flags = FLAG_ANY_INTERVAL ;
	    }
	  else
	    flags = 0 ;		/* require Position or Ends */
	}
      if (objp)
	{ bsGoto(item, 0);
	  *objp = item;
	}
      else
	bsDestroy(item);
      return flags;
    }
 
  if (count)
    *count = 0;	

  bsDestroy(item); /* mieg */
  return 0;
}

/**********************************/


/* This generates an array of all the objects on a map. It handles included 
   maps (one level only) and inherits_from maps. The output is a set of 
   (mapped object, map) pairs. The map is the map which directly contains the
   object, to be given as a hint to gMapGetMapObject.
   Note that the same object may not be mapped in more than one map, extra 
   mappings get ignored. The priority for this ensures that inheritance 
   works OK, and that objects mapped on a map and it's sub-map get the sub-map
   mapping ignored.
 */


Array gMapMapObjects(KEY chrom)
{ Array loci, results;
  int i;
  OBJ Chrom ;
  KEY key, locus, submap ;
  KEYPAIR *p;
  
  results = arrayCreate(300, KEYPAIR);
  loci = arrayCreate(200, BSunit);

  i = 256 ;
  while (i--) /* this is used to prevent infinite recursion in inheritance */
    lexClearClassStatus(i, CALCULSTATUS);

  while ((Chrom = bsCreate(chrom)))
  { 
    lexSetStatus(chrom, CALCULSTATUS) ;  /* prevent recursions */

    if (bsFindTag (Chrom,_Contains) && bsFlatten(Chrom, 2, loci))
      { for (i = 1; i<arrayMax(loci); i+=2)
	  { locus =  arr(loci, i, BSunit).k;
	    if (locus && !(CALCULSTATUS & lexGetStatus(locus)))
	      { lexSetStatus(locus, CALCULSTATUS);
		p = arrayp(results, arrayMax(results), KEYPAIR);
		p->obj = arr(loci, i, BSunit).k;
		p->map = chrom;
	      }
	  }
      }
    /* now deal with  included maps */
    /* note that these cannot inherit, and that only on level of 
       inclusion is allowed */
    
    if (bsGetKey(Chrom, str2tag("includes"), &submap))
      { do
	  { OBJ Subchrom = bsCreate(submap);
	    if (Subchrom &&
		bsFindTag (Subchrom,_Contains) && bsFlatten(Subchrom, 2, loci))
	      { for (i = 1; i<arrayMax(loci); i+=2)
		  { locus =  arr(loci, i, BSunit).k;
		    if (locus && !(CALCULSTATUS & lexGetStatus(locus)))
		      { lexSetStatus(locus, CALCULSTATUS);
			p = arrayp(results, arrayMax(results), KEYPAIR);
			p->obj = arr(loci, i, BSunit).k;
			p->map = submap;
		      }
		  }
	      }
	    bsDestroy(Subchrom);
	  }
      while (bsGetKey(Chrom, _bsDown, &submap));
      }
    
    /* now do the same in any inherited maps. */
   
    if (!bsGetKey(Chrom, _From_map, &key) || 
	(CALCULSTATUS & lexGetStatus(key)))
      { bsDestroy(Chrom);
	break;
      }
    
    bsDestroy(Chrom);
    chrom = key;
  }
  
  arrayDestroy(loci);
  return results;
}

KEY gMapKeyOnTag(KEY key, KEY tag)
{ OBJ obj;
  KEY ret;
  BOOL ok;

  if (!tag)
    return key;

  if (!(obj = bsCreate(key))) 
    return key;
  
  ok = bsGetKey(obj, tag, &ret);
  bsDestroy(obj);
  if (ok)
    return ret;
  else
    return key;
}


/****************** giface hook *****************/

void gMapDrawGIF (void)
{
  KEY key, view ;
  KEY _VView ;
  BOOL haveBounds = FALSE ;
  BOOL isWhole = FALSE ;
  BOOL isHideHeader = FALSE ;
  float x1, x2, hdiff ;
  char *word ;
  MAPCONTROL map ;

  lexaddkey ("View", &key, _VMainClasses) ;
  _VView = KEYKEY(key) ;
	      
  if (!freecheck("w"))
    goto usage ;
  
  if (!lexword2key((word=freeword()), &key, _VMap))
    { freeOutf ("// gif map error: map %s not known\n", word) ;
      return ;
    }

  view = 0 ;
  while (freestep ('-'))	/* options */
    if ((word = freeword()))
      {
	if (!strcmp (word, "view") && freecheck("w"))
	  {
	    if (!lexword2key ((word = freeword()), &view, _VView))
	      {
		freeOutf ("// gif map error: view %s not known\n", word);
		return ;
	      }
	  }
	else if (!strcmp (word, "coords") && freecheck("ff"))
	  {
	    freefloat (&x1) ; freefloat (&x2) ;
	    haveBounds = (x2 != x1) ;
	    if (!haveBounds)
	      {
		freeOutf ("// gif map error: coords given are the same (%g)\n", x1) ;
		return ;
	      }
	  }
	else if (!strcmp (word, "whole"))
	  isWhole = TRUE ;
	else if (!strcmp (word, "hideheader"))
	  isHideHeader = TRUE ;
	else
	  goto usage ;
      }
    else
      goto usage ;

  if (freecheck ("w"))
    goto usage ;

  { extern BOOL displayReportGif ;
    displayReportGif = FALSE ;
    display (key, view, GMAP) ;	/* the primary display */
    displayReportGif = TRUE ;
  }

  map = currentMapControl () ;

  if (isHideHeader)
    { map->control->hideHeader = TRUE ;
      map->control->topMargin = 1 ;
    }
  if (haveBounds)
    { map->centre = (x1 + x2) / 2.0 ;  
      map->mag = (map->control->graphHeight - 
		  map->control->topMargin - 2.0) / (x2 - x1) ;
    }
  else if (isWhole)
    {  BOOL isNeg = (map->mag < 0) ;
       map->centre = (map->max + map->min) / 2.0 ;  
       map->mag = (map->control->graphHeight - 
		   map->control->topMargin - 2.0) /
	 (map->max - map->min) ;
       if (isNeg) map->mag = -map->mag ;
    }

  if (haveBounds || isHideHeader || isWhole)
    controlDraw () ;

  hdiff = (map->control->graphHeight - 5) / (2.0 * map->mag) ;
  x1 = map->centre - hdiff ;
  x2 = map->centre + hdiff ;
  freeOutf ("// MAP  %s %g %g", freeprotect(name(key)), x1, x2) ;
  if (view) freeOutf (" view=%s", freeprotect(name(view))) ;
  freeOut ("\n") ;

  return ;

usage:
  freeOut ("// gif map error: usage: MAP map [-view view] [-coords x1 x2]\n") ;
}
 
/*******************************************************************/
/********************** end of file ********************************/
