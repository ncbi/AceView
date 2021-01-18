

/* Ace DB  statics shutdown */
void poorStaticsShutDown(void)
{
	void freeshutdown(void);
	void messAllocShutdown(void);

	freeshutdown();
	messAllocShutdown();
}











/*drawing wrapped text in xacembly */
static int drawBlockText(char * src,char * separ, int maxnum,int x, int y,int charrayLen)
{
  int i,k,iy=0,lastspace=0,lastpos=0;
  char dst[1024],ch;
  
  lastspace=charrayLen;
  lastpos=charrayLen;
  for(k=0,i=0;(ch=src[i])!=0;i++){
      if(ch==' ' || ch=='\n'){ 
		  lastspace = k; lastpos = i; 
	  }
      
		/* put endline and printf the piece */
		if(ch=='\n' || k>=charrayLen){
			i=lastpos;
			if(lastspace){
				dst[lastspace]=0;
				graphText(dst,x,y+iy);
			}
			k=0; /* and reset the desttination */
			iy++;
			lastspace=charrayLen;
			lastpos= i + charrayLen ;
			continue;
		}
      
		dst[k]=ch;k++;
		if (maxnum && i>= maxnum) 
		break;
    }
	if(k){
		dst[k]=0;
		graphText(dst,x,y+iy);
    }
	return iy;
}











/* external prodcedure calls in xacembly */
static void ficheEditorOut (FILE *f, void *vp)
{
}


static void ficheEditorIn(KEY gene)
{
	externalAsynchroneCommand ("ficheEditor.csh", fileName, assVoid (gene), ficheEditorOut) ;
}
































































#ifdef DONTOCOMPILE 

	/* build the colorMap and 3DMap */
	if( (ptr=argVal("-map")) && 
		((sscanf(ptr,"%dx%d",&dimX,&dimY))+1) && 
		(dimX && dimY) &&
		(mapClrArray=(clrDEF * )memNew(sizeof(clrDEF)*dimX*dimY)) &&
		(mapNumArray=(double * )memNew(sizeof(double)*dimX*dimY)) ){
		
		memset((void*)mapNumArray,0,sizeof(double)*dimX*dimY);
		/* build the map */
		graphMapTree(seqTree);

		/* output it as bitmap */
		clrSaveBmpImage((ptr=fileMakePath(".",treeFile,"bmp")),dimX,dimY,mapClrArray);memFree(ptr);
		/* output it as numerical 2D array of Z values */
		printMap((ptr=fileMakePath(".",treeFile,"dat")),dimX,dimY,mapNumArray);memFree(ptr);
	}
	memFree(mapClrArray);
	memFree(mapNumArray);
    
/*
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
_/
_/ Tree Mapping
_/
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/

static int mapSequence(char * subSeq,int hits,int level,tTREE * tTree)
{
	double	divisorX,hitsF;
	int		iX,iY,sx,sy,fx,fy,divisorY,ofset,ofsStpY;
	clrDEF	clr,r,g,b;

	/* if(hits<1)return 1; */
	/*hitsF=sqrt((1.*hits)/tTree->stat[level].maxHits); */
	hitsF=(1.*hits)/tTree->stat[level].averageHits; 
	/* hitsF=(hits*tTree->stat[level].inProbab/tTree->symbCnt); *//* probabilistic significance of the hits */
	

	/* determine coordinates on the map */
	sx=(int)(tTreeFixedPosition(tTree,subSeq,level,&divisorX)*dimX);
	divisorY=dimY/(tTree->nSubMax);
	sy=(level-1)*divisorY;
	fx=sx+(int)(ceil(divisorX*dimX));
	fy=sy+divisorY;
	if(fx==sx)fx++; /* at least on line */
	if(fy==sy)fy++; /* at least on column */
	
	/* determine the color to draw */
	clr=nucColor[(int)tTreeNuc(tTree,subSeq[level-1])];
	r=clrR(clr);g=clrG(clr);b=clrB(clr);
	r=(int)(r*hitsF);g=(int)(g*hitsF);b=(int)(b*hitsF);
	if(r>0xFF)r=0xFF;if(g>0xFF)g=0xFF;if(b>0xFF)b=0xFF;
	clr=clrRGB(r,g,b);clr=vtINTL4_TO_ORIGINALBYTEORDER(clr);
	
	ofsStpY=dimX-(fx-sx);
	for(ofset=sx+sy*dimX,iY=sy;iY<fy;iY++,ofset+=ofsStpY){
		for(iX=sx;iX<fx;iX++,ofset++){
			if(mapNumArray[ofset]>hitsF)continue;	/* this is smaller hit */
			mapNumArray[ofset]=hitsF; /* hitsF*/
			mapClrArray[ofset]=clr; 
		
		}
	}

	/*printMap("map.dat",dimX,dimY,mapNumArray); */
	return 1;
}

void graphMapTree(tTREE * tTree)
{
	tTreeEnumerate(tTree,0,0,(tTREEEnumFunction)mapSequence,(void*)tTree);
}

void	printMap(char * fileName,int dimX,int dimY,double * mapArray)
{
	int	iX,iY;
	FILE * fp;

	if(!(fp=fopen(fileName,"w")))
		return ;
	for(iX=0;iX<dimX;iX++){
		for(iY=0;iY<dimY;iY++){
			fprintf(fp,"%6.2f ",mapArray[iX+dimX*iY]);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
}



	    /* print frequency chart */
    if(argIs("-basis")){

		int		stp,i,startIt,endIt,lenIt,iP;
		double	possiB;
		
		fp=fopen((ptr=fileMakePath(".",fileName,"bas.dat")),"w");memFree(ptr);
		if(fp){	
			for(i=0;i<seqLen;i++){
				for(stp=2;(i%stp)!=(stp>>1)-1;stp<<=1); /* determine the divisor */
				
				startIt=i-(stp>>1);if(startIt<0)startIt=0;
				endIt=i+(stp>>1);if(endIt>seqLen)endIt=seqLen;
				lenIt=endIt-startIt+1;
				
				tKnot=tTreeGetKnot(seqTree,atgcSequence+startIt,lenIt);
				if(!tKnot)hits=0;
				else hits=tKnot->hits;
				
			
				for(possiB=1.,iP=0;iP<lenIt;iP++)possiB*=4;
				fprintf(fp,"%6d %-6d %6d %8.2f\n",i,stp,hits,(hits>1)  ? hits*possiB : 0 );
				
			}
			fclose(fp);
		}

	}

clrDEF	nucColor[]={(0x0000ff),	/* blue */
					(0x00ff00),	/* dark green */
					(0xff0000),	/* red */
					(0xffffff)	/* white - gray */
					};


clrDEF *	mapClrArray=0;
double *	mapNumArray=0;
int			dimX=480,dimY=24;

    
#endif
