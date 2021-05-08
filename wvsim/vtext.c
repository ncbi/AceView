#include <vstd.h>

/*
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
_/                                          _/
_/  Dictionary sorting functions            _/
_/                                          _/
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/


static vMEG * vtextSortDic=0;
static int vtextSortDirection=1;
static int _vtextSortByScore ( int * pi,int i, int j )
{
	int dif=vtextTreeWordScore(vtextSortDic,pi[j])-vtextTreeWordScore(vtextSortDic,pi[i]);
	if (dif) return vSi0(dif)*vtextSortDirection ;	
	return strcmp ((char *)vmegId(vtextSortDic,i), (char *)vmegId(vtextSortDic,j)) ;
}

static int _vtextSortByDicScore(vMEG * dic,int i, int j)
{
	int dif=vtextTreeWordScore(dic,j)-vtextTreeWordScore(dic,i);
	if (dif) return vSi0(dif)*vtextSortDirection;
	return strcmp ((char *)vmegId(dic,i), (char *)vmegId(dic,j)) ;
}

/*
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
_/                                          _/
_/  Text tree functions                     _/
_/                                          _/
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/

/* link two items if they were not linked before */
void vtextTreeLinkAttach(vMEX * ref,int  ix)
{
	int i,* pi=(int *)(ref->buf);
	
	for(i=0;i<ref->cnt;i++)if(pi[i]==ix)return; /* check if ix is attached ? */
	vmexAppend( ref, &ix, sizeof(ix));/* not yet - link them */
}

/* sorts all sublists in <dic> dictionary by their score in 
	other <vtextSortDic> dictionary .
	Sorts accending if sortType =1, descending if -1 and does not 
	sort if 0  */
void vtextTreeSort (vMEG * dic,int lnk,vMEG * dicSort,int sortType)
{
	int * stk, * ind,i,iWrd;
	vMEX * ref;

	/* allocate if needed */
	if(!dicSort || !(stk=messalloc(sizeof(int)*2*dicSort->cnt)) )return ;
	ind=stk+dicSort->cnt;
	vtextSortDic=dicSort;
	vtextSortDirection=sortType;

	for(iWrd=0;iWrd<dic->cnt;iWrd++){
		ref=vtextTreeWordRef(dic, iWrd, lnk );if(ref->cnt<=1)continue;
		vsortQuickCallbackIndex(0,ref->buf,ref->cnt,stk,ind,(vSORTFunction)_vtextSortByScore);

		/* update the sorting in the place */
		for(i=0;i<ref->cnt;i++)stk[i]=((int *)ref->buf)[ind[i]];
		for(i=0;i<ref->cnt;i++)((int *)ref->buf)[i]=stk[i];
	}

	messfree(stk); /* I created you, I shall destroy you! (God !) */
}


/*  creates dictionary-tree for sentences <dicStc> and <dicWrd> based on 
	text <in> wich iz DZ terminated set of sentneces. 
	Also sorts accending if sortType =1, descending if -1 and does not 
	sort if 0 */														 
void vtextTreeCreate (char * in,vMEG * dicStc, vMEG * dicWrd,int sortType,vCallbackLogical filterFunc)
{
	int iWrd, iStc, iPrv;
	char * stc, * wrd, * nxt,stcBuf[vTEXTMAXLINE],wrdBuf[vTEXTMAXLINE];
	vTEXTTREEWRD	lnk;
	#define	isendline(v_ch)  ( (v_ch)==';' || (v_ch)=='\r' || (v_ch)=='\n')

	vSet0(lnk);lnk.score=1;lnk.cnt=1;

	/* add the beginning and the end of the sentence as words */
	if(dicWrd->cnt==0){
		if( !vmegAdd(dicWrd,&iWrd,"__begin",&lnk,sizeof(lnk))  || 
			!vmegAdd(dicWrd,&iWrd,"__end",&lnk,sizeof(lnk)) ) 
			return ;
	}

	for( stc=in; stc; stc=vstrNextList(stc,1)){

		if(sscanf(stc,"%d",&lnk.score)){
			nxt=stc;
			if(!(stc=vstrSkipWords(stc,1,vSTRBLANK))){
				stc=nxt;continue;
			}
		}

		vstrCleanEnds(stc,stcBuf,0,vSTRBLANK,1);

		if(filterFunc && filterFunc((void *)stcBuf,0,dicStc,&lnk)) /* if filtered out */
			continue;
		
		if(!vmegAdd(dicStc,&iStc,stcBuf,&lnk,sizeof(lnk))) {
			vtextTreeWordScore(dicStc,iStc)+=lnk.score; /* if it was there before, just add it */
			vtextTreeWordCnt(dicStc,iStc)+=lnk.cnt; /* if it was there before, just add it */				
		}else {
			/* link the sentence with begin word */
			vtextTreeLinkAttach( vtextTreeWordRef(dicWrd,vTEXTTREEBEGINWRD,vTEXTTREELNK),iStc);
			vtextTreeLinkAttach( vtextTreeWordRef(dicStc,iStc,vTEXTTREELNK),vTEXTTREEBEGINWRD);
			/* link the sentence with end word */
			vtextTreeLinkAttach( vtextTreeWordRef(dicWrd,vTEXTTREEENDWRD,vTEXTTREELNK),iStc);
			vtextTreeLinkAttach( vtextTreeWordRef(dicStc,iStc,vTEXTTREELNK),vTEXTTREEENDWRD);
		}

		for(iPrv=vTEXTTREEBEGINWRD, wrd=stcBuf; wrd ; wrd=vstrSkipWords(wrd,1,vSTRBLANK), iPrv=iWrd){
				
			nxt=strpbrk(wrd,vSTRBLANK);if(nxt)*nxt=0;
			vstrCleanEnds(wrd,wrdBuf,0,vSTRBLANK,1);

			if(filterFunc && filterFunc((void *)wrdBuf,1,dicWrd,&lnk) ){ /* if filtered out */
				if(nxt)
				  *nxt=vSTRBLANK[0];
				iWrd=vTEXTTREEBEGINWRD;continue;
			}

			if( !vmegAdd(dicWrd,&iWrd,wrdBuf,&lnk,sizeof(lnk)) ) {
				vtextTreeWordScore(dicWrd,iWrd)+=lnk.score; /* if it was there before, just add it */
				vtextTreeWordCnt(dicWrd,iWrd)+=lnk.cnt; /* if it was there before, just add it */
			}

			/* link the sentence and the word uniquely */
			vtextTreeLinkAttach( vtextTreeWordRef(dicWrd,iWrd,vTEXTTREELNK),iStc);
			vtextTreeLinkAttach( vtextTreeWordRef(dicStc,iStc,vTEXTTREELNK),iWrd);
			/* attach neighbours */
			vtextTreeLinkAttach( vtextTreeWordRef(dicWrd,iWrd,vTEXTTREELNKPRV),iPrv);
			vtextTreeLinkAttach( vtextTreeWordRef(dicWrd,iPrv,vTEXTTREELNKNXT),iWrd);

			if(nxt)*nxt=vSTRBLANK[0];
		}

		if(iPrv!=vTEXTTREEBEGINWRD){ /* attach to the end sentence */
			vtextTreeLinkAttach( vtextTreeWordRef(dicWrd,vTEXTTREEENDWRD,vTEXTTREELNKPRV),iPrv);
			vtextTreeLinkAttach( vtextTreeWordRef(dicWrd,iPrv,vTEXTTREELNKNXT),vTEXTTREEENDWRD);
		}
	}
	
	for(iStc=0;iStc<dicStc->cnt;iStc++){ /* averaging the score per repeat */
		vtextTreeWordScore(dicStc,iStc)/=vtextTreeWordCnt(dicStc,iStc);
	}
	for(iWrd=0;iWrd<dicWrd->cnt;iWrd++){ /* averaging the score per repeat */
		vtextTreeWordScore(dicWrd,iWrd)/=vtextTreeWordCnt(dicWrd,iWrd);
	}

	if(sortType!=0){
		vtextTreeSort(dicWrd,vTEXTTREELNKNXT,dicWrd,sortType);
		vtextTreeSort(dicWrd,vTEXTTREELNKPRV,dicWrd,sortType);
		vtextTreeSort(dicStc,vTEXTTREELNK,dicWrd,sortType);
		vtextTreeSort(dicWrd,vTEXTTREELNK,dicStc,sortType);
	}
}

/* cleans all */
void vtextTreeEmpty (vMEG * dic)
{
	vMEX * sub;
	int	i,iFnd;

	for(iFnd=0;iFnd<dic->cnt;iFnd++){
		for(i=0;i<3;i++){
			sub = vtextTreeWordRef(dic,iFnd,i); 
			vmexEmpty(sub);
		}
	}
	vmegEmpty(dic);
}

/* prints the given <dic2,what> substructure of dictionary <dic1> tree */
void vtextTreePrint (vSTR * out, char * title,vMEG * dic1, vMEG * dic2,int what, int * indSort)
{
	int i1,i2,num,me;
	vMEX * ref;
	vSTR my;

	if(!out){vSet0(my);my.mode=vSTRECHO|vSTRDISCARD;out=&my;}

	if(title)vstrPrintf(out,"\n_/_/_/_/_/_/_/_/_/_/\n_/\n_/ %s\n_/\n\n",title);

	for(i1=0;i1<dic1->cnt;i1++){
		if(indSort)me=indSort[i1]; else me=i1;

		ref=vtextTreeWordRef(dic1,me,what);
		vstrPrintf(out,"<%-5d>  %8d %4d [%s]\n",me,vtextTreeWordScore(dic1,me),vtextTreeWordCnt(dic1,me), vmegId(dic1,me));
		
		if(!ref->cnt || !dic2)continue;
		for(i2=0;i2<ref->cnt;i2++){

			num=vmexArr(ref,i2,int);
			vstrPrintf(out,"    <%-5d> %8d %4d %s\n",i2,vtextTreeWordScore(dic2,num),vtextTreeWordCnt(dic2,num),vmegId(dic2,num));
		}
	}
}


/*
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
_/                                          _/
_/  Dictionary functions                    _/
_/                                          _/
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/

/* sorts dictionary by the score its own score  
	Sorts accending if sortType =1, descending if -1 and does not 
	sort if 0  */
void vtextDicSortIndex(vMEG * dic,int sortType, int * pInd, int maxInd)
{
	int * stk, * ind,i;

	/* allocate if needed */
	if(!dic || !dic->cnt || !(stk=messalloc(sizeof(int)*2*dic->cnt)) )return ;
	ind=stk+dic->cnt;
	vtextSortDirection=sortType;

	vsortQuickCallbackIndex(0,dic,dic->cnt,stk,ind,(vSORTFunction)_vtextSortByDicScore);
	for(i=0;i<dic->cnt && i < maxInd ;i++)pInd[i]=ind[i];

	messfree(stk); /* I created you, I shall destroy you! (God !) */
}


char * vtextDicSentenceCompose(vMEG * dicWrd,int * pWrdList, int maxwrd, char * pStc, int * pScore)
{
	char * stc,* wt;
	int iWrd;
	if(!(stc=pStc))return 0;

	if(pScore)*pScore=0;

	for ( iWrd=0 ; iWrd<maxwrd ; iWrd++) {
/*		if(pWrdList[iWrd]==vTEXTTREEBEGINWRD)wt="<";
		else if(pWrdList[iWrd]==vTEXTTREEENDWRD)wt=">";
		else */
		wt=(char *)vmegId(dicWrd,pWrdList[iWrd]);
 		stc+=sprintf(stc,"%s ",wt);
		if(pScore)*pScore+=vtextTreeWordScore(dicWrd,pWrdList[iWrd]);
	}
	if(*(stc-1)==' ')*(stc-1)=0;
	return pStc;
}


int vtextDicSentenceConsume(char * stc,vMEG * dicWrd,int * pCnt,int sizeDat)
{
	int iWrd,oriCnt;
	char * wrd,* nxt,wrdBuf[vTEXTMAXLINE];


	if(pCnt)*pCnt=0;
	oriCnt=dicWrd->cnt;
	for(wrd=stc; wrd && *wrd ; wrd=vstrSkipWords(wrd,1,vSTRBLANK)){
	
		nxt=strpbrk(wrd,vSTRBLANK);if(nxt)*nxt=0;
		vstrCleanEnds(wrd,wrdBuf,0,vSTRBLANK,1);

		vmegAdd(dicWrd,&iWrd,wrdBuf,0,sizeDat);

		if(nxt)*nxt=vSTRBLANK[0];
		if(pCnt)(*pCnt)++;
	}

	return dicWrd->cnt-oriCnt;
}



/*
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
_/                                          _/
_/  Sentence Comparison functions           _/
_/                                          _/
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/

int vtxtCompareSentences(char * seq1,char * seq2,vMEG * dicRes)
{
	int	i,j,id[2],iRes,len,maxLen;
	typedef struct {
		int i,j;
		int len;
		}_vtxtCMP;
	_vtxtCMP * tab;
	#define	Tab(_i,_j)	(tab[(_i)*maxLen+(_j)])

	maxLen=strlen(seq1);if((len=strlen(seq2))>maxLen)maxLen=len;

	/* allocate the memory */
	if(!(tab=(_vtxtCMP * )messalloc(sizeof(_vtxtCMP)*maxLen*maxLen)))
		return 0;
	memset(tab,0,sizeof(_vtxtCMP)*maxLen*maxLen);

	for(i=0;seq1[i];i++){
		for(j=0;seq2[j];j++){
			if(/*seq1[i]!=' ' && */seq1[i]==seq2[j]){
				Tab(i,j).len = 1 ;

				if(i>0 && j>0 && Tab(i-1,j-1).len ){
					len=Tab(i,j).len+=Tab(i-1,j-1).len;
					id[0]=Tab(i,j).i=Tab(i-1,j-1).i;
					id[1]=Tab(i,j).j=Tab(i-1,j-1).j;
					
					/* remember this hit in the dictionary */
					if(!vmegAddBin(dicRes,&iRes,id,sizeof(id),&len,sizeof(len) )) 
						* ((int *)vmegData(dicRes,iRes))=len;
				}
				else {
					Tab(i,j).i=i;
					Tab(i,j).j=j;
				}
			}
		}
	}
	messfree(tab);
	return dicRes->cnt;
}



