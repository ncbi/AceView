#include <vstd.h>


/*
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
_/                                          _/
_/  MISCELANEOUS FUNCTIONS                  _/
_/                                          _/
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/

/* returns bits limited hash for memory piece */
int vmemHashMem(void * mem, int len, int bits, int isDiff)
{
	int i ;
	unsigned int j, x = 0 ;
	int rotate = isDiff ? 21 : 13 ;
	int leftover = 8*sizeof(int) - rotate ;
	char * src=(char * )mem;

	/* preparing the hash */
	for(i=0;i<len;i++)
		x = src[i] ^ (( x >> leftover) | (x << rotate)) ; 

	/* compress down to given number of bits */
	for (j = x, i = bits ; i < sizeof(int) ; i += bits)j ^= (x >> i) ;
	j &= (1 << bits) - 1 ;

	if(isDiff)j|=1;
	return j ;
}

/*
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
_/                                          _/
_/  MEMBLOCK FUNCTIONS                      _/
_/                                          _/
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/

#define	memStrategyPiece		(2)
#define	memMinimalPiece			(1024)

/* 
	this function adds automatically reallocable piece of memory 
	add can be zero. In this case only the size of buffer will be changed 
	but nothing will be added to the buffer.
*/
int vmexReplace(vMEX * mem,int pozReplace,void * add,int sizeAdd,int sizeDel)
{   											     
    int     pos,sizeRequired,newSize;
    char *  newBfr;
    
    if(!mem || sizeAdd<0)
		return vINVIVAL;

	if(!mem->buf){ /* fix the size and position if the buffer was freed manually */
		mem->cnt=mem->curSize=mem->curPos=0;
		if(!mem->piece) mem->piece = memStrategyPiece;
		else if(mem->piece!=memStrategyPiece) /* always paragraph aligned number */
			mem->piece=(mem->piece+0x10)&0xFFFFFFF0;
		mem->buf=0;	
	}

	pos=mem->curPos ; /* for returning only */
	
	/* fix the position and size od replacement if out of */
	if(pozReplace<0)pozReplace=mem->curPos;
	if(pozReplace+sizeDel>mem->curPos)sizeDel=mem->curPos-pozReplace;

	sizeRequired=mem->curPos+sizeAdd-sizeDel;
	if(sizeRequired<=0){
		messfree(mem->buf);
		return vINVIVAL;
	}

	/* if memory of each entry should be aligned */
	if(mem->flags&vMEXALIGN){ 
		sizeRequired=(sizeRequired+0x10)&0xFFFFFFF0;
	}

	if(sizeRequired<memMinimalPiece)sizeRequired=memMinimalPiece;
    if(sizeRequired>mem->curSize) {

		newSize=mem->curSize;
		if(!newSize)newSize=memMinimalPiece;

		/* determine the new size */	
		while(newSize<sizeRequired) {

			if(mem->piece!=memStrategyPiece)  /* additive strategy */
				newSize+=mem->piece; 
			else 							 /* multiplicative strategy */
				newSize*=memStrategyPiece;
		}

		/* reallocate the buffer */	
		if( newSize && (newBfr=messalloc(newSize)) ) {
			
			/* copy the content before replacement pozition and after insertion pozition */
			if(pozReplace)
				memmove(newBfr,mem->buf,pozReplace); 
			if(mem->curPos-pozReplace)
				memmove(newBfr+pozReplace+sizeAdd,mem->buf+pozReplace+sizeDel,mem->curPos-pozReplace);
			
			messfree(mem->buf);
			mem->buf=newBfr;
			mem->curSize=newSize;
		}
		else return vINVIVAL;
    }
    
	if(sizeAdd){
		if(add){
	        memcpy(mem->buf+pozReplace,add,sizeAdd);
			mem->cnt++;
		}
		else if(mem->flags&vMEXSETZERO){
			memset(mem->buf+pozReplace,0,sizeAdd);
		}
	}
	if(sizeDel)
		mem->cnt--;

    mem->curPos+=sizeAdd-sizeDel;	
	if(mem->flags&vMEXALIGN){ 
		mem->curPos=(mem->curPos+0x10)&0xFFFFFFF0;
	}


    return pos;
}





/*
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
_/                                          _/
_/  MEMORY MAPPING FUNCTIONS                _/
_/                                          _/
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/

char * vmemMapFile(char * fileName,unsigned long * sizeMem, int * fileHandle)
{
	int		memFileHandle,zero[vFILEBLOCKSIZE];
	unsigned long fileLength,i,size;
	char	* mapAddress;

	/* open the map file */
	memFileHandle=open(fileName,O_CREAT|O_RDWR,S_IWUSR | S_IRUSR | S_IRGRP);
	if(	memFileHandle==-1 )return 0;
	
	size=*sizeMem;
	fileLength=lseek(memFileHandle,0,SEEK_END);
	if(fileLength<size){
		printf("extending memoryMap file %s from %lu up to %lu bytes...",fileName,fileLength,size);
                memset (zero, 0, sizeof(zero)) ;
		for(i=0;i<size-fileLength+1;i+=sizeof(zero))
			write(memFileHandle,&zero,sizeof(zero));
		printf("done\n");
	}else size=fileLength;
	*sizeMem=size;
	
	/* map it */
	mapAddress=mmap(0,size,PROT_READ|PROT_WRITE,MAP_SHARED,memFileHandle,0);
	if(!mapAddress){
		close(memFileHandle);
		return 0;
	}
	if(*fileHandle)*fileHandle=memFileHandle;
	return mapAddress;
}

vIN void vmemUnmapFile(char * mapAddress,int size,int fileHandle)
{
	if(mapAddress)munmap((void *)mapAddress,size);
	if(fileHandle!=-1)close(fileHandle);
}

                                           

/*
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
_/                                          _/
_/  MEMORY GRAFS FUNCTIONS                  _/
_/                                          _/
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/

#define	megExtraCollisionReducer 4
#define	megInitialDim 		64
#define	megInitialBit		6

/* find the slot of this id,lenId is in our dictionary. 
pLasthash is set to the last considered hash even if 
it was an empty slot. Returns the slot number. */
int _vmegFindSlot(vMEG * meg,int * pIndex,void * id,int lenId, int * pLastHash)
{
	int hash,dh=0,iFnd;
	vMEGITM * itm;
	char * idPtr;
	
	if(!meg || !id || !lenId || !(meg->dim) )
		return vFALSE;

	hash=vmemHashMem(id,lenId,meg->bits,0);
	if(pLastHash)*pLastHash=hash;

	while (vTRUE){ 
		if(!(iFnd=((vMEGITM *)meg->tbl.buf)[hash].hashTo)){ /* not found - hash points to nowhere */
			return vFALSE;
		}
		if(iFnd!=vMAXINUM){
			iFnd--;
			itm=&((vMEGITM *)meg->tbl.buf)[iFnd];

			/* compare if this is the item we are looking for */
		    if(	itm->idLen==lenId && 
				( (idPtr=itm->idOfs+meg->dat.buf) == id || /* check if meself (frequent)? */ 
				(!memcmp(idPtr,id,lenId)) ) 
				){
				if(pIndex)*pIndex=iFnd;
				return vTRUE;
			}
		}

		/* collision ! */
		if(!dh)dh = vmemHashMem(id,lenId,meg->bits,1);
	    hash += dh;if (hash >= meg->dim) hash-=meg->dim ;
		meg->collisions++;
		if(pLastHash)*pLastHash=hash;
	}

}			    


/* rehashes all the dictionary */
static void _vmegRehashAll(vMEG * meg)
{
	int 	i,hash;
	vMEGITM * itm;

	for(i=0;i<meg->cnt;i++)
		((vMEGITM *)meg->tbl.buf)[i].hashTo=0;
	meg->collisions=0;

	for(i=0;i<meg->cnt;i++){
		itm=&((vMEGITM *)meg->tbl.buf)[i];

		_vmegFindSlot(meg,0,itm->idOfs+meg->dat.buf,itm->idLen,&hash);
		((vMEGITM *)meg->tbl.buf)[hash].hashTo=i+1;
	}
}

/* adds the item into the dictionary */
int vmegAddBin(vMEG * meg,int * pIndex,void * id,int lenId, void * data,int lenData)
{
	int iFnd,obits,hash=0;
	vMEGITM * itm;

	if ( !meg || !id || !lenId || /* wrong data ? */
		(_vmegFindSlot(meg,pIndex, id, lenId,&hash) ) /* or this data are already known */
		){
		return vFALSE;
	}

	iFnd=meg->cnt;
	obits=meg->bits;

	/* determine the dimension of table and bitness; happens only once */
	if(!meg->dim){
		meg->tbl.flags|=vMEXSETZERO;
		meg->dat.flags|=vMEXSETZERO|vMEXALIGN;
		meg->dim=megInitialDim; meg->bits=megInitialBit;
	}
	/* reallocation is neccessary ? */
	while(meg->dim <= meg->cnt*megExtraCollisionReducer){
		meg->dim<<=1, meg->bits++;
	}

	/* in case if the size was changed : reserve enough element space for dictionary */
	if( obits!=meg->bits && (vINVIVAL==vmexExpand(&meg->tbl,sizeof(vMEGITM)*(meg->dim))) )
		return vFALSE;


	/* initialize the item in the table */
	itm=&((vMEGITM *)meg->tbl.buf)[meg->cnt];
	itm->idOfs=vmexAppend(&meg->dat,id,lenId);
	itm->idLen=lenId;
	if(lenData){itm->dataOfs=vmexAppend(&meg->dat,data,lenData);itm->dataLen=lenData;}
	meg->cnt++;
	
	if(meg->bits!=obits)_vmegRehashAll(meg);
	else {
		if(!hash)hash = vmemHashMem(id,lenId,meg->bits,0);/* if it wasn't hashed by vFind */
		((vMEGITM *)meg->tbl.buf)[hash].hashTo=iFnd+1;
	}

	if(pIndex)*pIndex=iFnd;

	return vTRUE;
}



/* finds the id in the dictionary */
vIN int vmegFindBin(vMEG * meg,int * pIndex, void * id,int lenId)
{
	return _vmegFindSlot(meg,pIndex,id,lenId,0);
}

/* empties the dictionary */
void vmegEmpty(vMEG * meg) 
{
	if(!meg)return;

	vmexEmpty(&meg->tbl);
	vmexEmpty(&meg->dat);
	meg->dim=meg->cnt=meg->bits=meg->collisions=0;
}

char * _vmegId(vMEG * meg,int iFnd){
	if(iFnd>=meg->cnt || iFnd<0)return "out of boundary";
	return (char *)vmegId(meg,iFnd);
} 
vMEGITM * _vmegItm(vMEG * meg,int iFnd){return vmegItm(meg,iFnd);}
char * _vmegData(vMEG *meg,int iFnd){return vmegData(meg,iFnd);}

