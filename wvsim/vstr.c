#include "regular.h"
#include <vstd.h>
#include <vtxt_.h>

/*
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
_/                                          _/
_/  BASIC STRING FUNCTIONS                  _/
_/                                          _/
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/


/* glues a list of strings given as arguments
    returns the destination address */
char * vstrGlueList(int cnt, char * firstPtr, ... )   /* CHECK : mem */
{
    int len=0,i;
    char * ptr = firstPtr,*bfr,*bfrPtr;
    va_list marker;

    va_start(marker, firstPtr );
    for(i=0;i<cnt;i++){
        if(ptr){len+=strlen(ptr);}
        ptr = va_arg(marker, char *);
    }
    va_end( marker );
    if(!len || !(bfr=messalloc(len+2)))return 0;

    bfrPtr=bfr;
    ptr=firstPtr;

    va_start(marker, firstPtr );
    for(i=0;i<cnt;i++){
        if(ptr){strcpy(bfrPtr,ptr);bfrPtr+=strlen(bfrPtr);}
        ptr = va_arg(marker, char *);
    }
    va_end( marker );

    return bfr;
}



/*
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
_/                                          _/
_/  DOUBLE ZERO TERMINATED STRING FUNCTIONS _/
_/                                          _/
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/

/* returns string number cnt in the string list where strings are
   separated by \0 symbol. The double "\0\0" is the end of list */
char * vstrNextList(char * cont,int cnt)
{
    if(!cont)return 0;
    while(cnt){
        while(*cont){cont++;}
        cont++;
        if(!(*cont))break;
        cnt--;
    }
    if(!(*cont) || cnt)return 0;
    return cont;
}

/* counts string list items */
int vstrCntList(char * cont)
{
    int cnt;

    if(!cont)return 0;
    for(cnt=0;(*cont);cnt++){
        while(*cont){cont++;}
        cont++;
    }
    return cnt;
}

/* compares if the curent src string is equal to one of strings in stringlist list
    returns -1 if not found, or the length of found string otherwise
	puts numfnd to the ordinal of found string if not zero */
int vstrIsOnList(char * src,char * list,int * numfnd,int isCaseInSensitive)
{
    int sc,i,cnt;

    if(!list[0]){if(numfnd)*numfnd=0;return 0;}

    for(cnt=0,sc=0;list[sc];cnt++) { /* check if current string is equal to one of column separators in the separator stringList */

      if(!isCaseInSensitive)for(i=0;list[sc] && src[i]==list[sc];i++,sc++); /* compare while equal */
	  else for(i=0;list[sc] && toupper((int)src[i])==toupper((int)list[sc]);i++,sc++); /* compare while equal case insensitive */

      if(!list[sc]){/* correposndence */
          if(numfnd)*numfnd=cnt; /* return the ordinal number of the found string  */
          return i;
      }
      while(list[sc])sc++; /* check next strin in stringlist for correspondence */
      sc++;
    }
    return -1;
}

/* double the zero termination */
vIN int vstr00(char * src)
{
	int len=strlen(src);

	src[len+1]=0;

	return len;
}


/*
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
_/                                          _/
_/  STRING SEARCH FUNCTIONS                 _/
_/                                          _/
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/
/*  looks for n-th occurence of one of given "find" stringlist before stopFind is reached
    src -source string
    find - start of the tag to look for (is a list ending with two zero symbols )
    stopFind - end of the tag to look for (is a list ending with two zero symbols )
                if stopFind=="" search will be performed until the end of the string
    occurence  - wich occurence to pick up
    returns the found position or zero if not found
*/
char * vstrSearchSubstring(char * src,char * find,int occurence,char * stopFind,int isCaseInSensitive)
{
    int i,iFound=-1,ifnd=0,pos;

    if(!src || !find)return 0;

   /* find the string find if that is before the string endTag */
    for(i=0;src[i];i++){

        /* scan if the current position corresponds to the find */
        pos=vstrIsOnList(src+i,find,0,isCaseInSensitive);
        if(pos!=-1){ /* found find let's see if that is the occurence we need */
            if(ifnd+1==occurence || occurence==-1)iFound=i;
            else ifnd++;
            if(occurence!=-1 && iFound!=-1)break;
        }

        /* scan if the current position corresponds to the endTag */
        if(stopFind && stopFind[0] && vstrIsOnList(src+i,stopFind,0,isCaseInSensitive)!=-1)
            break;
    }
    if(iFound==-1)return 0;
    return src+iFound;
}

/*  compares two strings str1,str2 until a symbol from the symblist string
	is reached.
	Returns  zero if not equal otherwise returns the length of equal part */
int vstrCompareUntilSymbol(char * str1,char * str2,char * symblist,int isCaseInSensitive)
{
    int   i=0,j;

    for(j=0;(*str1) && (*str2);j++){ /* until the end of the strings */

        if(symblist){
            for(i=0;symblist[i];i++){if(symblist[i]==(*str2))break;} /* look for symbol until the end of the symbol list */
            if(symblist[i])break; /* if found */
        }

        if(!isCaseInSensitive && ((*str1)-(*str2)))break;
        else if(((toupper((int)(*str1)))-(toupper((int)(*str2)))))break;

        str1++;str2++;
    }
    if(symblist)return symblist[i] ? 0 : j ;
	return j ;

}


/*  performs structured search
    char * stringStructuredSearch(char * src,       source
        char * begin,                               the stringlist for structure beginning tags
        char * end,                                 the stringlist for structure endinning tags
        int * pst                                   points to the beginning of the found block
        int * pfn                                   points to the end of the found block
        returns the beginning of the found block

		example:
		{
			begin
				start
					  what I am looking for
				finish
			end
		}
		stringStructuredSearch (src,"{\0begin\0start\0\0","}\0end\0finish\0\0",&st,&fn);
        */

char * vstrStructuredSearch(char * src,char * begin,char * end,int * pst,int * pfn)
{
    char * st,*fn,*bg,*ed,*fnp;

    st=src;fn=fnp=st+strlen(st);
    bg=begin;ed=end;
     while(*bg || *ed){
        if(*bg){
            st=strstr(st,bg);if(!st)return 0;
            if(st>=fnp)return 0;
            st+=strlen(bg);
        }
        if(*ed){
            fn=strstr(st,ed);if(!fn)return 0;
        }
        bg+=strlen(bg)+1;
        ed+=strlen(ed)+1;
        fnp=fn;
    }
    if(pst)*pst=st-src;
    if(pfn)*pfn=fn-src;

    return st;
}


/*  skips several words
    char * stringSkipWords(char * src,       source
    int num                             how many words to skip
    char * separators ,                 separators between words
    returns the position of the found word
*/
char * vstrSkipWords(char * src, int num,char * separators)
{
  int i;
  
  if(!separators)separators=vSTRBLANK;
  while(*src && strchr(separators,*src))src++; /* pass spaces */
  for(i=0;*src && i<num;i++){
    while(*src && !strchr(separators,*src))src++; /* pass 1st non spaces */
    while(*src && strchr(separators,*src))src++; /* pass spaces */
  }
  if(!(*src))return 0;
  return src;
}

/*
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
_/                                          _/
_/  STRING FIND+REPLACE FUNCTIONS           _/
_/                                          _/
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/


/*
    find symbol from given symbol set and replaces them
    char * stringFindReplaceSymbols(char * src,     source
            char * dst,                             destination where to copy matches
            char * find,                            symbols to find
            char * replacement,                     symbols will be replaced by this string
            int maxTags,                            how many times perform the search&replace
            int isMatch,                            1 if match is required 0 otherweise
            int isSkipMultiple)                     1 skip multiple repetitions0 otherwise
    returns the last treated position to continue from
*/
char * vstrFindReplaceSymbols(char * src,char * dst,int sizeDest,char * find, char * replacement,int maxTags,int isMatch,int isSkipMultiple)
{
    int i,k;        /* i= current symbol to read from src, k at destination */
    int sc;         /* index for searching markUp tags */
    int replaced=0; /* index for replacement string */
    int howmany=0;
    int ir; /* index for replacement string */

    if(!src || !src[0])return src;
    if(!maxTags)maxTags=vMAXINUM; /* default maxtags */
	if(!sizeDest)sizeDest=vMAXINUM;

    for(i=0,k=0;k<sizeDest-2 && src[i] && howmany<maxTags;i++){
        /* scan if the current position is in find string */
        for(sc=0;src[i] && find[sc] && src[i]!=find[sc];sc++);
        if((find[sc] && isMatch) || (!find[sc] && !isMatch)){if(replaced && isSkipMultiple)continue;
            if(!replacement){dst[k]=0;k++;}
            else for(ir=0;replacement[ir];ir++){dst[k]=replacement[ir];k++;}
            replaced=1;howmany++;
        }else {dst[k]=src[i];k++;replaced=0;}

    }dst[k]=0;dst[k+1]=0;

    return src+i;
}
/*
    find strings from given stringlist and replaces them
    char * stringFindReplaceStrings(char * src,     source
            char * dst,                             destination where to copy matches
            char * find,                            on of the strings in this stringlist will be found
            char * replacement,                     found will be replaced by this string
            int maxTags,                            how many times perform the search&replace
    returns the last treated position to continue from
*/
char * vstrFindReplaceStrings(char * src,char * dst,int sizeDest,char * find, char * replacement,int maxTags,int isCaseInSensitive)
{
    int i,k;        /* i= current symbol to read from src, k at destination */
    int howmany=0,fnum;
    int ir,pos; /* index for replacement string */
    char * rpl;

    if(!src || !src[0])return src;
    if(!maxTags)maxTags=vMAXINUM; /* default maxtags */
	if(!sizeDest)sizeDest=vMAXINUM;

    for(i=0,k=0;k<sizeDest-1 && src[i] && howmany<maxTags;i++){

        /* scan if the current position corresponds to the startTag */
        pos=vstrIsOnList(src+i,find,&fnum,isCaseInSensitive);
        if(pos!=-1){
            if(!replacement){dst[k]=0;k++;}
            else { /* find the replacer string fnum and replace */
                rpl=vstrNextList(replacement,fnum);if(!rpl)rpl=replacement;
                for(ir=0;rpl[ir];ir++){dst[k]=rpl[ir];k++;}
            }
            howmany++;i+=pos-1;
        }else {dst[k]=src[i];k++;}

    }dst[k]=0;

    return src+i;
}

/*  cleans markup in the text
    src -source string
    dst -destination string (can be the same as source)
    startTag - start of the tag to look for (is a list ending with two zero symbols )
    endTag - end of the tag to look for (is a list ending with two zero symbols )
    replacement - whatever is found between tags will be replaced by this
    maxTags - maximum number of tag pairs to process
    inside - = 1 if tags internal content is necessary to leave untouched 0 otherwise
    returns the last treated position to continue from
 */
char * vstrClnMarkup(char * src,char * dst,int sizeDest,char * startTag, char * endTag,char * replacement,int maxTags,int inside,int isCaseInSensitive)
{
    int i,k,pos;        /* i= current symbol to read from src, k at destination */
    int notreplaced=0;  /* index for replacement string */
    int inMark=0;   /* =1 inside of the markup tags; otherwise 1 */
    int iEnd=-1,iStart=-1;      /* where to change the mode */
    int howmany=0,isin=0;
    int ir; /* index for replacement string */

    if(!maxTags)maxTags=vMAXINUM; /* default maxtags */
	if(!sizeDest)sizeDest=vMAXINUM;

    for(i=0,k=0;k<sizeDest-2 && src[i] && howmany<maxTags;i++){
        /* scan if the current position corresponds to the startTag list */
        pos=vstrIsOnList(src+i,startTag,0,isCaseInSensitive);
        if(pos!=-1 && (!isin)){iStart=i;if(inside)iStart+=pos-1;if(iStart<0)iStart=0;}  /* if inside is asked through away the markup tag recognizers */

        /* scan if the current position corresponds to the endTag */
        pos=vstrIsOnList(src+i,endTag,0,isCaseInSensitive);
        if(pos>0 && isin){iEnd=i;isin=0;if(!inside)iEnd+=pos-1;} /* if inside is asked through away the markup tag recognizers */

        /* determine if current position is inside of the Markup tag or outside */
        if(i==iStart){inMark=1;isin=1;notreplaced=0;if(startTag[0])continue;}
        if(i==iEnd){inMark=0;notreplaced=0;howmany++;continue;}

        /* copy inside of markup if inside is asked or copy outside of markup if non-inside */
        if(inside==inMark){dst[k]=src[i];k++;}
        else if(!notreplaced){notreplaced=1;/* copy the replacement to the found position */
            if(!replacement){dst[k]=0;k++;}
            else for(ir=0;replacement[ir];ir++){dst[k]=replacement[ir];k++;}
        }

    }dst[k]=0;dst[k+1]=0;

    return src+i;
}

/*
  cleans start and end of string from the symbols of given set
  char * stringCleanEnds(char * src,            source  string
                char * dst,                     destination where to copy
                char * find,                    look for these symbols
                int isMatch)                    1 if match is requred 0 otherwise
    returns the length of resulting string
*/
int vstrCleanEnds(char * src,char * dst,int sizeDest,char * find,int isMatch)
{
    int i,k=0,sc;       /* i= current symbol to read from src, k at destination */

    if(!src || !src[0])return 0;
	if(!sizeDest)sizeDest=vMAXINUM;

    for(i=0;src[i];i++){
        /* scan if the current position is in find string */
        for(sc=0;src[i] && find[sc] && src[i]!=find[sc];sc++);
        if((find[sc] && !isMatch) || (!find[sc] && isMatch))break;
    }
    for(;k<sizeDest-2 && src[i];i++,k++){dst[k]=src[i];}dst[k]=0;

    for(;k>0;k--){
        for(sc=0;dst[k] && find[sc] && dst[k]!=find[sc];sc++);
        if((find[sc] && !isMatch) || (!find[sc] && isMatch))break;
    }

    dst[k+1]=0;
    return k;
}

/*  returns string from array of strings separated by some markup
    char * stringGetValueFromArray(char * src,      source
        char * dst,                                 destionation
        int nextStp,                                how many separators to pass
        char * nextSepar)                             separator stringlist>

 so i should do
 > while (cp = vstrSeparateSubstringFromString(untreat,bfr,0, i, ...)
 > { ...
 > i++
 > }
 >
 If you want to go over tokens , yes !
 Another way you could do is :
 for (cp=untreat; (cp = vstrSeparateSubstringFromString(cp,bfr,0, 0, ...))
 ; )
 In this case you move your pointer (cp) starting from untreat and look for
 zero-th token every time,
 
*/

char * vstrSeparateSubstringFromString(char * src,char * dst,int sizeDest,int nextStp,char * nextSepar, int isCaseInSensitive)
{
    int iNxt=0,i,k,pos;

    if(!src || !src[0])return 0;
    if(!sizeDest)sizeDest=vMAXINUM;
	dst[0]=0;

    for(i=0,k=0;k<sizeDest-2 && iNxt<=nextStp && src[i];i++){
       pos=vstrIsOnList(src+i,nextSepar,0,isCaseInSensitive);
       if(pos!=-1){i+=pos-1;iNxt++;continue;}  /* if current string is in the list */
       /* if we have finally reached what we needed */
       if(iNxt==nextStp){dst[k]=src[i];k++;}
       if(iNxt>nextStp)break;
    }
    dst[k]=0;dst[k+1]=0;
    return src+i;
}


/*  copies two strings str1,str2 until a symbol from the symblist string
	is reached.
	Returns the length of the copied buffer */
int vstrCopyUntilSymbol(char * strsrc,char * strdst,int sizeDest,char * symblist)
{
    int   i,j;
    char *   str1=strdst;

	if(!sizeDest)sizeDest=vMAXINUM;

    for(j=0;j<sizeDest && (*strsrc);j++){ /* untill the end of the strings */

        if(symblist){
            for(i=0;symblist[i];i++){if(symblist[i]==(*strsrc))break;} /*look for symbol until the end of the symbol list */
            if(symblist[i])break; /* if found */
        }

        *str1=*strsrc;
        str1++;strsrc++;
    }
    *str1=0;
    return j;
}


/*
	modifies the string <TextStr>[MaxLen] according to hungarian
	notation if <IsName>==vTRUE and CaseType==0.  Otherwise changes the case
	( CaseType = 1 -> Translate To UperCase CaseType = 2 -> Translate To LowerCase)
	and removes the internal and external blanks.
	returns the length of the new string

*/
int vstrHungarianText(char * TextStr,char * TextDest,int sizeDest,int CaseType,int IsName,int IsRemIntBlanks)
{
    char chR,chW;
    int iScan,iWrite;
    int SkipBlank=1;
    int UpCasechar=1;

        if(!TextStr)return vINVIVAL; /* check for error */
        if(!sizeDest)sizeDest=vMAXINUM;/* Correct MaxLen */

        for(iScan=0,iWrite=0;iWrite<sizeDest && TextStr[iScan]!=0;iScan++){   /* Scan the original string */

            chR=TextStr[iScan]; /* read the next character */
            chW=0;  /* remember that do not transfer anything as default */

            if(strchr(vSTRBLANK,chR)){   /* if blank character occured */
                if(!SkipBlank){ /* if blank must not be skipped */
                    chW=vSTRBLANK[0];    /* remember to write SPC character */
                    SkipBlank=1;    /* One blank character already written and skip other */
                }
            }
            else    {   /* if Not blank character */
                chW=chR;    /* remember to transfer not blank character */
                if(!IsRemIntBlanks)SkipBlank=0; /* if internal blanks mest be reminded set that do not skip next blank */
            }

            if(CaseType==vSTRCASEUP)chW=toupper((int)chW);  /* revert to that case is need */
            if(CaseType==vSTRCASELO)chW=tolower((int)chW);
            if(IsName && UpCasechar)chW=toupper((int)chW);  /* if Upper Case is need in the names */

            if(chW!=0)TextDest[iWrite++]=chW;    /* write transferring character */

            if(isalpha((int)chW))UpCasechar=0;  /* if ascii character do not uppercase the next character in the names */

            else UpCasechar=1; /* else the next must be uppercased in the names */
        }

        if(iWrite>0)if(strchr(vSTRBLANK,TextDest[iWrite-1]))iWrite--;     /* if last character is blank */
        TextDest[iWrite]=0;  /* finalizing the string */

	return iWrite;  /* return the length of the string */
}

int vstrCase(char * TextStr,char * TextDest,int sizeDest,int CaseType)
{
	int i;
	if(!sizeDest)sizeDest=vMAXINUM;/* Correct MaxLen */

	for(i=0;TextStr[i] && i<sizeDest;i++){
		if(CaseType==vSTRCASELO)TextDest[i]=tolower((int)TextStr[i]);
		else if(CaseType==vSTRCASEUP)TextDest[i]=toupper((int)TextStr[i]);
		else TextDest[i]=TextStr[i];
	}
	TextDest[i]=0;
	return i;
}

/*
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
_/                                          _/
_/  EXTENDED PRINTF AND SSCANF FUNCTIONS    _/
_/                                          _/
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/

int _vstrScanAfter(char * wb,char * scnFormat,char * KeyFormat,char * pVal,char * find,char * scanto)
{
    char *   pFormat;

    if(find)pFormat=strpbrk(KeyFormat,find);
    else pFormat=KeyFormat;
    if(pFormat){
        vstrCopyUntilSymbol((pFormat+1),wb,0,scanto);
    }else return vFALSE;

    if(scnFormat)if(!sscanf(wb,scnFormat,pVal))return vFALSE;

    return vTRUE;
}


int vstrExtendedSScanf(char * textScan,char * KeyFormat,void * pVal)
{
    int   iChoice,is,isScanned,iVal,iTmp,isNON;
    char   wb[vSTRMAXNAME],bFormat[vSTRMAXNAME];
    double   rVal,rTmp;
    char *   toEndSym="=><;",* pItem,* curT,* pFormat;

    vstrCopyUntilSymbol((KeyFormat+1),bFormat,0,toEndSym);pFormat=bFormat;
    /* is it undefined */
    isNON= vstrSearchSubstring(textScan,"default\0undefined\0\0",1,0,0) ? vTRUE : vFALSE ;


    switch(KeyFormat[0]){
    case 'i':iVal=*((int *)pVal);if(isNON){*((int *)pVal)=vINVIVAL;break;}
        if(!(*pFormat))pFormat="%i";
        isScanned=sscanf(textScan,pFormat,&iVal);

        /* if couldn't scan try default value */
        if(isScanned!=1 && _vstrScanAfter(wb,pFormat,KeyFormat,(char *)&iTmp,"=",toEndSym))iVal=iTmp;
        if(_vstrScanAfter(wb,pFormat,KeyFormat,(char *)&iTmp,">",toEndSym) && (iVal<iTmp))iVal=iTmp; /* minimum */
        if(_vstrScanAfter(wb,pFormat,KeyFormat,(char *)&iTmp,"<",toEndSym) && (iVal>iTmp))iVal=iTmp; /* maximum */
        *((int *)pVal)=iVal;
        break;
    case 'r':rVal=*((double *)pVal);if(isNON){*((double *)pVal)=vINVRVAL;break;}
        if(!(*pFormat))pFormat="%lf";
        isScanned=sscanf(textScan,pFormat,&rVal);

        /* if couldn't scan try default value */
        if(isScanned!=1 && _vstrScanAfter(wb,pFormat,KeyFormat,(char *)&rTmp,"=",toEndSym))rVal=rTmp;
        if(_vstrScanAfter(wb,pFormat,KeyFormat,(char *)&rTmp,">",toEndSym) && (rVal<rTmp))rVal=rTmp; /* minimum */
        if(_vstrScanAfter(wb,pFormat,KeyFormat,(char *)&rTmp,"<",toEndSym) && (rVal>rTmp))rVal=rTmp; /* maximum */
        *((double *)pVal)=rVal;
        break;
    case 't':
        if(!(*pFormat))pFormat="%s";
        if(textScan && (*textScan))strcpy((char * )pVal,textScan);
        else _vstrScanAfter((char * )pVal,0,KeyFormat,0,"=",toEndSym);
        break;
    case 'l':
        pItem=strpbrk(KeyFormat,"^;"); /* for logica variables if they aren't full defined define them */
        if(!pItem || (*pItem)!='^')KeyFormat="l%=0^FALSE^TRUE^OFF^ON^NO^YES^CANCEL^OK^0^1;";
    case 'e':iVal=*((int *)pVal);
        pItem=KeyFormat;
        isScanned=vFALSE;
		is=0;
        for(iChoice=0;;iChoice++){ /* scan all choices */
            pItem=strpbrk(pItem,"^;");if((*pItem)!='^')break; /* no more items */
            pItem++;

            is=vstrCompareUntilSymbol(textScan,pItem,"=^;",0);
            if(strchr(";^=",pItem[is])){isScanned=vTRUE;break;} /*if we have found it */
        }
        /* if couldn't find try default value otherwise get the value of the choice or it's number  */
        if(!isScanned && _vstrScanAfter(wb,"%i",KeyFormat,(char *)&iTmp,"=",";^"))iVal=iTmp;
        else if(pItem[is]=='=' && _vstrScanAfter(wb,"%i",pItem+is,(char *)&iTmp,0,";^"))iVal=iTmp;
        else iVal=iChoice;

        if(KeyFormat[0]=='l')iVal=(iVal%2); /* for logical variables only */
        *((int *)pVal)=iVal;
        break;
    case 'f':iVal=0;
        for(curT=textScan;curT && *curT;){ /* scan through all the pieces of the input string */
            pItem=KeyFormat;
            isScanned=vFALSE;
			is=0;
            for(iChoice=0;;iChoice++){ /* scan all choices */
                pItem=strpbrk(pItem,"|;");if((*pItem)!='|')break; /* no more items */
                pItem++;

                is=vstrCompareUntilSymbol(curT,pItem,"=|;",0);
                if(strchr(";|=",pItem[is])){isScanned=vTRUE;break;} /*if we have found it */
            }
            /* if couldn't find try default value otherwise get the value of the choice or it's number  */
            if(isScanned){
                if(pItem[is]=='=' && _vstrScanAfter(wb,"%lx",pItem+is,(char *)&iTmp,0,"|;"))iVal|=iTmp;
                else iVal|=(1<<iChoice);
            }
            curT=strpbrk(curT+is,"|;"); /* next piece */
            if(curT){curT++;
				while(*curT!=0 && (  (*curT==' ') || (*curT=='\t') || (*curT=='\r') || (*curT=='\r') )  )curT++;
			}
        }
        if(!iVal && _vstrScanAfter(wb,"%lx",KeyFormat,(char *)&iTmp,"=","|;"))iVal=iTmp;

        *((int *)pVal)=iVal;
        break;

    default:
        strcpy((char * )pVal,textScan);
        break;
    }


  return vTRUE;

}
/*
// printf to <RetBuf> string using extended KeyFormat to the variable pVal
*/
char *  vstrExtendedSPrintf(char * RetBuf,char * KeyFormat,void * pVal)
{
    int   iChoice,iVal,iTmp;
    char   wb[vSTRMAXNAME],bFormat[vSTRMAXNAME];
    double   rVal;
    char *   toEndSym="=><;",* pItem,* curT,* pFormat;

    vstrCopyUntilSymbol((KeyFormat+1),bFormat,0,toEndSym);pFormat=bFormat;

    switch(KeyFormat[0]){

    case 'i':iVal=*((int *)pVal);
        if(!(*pFormat))pFormat="%i";
        sprintf(RetBuf,pFormat,iVal);
        break;

	case 'r':rVal=*((double *)pVal);
        if(!(*pFormat))pFormat="%lf";
        sprintf(RetBuf,pFormat,rVal);
        break;

	case 't':
        if(!(*pFormat))pFormat="%s";
        sprintf(RetBuf,pFormat,(char * )pVal);
		break;

	case 'l':
        pItem=strpbrk(KeyFormat,"^;"); /* for logical variables if they aren't full defined define them */
        if(!pItem || (*pItem)!='^')KeyFormat="l%=0^FALSE^TRUE^OFF^ON^NO^YES^CANCEL^OK^0^1;";

	case 'e':iVal=*((int *)pVal);
        pItem=KeyFormat;
        for(iChoice=0;;iChoice++){ /* scan all choices */
            pItem=strpbrk(pItem,"^;");if((*pItem)!='^')break; /* no more items */
            pItem++;

            pFormat=strpbrk(pItem,"=^;");
            if(pFormat && pFormat[0]=='=' && _vstrScanAfter(wb,"%lx",pItem,(char *)&iTmp,0,"^;"));
            else iTmp=(iChoice);
            if(iTmp==iVal){vstrCopyUntilSymbol(pItem,RetBuf,0,"^;");break;}
        }
        if((*(pItem-1))!='^')sprintf(RetBuf,"%i",iVal);/* no match - print as a number */
        break;

    case 'f':iVal=*((int *)pVal);
        pItem=KeyFormat;curT=RetBuf;
        for(iChoice=0;;iChoice++){ /* scan all choices */
            pItem=strpbrk(pItem,"|;");if((*pItem)!='|')break; /* no more items */
            pItem++;

            pFormat=strpbrk(pItem,"=|;");
            if(pFormat && pFormat[0]=='=' && _vstrScanAfter(wb,"%lx",pFormat,(char *)&iTmp,0,"|;"));
            else iTmp=(1<<iChoice);
            if(iVal&iTmp){
                curT+=vstrCopyUntilSymbol(pItem,curT,0,"|;");
                *curT='|';curT++;
            }
        }
        if(*(curT-1)=='|')curT--;
        *curT=0;
        break;

    default:
        strcpy(RetBuf,(char * )pVal);
        break;
    }


  return RetBuf;
}

/*
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
_/                                          _/
_/  DYNAMIC SPRINTF FUNCTIONS               _/
_/                                          _/
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/

int vstrPrintfSizeOfArgList(char * formatDescription,va_list marker)
{
    int len=0,lenAdd,width,stopCountingWidth,isLong;
    char *fmt,*typeOf;
	void *pVal; 

	#define	textSIZEOFINT   	16	 /* should be width or strlen(sprintf(intnumber)) */
	#define	textSIZEOFREAL		32	 /* see above ... */
	#define	textSIZEOFADDRESS	16 /* should be (sizeof(void*)*2) each byte is two hexadecimal characters */
	#define textSIZEOFCHAR		1	/* 1 */
	#define textSIZEOFPERCENT	1	/* 1 */
	#define textSIZEOFDEFAULT	1	/* 1 */


	if(!formatDescription)return vNULL;

	for(fmt=formatDescription;*fmt;fmt++){
		if(*fmt=='%'){ /* scan for format specification */

			/* see the width of the variable, example: %10.3lf */
			width=0; stopCountingWidth=0;
/*			for(typeOf=fmt+1; *typeOf && *typeOf!='%' && !isalpha((int)(*typeOf)); typeOf++ ){*/
			for(typeOf=fmt+1; *typeOf && strchr("01234567890.+-",*typeOf) ; typeOf++ ){ 

				/* stop counting if dot is enountered after some numbers */
				if(*typeOf=='.')
					stopCountingWidth=1;
				/* accumulate the width */
				if(isdigit((int)(*typeOf)) && (!stopCountingWidth))
					width=width*10+(*typeOf)-'0';
			}

			if(*typeOf=='l'){isLong=1;typeOf++;} /* long format */
			else isLong=0;
			switch(tolower((int)*typeOf)){
				case 'i': case 'd': case 'x':
					if(isLong)va_arg(marker, int);
					else va_arg(marker, long int );
					lenAdd=textSIZEOFINT;
					break;
				case 'f': case 'g': case 'e':
					if(isLong)va_arg(marker, double );
					else va_arg(marker, double );
					lenAdd=textSIZEOFREAL;
					break;
				case 's':
					pVal=va_arg(marker, char * );
					lenAdd=strlen(pVal);
					break;
				case 'p':
					pVal=va_arg(marker, void * );
					lenAdd=textSIZEOFADDRESS;
					break;
				case 'c':
					va_arg(marker, int );
					lenAdd=textSIZEOFCHAR;
					break;
				case '%': /* %% or alike */
					lenAdd=textSIZEOFPERCENT;
					break;
				 default: 
				 	typeOf--;
					lenAdd=textSIZEOFDEFAULT;
				 	break;
			}
			if(width)lenAdd=vMax(lenAdd,width);
			fmt=typeOf;
			len+=lenAdd;
		}
		else
	        len++;
    }

    return len+1;
}

char * vstrDynaPrintfArgList(char * formatDescription, va_list ap)
{
  int	len;
  char * newBuf=0;
  va_list myap; 
  
#if ! defined(ALPHA)
  /* somehow va_copy is not found on my alpha, but this is a bug not to use it */
  va_copy (myap, ap) ;
  len = vstrPrintfSizeOfArgList(formatDescription, myap) ;
  va_end (myap) ;
#else
  messcrash ("vstr.c:va_copy not declare for ALPHA in vstr.c line 847") ;
  len = vstrPrintfSizeOfArgList(formatDescription, myap) ;
#endif  

  newBuf = messalloc(len+1) ;
  if (newBuf)
    {
      vsprintf(newBuf,formatDescription,ap) ;
      if (len < strlen(newBuf))
	messcrash ("vstr:sprintf size estimation failed \n estimated %d required=%d\n%s",len,strlen(newBuf),formatDescription);
    }
	
  return newBuf;
}

char * vstrDynaPrintf(char * formatDescription , ...)
{
	char * res;
	vCallVarg(res,vstrDynaPrintfArgList,formatDescription);
	return res;
}



/*
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
_/                                          _/
_/  Extensible string buffer Functions      _/
_/                                          _/
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/

int vstrPrintfArgList(vSTR * str,char * formatDescription,va_list marker)
{
    char * cont;
	int	ret=0,len;

    if(!str)return 0;

	if( (cont=vstrDynaPrintfArgList(formatDescription,marker)) &&
		(len=strlen(cont)) ){

	    if(!(str->mode&vSTRDISCARD)){ret=vmexAppend(&str->mem,cont,len+1);str->mem.curPos--;}
	    if(str->funcCallback)str->funcCallback(str,str->paramCallback,cont);
	    if(str->mode&vSTRECHO) printf ("%s",cont);
		if(str->outFile) vtxtFileAppend(str->outFile,cont);
		if(str->hFile)vtxtFileWrite(str->hFile,cont,len);
	}
	messfree(cont);

	return ret;
}

vIN int vstrPrintf(vSTR * str,char * formatDescription,...)
{
	int ret;
	vCallParaVarg(ret,vstrPrintfArgList,str,formatDescription);
	/*if((!ret) && (str)){ret=vstrLen(str);} */
	return ret;
}

int vstrPrintfWrapped(vSTR * str,char * src,char * separ,int charrayLen,int caseChar, int maxnum)
{
	int i,k,iy=0,lastspace=0,lastpos=0,anySpace=0;
	char dst[vTEXTMAXLINE],ch;

	/*
	lastspace=charrayLen;
	lastpos=charrayLen-1;
	*/
	lastspace=0;
	lastpos=0;

	if(!maxnum)maxnum=vMAXINUM;

	for(k=0,i=0;i<maxnum && (ch=src[i])!=0;i++){

		if((separ && strchr(separ,ch)) || strchr(vSTRENDLINE,ch) ){
			lastspace = k; lastpos = i;anySpace=1;
		}
		/* put endline and printf the piece */
		if(strchr(vSTRENDLINE,ch) || k>=charrayLen){
			if(anySpace)dst[lastspace]=0;
			else dst[k]=0;

			vstrPrintf(str,"%s\n",dst);
			
			if(!anySpace){ /* there were no places to break inside of the line */
				i--;lastpos=i;
			}else {
				i=lastpos;
			}
			
			lastspace=k=0; /* and reset the destination */
			iy++;
			anySpace=0;
			continue;
		}
        if(caseChar==vSTRCASEUP)ch=toupper((int)src[i]);
        else if(caseChar==vSTRCASELO)ch=tolower((int)src[i]);
        else ch=src[i];

		dst[k]=ch;k++;
    }
	if(k){
		dst[k]=0;
		vstrPrintf(str,"%s",dst);
	}
	return iy;
}




/*
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
_/                                          _/
_/  Regexp Functions                        _/
_/                                          _/
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/
/* in the string <src> looks for the regular expression  <rexpLine> puts the result into
	the array of pmatch[cntmatch> . Returns the first occurence or -1 if there is none ; */
	/*
int vstrSearchRegexp(char * src, char * rexpLine, int cntmatch, int * pmatch,int isCaseInSentitive)
{
	regex_t rx;
	int mymatch;

 	if(!cntmatch || !pmatch){cntmatch=1;pmatch=&mymatch;}

	if( regcomp (&rx, (char * )rexpLine, REG_EXTENDED | (isCaseInSentitive ? REG_ICASE : 0 ) )!=0)
	  return vFALSE;

	if(regexec (&rx, src, cntmatch, pmatch, 0) == 0){
	  cntFnd++;
	  foundIt=1;
	}
	return *pmatch;
}
*/
