#include    "../wac/ac.h"
#include    <vstd.h>
#include    <utils.h>
#include "dict.h"

/* seems to be used to recover the GI of th YK clones
   code by Vahan
*/

typedef struct tagINFO {
  char giBuf[128],acBuf[128];	
  int isUsed;
} INFO;


static void objOut(const char *nm,INFO * info)
{
  char * vers = "1" ;
  
  if ((vers = strrchr(info->acBuf,'.')))
    vers++;
  
  printf(	"Sequence \"%s\"\n"
		"Database Genbank \"%s\" \"%s\" %s\n\n"
		, nm, nm, info->acBuf, vers
		);
  
}

static int objFind (AC_OBJ obj, Array infoArray, DICT *dict)
{
  int iFnd = -1,pos;
  INFO * info;
  char nm[1024],buf[1024],*ptr;
  
  
  /* try to find it as it is */
  strcpy(nm,ac_name(obj));
  dictFind (dict, nm, &iFnd) ;
  
  
  /* try to find it after replacing stage infoo to Kohara letter k */
  if(iFnd==-1 && nm[0]=='y')
    {
      vstrFindReplaceStrings(nm,nm,0,"y1L"_0"y2L"_0"y3L"_0"y4L"_0"yd"_0"yc"_0,"yk"_00,0,1);
      dictFind (dict, nm, &iFnd) ;
      
      printf("// YK removed the stage infoo: %s\n",nm);
    }
  
  /* look for zerofied things */
  if (iFnd==-1 && !strncasecmp(nm,"yk",2) && (ptr=strpbrk(nm+2,"abcdefgh")) && ptr[2]=='.' )
    {
      pos=ptr-nm+1;
      strncpy(buf,nm,pos);
      strcpy(buf+pos+1,nm+pos);
      buf[pos]='0';
      dictFind (dict, nm, &iFnd) ;
      printf("// YK added zero in numbering: %s\n",buf);
  }
  
  /* look for REVERSE or FORWARD markings */
  if(iFnd==-1 && !strncasecmp(nm,"CE",2))
    {
      for(ptr=nm;isalpha((int)(*ptr));ptr++); /* skip the beginning of the name */
      for(;isdigit((int)(*ptr));ptr++);
      
      if(ace_upper((*ptr))=='F' || ace_upper((*ptr))=='R')
	*ptr=0;
      
      dictFind (dict, nm, &iFnd) ;
      printf("// CE removing Forward or Reverse in naming: %s\n",nm);
    }
  
  /* remove .p or adding .5 */
  if(iFnd==-1 && !strncasecmp(nm,"cm",2))
    {
      if((ptr=strchr(nm,'.')))*ptr=0;
      strcat(nm,".5");
      
       dictFind (dict, nm, &iFnd) ;
       printf("// CM cleaning the naming : %s\n",nm);
    }
  
  if(iFnd==-1)
    {
      printf("// failed to locate the accession infoo for %s\n\n", ac_name(obj));
      return 1;
    }
  
  
  info = arrayp (infoArray, iFnd, INFO);
  objOut (ac_name(obj), info) ;
  info->isUsed++;
  
  return 1;
}

int main(int argc, const char * argv[])
{   
  const char *par ;
  char *ptr;  
  AC_DB  db;
  INFO info;
  DICT *dict = dictCreate (10000) ;
  Array infoArray = arrayCreate (10000, INFO) ;
  char buf[16*1024],clBuf[1024];
  
  freeinit () ;
  
  if(argc<2)
    {
      printf ("\nUsage: cl2acc parameters < inputfile > outputfile\n"
	      "Input File has to be preliminary fetched from GenBank :\n"
	      "\t1. submit a query to nucleotide search :\n"
	      "\t     (EST[Keyword] AND \"caenorhabditis elegans\"[organism]) \n"
	      "\t   or \n"
	      "\t     \"caenorhabditis elegans\"[organism] gbdiv_est \n"
	      "\t2. Save the GI list \n"
	      "\t3. get the summaries for all GI's \n"
	      "\t     idfetch -dn -G GI_list_file_name_here -n -t 7 > result_file_here\n"
	      "\n"
	      "Parameters:\n"
	      "\t[-db host:port] /* to fetch accessions only for ests from this DB */ \n"
	      "\t  /* if this parameter is not specified all accessions \n"
	      "\t     from the input file will be dumped as ace file  */\n"
	      "\t[-q \"est_wildcard\"]  /* to fetch accessions only for these ests from DB\n"
	      "\t    /* deafult is '*' */ \n"
	      "\n"
	      "\n"
	      );
      exit (1) ;
    }
  
  while(fgets(buf,sizeof(buf)-2,stdin))
    {
      int ii ;
      INFO *infop ;

      memset (&info, 0, sizeof (INFO)) ;
      /* gi and accession */
      vstrSeparateSubstringFromString(buf,info.giBuf,sizeof(info.giBuf)-1,1,"|"_00,0);
      vstrSeparateSubstringFromString(buf,info.acBuf,sizeof(info.acBuf)-1,3,"|"_00,0);
      
      /* clone name and suffix 3' or 5' */
      if( !(ptr=strstr(buf,"clone")) || !(ptr=vstrSkipWords(ptr,1,0)) )continue;
      vstrCopyUntilSymbol(ptr,clBuf,sizeof(clBuf)-1," ,\n");
      if((ptr=vstrSkipWords(ptr,1,0)) && (ptr[0]=='3' || ptr[0]=='5') )
	sprintf(clBuf+strlen(clBuf),".%c",ptr[0]);
      
      /* put it to dictionary */
      dictAdd (dict, clBuf, &ii) ;
      infop = arrayp (infoArray, ii, INFO) ;
      *infop = info ;
  }
  
  if (getCmdLineOption (&argc, argv, "-db", &par))
    {
      const char *error = 0 ;
      const char *qq = "*" ;

      getCmdLineOption (&argc, argv, "-q", &qq) ;
      if((db = ac_open_db ( par , &error)))
	{
	  AC_ITER iter = ac_query_iter (db, 1, messprintf ("find est %s" , qq), 0, 0) ;
	  AC_OBJ obj ;
	  
	  while ((obj = ac_next_obj (iter)))
	    objFind (obj, infoArray, dict) ;

	  ac_db_close (db);
	}
    }
  else
    {
      int i;
      
      for(i=1; i <= dictMax (dict); i++)
	{
	  objOut (dictName (dict, i), arrayp (infoArray, i, INFO)) ;
	}
      
    }
  dictDestroy (dict) ;
  arrayDestroy (infoArray) ;

  return 0;
}
