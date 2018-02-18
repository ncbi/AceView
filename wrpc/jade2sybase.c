/* $Id: jade2sybase.c,v 1.5 2015/08/20 00:12:59 mieg Exp $ */

#include "regular.h"
#include "dict.h"
#include "freeout.h"
#include "call.h"

#define OK 200
#define GOODBYE 201
#define DEBUG 301
#define COMMENT 302
#define REDIRECT 303
#define PARAMERROR 401
#define SYNTAXERROR 402
#define UNIMPLEMENTED 403
#define TIMEOUTERROR 501
#define COMMERROR 502
#define MEMERROR 503
#define FILEERROR 504

static DICT *classDict = 0, *tableDict = 0 ;
static Array tableDef = 0, classDef = 0 ;
static BOOL debug = FALSE ;
static char* scriptName = 0 ;
static char* codeBase = 0 ;
static Stack currentList = 0 ;

BOOL isInteractive ;

static FREEOPT sybaseMenu[] =
{ {8, "Jade2Sybase"},
  {'?',"? : List of commands"},
  {'h', "HELP"},
  {'q', "QUIT"},
  { 't',"NTABLE"},
  {'Q',"QUERY"},
  { 'f',"FIND"},
  {'l', "LIST"},
  {'w', "SHOW"}
} ;

/* needed to link */
void tStatus(void) {return ;} 
char *linkDate = "" ;

void writeStatus (int status,char* message) {
  fprintf(stdout,"%d %s\r\n",status,message);
  fflush(stdout);
}

void printHelp () {
  char* helpText[] = {
    "- Commands:",
    "-    FIND   LIST NTABLE",
    "-    QUIT   HELP",
    "- You need to know the SQL to use this",
    "- service.  Sorry.",
    "End of HELP info",
    NULL
  };
  char** m = helpText;
  while (*m != NULL)
    writeStatus(COMMENT,*m++);
}


static Stack callMyScript(char *text)
{ char *cp ;
  FILE *f = callCdScriptPipe(0, scriptName,text) ;
  int level = freesetpipe(f,0) ;

  Stack s = stackCreate (2000) ;
  while (freecard(level))
    { cp = freepos() ;
      if (cp && *cp)
	pushText(s, freepos()) ;
    }
  return s ;
}

static int  isJade = 0 ;

static void sybaseFind(char *classe) 
{
  int nn = 0 ;
  char *cp = classe ;
  while (*cp)  { *cp = ace_lower(*cp) ; cp++ ;}
  messdump("sybaseFind %s\n", classe) ;
  if (isJade < 2) isJade = 0 ;
  if (dictFind(classDict, messprintf("%s",classe), &nn))
    {
      Stack s = callMyScript (arr(classDef,nn,char*)) ;
      stackCursor(s,0) ; nn = -1; /* -1 for the 200 message */
      while((cp = stackNextText(s))) if (cp && *cp && *cp != '.') nn++ ;
      freeOutf("%d %d objects\r\n",OK, nn) ;
      fflush(stdout) ;
      stackDestroy(currentList) ;
      currentList = s ;
      messdump("200 %d objects\r\n",nn) ;

    }
  else if ( !strcmp(classe,"jade"))   /* cheat */
    { 
      isJade = 1 ;messdump("200 1 objects\r\n") ;
      writeStatus(OK,"1 objects") ;
    }
  else
    { /* freeOutf("302 Unknown class %s\r\n",classe) ; */
    freeOutf("200 0 objects\r\n") ;
    messdump( "302 Unknown class %s\r\n",classe) ;
    }
}

static void sybaseList (void)
{
  char *cp ;
  int nn = 0 ;

  messdump("list: \n") ;

  freeOut("200 Results follow, terminated by \".\"\r\n") ;  
  if (stackExists(currentList))
    {
      Stack s = currentList ;
      stackCursor(s,0) ; 
      if ((cp = stackNextText(s))) /* jump the 200 message */
	while((cp = stackNextText(s)))
	  { freeOutf("%s\r\n",cp) ;
	  if (cp && *cp && *cp != '.' && nn++ < 12) messdump("list: %s\n", cp) ;
	  }
    }
  else if (isJade== 1)
    { nn = 1 ; freeOutf("?Jade?default?\n") ; } 

  freeOutf("\n.\r\n") ;

  messdump("listed %d objects\n", nn) ;
}

static void sybaseShow (void)
{
  freeOut("200 Results follow, terminated by \".\"\r\n") ;  
  if (isJade)
    { freeOutf("?Jade?default?\t?tag?Display?\t?txt?Map?\t?txt?jade.maps.Vmap?\r\n\n") ;
    messdump("?Jade?default?\t?tag?Display?\t?txt?Map?\t?txt?jade.maps.Vmap?\r\n\n") ;
    isJade = 2 ;
    }
  freeOut(".\r\n") ;
}

static void classList(void)
{ int i ;
  freeOut("200 Results follow, terminated by \".\"\r\n") ;  
  for (i=1 ; i < dictMax(classDict) ; i++)
    freeOutf ("?Class?%s?\r\n",dictName(classDict,i)) ;

  freeOutf ("\n\n.\n") ;
}


static void ask(int level)
{
  KEY option = 0 ;
  Stack s = 0 ;
  char *cp, *cq=0, *cr=0 ;
  int i, nn ;

  freeOutf("200 Sybaseserver ready and waiting\r\n") ;
  fflush(stdout) ;
  isInteractive = FALSE ;
  while (option = 0, TRUE)
  if (freelevelselect (level, &option, sybaseMenu))
    {
      cp = freepos() ;
      if (option)
	messdump("option %c ::%s\n",(char)option, cp ? cp : "") ;
     switch (option)
	{
	case (KEY)(-1): case 'q':   /* exit */
	  writeStatus (GOODBYE,"A bientot") ;
	  return ;
	  
	case '?':
	  freeOut ("Jade2Sybase list of commands:\n") ;
	  for (i = 1 ; i <= sybaseMenu[0].key ; i++)
	    { freeOut (messprintf("%s\n",sybaseMenu[i].text)) ;
	    }
	  break ;   
	  
	case 'h':
	  printHelp() ;
	  break ;
	case 't':
	  messdump( "Call table") ;
	  cp = freeword() ;
	  if (!cp) break ;
	  cq = strnew(cp, 0) ;
	  cp = cq ; while (*cp) { *cp = ace_lower(*cp) ; cp++ ;}
	  freenext() ;
	  messdump( "Call table %s\n",cq) ;
	  cp = freepos() ;
	  if (cp && *cp)
	    cr = strnew(cp,0) ;
	  else cr = 0 ;
	  if (!strcasecmp(cq,"Classes"))
	    classList() ;
	  else if (dictFind(tableDict, cq, &nn))
	    {
	      int level1 = freesettext(arr(tableDef,nn,char*),cr) ;
	      char *myquery ;
	      
	      freespecial("%%/\"\\\n") ;
	      freecard(level1) ;
	      myquery=freepos() ;
		messdump( "Call table %s :: %s  \n", 
		      cq, myquery) ;
		/*freeOut("200 Results follow, terminated by \".\"\r\n") ;  
		 */
	      s = callMyScript(myquery) ;
	      stackCursor(s,0) ; nn = 0 ;
	      while ((cp = stackNextText(s)))
		if (cp && *cp) { freeOutf ("%s\r\n",cp) ; nn++ ;}
	      stackDestroy (s) ;
	      messdump("  :: exported %d lines\n",nn) ; 
	      freeOut("\n\n.\r\n") ;  
	    }
	  else
	    {
	      freeOut("200 Results follow, terminated by \".\"\r\n") ;  
	      freeOut("\n\n.\r\n") ;  
	      messdump("Unknown table %s\r\n", cq) ;
	    }
	  messfree(cq) ; messfree (cr) ;
	  break ;
	case 'Q':
	  cp = freeword() ; freenext() ;
	  /* fall thru:: treat as Query Find .. */
	case 'f':
	  cp = freeword() ;
	  if (!cp) break ;
	  cq = strnew(cp, 0) ;
	  sybaseFind(cq) ;
	  messfree(cq) ;
	  break ;
	case 'l':
	  sybaseList() ; /*active list */
	  break;
	case 'w':
	  sybaseShow() ; /*active list */
	  break;
	default:
	  break ;
      }
      fflush(stdout) ;
    }
}
  
static void getConfig (int level)
{
  char *cp=0, *cq =0, *cr, *scriptPartialName = 0 ;
  int i, nn ;
  int state = 0 ;
 
  while (freecard(level))
    {
      cp = freeword() ;
      if (!cp) continue ;

      if (!strcmp(cp,"CODE_BASE"))
	{
	  state = 0 ;
	  if ((cp = freeword()))
	    codeBase = strnew (cp,0) ;
	  continue ;
	}
	 
      if (!strcmp(cp,"EXTERNAL_CALL"))
	{
	  state = 0 ;
	  if ((cp = freeword()))
	    scriptPartialName = strnew (cp,0) ;
	  else
	    messcrash 
("No program name specified on the EXTERNAL_CALL line of the configuration file") ;
	  continue ;
	}
	 
	 
      else if (!strcmp(cp,"CLASS_LIST"))
	{ state = 1 ;
	  continue ;
	}

	 
      else if (!strcmp(cp,"TABLE_LIST"))
	{ state = 2 ;
	  continue ;
	}
   
      else if (!strcmp(cp,"IGNORE"))
	{ state = -1 ;
	  continue ;
	}
   

      cq = strnew(cp,0) ;
      freenext() ;
      cp = freepos() ;
      if (cp && *cp)
	{
	  cr = codeBase ? strnew(messprintf("%s/%s",codeBase,cp), 0) : strnew(cp,0) ;
	  cp = cq ; while (*cp) { *cp = ace_lower(*cp) ; cp++ ;}
	  switch (state)
	    {
	    case 1: 
	      dictAdd(classDict, cq, &nn) ;
	      array(classDef,nn,char*) = cr ;
	      break ;
	    case 2:
	      dictAdd(tableDict, cq, &nn) ;
	      array(tableDef,nn,char*) = cr ;
	      break ;
	    default:
	      break ;
	    }
	}
      messfree (cq) ;
    }
  if (debug)
    for (i=1 ; i <= dictMax(tableDict); i++)
      messdump("Table %s,  sql: %s\n",dictName(tableDict,i), arr(tableDef,i,char*)) ;
  if (codeBase)
    scriptName = strnew(messprintf("%s/%s",codeBase, scriptPartialName), 0) ;
  else
    scriptName = scriptPartialName ? strnew(scriptPartialName, 0) : 0 ;
}

int main(int argc, const char **args)
{
  int level ;
  FILE *f = 0 ;

  freeinit() ;
  if (argc<2)
    messcrash(
   "//Usage: jade2sybase <config.file> \nYou said %s %s %s %s\n",
   args[0], argc>1 ? args[1]:"",argc>2 ? args[2]:"",argc>3 ? args[3]:"") ;

  tableDef = arrayCreate (20,char*) ;
  classDef = arrayCreate (20,char*) ;
  classDict = dictCreate (50) ;
  tableDict = dictCreate (50) ;


  f = fopen(args[1],"r") ;
  if (!f) messcrash ("Cannot open the configuration file %s\n",args[1]) ;
  level = freesetfile(f,0) ;
  freespecial("/\"\\\n") ;
  getConfig(level) ;

  if (!scriptName || !filName(scriptName,0,"x"))
    messcrash("Cannot locate script %s indicated in the configuration file  %s\n",
	      scriptName, args[1]) ;
  freeinit() ;
  freeOutInit() ;
  freeOutSetFile(stdout) ;

  level = freesetfile(stdin,0) ;
  freespecial ("/\"\\\n") ;
  ask(level) ;
  return 0 ;
}
