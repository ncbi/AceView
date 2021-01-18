/* $Id: prefsubs.c,v 1.3 2017/05/24 16:33:26 mieg Exp $ */
/*  Last edited: Nov 19 15:53 1998 (fw) */

#include "acedb.h"
#include "lex.h"
#include "pick.h"
#include "systags.h"
#include "sysclass.h"
#include "pref.h"
#include "session.h"
#include "bs.h"
#include "query.h"


#include <ctype.h>

#define BOOLEAN 1
#define FLOAT 2
#define INT 3
#define COLOUR 4

#ifndef NON_GRAPHIC
#include "display.h"
#include "main.h"
#else

#include "colours.h"		/* for enum Colour {..} */

#endif

/*
il-27/01/97

 Examples of preferences file lines
  
_NAME -T       |set Name to TRUE
_NAME1 -F      |set Name1 to FALSE
_NAME2 -I 2    |set Integer Name2 to 2.
_NAME3 -R 5.5  |set Float Name3 to 5.5.
_NAME4 -C 1    |set colour NAME4 to WHITE (by number)
_NAME4 -C WHITE    |set colour NAME4 to WHITE (by name)


where "_" must be the first character on the line followed by the preferences name
then (tab or space delimited) a "-" the type and for int,floats and colours the value.
Any other lines will be ignored.
To use this later on in the program use prefValue(NAME) which will return the BOOL value.
prefInt(NAME) will return the value for NAME etc.
*/


void editPreferences(void);
static void prefRead(FILE *prefFile);

typedef struct{
  char name[32];
  int type;
  BOOL display;
  union value_tag{
    BOOL bval;
    int ival; /* used for int and colours */
    float fval;
  } value;
}PREFLIST;

struct configLocals {
  union value_tag2{
    BOOL bval;
    int ival; /* used for int and colours */
    float fval;
  } value[100];
} ;



Array prefList = 0;
#ifndef NON_GRAPHIC
static struct configLocals *cf;
Graph prefsGraph=0;


static void savePrefs(void) {
/* if acedb then save the file as the system defaults */
/* else the file should be under the user name */
  char filename[FIL_BUFFER_SIZE];
  FILE * prefFile=0;
  int i=0;
  PREFLIST *item;

  if(getenv("HOME")){
    strcpy(filename,getenv("HOME"));
    strcat(filename,"/.acedbrc");
    prefFile = filopen(filename,"","w");
  }
  if (!prefFile) return ;
  fprintf(prefFile,"This file contains the preferences for ACEDB.\n");
  fprintf(prefFile,"Do not edit by hand. Use the preferences from the main menu.\n");
  
  for(i=0;i<arrayMax(prefList);i++){
    item = arrp(prefList,i,PREFLIST);
    if(item->type == BOOLEAN){
      if(cf)
	item->value.bval = cf->value[i].bval;
      fprintf(prefFile,"_%s -%s\n",item->name,(item->value.bval ? "TRUE":"FALSE"));
    }
    else if(item->type  == INT){
      if(cf)
	item->value.ival = cf->value[i].ival;
      fprintf(prefFile,"_%s -I %d\n",item->name,item->value.ival);
    }
    else if(item->type  == COLOUR){
      if(cf)
	item->value.ival = cf->value[i].ival;
      fprintf(prefFile,"_%s -C %d\n",item->name,item->value.ival);
    }
    else if(item->type  == FLOAT){
      if(cf)
	item->value.fval = cf->value[i].fval;
      fprintf(prefFile,"_%s -R %f\n",item->name,item->value.fval);
    }
  }
  filclose(prefFile);
}


static void saveHomePrefs(){
  savePrefs();
}
#endif
static void makePrefDefaults(){
/* This creates the new defaults if none exist. You will have to add here any new */
/* preferences that you want to add. Remeber to add them to the .acedbrc file, */
/* the $ACEDB/wspec/acedbrc or remove the two completely. To have them take effect. */

  PREFLIST item;

  prefList  = arrayReCreate (prefList, 10, PREFLIST) ; 
  
  strcpy(item.name,"LARGER_FONTS");
  item.type = BOOLEAN; 
  item.value.bval = FALSE ;
  item.display = TRUE;
  array(prefList,arrayMax(prefList),PREFLIST) = item;

  strcpy(item.name,"OLD_STYLE_MAIN_WINDOW");
  item.type = BOOLEAN; 
  item.value.bval = FALSE ;
  item.display = TRUE;
  array(prefList,arrayMax(prefList),PREFLIST) = item;

  strcpy(item.name,"HORIZONTAL_TREE");
  item.type = BOOLEAN; 
  item.value.bval = TRUE;
  item.display = TRUE;
  array(prefList,arrayMax(prefList),PREFLIST) = item;

  strcpy(item.name,"NO_MESSAGE_WHEN_DISPLAY_BLOCK");
  item.type = BOOLEAN; 
  item.value.bval = FALSE;
  item.display = TRUE;
  array(prefList,arrayMax(prefList),PREFLIST) = item;

  strcpy(item.name,"SHOW_EMPTY_OBJECTS");
  item.type = BOOLEAN; 
  item.value.bval = FALSE;
  item.display = TRUE;
  array(prefList,arrayMax(prefList),PREFLIST) = item;

  strcpy(item.name,"TAG_COLOUR_IN_TREE_DISPLAY");
  item.type = COLOUR; 
  item.value.ival = BROWN ;
  item.display = TRUE;
  array(prefList,arrayMax(prefList),PREFLIST) = item;

  strcpy(item.name,"HIDE_NON_GOLD_TG");
  item.type =  BOOLEAN;
  item.value.ival = FALSE ;
  item.display = TRUE;
  array(prefList,arrayMax(prefList),PREFLIST) = item;

  return;
}

static void copySystemPrefs()
{
 /* Read in the system defaults if they exist. Else create a new set of defaults */
  FILE * prefFile = NULL;
  char *filename1=0;
  
  filename1 = filName("wspec/acedbrc","","r");
  if(filename1){
    prefFile = filopen(filename1,"","r");
    prefRead(prefFile);
  }
  else
    makePrefDefaults();
}

static FILE* openPrefFile(BOOL new){
  char filename[DIR_BUFFER_SIZE];
  char dirname[DIR_BUFFER_SIZE];
  FILE * prefFile = NULL;
  
  if(!new)
  {  if(getenv("HOME"))
       strcpy(filename,getenv("HOME"));
     strcat(filename,"/.acedbrc");

     if(!(access (filename, R_OK)))
       prefFile = filopen(filename,"","r");  
  }
  else
    { if(getenv("HOME"))
        strcpy(dirname,getenv("HOME"));
      prefFile = filqueryopen(dirname,filename,"","r","Choose Preferences file to load");
    }

  return prefFile;
}
static void prefAdd(PREFLIST *item){
  BOOL found =FALSE;
  int i;

  for(i=0;i<arrayMax(prefList);i++)
    { if(!strcmp(item->name,arrp(prefList,i,PREFLIST)->name)) /* is it the same name */
      {  if(arrp(prefList,i,PREFLIST)->type != item->type)
	{ messout("ERROR: Prefile reading in the preference %s but does not have the same type defintion as that in the master copy. This will be ignored.",item->name);
	  return;
	}
        if(arrp(prefList,i,PREFLIST)->type == BOOLEAN)
	  arrp(prefList,i,PREFLIST)->value.bval = item->value.bval;
	else if(arrp(prefList,i,PREFLIST)->type == INT ||arrp(prefList,i,PREFLIST)->type == COLOUR )
	  arrp(prefList,i,PREFLIST)->value.ival = item->value.ival; 
        else if(arrp(prefList,i,PREFLIST)->type == FLOAT)
	  arrp(prefList,i,PREFLIST)->value.fval = item->value.fval;
	else
	  { messout("ERROR: Prefile reading in the preference %s with invalid type, ignoring.",item->name);
	    return;
	  }
	arrp(prefList,i,PREFLIST)->display = TRUE;
	found = TRUE;
      }
    }
  if(!found){
    item->display = FALSE;
    array(prefList,arrayMax(prefList),PREFLIST) = *item;
  }
  return;
}
  

static void prefRead(FILE *prefFile){
  PREFLIST item;
  int i,j;  
  char junk[255],*cp;
  KEY colKey = 0 ;

  makePrefDefaults();
  while(fgets(junk,255,prefFile) != NULL){                 /* read in the preferences */
    if(junk[0] == '_'){
      i =1; j=0;
      while((isalpha((int)junk[i]) || junk[i] == '_') /*&&(junk[i]!=' ') &&junk[i]!='\t'*/ && j<30)
	item.name[j++]= junk[i++];
      item.name[j++] = '\0';
      if(j == 1)
	messout("Error reading name after _ in %s .acdbrc \n",junk);
      else{
	while(junk[i] != '-' && junk[i]!='\n')
	  i++;
	if(junk[i]=='-'){
	  if(junk[i+1] == 't' || junk[i+1] == 'T'){
	    item.value.bval = TRUE;
	    item.type = BOOLEAN;
	    prefAdd(&item);	      
	  }
	  else if(junk[i+1] == 'f' || junk[i+1] == 'F'){
	    item.value.bval = FALSE;
	    item.type = BOOLEAN;
	    prefAdd(&item);
	  }
	  else if(junk[i+1] == 'I' || junk[i+1] == 'i'){
	    cp = &junk[i+2];
	    if(sscanf(cp,"%d",&item.value.ival)){
	      item.type=INT;
	      prefAdd(&item);
	    }
	    else
	      messout("ERROR reading integer after %s",item.name);
	  }
	  else if(junk[i+1] == 'C' || junk[i+1] == 'c'){
	    cp = &junk[i+2];
	    if(sscanf(cp,"%d",&item.value.ival)){
	      item.type=COLOUR;
	      prefAdd(&item);
	    }
	    else if (lexword2key(cp,&colKey,0) && 
		     colKey >= _WHITE && 
		     colKey < _WHITE + 32){
	      item.type=COLOUR;
	      item.value.ival = WHITE + colKey - _WHITE ; /* tags to color enum */
	      prefAdd(&item);
	    }	      
	    else
	      messout("ERROR reading colour after %s",item.name);
	  }
	  else if(junk[i+1] == 'R' || junk[i+1] == 'r'){
	    cp = &junk[i+2];
	    if(sscanf(cp,"%f",&item.value.fval)){
	      item.type=FLOAT;
	      prefAdd(&item);
	    }
	    else
	      messout("ERROR reading float after %s",item.name);
	  }
	  else{
	    messout("Error reading value for %s.",item.name);
	  }
	}
      }
    }
  }
  filclose(prefFile);
}

#ifndef NON_GRAPHIC
void prefReRead(void)
{
  FILE * prefFile=NULL;

  prefFile = openPrefFile(TRUE);
  if(!prefFile)
    return;
  
  prefRead(prefFile);
  editPreferences();
}

void reReadPrefs()
{
  prefReRead() ;
}
#endif

void prefInit()    
{
  FILE * prefFile=NULL;

  prefFile = openPrefFile(FALSE);    /* open ~/.acedbrc for preferences */
  if(!prefFile){                /* if none create some */
    copySystemPrefs();
  }
  else
    prefRead(prefFile);
}

BOOL prefValue(char * name)
{
/* Returns the value for the char * name if it exists in the preferences. */
/* If it does not exist then returns the default FALSE. */
  int i;
  PREFLIST *item;

  if (!name || !*name || !arrayExists(prefList)) /* should not need this ??? */
    return FALSE ;
  for(i=0;i<arrayMax(prefList);i++){
    item = arrp(prefList,i,PREFLIST);
    if(strstr(name,item->name)){
      if(item->type == BOOLEAN)
	return item->value.bval;
      else{
	messout("%s not of type boolean",item->name);
	return FALSE;
      }
    }
  }
  return FALSE;   /* Default is set to false */
}
float prefFloat(char * name)
{
/* Returns the value for the char * name if it exists in the preferences. */
/* If it does not exist then returns the default 0.0. */
  int i;
  PREFLIST *item;

  if (!name || !*name || !arrayExists(prefList)) /* should not need this ??? */
    return 0.0 ;
  for(i=0;i<arrayMax(prefList);i++){
    item = arrp(prefList,i,PREFLIST);
    if(strstr(name,item->name)){
      if(item->type == FLOAT)
	return item->value.fval;
      else{
	messout("%s not of type float",item->name);
	return 0.0;
      }
    }
  }
  return 0.0;   /* Default is set to 0.0 */
}

int prefInt(char * name)
{
/* Returns the value for the char * name if it exists in the preferences. */
/* If it does not exist then returns the default FALSE. */
  int i;
  PREFLIST *item;

  if (!name || !*name || !arrayExists(prefList)) /* should not need this ??? */
    return 0 ;
  for(i=0;i<arrayMax(prefList);i++){
    item = arrp(prefList,i,PREFLIST);
    if(strstr(name,item->name)){
      if(item->type == INT || item->type == COLOUR) /* should really use prefColour */
	return item->value.ival;
      else{
	messout("%s not of type integer",item->name);
	return 0;
      }
    }
  }
  return 0;   /* Default is set to 0 */
}

int prefColour (char * name)
{
/* Returns the value for the char * name if it exists in the preferences. */
/* If it does not exist then returns the default FALSE. */
  int i;
  PREFLIST *item;

  if (!name || !*name || !arrayExists(prefList)) /* should not need this ??? */
    return 0 ;
  for(i=0;i<arrayMax(prefList);i++){
    item = arrp(prefList,i,PREFLIST);
    if(strstr(name,item->name)){
      if(item->type == COLOUR)
	return item->value.ival;
      else{
	messout("%s not of type colour",item->name);
	return 0;
      }
    }
  }
  return 0;   /* Default is set to 0 */
}

#ifndef NON_GRAPHIC
static void cancelPrefs(){
  /* remove cf struct and kill graph*/
  messfree(cf);
  graphDestroy();
}
static void donePrefs(){
  int i;
  PREFLIST *item;

  /* copy new values */
  for(i=0;i<arrayMax(prefList);i++){
    item = arrp(prefList,i,PREFLIST);
    if(item->type == BOOLEAN)
      item->value.bval = cf->value[i].bval;
    else if(item->type == INT || item->type == COLOUR)
      item->value.ival = cf->value[i].ival;
    else if(item->type == FLOAT)
      item->value.fval = cf->value[i].fval;
  }
  /* remove cf struct and kill graph */
  messfree(cf);
  graphDestroy();
  pickDraw();
}

static MENUOPT prefmenu[] = {
  {cancelPrefs, "Quit"},
  {graphPrint, "Print"},  
  {0, 0} };

void editPreferences(void)
{
/* Enables the user to set BOOL values to the preferences vis radio buttons */
  int i,line =3;
  PREFLIST *item;

  if (graphActivate (prefsGraph))
    { graphPop () ;
      graphClear();
    }
  else
    { prefsGraph = graphCreate (TEXT_SCROLL, "Preferences Menu", 
                                0, 0, 0.4, 0.4) ;
      graphMenu (prefmenu);
      cf = (struct configLocals *) messalloc(sizeof(struct configLocals));  
    }
  for(i=0;i<arrayMax(prefList);i++){
    item = arrp(prefList,i,PREFLIST);
    if(item->display)
      {
	if(item->type == BOOLEAN){
	  cf->value[i].bval = item->value.bval;
	  graphToggleEditor(item->name,&cf->value[i].bval,4.0,line++);
	}
	else if(item->type == INT){
	  cf->value[i].ival = item->value.ival;
	  graphIntEditor(item->name,&cf->value[i].ival,4.0,line++,0);
	}
	else if(item->type == COLOUR){
	  cf->value[i].ival = item->value.ival;
	  graphColourEditor(item->name," ",&cf->value[i].ival,4.0,line++);
	}
	else if(item->type == FLOAT){
	  cf->value[i].fval = item->value.fval;
	  graphFloatEditor(item->name,&cf->value[i].fval,4.0,line++,0);
	}
	line++;
      }    
  }
  graphButton("Quit",cancelPrefs,5.0,line);
  graphButton("Apply",donePrefs,11.0,line);
  graphButton("Save",saveHomePrefs,18.0,line);
  graphButton("Load",reReadPrefs,24.0,line);

  graphRedraw();
}
#endif 
 
