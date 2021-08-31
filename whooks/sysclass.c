/*  File: sysclass.c
 *  Author: Danielle et Jean Thierry-Mieg (mieg@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1994
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmba.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
      Definition of the system models, tags and classes needed by the acedb kernel
      inparticular during the bootstrap phase.
 
      DO NOT EDIT without consulting the authors.
 * Exported functions:
 * HISTORY:
 * Last edited: Oct 22 13:27 1998 (edgrif)
 * * Oct 22 13:26 1998 (edgrif): Move dec. of hardinit to sysclass.h
 * * Nov 17 17:58 1994 (mieg): incorporated sysoptions and graphTypes
 * Created: Wed Sep  7 16:10:30 1994 (mieg)
 *-------------------------------------------------------------------
 */

/* @(#)sysclass.c	1.35 1/29/97 */

#include "acedb.h"
#include "lex.h"
#include "pick.h"
#include "session.h"
#include "dbpath.h"

#include "whooks/sysclass.h"
#include "whooks/systags.h"

#include "disptype.h"

/********************************************************/
/********************************************************/

/* The sysmodel defined as a static string in this way
   no longer need to be read at run time from an external file
*/

Stack sysModels (void)
{ Stack s = stackCreate(1000) ;

  pushText(s, "// Acedb system models, please do not edit. ") ;
  pushText(s, " ") ;
  pushText(s, "?Model ANY   // must come first ") ;
  pushText(s, " ") ;
  pushText(s, "?Text   Quoted_in ANY ") ;
  pushText(s, " ") ;
  pushText(s, "?Comment  Quoted_in ANY ") ;
  pushText(s, " ") ;
  /*   2014_04_07   removed class keyword from the automatic system
    pushText(s, "?Keyword Quoted_in ANY ") ;
    pushText(s, " ") ;
  */
  pushText(s, "?Class   Compiled_as UNIQUE Text UNIQUE Int    // Consistency check  ") ;
  pushText(s, "         Belongs_to_class UNIQUE Int ") ;
  pushText(s, "         Is_a_subclass_of ?Class XREF Is_a_superclass_of UNIQUE Int // the mask of that subclass ") ;
  pushText(s, "         Is_a_superclass_of ?Class XREF Is_a_subclass_of UNIQUE Text UNIQUE Int //Query and mask ") ;
  pushText(s, "         Composite_of ?Class");
  pushText(s, "         Uses_tags ?Tag XREF Appears_in_class ") ;
  pushText(s, "         Visibility UNIQUE Buried ") ;
  pushText(s, "                           Hidden ") ;
  pushText(s, "                           Visible ") ;
  pushText(s, "         ID UNIQUE ID_local ID_template UNIQUE Text") ;
  pushText(s, "                            ID_nametype Text #ID_nametype_qualifiers") ;
  pushText(s, "                   ID_remote UNIQUE Text UNIQUE Int  // host and port") ;
  pushText(s, "         Mask   UNIQUE Int ") ;
  pushText(s, "         Filter UNIQUE Text ") ;
  pushText(s, "         Constraints Text ") ;
  pushText(s, "         Index ?Tag") ;
  pushText(s, "         Web_location UNIQUE ?WWW_server") ;
  pushText(s, " ") ;
  pushText(s, "#ID_nametype_qualifiers Is_unique") ;
  pushText(s, "                        Is_primary") ;
  pushText(s, "                        Users_with_write_access Text") ;
  pushText(s, " ") ;
  pushText(s, "#ID_history_action Event Created UNIQUE Text UNIQUE Text UNIQUE Text") ;
  pushText(s, "                         Killed") ;
  pushText(s, "                         Resurrected") ;
  pushText(s, "                         Merged_into UNIQUE Text") ;
  pushText(s, "                         Acquires_merge UNIQUE Text") ;
  pushText(s, "                         Split_from UNIQUE Text UNIQUE Text UNIQUE Text") ;
  pushText(s, "                         Split_into UNIQUE Text") ;
  pushText(s, "                         Add_name UNIQUE Text UNIQUE Text UNIQUE Text") ;
  pushText(s, "                         Remove_name UNIQUE Text UNIQUE Text UNIQUE Text") ;
  pushText(s, "                         Imported") ;
  pushText(s, "  // the three Texts are <id> <type> <name>") ;
  pushText(s, "  // Merged_into, Acquires_merge, Split_into, Split_from use <id>") ;
  pushText(s, "  // Add_name, Remove_name use <type> <name>") ;
  pushText(s, "  // if <type> is unique then Add_name replaces, as normal in acedb") ;
  pushText(s, " ") ;
  pushText(s, "?Tag     Parent_tag UNIQUE Text UNIQUE Int   // Name and value of Tag for Consistency check  ") ;
  pushText(s, "         Related_tags ?Tag XREF Related_tags ") ;
  pushText(s, "         Used_for Text Text  // Database_name semantics ") ;
  pushText(s, "         Appears_in_class ?Class XREF Uses_tags ") ;
  pushText(s, "         Appears_in_source_code ?SourceCode XREF Uses_tags ") ;
  pushText(s, "          ") ;
  pushText(s, "?SourceCode SourceFile ?LongText ") ;
  pushText(s, "	    Uses_tags ?Tag XREF Appears_in_source_code ") ;
  pushText(s, "	    Uses_class ?Class ") ;
  pushText(s, "            Includes ?Include XREF Included_by ") ;
  pushText(s, " ") ;
  pushText(s, "?Include SourceFile ?LongText ") ;
  pushText(s, "         Included_by ?SourceCode XREF Includes ") ;
  pushText(s, " ") ;
  pushText(s, "?Session        Session UNIQUE Int  ") ;
  pushText(s, "		Permanent_session     // If set, prevents automatic destruction ") ;
  pushText(s, "                Date   UNIQUE Text  ") ;
  pushText(s, "                User   UNIQUE Text  ") ;
  pushText(s, "                Session_Title  UNIQUE Text  ") ;
  pushText(s, "                CodeRelease UNIQUE Int UNIQUE Int  ") ;
  pushText(s, "                DataRelease UNIQUE Int UNIQUE Int  ") ;
  pushText(s, "		Created_from UNIQUE ?Session  ") ;
  pushText(s, "		Up_linked_to UNIQUE ?Session  ") ;
  pushText(s, "		Destroyed_by_session UNIQUE Int  ") ;
  pushText(s, "		GlobalLex   UNIQUE Int  ") ;
  pushText(s, "		SessionLex  UNIQUE Int UNIQUE Int UNIQUE Int UNIQUE Int  ") ;
  pushText(s, "		VocLex  UNIQUE Int UNIQUE Int UNIQUE Int UNIQUE Int  ") ;
  pushText(s, "		GlobalBat UNIQUE Int UNIQUE Int  // address # blocks used ") ;
  pushText(s, "		BatPlus  UNIQUE ?Bat  UNIQUE Int // # block set ") ;
  pushText(s, "		BatMinus UNIQUE ?Bat  UNIQUE Int // # blocks freed ") ;
  pushText(s, "		IndexVersion UNIQUE Int // automatic indexing ") ;
  pushText(s, " ") ;
  pushText(s, "?Display  Visibility // not yet used, but i need a model ") ;
  pushText(s, "                     // i will later turn wspec/displays.wrm into an ") ;
  pushText(s, "                     // ace file ") ;
  pushText(s, " ") ;
  pushText(s, "?Table  Title UNIQUE Text") ;
  pushText(s, "        Comment Text") ;
  pushText(s, "        Precompute ?TableResult Int // Int is session number") ;
  pushText(s, "        Sortcolumn UNIQUE Int") ;
  pushText(s, "        Parameters UNIQUE Text") ; /* default set of parameters */
  pushText(s, "        Colonne Int #Table_definition") ;
  pushText(s, " ") ;
  pushText(s, "#Table_definition  Subtitle UNIQUE Text") ;
  pushText(s, "                   Legend UNIQUE Text") ;
  pushText(s, "                   Origin UNIQUE From UNIQUE Int") ;
  pushText(s, "                                 Right_of UNIQUE Int") ;
  pushText(s, "                                 Copy UNIQUE Int") ;
  pushText(s, "                                 Calcul UNIQUE Int") ;
  pushText(s, "                   Tag UNIQUE Text") ;
  pushText(s, "                   DNA UNIQUE Text UNIQUE Text") ;
  pushText(s, "                   Visibility UNIQUE Hidden") ;
  pushText(s, "                                     Visible") ;
  pushText(s, "                   Width UNIQUE Int") ;
  pushText(s, "                   Presence UNIQUE Optional") ;
  pushText(s, "                                   Mandatory") ;
  pushText(s, "                                   Null") ;
  pushText(s, "                   Condition UNIQUE Text") ;
  pushText(s, "                   Type UNIQUE A_Class UNIQUE Text") ;
  pushText(s, "                               An_Int") ;   /* cannot use Int here */
  pushText(s, "                               A_Float") ;
  pushText(s, "                               A_Text") ;
  pushText(s, "                               A_Date") ;
  pushText(s, "                               Boolean"); 
  pushText(s, "                               Show_Tag") ;
  pushText(s, "                               Next_Tag") ;
  pushText(s, "                               Next_Key UNIQUE Text") ;
  pushText(s, "                   Extract UNIQUE All" ) ;
  pushText(s, "                                  MultiData") ;
  pushText(s, "                                  Min") ;
  pushText(s, "                                  Max") ;
  pushText(s, "                                  Average") ;
  pushText(s, "                                  Variance") ;
  pushText(s, "                                  Sum") ;
  pushText(s, "                                  Compute") ;
  pushText(s, "                                  Count") ;
  pushText(s, " ") ;
  pushText(s, "?Colour UNIQUE WHITE ") ;
  pushText(s, "               BLACK ") ;
  pushText(s, "               LIGHTGRAY ") ;
  pushText(s, "               DARKGRAY ") ;
  pushText(s, "               RED ") ;
  pushText(s, "               GREEN ") ;
  pushText(s, "               BLUE	 ") ;
  pushText(s, "               YELLOW ") ;
  pushText(s, "               CYAN ") ;
  pushText(s, "               MAGENTA ") ;
  pushText(s, "               LIGHTRED ") ;
  pushText(s, "               LIGHTGREEN ") ;
  pushText(s, "               LIGHTBLUE ") ;
  pushText(s, "               DARKRED ") ;
  pushText(s, "               DARKGREEN ") ;
  pushText(s, "               DARKBLUE ") ;
  pushText(s, "               PALERED ") ;
  pushText(s, "               PALEGREEN ") ;
  pushText(s, "               PALEBLUE ") ;
  pushText(s, "               PALEYELLOW ") ;
  pushText(s, "               PALECYAN ") ;
  pushText(s, "               PALEMAGENTA ") ;
  pushText(s, "               BROWN ") ;
  pushText(s, "               ORANGE ") ;
  pushText(s, "               PALEORANGE ") ;
  pushText(s, "               PURPLE ") ;
  pushText(s, "               VIOLET ") ;
  pushText(s, "               PALEVIOLET ") ;
  pushText(s, "               GRAY ") ;
  pushText(s, "               PALEGRAY ") ;
  pushText(s, "               CERISE ") ;
  pushText(s, "               MIDBLUE ") ;
  pushText(s, "               LIGHTMAGENTA") ;
  pushText(s, "               LIGHTCYAN") ;
  pushText(s, "               DARKVIOLET") ;
  pushText(s, "               LAVANDER") ;
  pushText(s, "               BLUE1") ;
  pushText(s, "               BLUE2") ;
  pushText(s, "               BLUE3") ;
  pushText(s, "               BLUE4") ;
  pushText(s, "               BLUE5") ;
  pushText(s, "               BLUE6") ;
  pushText(s, "               BLUE7") ;
  pushText(s, "               BLUE8") ;
  pushText(s, "               GREEN1") ;
  pushText(s, "               GREEN2") ;
  pushText(s, "               GREEN3") ;
  pushText(s, "               GREEN4") ;
  pushText(s, "               GREEN5") ;
  pushText(s, "               GREEN6") ;
  pushText(s, "               GREEN7") ;
  pushText(s, "               GREEN8") ;
  pushText(s, "               RED1") ;
  pushText(s, "               RED2") ;
  pushText(s, "               RED3") ;
  pushText(s, "               RED4") ;
  pushText(s, "               RED5") ;
  pushText(s, "               RED6") ;
  pushText(s, "               RED7") ;
  pushText(s, "               RED8") ;
  pushText(s, "               LIGHTVIOLET ") ;
  pushText(s, "               DARKCYAN ") ;
  pushText(s, "               LIGHTORANGE ") ;
  pushText(s, " ") ;
  pushText(s, "?UserSession   Session UNIQUE ?Session ") ;
  pushText(s, "               Start   UNIQUE DateType ") ;
  pushText(s, "               Finish  UNIQUE DateType ") ;
  pushText(s, "               User    UNIQUE Text  // use text") ;
  pushText(s, "               Keysets New UNIQUE ?KeySet") ;
  pushText(s, "                       Touched UNIQUE ?KeySet") ;
  pushText(s, "                       Deleted UNIQUE ?KeySet") ;
  pushText(s, "               Program UNIQUE Client") ;
  pushText(s, "                              Server") ;
  pushText(s, "                              tace") ;
  pushText(s, "                              xace") ;
  pushText(s, "               Remark ?Text") ;
  pushText(s, " ") ;
  /*  // Jade Display  Class, jade-class, jade-table, ace-table
      //      StandAloneDisplay Text
      */
  pushText(s, "?Jade  Display Text Text Text UNIQUE Text") ; 
  pushText(s, "       StandAloneDisplay Text") ;
  pushText(s, " ") ;
  pushText(s, "?View  Type") ; /*ensures a non null model */
  pushText(s, " ") ;
  pushText(s, "// End of the system models ") ;
  
  return s ;
}

/********************************************************/
/********************************************************/

 /* default MUST be hidden, so that wspec can fully ignore those classes */

void sysClassOptions(void) 
{ 
  extern BOOL READING_MODELS ;
  BOOL old = READING_MODELS ;
  READING_MODELS = TRUE ;

  /*
  * pickSetClassOption sets various attributes of the class.  Here, we
  * set several special classes.
  *
  * name is the class name.
  * visibility is whether you see it on the GUI in xace - h=hidden, v=visible
  * protected means only the kernel can manipulate this class
  * private is private to kernel, do not dump, invisibleKey 
  * case_sensitive is whether the object names in this class are case sensitive
  * update_names means that we set the object name every time it gets edited.  This
  *    means you can change the case of the object name, even though the case-insensitive
  *    name remains the same.  Has no effect if case_sensitive is true
  */

                   /*                 visibility   protected   case_sensitive        */
                   /*  name        type   |   display |    private  |   update_names */
	 	   /*  |             |    |    |      |     |       |     |          */
  pickSetClassOption("System",      'A', 'h', ZERO, TRUE , TRUE , FALSE, FALSE) ;
  pickSetClassOption("Global",      'A', 'h', ZERO, TRUE , TRUE , FALSE, FALSE) ;
  pickSetClassOption("MainClasses", 'A', 'h', ZERO, TRUE , TRUE , FALSE, FALSE) ;
  pickSetClassOption("Session",     'B', 'H', TREE, TRUE , FALSE, FALSE, FALSE) ;
  pickSetClassOption("Voc",         'A', 'h', ZERO, TRUE , TRUE , FALSE, FALSE) ;
  pickSetClassOption("Bat",         'A', 'h', ZERO, TRUE , TRUE , FALSE, FALSE) ;

  pickSetClassOption("KeySet",      'A', 'v', DtKeySet, 
                                                    FALSE, FALSE, FALSE, FALSE) ;
  pickSetClassOption("Model",       'B', 'v', TREE, TRUE , FALSE, FALSE, FALSE) ; 

  /* up to _VModel, the model cannot 
  * use #constructs 
  */

  pickSetClassOption("Calcul",      'A', 'h', ZERO, TRUE , TRUE , FALSE, FALSE) ;
  pickSetClassOption("Colour",      'B', 'h', ZERO, TRUE , FALSE, FALSE, FALSE) ;

  /* if Storage_type is X, auto 
  * XREF is enabled 
  */

  pickSetClassOption("Text",        'X', 'h', TREE, FALSE, FALSE, FALSE, TRUE ) ;
  pickSetClassOption("Comment",     'X', 'h', TREE, TRUE , FALSE, FALSE, TRUE ) ;
  /*   pickSetClassOption("Keyword",     'B', 'h', TREE, FALSE, FALSE, TRUE , FALSE) ; */
  pickSetClassOption("UserSession", 'B', 'h', TREE, TRUE , FALSE, FALSE, FALSE) ;

  pickSetClassOption("LongText",    'A', 'h', DtLongText, 
                                                    FALSE, FALSE, FALSE, FALSE) ;

  pickSetClassOption("Class",       'B', 'H', TREE, TRUE , FALSE, FALSE, FALSE) ; /* since needed to check consistency  */
  pickSetClassOption("Table",       'B', 'h', TREE, FALSE, FALSE, FALSE, FALSE) ; 
  pickSetClassOption("TableResult", 'A', 'h', ZERO, FALSE, FALSE, FALSE, FALSE) ; 
  pickSetClassOption("Jade",        'B', 'h', TREE, FALSE, FALSE, FALSE, FALSE) ; 
  pickSetClassOption("View",        'B', 'h', TREE, FALSE, FALSE, FALSE, FALSE) ; 
  pickSetClassOption("Binary",      'A', 'h', ZERO, FALSE, FALSE, FALSE, FALSE) ;

  pickSetClassOption("SourceCode",  'B', 'h', TREE, FALSE, FALSE, FALSE, FALSE) ;
  pickSetClassOption("Include",     'B', 'h', TREE, FALSE, FALSE, FALSE, FALSE) ;

  pickSetClassOption("Query",       'A', 'h', ZERO, TRUE , FALSE, FALSE, FALSE) ;
  pickSetClassOption("Constraint",  'A', 'h', ZERO, TRUE , FALSE, FALSE, FALSE) ;
  pickSetClassOption("Display",     'B', 'h', ZERO, TRUE , FALSE, FALSE, FALSE) ;
  
  pickSetClassOption("DNA",         'A', 'h', FMAP, FALSE, FALSE, FALSE, FALSE) ;
  pickSetClassOption("Peptide",     'A', 'h', PEPMAP, 
                                                    FALSE, FALSE, FALSE, FALSE) ;

  pickSetClassOption("MatchTable",  'A', 'h', ZERO, FALSE, FALSE, FALSE, FALSE) ;
  pickSetClassOption("BaseCall",    'A', 'h', ZERO, FALSE, FALSE, FALSE, FALSE) ;
  pickSetClassOption("cMap",        'A', 'h', CMAP, FALSE, FALSE, FALSE, FALSE) ;
  pickSetClassOption("gMap",        'A', 'h', PMAP, FALSE, FALSE, FALSE, FALSE) ;
  pickSetClassOption("pMap",        'A', 'h', PMAP, FALSE, FALSE, FALSE, FALSE) ;
  pickSetClassOption("vMap",        'A', 'h', VMAP, FALSE, FALSE, FALSE, FALSE) ;
  pickSetClassOption("FicheView",   'A', 'h', ZERO, FALSE, FALSE, FALSE, FALSE) ;

  READING_MODELS = old ;
}

/********************************************************/
/********************************************************/

int
  _VBat, _VKeySet, _VCalcul,
  _VClass,  /* the class hierarchy and class aliases */
  _VModel,  /* the model of all the main classes and the subtypes */
  _VText, _VLongText, _VImage, _VTable, _VTableResult, _VJade,
  _VComment, _VUserSession, _VView,
  _VQuery, _VConstraint, _VFicheView, _VBinary ;

extern 
  void lexReHashClass(int classe) ;
     
extern void lexHardDefineKey (int table, KEY key,  char *cp) ;


void sysClassInit (void)
{ KEY key ;
  extern BOOL READING_MODELS ;
  BOOL old = READING_MODELS ;
  READING_MODELS = TRUE ;

     /* note that i cannot register the zeroth class _VSystem */
  lexHardDefineKey (_VMainClasses, _VGlobal, "Global") ;
  lexHardDefineKey (_VMainClasses, _VSession , "Session") ;
  lexHardDefineKey (_VMainClasses, _VVoc, "Voc") ;

  lexHardDefineKey (_VMainClasses, _VDisplay, "Display") ;
  lexHardDefineKey (_VMainClasses, _VMainClasses, "MainClasses") ;

  lexReHashClass(_VMainClasses) ;  /* needed following the hardDefine calls */ 

  lexaddkey ("Bat", &key, _VMainClasses) ; _VBat = KEYKEY(key) ;
  lexaddkey ("KeySet", &key, _VMainClasses) ; _VKeySet = KEYKEY(key) ;
  lexaddkey ("Calcul", &key, _VMainClasses) ; _VCalcul = KEYKEY(key) ;

  lexaddkey ("Class", &key, _VMainClasses) ; _VClass = KEYKEY(key) ;

  lexaddkey ("Model", &key, _VMainClasses) ; _VModel = KEYKEY(key) ;

  lexaddkey ("Text", &key, _VMainClasses) ; _VText = KEYKEY(key) ;
  lexaddkey ("LongText", &key, _VMainClasses) ; _VLongText = KEYKEY(key) ;
  lexaddkey ("Image", &key, _VMainClasses) ; _VImage = KEYKEY(key) ;
  lexaddkey ("Table", &key, _VMainClasses) ; _VTable = KEYKEY(key) ;
  lexaddkey ("TableResult", &key, _VMainClasses) ; _VTableResult = KEYKEY(key) ;
  lexaddkey ("Jade", &key, _VMainClasses) ; _VJade = KEYKEY(key) ;
  lexaddkey ("View", &key, _VMainClasses) ; _VView = KEYKEY(key) ;
  lexaddkey ("FicheView", &key, _VMainClasses) ; _VFicheView = KEYKEY(key) ;

  lexaddkey ("Comment", &key, _VMainClasses) ; _VComment = KEYKEY(key) ;
  lexaddkey ("UserSession", &key, _VMainClasses) ; _VUserSession = KEYKEY(key) ;

  lexaddkey ("Query", &key, _VMainClasses) ; _VQuery = KEYKEY(key) ;
  lexaddkey ("Constraint", &key, _VMainClasses) ; _VConstraint = KEYKEY(key) ;
  lexaddkey ("Binary", &key, _VMainClasses) ; _VBinary = KEYKEY(key) ;

  READING_MODELS = old ;
}

/******************************************************************/

/* the systags are #defined in systags.wrm
   here we simply initialise the system lexique
   in a garanteed compatible way
*/

void lexDefineSystemTags (void) 
{ KEY key ;

  lexHardDefineKey(0,_Text,"Text") ; 
  lexHardDefineKey(0,_AddressText,"AddressText") ; 
  lexHardDefineKey(0,_Greek,"Greek") ; 
  lexHardDefineKey(0,_Russian,"Russian") ; 
  lexHardDefineKey(0,_Text1,"Text1") ; 
  lexHardDefineKey(0,_Text2,"Text2") ; 
  lexHardDefineKey(0,_Text3,"Text3") ; 
  lexHardDefineKey(0,__DNA1,"_DNA1") ; 
  lexHardDefineKey(0,__DNA2,"_DNA2") ; 
  lexHardDefineKey(0,__DNA3,"_DNA3") ; 
  lexHardDefineKey(0,__DNA4,"_DNA4") ; 
  lexHardDefineKey(0,__DNA5,"_DNA5") ; 
  lexHardDefineKey(0,_LongFloat,"LongFloat") ; 
  lexHardDefineKey(0,_LongInt,"LongInt") ; 
  lexHardDefineKey(0,__RNA1,"_RNA1") ; 
  lexHardDefineKey(0,__RNA2,"_RNA2") ; 
  lexHardDefineKey(0,__RNA3,"_RNA3") ; 
  lexHardDefineKey(0,__RNA4,"_RNA4") ; 
  lexHardDefineKey(0,__Protein1,"_Protein1") ; 
  lexHardDefineKey(0,__Protein2,"_Protein2") ; 
  lexHardDefineKey(0,__Protein3,"_Protein3") ; 
  lexHardDefineKey(0,_NextC,"NextC") ; 
  lexHardDefineKey(0,_LastC,"LastC") ; 
  lexHardDefineKey(0,_Int,"Int") ; 
  lexHardDefineKey(0,_Unsigned,"Unsigned") ; 
  lexHardDefineKey(0,_zLong,"zLong") ; 
  lexHardDefineKey(0,_zLong_Unsigned,"zLong_Unsigned") ; 
  lexHardDefineKey(0,_Float,"Float") ; 
  lexHardDefineKey(0,_DateType,"DateType") ; 
  lexHardDefineKey(0,_continuationKey,"continuationKey") ; 
  lexHardDefineKey(0,_LastN,"LastN") ; 

  lexHardDefineKey(0,_bsHere,"bsHere") ; 
  lexHardDefineKey(0,_bsRight,"bsRight") ; 
  lexHardDefineKey(0,_bsDown,"bsDown") ; 

  lexHardDefineKey(0,_UNIQUE,"UNIQUE") ; 
  lexHardDefineKey(0,_XREF,"XREF") ; 
  lexHardDefineKey(0,_ANY,"ANY") ; 
  lexHardDefineKey(0,_FREE,"FREE") ; 
  lexHardDefineKey(0,_REPEAT,"REPEAT") ; 
  lexHardDefineKey(0,_COORD,"COORD") ; 
  lexHardDefineKey(0,_SORTED,"SORTED") ; 

  lexHardDefineKey(0,_Date,"Date") ; 
  lexHardDefineKey(0,_User,"User") ; 
  lexHardDefineKey(0,_Session,"Session") ; 
  lexHardDefineKey(0,_BatPlus,"BatPlus") ; 
  lexHardDefineKey(0,_BatMinus,"BatMinus") ; 
  lexHardDefineKey(0,_GlobalBat,"GlobalBat") ; 
  lexHardDefineKey(0,_Session_Title,"Session_Title") ; 
  lexHardDefineKey(0,_SessionLex,"SessionLex") ; 
  lexHardDefineKey(0,_VocLex,"VocLex") ; 
  lexHardDefineKey(0,_GlobalLex,"GlobalLex") ; 
  lexHardDefineKey(0,_Quoted_in,"Quoted_in") ; 
  lexHardDefineKey(0,_CodeRelease,"CodeRelease") ; 
  lexHardDefineKey(0,_DataRelease,"DataRelease") ; 
  lexHardDefineKey(0,_Created_from,"Created_from") ; 
  lexHardDefineKey(0,_Image,"Image") ; 
  lexHardDefineKey(0,_Pick_me_to_call,"Pick_me_to_call") ; 
  lexHardDefineKey(0,_File,"File") ; 
  lexHardDefineKey(0,_Non_graphic,"Non_graphic") ; 
  lexHardDefineKey(0,_Centre,"Centre") ; 
  lexHardDefineKey(0,_Destroyed_by_session,"Destroyed_by_session") ; 
  lexHardDefineKey(0,_Up_linked_to,"Up_linked_to") ; 
  lexHardDefineKey(0,_Permanent_session,"Permanent_session") ; 
  lexHardDefineKey(0,_Secret,"Secret") ; 
  lexHardDefineKey(0,_This_session,"This_session") ; 
  lexHardDefineKey(0,_Related_tags,"Related_tags") ; 
  lexHardDefineKey(0,_Used_for,"Used_for") ; 
  lexHardDefineKey(0,_Appears_in_source_code,"Appears_in_source_code") ; 
  lexHardDefineKey(0,_Uses_tags,"Uses_tags") ; 
  lexHardDefineKey(0,_SourceFile,"SourceFile") ; 
  lexHardDefineKey(0,_Includes,"Includes") ; 
  lexHardDefineKey(0,_Included_by,"Included_by") ; 
  lexHardDefineKey(0,_Parent_tag,"Parent_tag") ; 
  lexHardDefineKey(0,_Appears_in_class,"Appears_in_class") ; 
  lexHardDefineKey(0,_Compiled_as,"Compiled_as") ; 
  lexHardDefineKey(0,_Is_a_subclass_of,"Is_a_subclass_of") ; 
  lexHardDefineKey(0,_Is_a_superclass_of,"Is_a_superclass_of") ; 
  lexHardDefineKey(0,_Uses_class,"Uses_class") ; 
  lexHardDefineKey(0,_Visible,"Visible") ; 
  lexHardDefineKey(0,_Hidden,"Hidden") ; 
  lexHardDefineKey(0,_Mask,"Mask") ; 
  lexHardDefineKey(0,_Belongs_to_class,"Belongs_to_class") ; 
  lexHardDefineKey(0,_Filter,"Filter") ; 
  lexHardDefineKey(0,_Visibility,"Visibility") ; 
  lexHardDefineKey(0,_Start,"Start") ; 
  lexHardDefineKey(0,_Finish,"Finish") ; 
  lexHardDefineKey(0,_PROTECTED,"PROTECTED") ; 
  lexHardDefineKey(0,_Constraints,"Constraints") ; 
  lexHardDefineKey(0,_File_name,"File_name") ; 
  lexHardDefineKey(0,_IndexVersion,"IndexVersion") ; 

  lexHardDefineKey(0,_WHITE,"WHITE") ; 
  lexHardDefineKey(0,_BLACK,"BLACK") ; 
  lexHardDefineKey(0,_LIGHTGRAY,"LIGHTGRAY") ; 
  lexHardDefineKey(0,_DARKGRAY,"DARKGRAY") ; 
  lexHardDefineKey(0,_RED,"RED") ; 
  lexHardDefineKey(0,_GREEN,"GREEN") ; 
  lexHardDefineKey(0,_BLUE,"BLUE") ; 
  lexHardDefineKey(0,_YELLOW,"YELLOW") ; 
  lexHardDefineKey(0,_CYAN,"CYAN") ; 
  lexHardDefineKey(0,_MAGENTA,"MAGENTA") ; 
  lexHardDefineKey(0,_LIGHTRED,"LIGHTRED") ; 
  lexHardDefineKey(0,_LIGHTGREEN,"LIGHTGREEN") ; 
  lexHardDefineKey(0,_LIGHTBLUE,"LIGHTBLUE") ; 
  lexHardDefineKey(0,_DARKRED,"DARKRED") ; 
  lexHardDefineKey(0,_DARKGREEN,"DARKGREEN") ; 
  lexHardDefineKey(0,_DARKBLUE,"DARKBLUE") ;
  lexHardDefineKey(0,_PALERED,"PALERED") ;
  lexHardDefineKey(0,_PALEGREEN,"PALEGREEN") ;
  lexHardDefineKey(0,_PALEBLUE,"PALEBLUE") ;
  lexHardDefineKey(0,_PALEYELLOW,"PALEYELLOW") ;
  lexHardDefineKey(0,_PALECYAN,"PALECYAN") ;
  lexHardDefineKey(0,_PALEMAGENTA,"PALEMAGENTA") ;
  lexHardDefineKey(0,_BROWN,"BROWN") ;
  lexHardDefineKey(0,_ORANGE,"ORANGE") ;
  lexHardDefineKey(0,_PALEORANGE,"PALEORANGE") ;
  lexHardDefineKey(0,_PURPLE,"PURPLE") ;
  lexHardDefineKey(0,_VIOLET,"VIOLET") ;
  lexHardDefineKey(0,_PALEVIOLET,"PALEVIOLET") ;
  lexHardDefineKey(0,_GRAY,"GRAY") ;
  lexHardDefineKey(0,_PALEGRAY,"PALEGRAY") ;
  lexHardDefineKey(0,_CERISE,"CERISE") ;
  lexHardDefineKey(0,_MIDBLUE,"MIDBLUE") ;
  lexHardDefineKey(0,_LIGHTMAGENTA,"LIGHTMAGENTA") ;
  lexHardDefineKey(0,_LIGHTCYAN,"LIGHTCYAN") ;
  lexHardDefineKey(0,_DARKVIOLET,"DARKVIOLET") ;
  lexHardDefineKey(0,_LAVANDER,"LAVANDER") ;
  lexHardDefineKey(0,_BLUE1,"BLUE1") ;
  lexHardDefineKey(0,_BLUE2,"BLUE2") ;
  lexHardDefineKey(0,_BLUE3,"BLUE3") ;
  lexHardDefineKey(0,_BLUE4,"BLUE4") ;
  lexHardDefineKey(0,_BLUE5,"BLUE5") ;
  lexHardDefineKey(0,_BLUE6,"BLUE6") ;
  lexHardDefineKey(0,_BLUE7,"BLUE7") ;
  lexHardDefineKey(0,_BLUE8,"BLUE8") ;
  lexHardDefineKey(0,_GREEN1,"GREEN1") ;
  lexHardDefineKey(0,_GREEN2,"GREEN2") ;
  lexHardDefineKey(0,_GREEN3,"GREEN3") ;
  lexHardDefineKey(0,_GREEN4,"GREEN4") ;
  lexHardDefineKey(0,_GREEN5,"GREEN5") ;
  lexHardDefineKey(0,_GREEN6,"GREEN6") ;
  lexHardDefineKey(0,_GREEN7,"GREEN7") ;
  lexHardDefineKey(0,_GREEN8,"GREEN8") ;
  lexHardDefineKey(0,_RED1,"RED1") ;
  lexHardDefineKey(0,_RED2,"RED2") ;
  lexHardDefineKey(0,_RED3,"RED3") ;
  lexHardDefineKey(0,_RED4,"RED4") ;
  lexHardDefineKey(0,_RED5,"RED5") ;
  lexHardDefineKey(0,_RED6,"RED6") ;
  lexHardDefineKey(0,_RED7,"RED7") ;
  lexHardDefineKey(0,_RED8,"RED8") ;
  lexHardDefineKey(0,_LIGHTVIOLET,"LIGHTVIOLET") ;
  lexHardDefineKey(0,_DARKCYAN,"DARKCYAN") ;
  lexHardDefineKey(0,_LIGHTORANGE,"LIGHTORANGE") ;

  lexHardDefineKey(0,_LastSystemTag,"LastSystemTag") ;

  lexReHashClass(0) ;  /* needed following the hardDefine calls */
  
  lexaddkey("_Global_Bat",&key,1) ;
  lexaddkey("_lexi1",&key,1) ;
  lexaddkey("_lexa1",&key,1) ;
  lexaddkey("_voc1",&key,1) ;
  lexaddkey("_lexi2",&key,1) ;
  lexaddkey("_lexa2",&key,1) ;
  lexaddkey("_voc2",&key,1) ;
  lexaddkey("_lexi3",&key,1) ;
  lexaddkey("_lexa3",&key,1) ;
  lexaddkey("_voc3",&key,1) ;
  lexaddkey("_batPlus",&key,1) ;
  lexaddkey("_batMinus",&key,1) ;
  lexaddkey("_oldPlus",&key,1) ;
  lexaddkey("_oldMinus",&key,1) ;
  lexaddkey("_lexh1",&key,1) ;
  lexaddkey("_lexh2",&key,1) ;
  lexaddkey("_lexh3",&key,1) ;
  lexaddkey("_lext2",&key,1) ;
  lexaddkey("_lext3",&key,1) ;
  lexaddkey("_lext1",&key,1) ;
  lexaddkey("_superkey",&key,1) ;
}

/********************************************************************************/
/* classHardInit, a fossible from acedb1, removed */
/********************************************************************************/
/********************************************************************************/
 
 
