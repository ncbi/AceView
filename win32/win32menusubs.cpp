/*  File: win32menusubs.cpp
 *  Author: Richard M. Bruskiewich (rbrusk@octogene.medgen.ubc.ca)
 *          inspired by xtsub.c code designed by Richard Durbin (rd@mrc-lmb.cam.ac.uk)
 *          and Christopher Lee
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1992
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: Menu for Win32 version of graph package
 * Exported functions: graphBoxMenu, graphFreeMenu, graphNewBoxMenu
 *						graphDevFinish, graphMenuCleanup
 * HISTORY:
 * Last edited: Jun 11 00:45 1996 (rbrusk): normal looking freemenus
 *		-	free menus made to look like other menus
 * * Dec 27 01:15 1995 (rmb): rewritten cleaner? 
 * Recreated: Jun  27 20:33 1995 (rmb)
 *-------------------------------------------------------------------
 */

/* $Id: win32menusubs.cpp,v 1.1.1.1 2002/07/19 20:23:25 sienkiew Exp $ */

#include "stdafx.h"

extern "C"
{
#include "regular.h"
#include "menu_.h"
}

#include "winace.h"
#include "cgraph.h" // gDev instantiation

#if defined(VERBOSE_DEBUG)
#undef VERBOSE_DEBUG
#endif

#ifdef _DEBUG
#undef THIS_FILE
static char BASED_CODE THIS_FILE[] = __FILE__;
#endif

#define new DEBUG_NEW

/******* menus ********/

extern "C" int menuBox ;


// WIN32 specific: dynamically assigns a new WM_COMMAND ID number for menus
int nextCmdID()
{
	static UINT nextID = _APS_NEXT_COMMAND_VALUE ; // pick up where AppStudio left off
	return nextID++ ;
}

void devMenuDestroy(void)
{
// ACEDB likes to recycle menus, so I won't destroy them
// until the program exits; see graphMenuFinish() ;
}

static void destroyNormalMenu(void *menu) 
{
	int i ;
	NORMALMENUBITS *normalBits = (NORMALMENUBITS *)menu ;

	delete(normalBits->menuObj) ;
	delete(normalBits->menuMap) ;
	for (i=0; (normalBits->opts)[i].f; i++)
		messfree((normalBits->opts)[i].text) ;
	messfree(normalBits->opts) ;
	messfree(normalBits) ;
}

static void destroyFreeMenu(void *menu) 
{
	unsigned int i ;
	FREEMENUBITS *freeBits = (FREEMENUBITS *)menu ;

	delete(freeBits->menuObj) ;
	delete(freeBits->menuMap) ;
  	for (i=0; i<= freeBits->opts->key; i++)
		messfree((freeBits->opts)[i].text) ;
	messfree(freeBits->opts) ;
	messfree(freeBits) ;
}


static void destroyMenuObj( CMenu *menuObj ) 
{
	int pos = menuObj->GetMenuItemCount() - 1 ;
	do
	{
		if( menuObj->GetSubMenu(pos) != NULL ) // if submenu exists
			menuObj->RemoveMenu(pos,MF_BYPOSITION);	// then remove it
		// Note: the submenu object will be deleted independently later?
	} 	while( pos-- ) ; // Until item 0 has been processed
	delete(menuObj) ;
}

// This routine should be invoked independently for all 
// distinct main and sub NEWMENUBITS records in newMenus
static void destroyNewMenu(void *menu) 
{
	NEWMENUBITS *newBits = (NEWMENUBITS *)menu ;
	destroyMenuObj( newBits->menuObj ) ;

	// MENU memory just release from handle0 now
	// menuDestroy(newBits->m) ;

	if( newBits->isaMainMenu ) delete newBits->menuMap ;
	messfree(newBits) ;
}

static void dumpNormalMenu(CDumpContext& dc, void *menu) 
{
	int i ;
	NORMALMENUBITS *normalBits = (NORMALMENUBITS *)menu ;

	dc 	<< 	"\n*** Start of Normal Menu Item Dump ***\n"
		<<	"\n\tMenu item ptr == " << menu 
		<<	"\n\n\tCMenu:\n" ;
	normalBits->menuObj->Dump(dc) ;

	dc 	<< "\n\tMenu Map:\n" ;
	normalBits->menuMap->Dump(dc) ;

	dc 	<< "\n\tMENUOPTs:\n" ;
	for (i=0; (normalBits->opts)[i].f; i++)
		dc 	<< "\t\tText == " << (normalBits->opts)[i].text
			<< ", Function Ptr == " << (normalBits->opts)[i].f << "\n";
	dc 	<< "\n\n*** End of Normal Menu Item Dump ***\n" ;
}
static void dumpFreeMenu(CDumpContext& dc, void *menu) 
{
	unsigned int i ;
	FREEMENUBITS *freeBits = (FREEMENUBITS *)menu ;

	dc 	<< 	"\n*** Start of Free Menu Item Dump ***\n"
		<<	"\n\tMenu item ptr == " << menu 
		<<	"\n\n\tCFreeMenu:\n" ;
	freeBits->menuObj->Dump(dc) ;

	dc 	<< "\n\tMenu Map:\n" ;
	freeBits->menuMap->Dump(dc) ;

	dc 	<< "\n\tFREEOPTs:\n" ;
  	for (i=0; i<= freeBits->opts->key; i++)
		dc 	<< "\t\tKey == " << (freeBits->opts)[i].key
			<< ", Text == " << (freeBits->opts)[i].text << "\n";
	dc 	<< "\n\n*** End of Free Menu Item Dump ***\n" ;
}

static void dumpNewMenuItems(CDumpContext& dc, MENU menu )
{
	dc 	<< "\n\t\tMenu Title == " << menu->title << "\n" ;
	for (MENUITEM item = menu->items ; item ; item = item->down)
	{
		dc 	<< "\n\t\tMenu Item:"
			<< "\n\t\t\tLabel ==\t\t" 	<< item->label << "\n" 
			<< "\n\t\t\tFunction ==\t\t"<< item->func << "\n"
			<< "\n\t\t\tFlags ==\t\t" 	<< item->flags << "\n"
			<< "\n\t\t\tCall ==\t\t\t" 	<< item->call << "\n"
			<< "\n\t\t\tValue ==\t\t" 	<< item->value << "\n"
			<< "\n\t\t\tPtr ==\t\t\t" 	<< item->ptr << "\n\n" ;

		if(item->submenu)
		{
			dc 	<< "\n\t\tSUBMENU:\n" ;
			dumpNewMenuItems(dc, item->submenu ) ;
		}
	}
} 

static void dumpNewMenu(CDumpContext& dc, void *menu) 
{
	NEWMENUBITS *newBits = (NEWMENUBITS *)menu ;

	dc 	<<	"\n*** Start of New Menu Item Dump ***\n"
		<<	"\n\tMenu item ptr ==" << menu 
		<<	"\n\n\tCMenu:\n" ;
	newBits->menuObj->Dump(dc) ;

	dc 	<< "\n\tMenu Map:\n" ;
	newBits->menuMap->Dump(dc) ;

	dc 	<< "\n\tMENU:\n" ;
	dumpNewMenuItems(dc, newBits->m ) ;

	dc 	<< "\n\n*** End of New Menu Item Dump ***\n" ;
}

static CMultiAss	*normalMenus	= NULL,
					*freeMenus		= NULL,
					*newMenus		= NULL ;

void graphMenuInit()
{
	normalMenus = new CMultiAss(destroyNormalMenu,dumpNormalMenu) ;
	freeMenus = new CMultiAss(destroyFreeMenu,dumpFreeMenu) ;
	newMenus = new CMultiAss(destroyNewMenu,dumpNewMenu) ;
}

void graphMenuFinish()
{
	delete normalMenus ; normalMenus = NULL ;
	delete freeMenus ; freeMenus = NULL ;
	delete newMenus ; newMenus = NULL ;
}

static Box getGraphBox(int k)
{
	ASSERT( k >= 0 ) ;
	Box box = gBoxGet (k) ;
  
#if !defined(NEW_WIN32_GRAPHS)  
  	if( k > BOXVIEWS.GetUpperBound() || BOXVIEW( k ) == NULL)
  	{
  		CGraphBox *newBox = new CGraphBox( VIEWPTR, k, box )	;
		TRY
		{	
   			BOXVIEWS.SetAtGrow( k, newBox ) ;
		}
		CATCH(CMemoryException, e )
		{
			// I should do something intelligent here?
			delete newBox ; // like this?
		}
		END_CATCH
  	}
#else // defined(NEW_WIN32_GRAPHS)
	BOOL NeedNewBox = FALSE ;
	BOXVIEWARRAY *bva ;
	USESUBDEVPTR(bva,GetBoxViews,TRUE )  // gets a BOXVIEWARRAY ptr for current gSubDev
	CGraphBox *gb ;

	if( k > bva->GetUpperBound() )
		NeedNewBox = TRUE ;
	else
	{
		USESUBDEVPTR(gb,BoxView,k) // gets a CGraphBox ptr for box 'k' in current gSubDev
		if(gb == NULL) NeedNewBox = TRUE ;
	}
	if( NeedNewBox )
  	{
  		gb = new CGraphBox( VIEWPTR, k, box )	;
		TRY
		{	
   			bva->SetAtGrow( k, gb ) ;
		}
		CATCH(CMemoryException, e )
		{
			// I should do something intelligent here?
			delete gb ; // like this?
		}
		END_CATCH
  	}
#endif

	return box ;
}		

static void registerMenu (int k, BoxMenuType bType, void *menu)
{
	Box box = getGraphBox (k) ;

#if defined(VERBOSE_DEBUG)
	TRACE("Registering a %s menu (CObject* &%p) for box #%d of graph #%d\n",
			(const char *)MenuType(bType),menu,k,gActive->id) ;
#endif

#if !defined(NEW_WIN32_GRAPHS)  
	BOXVIEW(k)->SetMenu( bType, menu ) ;
#else // defined(NEW_WIN32_GRAPHS)
	CGraphBox *gb ;
	USESUBDEVPTR(gb,BoxView,k) // gets a CGraphBox ptr for box 'k' in current gSubDev
	gb->SetMenu( bType, menu ) ;
#endif

	box->flag |= GRAPH_BOX_MENU_FLAG ;
}

static void unregisterMenu (int k)
{
	Box box = getGraphBox (k) ;

#if defined(VERBOSE_DEBUG)
	TRACE("Unregistering menu for box #%d of graph #%d\n",k,gActive->id) ;
#endif

#if !defined(NEW_WIN32_GRAPHS)  
	BOXVIEW(k)->SetMenu( NO_MENU ) ;
#else // defined(NEW_WIN32_GRAPHS)
	CGraphBox *gb ;
	USESUBDEVPTR(gb,BoxView,k) // gets a CGraphBox ptr for box 'k' in current gSubDev
	gb->SetMenu( NO_MENU ) ;
#endif

	box->flag &= ~GRAPH_BOX_MENU_FLAG ;
}

static void SetMenuItem( MENUMAP *map, WORD menuItemID, void *callBack )
{
	TRY
	{
		map->SetAt( menuItemID, callBack );
	}
	CATCH(CMemoryException,e)
	{
		// Should do something intelligent here
		messcrash("Memory exception in GraphBox::SetMenuItem()") ;
	}
	END_CATCH
}

static void *GetMenuItem( MENUMAP *map, WORD menuItemID ) 
{
	void *callBack ;
	if( !map->Lookup( menuItemID, callBack ) )
		return NULL ;
	else return callBack ;
}

void graphBoxMenu (int box, MENUOPT *options)
{
  NORMALMENUBITS *menu;
  MENUOPT *o;			 
  int i;

  if (!gDev ) return ;

  if (!options)
    { unregisterMenu (box) ;
      return ;
    }

  /* If we find an existing menu here, reuse it */
  for( 	menu = (NORMALMENUBITS *)(normalMenus->GetFirstItem( options->f )) ;
  		menu ;
		menu = (NORMALMENUBITS *)(normalMenus->GetNextItem()) )
	{
	    o = menu->opts;
		for (i=0; o[i].f && options[i].f; i++)
		 if ( o[i].f != options[i].f ||
		      strcmp(o[i].text, options[i].text) )
		     	break ;
		if (!o[i].f && !options[i].f)
		{  
		   registerMenu (box,  NORMAL_MENU, (void *)menu);
		   return;
		}
	}

  /*  otherwise, make a new one here. */

  for (i=0; options[i].f; i++);
  o = (MENUOPT *) messalloc((i+1)*sizeof(MENUOPT));

  for (i=0; options[i].f; i++)
    { o[i].f = options[i].f;
      o[i].text = (char *) messalloc(1+strlen(options[i].text));
      strcpy(o[i].text, options[i].text);
    }
  o[i].f =  0;
  o[i].text = 0;
  menu = (NORMALMENUBITS *)messalloc(sizeof(NORMALMENUBITS));
  menu->opts = o;
  
   
  //menu->widg = XtCreatePopupShell ("menu", simpleMenuWidgetClass, 
	//			   root_widget, args, 0) ;

  menu->menuMap = new MENUMAP ;
  menu->menuObj = new CMenu ;
  menu->menuObj->CreatePopupMenu() ;

  normalMenus->InsertItem( options->f, menu ) ;

  for(i=0; o[i].f; i++)
  {
    // entry = XtCreateManagedWidget (o[i].text,
	//			     smeBSBObjectClass,
	//			     menu->widg, NULL, 0) ;
    //  XtAddCallback (entry, XtNcallback,
	//	     oldMenuSelect, (XtPointer) o[i].f) ;

	UINT menuItemID = nextCmdID() ;
	SetMenuItem( menu->menuMap, (WORD) menuItemID, (void *)o[i].f ) ;
	menu->menuObj->AppendMenu( MF_ENABLED, menuItemID, o[i].text ) ;

	// WIN32 menu item callback mechanism implemented with
	// OnMenuSelect() override in CGraphView
  }
  registerMenu (box,  NORMAL_MENU, (void *)menu) ;
}

void graphBoxFreeMenu (int box, FreeMenuFunction proc, FREEOPT *options)
{
	FREEMENUBITS *menu ;
 	FREEOPT *o;
	unsigned int i;

	if (!gDev) return ;

	if (!proc || !options)
	{
		unregisterMenu (box) ;
		return ;
	}

  for( 	menu = (FREEMENUBITS *)(freeMenus->GetFirstItem( proc ) ) ;
  		menu ;
		menu = (FREEMENUBITS *)(freeMenus->GetNextItem()) )
	{
		o = menu->opts;
		if (menu->proc != proc) goto notfound;
		if (o->key != options->key ) goto notfound;
		for (i=1; i <= o->key; i++)
		{
			if (o[i].key != options[i].key) goto notfound;
			if ( 0 != strcmp(o[i].text, options[i].text)) goto notfound;
		}
		/* Found an existing menu here, reuse it */
		registerMenu (box, FREE_MENU, (void *)menu);
		return;
	notfound:
		continue;
    }

	/*  make a new one here. */
	o = (FREEOPT *)messalloc(((options->key)+1)*sizeof(FREEOPT));

	for (i=0; i<= options->key; i++)
	{
		o[i].key = options[i].key;
		o[i].text = (char *) messalloc(1+strlen(options[i].text));
		strcpy(o[i].text, options[i].text);
	}
	menu = (FREEMENUBITS *)messalloc(sizeof(FREEMENUBITS));
	menu->opts = o;
	menu->proc = proc;

	//menu->widg = XtCreatePopupShell ("menu", simpleMenuWidgetClass, 
	//				   root_widget, args, 0) ;

	menu->menuMap = new MENUMAP ;
	menu->menuObj  = new CMenu ;
	menu->menuObj->CreatePopupMenu() ;
 	
    freeMenus->InsertItem( proc, menu ) ;

	for (i=1; i <= o->key; i++)
	{
		UINT menuItemID = nextCmdID() ;
		SetMenuItem( menu->menuMap, (WORD) menuItemID, (void *)o[i].key ) ;
		menu->menuObj->AppendMenu( MF_ENABLED, menuItemID, o[i].text ) ;
	}
	registerMenu (box, FREE_MENU, (void *)menu) ;
}


/******************************* new menus *************************/

//********************************************************************
// createNewMenuItems(): creates the MENUITEM's of the menu are the
// contents of a standalone popup menu or of the popup associatedwith a 
// top level menubar;  both are created as popups here in this function.
//********************************************************
static NEWMENUBITS *createNewMenu( MENU menu, MENUMAP *menuMap = NULL ) ;  // Forward declaration 

static CMenu *createNewMenuItems( MENUITEM menuItems, MENUMAP *menuMap )
{
	CMenu *menuObj = new CMenu ; // the popup menu
	menuObj->CreatePopupMenu() ;

	for (MENUITEM item = menuItems ; item ; item = item->down)
	{
		unsigned int flags = item->flags ;
		if (flags & MENUFLAG_HIDE) // Don't show hidden menu items 
			continue ;
		if (flags & MENUFLAG_SPACER) // Use MF_SEPARATORs for menu spacer items
		{
			menuObj->AppendMenu( MF_SEPARATOR ) ;
			continue ;
		}

		UINT winFlags = MF_STRING ;
		if (flags & MENUFLAG_DISABLED) 
			winFlags |= MF_GRAYED ;  // grayed and disabled item?
		else
			winFlags |= MF_ENABLED ;

		if(item->submenu)
		{
			winFlags |= MF_POPUP ;
			NEWMENUBITS *submenuBits = createNewMenu( item->submenu, menuMap ) ;
			menuObj->AppendMenu( winFlags,(UINT)(submenuBits->menuObj->GetSafeHmenu()), item->label ) ;
		}								   
		else
		{
			// Record callback for menu item even if it is to be MF_GRAYED...
			// This permits later dynamic enabling of the item if desired
			UINT menuItemID = nextCmdID() ;
			SetMenuItem( menuMap, (WORD) menuItemID, (void *)item ) ;
			menuObj->AppendMenu( winFlags, menuItemID, item->label ) ;
		}
	}
	return menuObj ;
} 

//*************************************************************************
// createNewMenu(): creates a CMenu menuObj patterned after the MENU.
// menu->title is generally ignored within this function but may be
// used in the calling routine for labelling conventional WIN32 menuBars?
// This function is indirectly called recursively for MENUITEM submenus.
//*************************************************************************
static NEWMENUBITS *createNewMenu( MENU theMenu, MENUMAP *menuMap ) 
{
	NEWMENUBITS *menuBits ;
	MENUITEM item, item2 ;

	// Can an existing menu be found for reuse?
    for(menuBits = (NEWMENUBITS *)(newMenus->GetFirstItem( theMenu ) ) ;
  		menuBits ;
		menuBits = (NEWMENUBITS *)(newMenus->GetNextItem()) )
	{
		if (strcmp (theMenu->title, menuBits->m->title))
			continue ;
		for ( item = theMenu->items, item2 = menuBits->m->items ; 
   			  item && item2 ;
   			  item = item->down, item2 = item2->down )
			if ( item->label != item2->label ||
				 item->func != item2->func ||
				 item->flags != item2->flags ||
				 item->call != item2->call ||
				 item->value != item2->value ||
				 item->ptr != item2->ptr ||
				 item->submenu != item2->submenu )
  				break ;
		if (!item && !item2)	/* found an existing menu - reuse it */
			return menuBits;
	}

	// Otherwise, create one...
	menuBits = (NEWMENUBITS *) messalloc (sizeof (NEWMENUBITS)) ;
	menuBits->m = menuCopy (theMenu) ;

	if(menuMap == NULL)  // isaMainMenu?
	{
		menuBits->menuMap = new MENUMAP ; // then, create a new menuMap...
		menuBits->isaMainMenu = TRUE ;  // flag for menuMap deletion at menu destruction
	}
	else // isaSubMenu?
	{
		menuBits->menuMap = menuMap ; // then, use the existing menuMap...
		menuBits->isaMainMenu = FALSE ; // and don't touch the menuMap later...
	}
	menuBits->menuObj = createNewMenuItems( menuBits->m->items, menuBits->menuMap) ;

	newMenus->InsertItem( theMenu, menuBits ) ; // Keep track of the menu for reuse
	
	return menuBits ;
}

void graphNewBoxMenu (int box, MENU newMenu)
{
	NEWMENUBITS *menuBits ;
	if (!gDev) return ;

	if (!newMenu) { unregisterMenu (box) ; return ; }

	// Get a new menu (may be pre-existing) 
	menuBits = createNewMenu( newMenu /*, NULL menuMap => isaMainMenu */ ) ;

	// then, register it with the current box
	registerMenu (box, NEWBOX_MENU,(void *)menuBits) ;
}

extern UINT selectedMenuItem ;	// see CGraphView.cpp

void displayMenu( CPoint menuPos )
{
	ASSERT( menuBox >= 0 );
	getGraphBox(menuBox) ;  // just-in-case CGraphBox is not yet defined
							// this call may create a NULL box 

	// Convert to absolute screen coordinates...
	VIEWPTR->ClientToScreen( &menuPos ) ;

#if !defined(NEW_WIN32_GRAPHS)  
	BoxMenuType bmt = BOXVIEW( menuBox )->GetMenuType() ;
	
	if( bmt == NO_MENU ) return ;

	void *pBoxMenu = BOXVIEW( menuBox )->GetMenu() ;
#else // defined(NEW_WIN32_GRAPHS)
	CGraphBox *gb ;
	USESUBDEVPTR(gb,BoxView,menuBox) // gets a CGraphBox ptr for box 'k' in current gSubDev

	BoxMenuType bmt = gb->GetMenuType() ;
	
	if( bmt == NO_MENU ) return ;

	void *pBoxMenu = gb->GetMenu() ;
#endif

	if(pBoxMenu == NULL)
	{
		// Technically, executing this ASSERT is a bug because displayMenu()
		// should not be called if GRAPH_BOX_MENU_FLAG is not set
		// and the BoxMenuType should be NO_MENU
		ASSERT(FALSE) ;
		return ;
	}

	switch( bmt )
	{
		case FREE_MENU:
			{
				FREEMENUBITS *pFreeMenu = (FREEMENUBITS *)pBoxMenu ;
	
				if( !( pFreeMenu->menuObj->TrackPopupMenu( TPM_CENTERALIGN | TPM_RIGHTBUTTON,
										  menuPos.x, menuPos.y, VIEWPTR , NULL ) ) )
				{
					ASSERT(FALSE) ;  // if could not succeed in opening popup
					return ;
				}

				if( selectedMenuItem )
				{
					KEY arg = (KEY)GetMenuItem(pFreeMenu->menuMap, (WORD) selectedMenuItem ) ;

					if( pFreeMenu->proc != NULL )
						(pFreeMenu->proc)( arg, menuBox ) ;
				}
			}
			break ;

		case NORMAL_MENU:
			{
				NORMALMENUBITS *pNormalMenu = (NORMALMENUBITS *)pBoxMenu ;
	
				if( !( pNormalMenu->menuObj->TrackPopupMenu( TPM_CENTERALIGN | TPM_RIGHTBUTTON,
										  menuPos.x, menuPos.y, VIEWPTR , NULL ) ) )
				{
					ASSERT(FALSE) ;  // if could not succeed in opening popup
					return ;
				}
	
				if( selectedMenuItem )
				{
					MenuFunction fn = (MenuFunction)
						GetMenuItem(pNormalMenu->menuMap, (WORD) selectedMenuItem ) ;

					if( fn != NULL )
						fn(menuBox) ;  // callback function
				}
			}
			break ;

		case NEWBOX_MENU:
			{
				NEWMENUBITS *pNewMenu = (NEWMENUBITS *)pBoxMenu ;

				if( !( pNewMenu->menuObj->TrackPopupMenu( TPM_CENTERALIGN | TPM_RIGHTBUTTON,
										  menuPos.x, menuPos.y, VIEWPTR , NULL ) ) )
				{
					ASSERT(FALSE) ;  // if could not succeed in opening popup
					return ;
				}

				if( selectedMenuItem )
				{
					MENUITEM item = (MENUITEM)
						GetMenuItem( pNewMenu->menuMap, (WORD) selectedMenuItem ) ;

					if( item != NULL )
						(item->func)(item) ;  // callback function?
				}
			}
			break ;
		default: 
			ASSERT(FALSE) ;  // bug... invalid menu type
	}

}

//******* Right button menu pick adapted from gLeftDown() in graphsub.c
//        invoked in the OnRButtonDown() event handler

void gRightDown (CPoint point, CPoint origin ) 
{
	Box box ;
	float ux, uy ;
	// ACEDB graph coordinate in virtual viewport
	CPoint acePoint = point ;
	acePoint.Offset(origin) ;

	// Determine which box captures the menu invocation in current gActive graph
	ux = (float)XtoUabs( acePoint.x ) ;
	uy = (float)YtoUabs( acePoint.y ) ;
	for (menuBox = gActive->nbox ; menuBox-- ;)
	{
		box = gBoxGet (menuBox) ;
		if ( menuBox && !(box->flag & GRAPH_BOX_MENU_FLAG ))
			continue ;
		if (ux >= box->x1 && ux <= box->x2 && 
			uy >= box->y1 && uy <= box->y2)
			break ;
	}

	/* no menu found or box outside even the whole drawing area! */
	if (menuBox == -1) 	
		return ;

	// At this point, you should have the identity of the box
	// containing a menu, hence capturing the right button ?

	displayMenu( point ) ; // at device window relative point (x,y)

	//	*** Don't know whether any pass through is required to the
	//      ACEDB graph RIGHT_DOWN event callback ***
	//
	// if (gActive->func[RIGHT_DOWN])
	//	{
	//		mfn = (MouseFunc)(gActive->func[RIGHT_DOWN]);
	//		(*mfn)((double)x,(double)y) ;
	//	}
}


/************ end of file *************/
 
 
