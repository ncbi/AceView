// $Id: cmultiass.cpp,v 1.1.1.1 2002/07/19 20:23:25 sienkiew Exp $
/* CMultiAss.cpp : implementation file
 *  Author: Richard Bruskiewich
 *  Copyright (C) R. Bruskiewich, J Thierry-Mieg and R Durbin, 1996
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmba.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: WIN32 emulation of Multi-Insert Associator
 *				Constructor may be supplied with Destructor and
 *				Dump callback functions for user cleanup/destruction
 *				of CObjects contained in the Associator
 * Exported:  class CMultiAss
 * HISTORY:
 * Last edited:
 * Created: Dec 30 11:56:08 1995 (rmb)
 *-------------------------------------------------------------------
 */

#include "stdafx.h"
#include "win32.h"
#include "WinAce.h"
#include "CGraph.h"

#ifdef _DEBUG
#undef THIS_FILE
static char BASED_CODE THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CMultiAss

IMPLEMENT_DYNAMIC(CMultiAss, CObject)

#define new DEBUG_NEW

/////////////////////////////////////////////////////////////////////////////
// CMultiAss constructor

CMultiAss::CMultiAss(	void (*destroyItemFunc)(void *),
						void (*dumpItemFunc)(CDumpContext&, void *) )
{
	m_currentKey = NULL ;
	m_currentBucket = NULL ;
	m_currentPosition = m_nextPosition = NULL ;
	m_destroyItemFunc = destroyItemFunc ;
	m_dumpItemFunc = dumpItemFunc ;
}

/////////////////////////////////////////////////////////////////////////////
// CMultiAss destructor

CMultiAss::~CMultiAss()
{
	POSITION mapPos ;
	Bucket *bucket ;
	void * key ;
	void * item ;

	mapPos = m_map.GetStartPosition() ;
	while( mapPos != NULL )
	{
		m_map.GetNextAssoc(mapPos,key,bucket) ;
		m_nextPosition = bucket->GetHeadPosition() ;
		while( (m_currentPosition = m_nextPosition) != NULL )
		{
			bucket->GetNext(m_nextPosition) ;
			item = bucket->GetAt(m_currentPosition) ;
			bucket->RemoveAt(m_currentPosition) ;
			if(m_destroyItemFunc) (m_destroyItemFunc)(item) ;
			// else, user is responsible for destroying items! 
		}
		delete bucket ;
	}
}

// Saves menu objects in a map/bucket pool
// based upon a many-to-one key mapping
void CMultiAss::InsertItem( void *key, void *item )
{
  if ( !m_map.Lookup( key, m_currentBucket ) )
  {													  
	m_currentBucket	= new Bucket ; 
	TRY
	{	
		m_map.SetAt( key, m_currentBucket ) ;
	}
	CATCH(CMemoryException, e )
	{
		delete m_currentBucket ;
		messcrash("Out of Memory Exception during Heap creation?") ;
	}
	END_CATCH
  }
  TRY
  {	
	m_currentBucket->AddHead( item ) ;
  }
  CATCH(CMemoryException, e )
  {
 	messcrash("Out of Memory Exception during Heap Bucket addition?") ;
  }
  END_CATCH

}

void *CMultiAss::GetFirstItem(void *key)
{
	m_currentKey = key;
	if ( m_map.Lookup( m_currentKey, m_currentBucket ) )
	{
		m_nextPosition = m_currentPosition = m_currentBucket->GetHeadPosition() ;
		return( m_currentBucket->GetNext( m_nextPosition ) ) ;
	}
	else return NULL ;
}

void *CMultiAss::GetNextItem()
{
	if( (m_currentPosition = m_nextPosition) == NULL)
		return NULL ;
	else
		return( m_currentBucket->GetNext( m_nextPosition ) ) ;
}

// If keyitem == NULL, then first item found under key is returned
// Returns item (which is NULL, if not found)
void *CMultiAss::FindItem(void *key, void *keyItem)
{
	// First, find the item
	void *item = GetFirstItem(key) ;

	// If keyItem is defined, the item must match it as well
	if ( item && keyItem )
		while( item && item != keyItem )
			item = GetNextItem() ;

	return item ; // Could be NULL sometimes
} 

// If keyitem == NULL, then first item under key is removed;
// returns item deleted, if any item under key was found and deleted from CMultiAss
// otherwise, returns NULL
void *CMultiAss::RemoveItem(void *key, void *keyItem)
{
	void *item ;
	// If no item found, then fail
	// FindItem() should set m_currentPosition and m_NextPosition appropriately
	if(!( item = FindItem(key, keyItem) ) ) return NULL ;

	// otherwise, remove the item at m_currentPosition in the m_currentBucket
	m_currentBucket->RemoveAt(m_currentPosition) ;
	m_currentPosition = m_nextPosition ;  // Just a precaution...

	// destroy empty buckets...
	if( m_currentBucket->IsEmpty() )
	{
		delete m_currentBucket ;
		m_map.RemoveKey(key) ;
		m_currentKey = NULL ;
	}
	return item ;
}

// A RemoveItem() which uses m_destroyItemFunc to destroy the item
// return TRUE,  if any item under key was found and destroyed
// return FALSE, if no m_destroyItemFunc() is defined or if removal failed 
BOOL CMultiAss::DestroyItem(void *key, void *keyItem)
{	
	void *item ;
	if( !m_destroyItemFunc ||
		!(item = RemoveItem(key,keyItem) ) )
		return FALSE ;
	else
	{
		(m_destroyItemFunc)(item) ;
		return TRUE ;
	}
}

/////////////////////////////////////////////////////////////////////////////
// CGraph diagnostics

#ifdef _DEBUG
void CMultiAss::AssertValid() const
{
	CObject::AssertValid();
}

void CMultiAss::Dump(CDumpContext& dc) const
{
	Bucket *bucket ;
	void * item ;

	// touching m_currentKey etc. would be an unwanted side effect?
	POSITION hpos, bpos ;
	void * key ;

	dc	<< "\n\t" ;
	CObject::Dump(dc);

	dc	<< "\nCMultiAss specific data members:\n\n"
		<< "\tCurrentKey == " << m_currentKey 
		<< "\n\tCurrent Bucket Position ==\t"	<< m_currentPosition 
		<< "\n\tNext Bucket Position ==\t\t"	<< m_nextPosition 
		<< "\n\n\tData Map:\n" ;

	// Now dump the CMultiAss specific data? 
	hpos = m_map.GetStartPosition() ;
	while( hpos != NULL )
	{
		m_map.GetNextAssoc(hpos,key,bucket) ;
		dc	<< "\n\tAt Bin POSITION == " << hpos 
			<< ",\tKey == " << key 
			<< "\n\n\tBucket Dump:\n" ;

		bpos = bucket->GetHeadPosition() ;
		while( bpos != NULL )
		{
			dc	<< "\n\t\tAt Bucket POSITION == " << bpos 
				<< ",\tItem:\n"  ;

			item = bucket->GetNext(bpos) ;
			if(!m_dumpItemFunc)
				dc	<< item << "\n" ;
			else
				(m_dumpItemFunc)( dc, item) ; 
		}
	}
}
#endif //_DEBUG
 
