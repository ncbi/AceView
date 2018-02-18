/* $Id: cmultiass.h,v 1.1.1.1 2002/07/19 20:23:26 sienkiew Exp $ */
/* CMultiAss.h : header file
 *  Author: Richard Bruskiewich
 *  Copyright (C) R. Bruskiewich, J Thierry-Mieg and R Durbin, 1996
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmba.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:   WIN32 emulation of Multi-Insert Associator
 * Exported:  class CMultiAss
 * HISTORY:
 * Last edited:
 * Created: Dec 30 11:56:08 1995 (rmb)
 *-------------------------------------------------------------------
 */
////////////////////////////////////////////////////////////
// WIN32 specific implementation of a multi-insertion
// associative memory heap
class CMultiAss : public CObject
{
public:

// Constructor
//	"destroyItemFunc" is an optional user specified item deletion function  used in destructor
//	"dumpItemFunc" is an optional user specified item dc << dump function used in CMultiAss::Dump()
	CMultiAss( 	void (*destroyItemFunc)(void *) = NULL,
				void (*dumpItemFunc)( CDumpContext& , void *) = NULL ) ;

protected:
	DECLARE_DYNAMIC(CMultiAss)

private:
	typedef CPtrList Bucket ;
	CTypedPtrMap<CMapPtrToPtr,void *,Bucket *> m_map ;
	void *m_currentKey ;
	Bucket *m_currentBucket ;

	// POSITIONs in the m_currentBucket
	POSITION m_currentPosition, m_nextPosition ;

	void (*m_destroyItemFunc)( void *item ) ;
	void (*m_dumpItemFunc)( CDumpContext& dc, void *item ) ;

public:
// Operations
	void  InsertItem(void *key, void *item) ;

	// If keyItem == NULL, then first instance of key is removed
	// return item deleted if any item under key was found and deleted from CMultiAss
	// otherwise, returns NULL
	void *RemoveItem(void *key, void *keyItem = NULL) ;

	// A RemoveItem() which uses m_destroyItemFunc to destroy the item
	// return TRUE if any item under key was found and destroyed; else FALSE
	BOOL DestroyItem(void *key, void *keyItem = NULL) ;
	
	void *GetFirstItem(void *key) ;
	void *GetNextItem() ;
	void *FindItem(void *key, void *keyItem = NULL) ;

// Implementation
	virtual ~CMultiAss();

#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

};
 
