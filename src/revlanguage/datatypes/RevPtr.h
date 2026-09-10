/**
 * @file RevPtr.h
 *
 * SuperSmartPointer:
 *
 * From Chapter 25 of Professional C++ 
 * by Nicholas A. Solter and Scott J. Kleper 
 * (Wrox, 2005) 
 *
 * @brief Declaration of ReferenceCount
 *
 * (c) Copyright 2009-
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 */

#ifndef RevPtr_H
#define RevPtr_H

#include <string>

namespace RevLanguage {

template <typename ptrType>
class RevPtr {
    
public:
    RevPtr(ptrType* inPtr = NULL);
    ~RevPtr(void);
    RevPtr(const RevPtr<ptrType>& src);
    RevPtr&                                     operator=(const RevPtr<ptrType>& rhs);
    ptrType&                                    operator*(void) const;
    ptrType*                                    operator->(void) const;
    operator                                    ptrType*(void) const { return mPtr; }
    static bool                                 decrementCountForAddress(ptrType* qPtr);
    static void                                 incrementCountForAddress(ptrType* qPtr);

                                                operator RevPtr<const ptrType> () const;
    
private:
    //size_t                                      countForAddress(const ptrType* qPtr);
    void                                        finalizePointer(void);
    void                                        initializePointer(ptrType* inPtr);
    
    ptrType*                                    mPtr;
};
    
}


#include "RbConstants.h"

#include <iostream>


template <typename ptrType>
RevLanguage::RevPtr<ptrType>::RevPtr(ptrType* inPtr) : mPtr( NULL )
{
    
    initializePointer(inPtr);
}

template <typename ptrType>
RevLanguage::RevPtr<ptrType>::RevPtr(const RevPtr<ptrType>& src) : mPtr( NULL )
{
    
    initializePointer(src.mPtr);
}

template <typename ptrType>
RevLanguage::RevPtr<ptrType>::operator RevPtr<const ptrType> () const
{
    return RevPtr<const ptrType>(mPtr);
}

template <typename ptrType>
RevLanguage::RevPtr<ptrType>::~RevPtr(void)
{
    
    finalizePointer();
}

template <typename ptrType>
RevLanguage::RevPtr<ptrType>& RevLanguage::RevPtr<ptrType>::operator=(const RevPtr<ptrType>& rhs)
{
    
    if (this == &rhs)
    {
        return (*this);
    }
    finalizePointer();
    initializePointer(rhs.mPtr);
    return (*this);
}


template <typename ptrType>
void RevLanguage::RevPtr<ptrType>::initializePointer(ptrType* inPtr) {
    
    mPtr = inPtr;
    
    incrementCountForAddress(mPtr);
}

template <typename ptrType>
void RevLanguage::RevPtr<ptrType>::finalizePointer(void)
{
    
    if ( decrementCountForAddress(mPtr) )
    {
        delete mPtr;
    }
    
}

template <typename ptrType>
ptrType* RevLanguage::RevPtr<ptrType>::operator->(void) const
{
    
    return (mPtr);
}

template <typename ptrType>
ptrType& RevLanguage::RevPtr<ptrType>::operator*(void) const
{
    
    if ( mPtr == NULL )
    {
        std::cerr << "Major BUG: Dereferencing NULL pointer!!!" << std::endl;
    }
    
    return (*mPtr);
}


/*template <typename ptrType>
size_t RevLanguage::RevPtr<ptrType>::countForAddress(const ptrType* qPtr)
{
    
    // check if we got the NULL pointer
    if (qPtr == 0)
    {
        // size_t cannot be nan (WP)
        return RbConstants::Size_t::nan;
    }
    
    return qPtr->getReferenceCount();
}*/


template <typename ptrType>
void RevLanguage::RevPtr<ptrType>::incrementCountForAddress(ptrType* qPtr)
{
    
    // check if we got the NULL pointer
    if (qPtr == 0) return;
    
    qPtr->incrementReferenceCount();
}

template <typename ptrType>
bool RevLanguage::RevPtr<ptrType>::decrementCountForAddress(ptrType* qPtr)
{
    
    // check if we got the NULL pointer
    if (qPtr == 0) return false;
    
    size_t refCount = qPtr->decrementReferenceCount();
    
    return refCount == 0;
}


#endif

