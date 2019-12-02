#ifndef RbOrderedSet_H
#define RbOrderedSet_H

#include "StringUtilities.h"


#include <vector>
#include <iostream>
#include <sstream>

namespace RevBayesCore {
    
    template <class valueType>
    class RbOrderedSet {
        
    public:
        typedef typename std::vector<valueType>::iterator          iterator;
        typedef typename std::vector<valueType>::const_iterator    const_iterator;
        
        // constructor(s)
        RbOrderedSet();
        virtual                                            ~RbOrderedSet(void);
        
        // public member functions
        iterator                                            begin(void);
        const_iterator                                      begin(void) const;
        void                                                clear(void);
        iterator                                            end(void);
        const_iterator                                      end(void) const;
        void                                                erase(const valueType &v);
        iterator                                            find(const valueType &v);
        void                                                insert(const valueType &v);
        size_t                                              size(void) const;

//        RbOrderedSet<valueType>*                            clone(void) const;                                                                      //!< Create an independent clone
//        void                                                printElement(std::ostream &o, size_t i, std::string sep="\t", int l=-1, bool left=true) const;                                          //!< Print the i-th element
        
//        void                                                sort(void);
        
    private:
        
        std::vector<valueType>                              values;
        
    };

}


#include <algorithm>


template <class valueType>
RevBayesCore::RbOrderedSet<valueType>::RbOrderedSet()
{
    
}




template <class valueType>
RevBayesCore::RbOrderedSet<valueType>::~RbOrderedSet( void )
{
    
}


//template <class valueType>
//RevBayesCore::RbOrderedSet<valueType>* RevBayesCore::RbOrderedSet<valueType>::clone(void) const {
//    
//    return new RbOrderedSet<valueType>( *this );
//}



template <class valueType>
typename RevBayesCore::RbOrderedSet<valueType>::const_iterator RevBayesCore::RbOrderedSet<valueType>::begin(void) const
{

    return values.begin();
}


template <class valueType>
typename RevBayesCore::RbOrderedSet<valueType>::iterator RevBayesCore::RbOrderedSet<valueType>::begin(void)
{
    
    return values.begin();
}



template <class valueType>
void RevBayesCore::RbOrderedSet<valueType>::clear( void )
{
    
    return values.clear();
}


template <class valueType>
typename RevBayesCore::RbOrderedSet<valueType>::const_iterator RevBayesCore::RbOrderedSet<valueType>::end(void) const
{
    
    return values.end();
}


template <class valueType>
typename RevBayesCore::RbOrderedSet<valueType>::iterator RevBayesCore::RbOrderedSet<valueType>::end(void)
{
    
    return values.end();
}


template <class valueType>
void RevBayesCore::RbOrderedSet<valueType>::erase(const valueType &v)
{
    typename std::vector<valueType>::iterator pos = std::find(values.begin(), values.end(), v);
    if ( pos != values.end() )
    {
        values.erase( pos );
    }
    
}




template <class valueType>
typename RevBayesCore::RbOrderedSet<valueType>::iterator RevBayesCore::RbOrderedSet<valueType>::find(const valueType &v)
{
    
    typename std::vector<valueType>::iterator pos = std::find(values.begin(), values.end(), v);
    
    return pos;
}


template <class valueType>
void RevBayesCore::RbOrderedSet<valueType>::insert(const valueType &v)
{
    typename std::vector<valueType>::iterator pos = std::find(values.begin(), values.end(), v);
    if ( pos == values.end() )
    {
        values.push_back( v );
    }
    
}


template <class valueType>
size_t RevBayesCore::RbOrderedSet<valueType>::size( void ) const
{
    
    return values.size();
}



#endif

