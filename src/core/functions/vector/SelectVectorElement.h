#ifndef SelectVectorElement_H
#define SelectVectorElement_H

#include "RbVector.h"
#include "TypedFunction.h"
#include "RbException.h"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace RevBayesCore {

    template <class valueType>
    class SelectVectorElement : public TypedFunction<valueType> {
        
    public:
        SelectVectorElement(const TypedDagNode<std::int64_t>* index, const std::vector<const TypedDagNode<valueType>*>& elements);
        virtual                                            ~SelectVectorElement(void);                                       //!< Virtual destructor
        
        SelectVectorElement*                                clone(void) const;                                               //!< Create an independent clone
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode* oldP, const DagNode* newP); //!< Implementation of swapping parameters
        
    private:
        
        // Members
        const TypedDagNode<std::int64_t>*                   idx;
        std::vector<const TypedDagNode<valueType>*>         elems;
    };
}


template <class valueType>
RevBayesCore::SelectVectorElement<valueType>::SelectVectorElement(const TypedDagNode<std::int64_t>* index, const std::vector<const TypedDagNode<valueType>*>& elements) : TypedFunction<valueType>( new valueType() ),
    idx( index ),
    elems( elements )
{
    if ( idx == nullptr )
    {
        throw RbException( "SelectVectorElement: index must not be NULL." );
    }
        
    if ( elems.empty() )
    {
        throw RbException( "SelectVectorElement: elements must not be empty." );
    }

    // DAG wiring: depend on the index and *all* candidate elements. This keeps the DAG structure static while allowing the output to update
    // when the selected index value changes.
    this->addParameter( idx );

    for ( const auto* e : elems )
    {
        this->addParameter( e );
    }

    update();
}


template <class valueType>
RevBayesCore::SelectVectorElement<valueType>::~SelectVectorElement( void )
{
    // We don't delete parameters because they are owned by the model/DAG.
}


template <class valueType>
RevBayesCore::SelectVectorElement<valueType>* SelectVectorElement<valueType>::clone(void) const
{
    return new SelectVectorElement<valueType>( *this );
}


template <class valueType>
void RevBayesCore::SelectVectorElement<valueType>::update(void)
{
    // Natural indices in Rev scripts are 1-based; elems is stored 0-based.
    const std::int64_t k = idx->getValue();
    if ( k < 1 || k > static_cast<std::int64_t>( elems.size() ) )
    {
        throw RbException("Index out of range. Index = " + std::to_string( k ) + ", vector size = " + std::to_string( elems.size() ) + "." );
    }

    const std::size_t sel = static_cast<std::size_t>( k - 1 );
    *this->value = elems[ sel ]->getValue();
}


template <class valueType>
void RevBayesCore::SelectVectorElement<valueType>::swapParameterInternal(const DagNode* oldP, const DagNode* newP)
{
    if ( oldP == idx )
    {
        idx = static_cast<const TypedDagNode<std::int64_t>*>( newP );
    }

    for ( size_t i = 0; i < elems.size(); ++i )
    {
        if ( oldP == elems[i] )
        {
            elems[i] = static_cast<const TypedDagNode<valueType>*>( newP );
        }
    }
}


#endif
