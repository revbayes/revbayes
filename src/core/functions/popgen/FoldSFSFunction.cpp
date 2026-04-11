#include "FoldSFSFunction.h"

#include "RbException.h"
#include "TypedDagNode.h"
#include "Cloner.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;


FoldSFSFunction::FoldSFSFunction(const TypedDagNode< RbVector<std::int64_t> >* u) :
    TypedFunction< RbVector<std::int64_t> >( new RbVector<std::int64_t>() ),
    unfolded( u )
{
    addParameter( unfolded );
    update();
}


FoldSFSFunction::~FoldSFSFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


FoldSFSFunction* FoldSFSFunction::clone( void ) const
{
    return new FoldSFSFunction( *this );
}


void FoldSFSFunction::update( void )
{
    const RbVector<std::int64_t>& uf = unfolded->getValue();
    size_t len = uf.size();

    if ( len < 2 )
    {
        throw RbException() << "fnFoldSFS: site frequency spectrum must have at least 2 entries (got " << len << ").";
    }

    size_t n = len - 1;                   // number of individuals
    size_t folded_len = n / 2 + 1;        // floor(n/2) + 1

    value->clear();
    value->resize( folded_len, 0.0 );

    for ( size_t i = 0; i <= n / 2; ++i )
    {
        if ( i != n - i )
        {
            (*value)[i] = uf[i] + uf[n - i];
        }
        else
        {
            (*value)[i] = uf[i];
        }
    }
}


void FoldSFSFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    if ( oldP == unfolded )
    {
        unfolded = static_cast< const TypedDagNode< RbVector<std::int64_t> >* >( newP );
    }
}
