#include "OrderedElementSlideProposal.h"

#include <cstddef>
#include <ostream>
#include <vector>

#include "MatrixReal.h"
#include "OrderedVector.h"
#include "OrderedVectorView.h"
#include "Proposal.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbException.h"
#include "StochasticNode.h"

using namespace RevBayesCore;


OrderedElementSlideProposal::OrderedElementSlideProposal( StochasticNode<MatrixReal> *n ) : Proposal(),
    variable( n ),
    stored_vector( 0 ),
    stored_entry( 0 ),
    stored_value( 0.0 ),
    failed( false )
{
    addNode( variable );

    if ( dynamic_cast<OrderedVectorView*>( &variable->getValue() ) == NULL )
    {
        throw RbException("mvOrderedElementSlide needs a value whose parts are ordered vectors.");
    }
}


void OrderedElementSlideProposal::cleanProposal( void )
{
}


OrderedElementSlideProposal* OrderedElementSlideProposal::clone( void ) const
{
    return new OrderedElementSlideProposal( *this );
}


const std::string& OrderedElementSlideProposal::getProposalName( void ) const
{
    static std::string name = "OrderedElementSlide";

    return name;
}


double OrderedElementSlideProposal::getProposalTuningParameter( void ) const
{
    // the neighbours set the window, so there is nothing to tune
    return RbConstants::Double::nan;
}


double OrderedElementSlideProposal::doProposal( void )
{
    RandomNumberGenerator* rng = GLOBAL_RNG;

    OrderedVectorView& value = dynamic_cast<OrderedVectorView&>( variable->getValue() );

    failed = true;

    size_t n = value.numOrderedVectors();
    if ( n == 0 ) return RbConstants::Double::neginf;

    // start from a random vector and take the first with somewhere to move, so a value whose
    // vectors are mostly pinned still proposes rather than failing at a rate we would have to tune
    size_t offset = size_t( rng->uniform01() * n );

    for (size_t s = 0; s < n; ++s)
    {
        size_t k = ( offset + s ) % n;
        OrderedVector v = value.getOrderedVector( k );

        std::vector<size_t> movable;
        for (size_t j = 0; j < v.size(); ++j)
        {
            if ( v.isFree(j) ) movable.push_back( j );
        }
        if ( movable.empty() ) continue;

        // an entry already outside its neighbours means something upstream is wrong; refuse to
        // act rather than quietly sort it, which would be a one-way correction, not a proposal
        if ( v.isOrdered() == false ) return RbConstants::Double::neginf;

        size_t j = movable[ size_t( rng->uniform01() * movable.size() ) ];

        double lo = v.lowerBound(j);
        double hi = v.upperBound(j);
        if ( hi <= lo ) return RbConstants::Double::neginf;

        stored_vector = k;
        stored_entry  = j;
        stored_value  = v[j];
        failed        = false;

        v.set( j, lo + rng->uniform01() * ( hi - lo ) );

        variable->addTouchedElementIndex( k * v.size() + j );

        // uniform on an interval the reverse move draws from identically
        return 0.0;
    }

    return RbConstants::Double::neginf;
}


void OrderedElementSlideProposal::prepareProposal( void )
{
}


void OrderedElementSlideProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
}


void OrderedElementSlideProposal::undoProposal( void )
{
    if ( failed == true ) return;

    OrderedVectorView& value = dynamic_cast<OrderedVectorView&>( variable->getValue() );
    OrderedVector v = value.getOrderedVector( stored_vector );
    v.set( stored_entry, stored_value );
}


void OrderedElementSlideProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    variable = static_cast<StochasticNode<MatrixReal>* >(newN);
}


void OrderedElementSlideProposal::setProposalTuningParameter(double tp)
{
}


void OrderedElementSlideProposal::tune( double rate )
{
}
