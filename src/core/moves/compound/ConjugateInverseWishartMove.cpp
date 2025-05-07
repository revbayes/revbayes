#include <cstddef>
#include <iosfwd>
#include <vector>

#include "DistributionInverseWishart.h"
#include "InverseWishartDistribution.h"
#include "MultivariateNormalDistribution.h"
#include "RandomNumberFactory.h"
#include "ConjugateInverseWishartMove.h"
#include "AbstractGibbsMove.h"
#include "Cloneable.h"
#include "DagNode.h"
#include "Distribution.h"
#include "MatrixReal.h"
#include "RbException.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "StochasticNode.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore { class RandomNumberGenerator; }

using namespace RevBayesCore;

ConjugateInverseWishartMove::ConjugateInverseWishartMove(StochasticNode<MatrixReal>* s, double w) : AbstractGibbsMove(w),
    sigma( s )
{
    InverseWishartDistribution* dist = dynamic_cast<InverseWishartDistribution *>( &sigma->getDistribution() );
    if (dist == NULL)
    {
        throw(RbException("ConjugateInverseWishart move only works for parameters with an Inverse Wishart prior distribution"));
    }

    addNode( sigma );
    
}


/**
 * Basic destructor doing nothing.
 */
ConjugateInverseWishartMove::~ConjugateInverseWishartMove( void )
{
    
    // nothing to do
    // the base class does the work
    
}




/* Clone object */
ConjugateInverseWishartMove* ConjugateInverseWishartMove::clone( void ) const
{
    
    return new ConjugateInverseWishartMove( *this );
}



const std::string& ConjugateInverseWishartMove::getMoveName( void ) const
{
    
    static std::string name = "ConjugateInverseWishart";
    
    return name;
}


/** Perform the move */
void ConjugateInverseWishartMove::performGibbsMove( void )
{
    InverseWishartDistribution* inverse_wishart = dynamic_cast<InverseWishartDistribution *>( &sigma->getDistribution() );

    int df = inverse_wishart->getDF();
    
    const TypedDagNode<MatrixReal>* sigma0             = inverse_wishart->getSigma0();
    const TypedDagNode<RbVector<double> >* kappaVector = inverse_wishart->getKappaVector();
    const TypedDagNode<double>* kappa                  = inverse_wishart->getKappa();

    // Get random number generator
    RandomNumberGenerator* rng = GLOBAL_RNG;
    
    size_t dim = sigma->getValue().getDim();
    
    MatrixReal A = MatrixReal(dim);

    size_t num_children = 0;
    
    const std::vector<DagNode*>& c = sigma->getChildren();
    for (std::vector<DagNode*>::const_iterator it = c.begin(); it != c.end(); ++it)
    {
        if ( (*it)->isStochastic() )
        {
            MultivariateNormalDistribution* dist = dynamic_cast<MultivariateNormalDistribution *>( &(*it)->getDistribution() );
            if ( dist != NULL )
            {
                A += dist->computeContrasts();
                num_children++;
            }
        }
    }
    
    if ( num_children == 0 )
    {
        throw(RbException("Inverse Wishart parameter for ConjugateInverseWishart move has no Multivariate Normal children"));
    }

    if (sigma0 != NULL)
    {
        A += sigma0->getValue();
    }
    else if ( kappaVector != NULL)
    {
        for (size_t i=0; i<dim; ++i)
        {
            A[i][i] += kappaVector->getValue()[i];
        }
    }
    else
    {
        for (size_t i=0; i<dim; ++i)
        {
            A[i][i] += kappa->getValue();
        }
    }
    
    // resample sigma based on new sufficient statistics
    sigma->setValue( RbStatistics::InverseWishart::rv(A, num_children + df, *rng).clone() );

    sigma->touch();
    sigma->keep();
}


void ConjugateInverseWishartMove::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    if (oldN == sigma)
    {
        sigma = static_cast<StochasticNode<MatrixReal>* >(newN) ;
    }
    
}


