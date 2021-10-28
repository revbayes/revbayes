#include "SortedDirichletDistribution.h"
#include "DistributionDirichlet.h"
#include "RandomNumberFactory.h"
#include "RbVector.h"
#include "TypedDagNode.h"

#include <algorithm>


namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/** SortedDirichletDistribution Constructor
 * @param a vector of concentration parameters
 */
SortedDirichletDistribution::SortedDirichletDistribution( const TypedDagNode< RbVector<double> > *a )
    : TypedDistribution< Simplex >( new Simplex() ), alpha( a )
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( alpha );

    // sort the value
    *value = RbStatistics::Dirichlet::rv(alpha->getValue(), *GLOBAL_RNG);
    value->sort(false);
}


/** SortedDirichletDistribution Destructor */
SortedDirichletDistribution::~SortedDirichletDistribution( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


/** Clone the SortedDirichletDistribution */
SortedDirichletDistribution* SortedDirichletDistribution::clone( void ) const
{
    return new SortedDirichletDistribution( *this );
}


/** Compute log probability */
double SortedDirichletDistribution::computeLnProbability( void )
{

    double cur = value->operator[](1);
    for (size_t i = 2; i < value->size(); ++i )
    {
        if ( value->operator[](i) > cur )
        {
            return RbConstants::Double::neginf;
        }
        cur = value->operator[](i);
    }

    return RbStatistics::Dirichlet::lnPdf(alpha->getValue(), *value);
}


/** Redraw value from distribution */
void SortedDirichletDistribution::redrawValue( void )
{
    *value = RbStatistics::Dirichlet::rv(alpha->getValue(), *GLOBAL_RNG);
    value->sort(false);

}


/** Swap a parameter of the distribution
 * @param oldP old value
 * @param newP new value
 */
void SortedDirichletDistribution::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    if (oldP == alpha)
    {
        alpha = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
}
