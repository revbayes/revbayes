#include "DirichletDistributionSorted.h"

#include "DistributionDirichlet.h"
#include "RandomNumberFactory.h"
#include "RbVector.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/** DirichletDistributionSorted Constructor
 * @param a vector of concentration parameters
 */
DirichletDistributionSorted::DirichletDistributionSorted(const TypedDagNode< RbVector<double> > *a) : TypedDistribution< Simplex >( new Simplex() ),
    alpha( a )
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( alpha );

    // TODO - actually sort here later
    *value = RbStatistics::Dirichlet::rv(alpha->getValue(), *GLOBAL_RNG);
    value->sort(false);

    //std::stringstream ss;
    //for (size_t i = 0; i < value->size(); ++i)
    //  value->printElement(ss,i,"\t",-1,true);
    //std::cout << ss.str();
}


DirichletDistributionSorted::~DirichletDistributionSorted( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



DirichletDistributionSorted* DirichletDistributionSorted::clone( void ) const
{
    return new DirichletDistributionSorted( *this );
}


double DirichletDistributionSorted::computeLnProbability( void )
{
  //this->operator[](i)
  for (size_t i = 1; i < value->size(); ++i )
    if (value->operator[](i-1) < value->operator[](i)) { return RbConstants::Double::neginf; }

    return RbStatistics::Dirichlet::lnPdf(alpha->getValue(), *value);
}


void DirichletDistributionSorted::redrawValue( void )
{
    *value = RbStatistics::Dirichlet::rv(alpha->getValue(), *GLOBAL_RNG);
    value->sort(false);
}

/** Swap a parameter of the distribution */
void DirichletDistributionSorted::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{

    if (oldP == alpha)
    {
        alpha = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
}
