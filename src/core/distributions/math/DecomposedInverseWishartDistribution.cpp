#include "DecomposedInverseWishartDistribution.h"

#include <cstddef>
#include <cstdint>

#include "RandomNumberFactory.h"
#include "DistributionDecomposedInverseWishart.h"
#include "Cloneable.h"
#include "RbException.h"
#include "RbVector.h"

namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { class RandomNumberGenerator; }

using namespace RevBayesCore;

DecomposedInverseWishartDistribution::DecomposedInverseWishartDistribution(const TypedDagNode<MatrixReal> *insigma0, const TypedDagNode<std::int64_t>* indf)  :
TypedDistribution<RevBayesCore::MatrixReal>(new MatrixReal(insigma0->getValue().getDim())),
sigma0(insigma0),
kappaVector(NULL),
kappa(NULL),
df(indf),
dim( NULL )
{
    
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( insigma0 );
    addParameter( df );
    
    redrawValue();
}

DecomposedInverseWishartDistribution::DecomposedInverseWishartDistribution(const TypedDagNode<RbVector<double> > *inkappaVector, const TypedDagNode<std::int64_t>* indf)  :
TypedDistribution<RevBayesCore::MatrixReal>(new MatrixReal( inkappaVector->getValue().size()) ),
    sigma0(NULL),
    kappaVector(inkappaVector),
    kappa(NULL),
    df(indf),
    dim( NULL ) {
    
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( kappaVector );
    addParameter( df );
    
    redrawValue();
}

DecomposedInverseWishartDistribution::DecomposedInverseWishartDistribution(const TypedDagNode<std::int64_t>* indim, const TypedDagNode<double> *inkappa, const TypedDagNode<std::int64_t>* indf)  :
TypedDistribution<RevBayesCore::MatrixReal>(new MatrixReal( size_t(indim->getValue()) )),
    sigma0(NULL),
    kappaVector(NULL),
    kappa(inkappa),
    df(indf),
    dim(indim) {

    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( inkappa );
    addParameter( indf );
    addParameter( indim );
    redrawValue();
}

DecomposedInverseWishartDistribution* DecomposedInverseWishartDistribution::clone(void) const
{

    return new DecomposedInverseWishartDistribution(*this);
}


void DecomposedInverseWishartDistribution::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{

    if (oldP == sigma0)
    {
        sigma0 = static_cast<const TypedDagNode<MatrixReal>* >( newP );
    }
    if (oldP == kappaVector)
    {
        kappaVector = static_cast<const TypedDagNode<RbVector<double> >* >(newP);
    }
    if (oldP == kappa)
    {
        kappa = static_cast<const TypedDagNode<double>* >(newP);
    }
    if (oldP == dim)
    {
        dim = static_cast<const TypedDagNode<std::int64_t>* >(newP);
    }
    if (oldP == df)
    {
        df = static_cast<const TypedDagNode<std::int64_t>* >(newP);
    }
}


double DecomposedInverseWishartDistribution::computeLnProbability(void)
{
    
    double ret = RbStatistics::DecomposedInverseWishart::lnPdf(df->getValue(), getValue());
    return ret;
}

void DecomposedInverseWishartDistribution::redrawValue(void)
{

    RandomNumberGenerator* rng = GLOBAL_RNG;
    if ( kappa != NULL )
    {
        setValue( RbStatistics::DecomposedInverseWishart::rv(getValue().getDim(), df->getValue(), *rng).clone() );
    }
    else
    {
        throw RbException("error in decomposed inverse wishart: no parameter\n");
    }

}
