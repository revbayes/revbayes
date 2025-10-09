#include "WishartDistribution.h"

#include <cstdlib>
#include <iostream>

#include "RandomNumberFactory.h"
#include "DistributionWishart.h"
#include "Cloneable.h"

namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { class RandomNumberGenerator; }

using namespace RevBayesCore;

/**
 * Default Constructor for the Wishart Distribution
 * @param inomega0 A scale matrix of positive real numbers
 * @param indf a positive std::int64_t number for the degrees of freedom
 *
 */
WishartDistribution::WishartDistribution(const TypedDagNode<MatrixReal> *inomega0, const TypedDagNode<std::int64_t>* indf)  :
TypedDistribution<RevBayesCore::MatrixReal>(new MatrixReal(inomega0->getValue().getDim())),
    omega0(inomega0),
    kappa(NULL),
    df(indf),
    dim( NULL )
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( inomega0 );
    addParameter( df );
    
    redrawValue();
}


/**
 * Constructor for the Wishart Distribution
 * @param inkappa a value for all diagonal elements of the scaling matrix
 * @param indim the number of dimensions for the scaling matrix
 * @param indf a positive std::int64_t number for the degrees of freedom
 *
 */
WishartDistribution::WishartDistribution(const TypedDagNode<std::int64_t>* indim, const TypedDagNode<double> *inkappa, const TypedDagNode<std::int64_t>* indf)  :
TypedDistribution<RevBayesCore::MatrixReal>(new MatrixReal( size_t(indim->getValue()) )),
    omega0(NULL),
    kappa(inkappa),
    df(indf),
    dim(indim)
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( kappa );
    addParameter( df );
    addParameter( dim );
    
    redrawValue();
}

WishartDistribution* WishartDistribution::clone(void) const
{

    return new WishartDistribution(*this);
}


/** Swap a parameter of the distribution */
void WishartDistribution::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == omega0)
    {
        std::cerr << "omega0??\n";
        exit(1);
//        omega0 = static_cast<const TypedDagNode<MatrixReal>* >( newP );
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


double WishartDistribution::computeLnProbability(void)
{
    
    double ret = 0;
    
    if ( omega0 != NULL )
    {
        ret = RbStatistics::Wishart::lnPdf(omega0->getValue(),df->getValue(),getValue());
    }
    else
    {
        ret = RbStatistics::Wishart::lnPdf(kappa->getValue(),df->getValue(),getValue());        
    }

    return ret;
}

void WishartDistribution::redrawValue(void)
{

    RandomNumberGenerator* rng = GLOBAL_RNG;

    if ( omega0 != NULL )
    {
        getValue() = RbStatistics::Wishart::rv(omega0->getValue(),df->getValue(), *rng);
    }
    else
    {
        getValue() = RbStatistics::Wishart::rv(kappa->getValue(),getValue().getDim(),df->getValue(), *rng);        
    }

}
