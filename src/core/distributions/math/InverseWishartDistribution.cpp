

#include "InverseWishartDistribution.h"

#include <cstddef>
#include <cstdint>

#include "RandomNumberFactory.h"
#include "DistributionInverseWishart.h"
#include "Cloneable.h"
#include "RbException.h"
#include "RbVector.h"

namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { class RandomNumberGenerator; }

using namespace RevBayesCore;


/**
 * Default Constructor for the Inverse Wishart Distribution
 * @param insigma0 A scale matrix of positive real numbers
 * @param indf a positive std::int64_t number for the degrees of freedom
 *
 */
InverseWishartDistribution::InverseWishartDistribution(const TypedDagNode<MatrixReal> *insigma0, const TypedDagNode<std::int64_t>* indf)  :
TypedDistribution<RevBayesCore::MatrixReal>(new MatrixReal(insigma0->getValue().getDim())),
sigma0(insigma0),
kappaVector(NULL),
kappa(NULL),
df(indf),
dim( NULL )  {
    
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( insigma0 );
    addParameter( df );
    
    redrawValue();
}


/**
 * Constructor for the Inverse Wishart Distribution
 * @param inkappaVector A vector for the diagonal of the scaling matrix
 * @param indf a positive std::int64_t number for the degrees of freedom
 *
 *@note For this parameterization the scaling matrix is calculated as:
 *@note sigma0 = Diagonal(kappaVector)
 *
 *
 */
InverseWishartDistribution::InverseWishartDistribution(const TypedDagNode<RbVector<double> > *inkappaVector, const TypedDagNode<std::int64_t>* indf)  :
TypedDistribution<RevBayesCore::MatrixReal>(new MatrixReal( inkappaVector->getValue().size()) ),
    sigma0(NULL),
    kappaVector(inkappaVector),
    kappa(NULL),
    df(indf),
    dim( NULL )
{
    
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( kappaVector );
    addParameter( df );
    
    redrawValue();
}

/**
 * Constructor for the Inverse Wishart Distribution
 * @param inkappa A value for the diagonal of the scaling matrix
 * @param indim The number of dimensions on the scaling matrix
 * @param indf a positive std::int64_t number for the degrees of freedom
 *
 *
 *@note For this parameterization the scaling matrix is calculated as:
 *@note sigma0=kappa * Identitymatrix
 *@note Where the identity matrix has indim dimensions
 *
 *
 */
InverseWishartDistribution::InverseWishartDistribution(const TypedDagNode<std::int64_t>* indim, const TypedDagNode<double> *inkappa, const TypedDagNode<std::int64_t>* indf)  :
TypedDistribution<RevBayesCore::MatrixReal>(new MatrixReal( size_t(indim->getValue()) )),
    sigma0(NULL),
    kappaVector(NULL),
    kappa(inkappa),
    df(indf),
    dim(indim)
{
        
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( inkappa );
    addParameter( indf );
    addParameter( indim );
    
    redrawValue();
}

InverseWishartDistribution* InverseWishartDistribution::clone(void) const   {

    return new InverseWishartDistribution(*this);
}


void InverseWishartDistribution::swapParameterInternal(const DagNode *oldP, const DagNode *newP) {
    
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


double InverseWishartDistribution::computeLnProbability(void)
{
    
    double ret = 0;
    
    if ( sigma0 != NULL )
    {
        ret = RbStatistics::InverseWishart::lnPdf(sigma0->getValue(),df->getValue(),getValue());
    }
    else if ( kappaVector != NULL )
    {
        ret = RbStatistics::InverseWishart::lnPdf(kappaVector->getValue(),df->getValue(),getValue());        
    }
    else if ( kappa != NULL )
    {
        ret = RbStatistics::InverseWishart::lnPdf(kappa->getValue(),df->getValue(),getValue());        
    }
    else
    {
        throw RbException("error in inverse wishart: no parameter.");
    }
    
    return ret;
}

void InverseWishartDistribution::redrawValue(void)
{

    RandomNumberGenerator* rng = GLOBAL_RNG;

    if ( sigma0 != NULL )
    {
        setValue( RbStatistics::InverseWishart::rv(sigma0->getValue(),df->getValue(), *rng).clone() );
    }
    else if ( kappaVector != NULL )
    {
        setValue( RbStatistics::InverseWishart::rv(kappaVector->getValue(),df->getValue(), *rng).clone() );
    }
    else if ( kappa != NULL )
    {
        setValue( RbStatistics::InverseWishart::rv(kappa->getValue(),getValue().getDim(),df->getValue(), *rng).clone() );
    }
    else
    {
        throw RbException("error in inverse wishart: no parameter\n");
    }

}
