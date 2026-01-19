#include "BirthDeathJeffreysDistribution.h"
#include "BD_FIM_2x2.hpp"

#include <cmath>
#include <string>

#include "RandomNumberFactory.h"
#include "RbConstants.h"
#include "RbException.h"
#include "StochasticNode.h"
#include "Cloner.h"
#include "RbVectorImpl.h"
#include "RevVariable.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { class RandomNumberGenerator; }

using namespace RevBayesCore;

/* Multivariate Normal Distribution Constructor
 *
 * @param m A vector of doubles for the location parameters
 * @param cov a matrix of reals for that represents the covariance matrix
 */

BirthDeathJeffreysDistribution::BirthDeathJeffreysDistribution(const TypedDagNode<double> *ra,
		const std::string &cdt,
        bool uo,
		const TypedDagNode<double> *r,
		double l) :
    TypedDistribution< RbVector<double> >( new RbVector<double>() ),
    process_age( ra ),
    use_origin( uo ),
    condition( cdt ),
	rho( r ),
	limit( l )
{

    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( process_age );
    addParameter( rho );

    redrawValue();
}


BirthDeathJeffreysDistribution::~BirthDeathJeffreysDistribution( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



BirthDeathJeffreysDistribution* BirthDeathJeffreysDistribution::clone( void ) const
{
    return new BirthDeathJeffreysDistribution( *this );
}

double BirthDeathJeffreysDistribution::computeLnProbability( void )
{

	// get variables
	double lambda = (*value)[0];
	double mu     = (*value)[1];
	double p      = rho->getValue();
	double age    = process_age->getValue();

	if (lambda > limit || mu > limit) {
		return RbConstants::Double::neginf;
	}

	// compute the 2x2 FIM
	std::array<std::array<double, 2>, 2> fim = BDGrads::TwoByTwo::BD_FIM(lambda, mu, p, age);
	double a = fim[0][0];
	double b = fim[0][1];
	double c = fim[1][1];
	if (!use_origin) {
		a *= 2.0;
		b *= 2.0;
		c *= 2.0;
	}
	double s = a - (b * b) / c;  // compute Schur complement
	if (s <= 0.0) {
		// numerically singular
		return RbConstants::Double::neginf;
	}

	double log_det = std::log(c) + std::log(s);
	double ln_p    = 0.5 * log_det;

	return ln_p;

}



void BirthDeathJeffreysDistribution::redrawValue( void )
{

	// base initial value on age of process
	double age    = process_age->getValue();

	// initialize with some fixed starting values
	*value = RbVector<double>(2);
//	(*value)[0] = limit * 0.5;
//	(*value)[1] = limit * 0.25;
	(*value)[0] = 1.0 / age;
	(*value)[1] = 0.5 / age;

}

/** Swap a parameter of the distribution */
void BirthDeathJeffreysDistribution::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == process_age)
    {
    	process_age = static_cast<const TypedDagNode< double >* >( newP );
    }
    if (oldP == rho)
    {
    	rho = static_cast<const TypedDagNode< double >* >( newP );
    }
    
}


