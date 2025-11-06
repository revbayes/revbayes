#include "DiscretizeLognormalQuadratureFunction.h"

#include <cstddef>
#include <cmath>
#include <map>
#include <utility>
#include <vector>

#include "EigenSystem.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }



RevBayesCore::DiscretizeLognormalQuadratureFunction::DiscretizeLognormalQuadratureFunction(const TypedDagNode<double> *m, const TypedDagNode<double> *s, const TypedDagNode<std::int64_t> *nc) : TypedFunction< MatrixReal >( new MatrixReal(2, nc->getValue(), 1.0) ),
    mean( m ),
    sd( s ),
    numCats(nc)
{
    
    addParameter( mean );
    addParameter( sd );
    addParameter( numCats );
        
}



RevBayesCore::DiscretizeLognormalQuadratureFunction* RevBayesCore::DiscretizeLognormalQuadratureFunction::clone( void ) const
{
    
    return new DiscretizeLognormalQuadratureFunction(*this);
}


void RevBayesCore::DiscretizeLognormalQuadratureFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP) {
    
    if (oldP == mean) 
    {
        mean = static_cast<const TypedDagNode<double>* >( newP );
    }
    
    if (oldP == sd)
    {
        sd = static_cast<const TypedDagNode<double>* >( newP );
    }
    
    if (oldP == numCats)
    {
        numCats = static_cast<const TypedDagNode<std::int64_t>* >( newP );
    }
    
}

void RevBayesCore::DiscretizeLognormalQuadratureFunction::update( void )
{
    
    const double mu = mean->getValue();
    const double sigma = sd->getValue();
    const size_t nCats = static_cast<size_t>(numCats->getValue());
    
    // Fernandes and Atchley, 2006
    // calculate the recursion coefficients for the Lognormal distribution
    double lambda = exp(mu);
    double lambda2 = lambda * lambda;
    double zeta = exp(sigma * sigma);
    std::vector<double> a_prime(nCats, 1.0);
    std::vector<double> b_prime(nCats, 1.0);
    
    for (int i = 0; i < nCats; ++i)
    {
        a_prime[i] = pow(zeta, (2*i-1)/2.0) * (pow(zeta, i) * (zeta + 1.0) - 1.0) * lambda;
        b_prime[i] = pow(zeta, 3*i-2) * (pow(zeta, i) - 1.0) * lambda2;
    }
    b_prime[0] = 1.0;
    
    // fill in the Jacobi matrix with the recursion coefficients
    MatrixReal Jacobi(nCats);
    Jacobi[0][0] = a_prime[0];
    
    for (size_t i = 1; i < nCats; ++i)
    {
        Jacobi[i][i] = a_prime[i];
        Jacobi[i][i-1] = sqrt(b_prime[i]);
        Jacobi[i-1][i] = Jacobi[i][i-1];
    }
    
    // compute the eigenvalues and eigenvectors for the Jacobi matrix
    const EigenSystem& theEigenSystem = Jacobi.getEigenSystem();
    const std::vector<double>& eigenValue = theEigenSystem.getRealEigenvalues();
    const MatrixReal& eigenVector = theEigenSystem.getEigenvectors();
    
    // compute rates and weights according to the eigenvalues and eigenvectors, and sort them by rates in ascending order
    std::map<double, double> weights_map;
    for (size_t j=0; j<nCats; j++)
    {
        weights_map.insert(std::pair<double, double> (eigenValue[j], b_prime[0] * pow(eigenVector[0][j], 2)));
    }
    
    size_t j = 0;
    for ( std::map<double, double>::const_iterator it = weights_map.begin(); it != weights_map.end(); ++it)
    {
        (*value)[0][j] = it->first;
        (*value)[1][j] = it->second;
        j++;
    }
    
}
