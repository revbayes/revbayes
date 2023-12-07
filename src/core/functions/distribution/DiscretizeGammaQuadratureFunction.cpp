#include "DiscretizeGammaQuadratureFunction.h"

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



RevBayesCore::DiscretizeGammaQuadratureFunction::DiscretizeGammaQuadratureFunction(const TypedDagNode<double> *s, const TypedDagNode<double> *r, const TypedDagNode<long> *nc) : TypedFunction< MatrixReal >( new MatrixReal(2, nc->getValue(), 1.0) ),
    shape( s ),
    rate( r ),
    numCats(nc)
{
    
    addParameter( shape );
    addParameter( rate );
    addParameter( numCats );
        
}



RevBayesCore::DiscretizeGammaQuadratureFunction* RevBayesCore::DiscretizeGammaQuadratureFunction::clone( void ) const
{
    
    return new DiscretizeGammaQuadratureFunction(*this);
}


void RevBayesCore::DiscretizeGammaQuadratureFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP) {
    
    if (oldP == shape) 
    {
        shape = static_cast<const TypedDagNode<double>* >( newP );
    }
    
    if (oldP == rate)
    {
        rate = static_cast<const TypedDagNode<double>* >( newP );
    }
    
    if (oldP == numCats)
    {
        numCats = static_cast<const TypedDagNode<long>* >( newP );
    }
    
}

void RevBayesCore::DiscretizeGammaQuadratureFunction::update( void )
{
    
    const double a = shape->getValue();
    const double b = rate->getValue();
    const size_t nCats = static_cast<size_t>(numCats->getValue());
    
    // Fernandes and Atchley, 2006
    // calculate the recursion coefficients for the Gamma distribution
    double lambda = 1.0/b;
    double lambda2 = lambda * lambda;
    std::vector<double> a_prime(nCats, 1.0);
    std::vector<double> b_prime(nCats, 1.0);
    
    for (int i = 0; i < nCats; ++i)
    {
        a_prime[i] = (a + 2 * i) * lambda;
        b_prime[i] = (i * (a + i - 1)) * lambda2;
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
