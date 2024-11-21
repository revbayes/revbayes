#include "DiscretizeBetaQuadratureFunction.h"

#include <cstddef>
#include <cstdint>
#include <cmath>
#include <map>
#include <utility>
#include <vector>

#include "EigenSystem.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }



RevBayesCore::DiscretizeBetaQuadratureFunction::DiscretizeBetaQuadratureFunction(const TypedDagNode<double> *a, const TypedDagNode<double> *b, const TypedDagNode<std::int64_t> *nc) : TypedFunction< MatrixReal >( new MatrixReal(2, nc->getValue(), 1.0) ),
alpha( a ),
beta( b ),
numCats(nc)
{
    
    addParameter( alpha );
    addParameter( beta );
    addParameter( numCats );
    
}



RevBayesCore::DiscretizeBetaQuadratureFunction* RevBayesCore::DiscretizeBetaQuadratureFunction::clone( void ) const
{
    
    return new DiscretizeBetaQuadratureFunction(*this);
}


void RevBayesCore::DiscretizeBetaQuadratureFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP) {
    
    if (oldP == alpha)
    {
        alpha = static_cast<const TypedDagNode<double>* >( newP );
    }
    
    if (oldP == beta)
    {
        beta = static_cast<const TypedDagNode<double>* >( newP );
    }
    
    if (oldP == numCats)
    {
        numCats = static_cast<const TypedDagNode<std::int64_t>* >( newP );
    }
    
}

void RevBayesCore::DiscretizeBetaQuadratureFunction::update( void )
{
    
    double a = alpha->getValue();
    double b = beta->getValue();
    const size_t nCats = static_cast<size_t>(numCats->getValue());
    
    // Fernandes and Atchley, 2006
    // calculate the recursion coefficients for the Beta distribution
    double gamma = a + b;
    std::vector<double> a_prime(nCats, 1.0);
    std::vector<double> b_prime(nCats, 1.0);
    
    for (int i = 0; i < nCats; ++i) {
        if (i == 0 && a == 1.0 && b == 1.0) {
            a_prime[i] = 0.5;
        } else {
            a_prime[i] = (a*gamma+(2*i-2)*a+2*i*b+i*(2*i-2))/((gamma+2*i)*(gamma+2*i-2));
        }
        b_prime[i] = (i*(gamma+i-2)*(a+i-1)*(b+i-1))/((gamma+2*i-1)*pow((gamma+2*i-2),2)*(gamma+2*i-3));
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
