#ifndef DirichletDistribution_H
#define DirichletDistribution_H

#include "Simplex.h"
#include "TypedDistribution.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class RbVector;
template <class valueType> class TypedDagNode;
    
/**
* @brief Dirichlet distribution class.
*
* The Dirichlet probability distribution is a multivariate distribution for random variables @f${\mathbf x} = (x_1, x_2, \ldots, x_k)@f$ which form a @f$k-1@f$ dimension simplex
* (that is, @f$0 \leq x_i \leq 1@f$ and @f$ \sum_{i=1}^k x_i = 1@f$). The probability density for the Dirichlet distribution is
* @f[ f({\mathbf x} | \alpha) = {\Gamma(\alpha_0) \over \prod_{i=1}^k \Gamma(\alpha_i)} \prod_{i=1}^k x_i^{\alpha_i-1}@f]
* where @f$ \alpha = (\alpha_1, \alpha_2,\ldots,\alpha_k)@f$ are the concentration parameters (@f$ \alpha_i \geq 0@f$) and @f$\alpha_0 = \sum_{i=1}^k \alpha_i@f$.
* A Dirichlet distribution with all @f$\alpha_i = 1@f$ gives equal probabiity density to all combinations of the random variable and is called a "Flat Dirichlet" distribution.
*
*/
    class DirichletDistribution : public TypedDistribution< Simplex > {
        
        public:
                                                        DirichletDistribution(const TypedDagNode< RbVector<double> > *l);
            virtual                                    ~DirichletDistribution(void);                                                //!< Virtual destructor
            DirichletDistribution*                      clone(void) const;                                                          //!< Create an independent clone
            double                                      computeLnProbability(void);                                                 //!< Natural log of the probability density
            void                                        redrawValue(void);
            
        protected:
            void                                        swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Swap a parameter
            
        private:
            const TypedDagNode< RbVector<double> >*     alpha;
    };
}

#endif
