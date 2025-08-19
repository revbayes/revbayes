/**
 * @brief
 * This file contains the declaration of the Multivariate Normally distributed random variable class.
 * This class is derived from the stochastic node and each instance will represent a random variable
 * from a normal distribution in the model graph.
 * The Multivariate Normal distribution represents a family of distributions
 * defined on the real numbers. The Multivariate Normal distribution has a varying number of parameters depending on the dimensionality of the distribution:
 *
 *@param Mu A location parameter for each dimension
 *@param Sigma a variance-covariance matrix
 */



#ifndef MultivariateNormalDistribution_H
#define MultivariateNormalDistribution_H

#include <cstddef>
#include <vector>
#include <iosfwd>

#include "RbVector.h"
#include "TypedDistribution.h"
#include "MatrixReal.h"
#include "RevPtr.h"

namespace RevLanguage { class RevVariable; }

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;
    
    class MultivariateNormalDistribution : public TypedDistribution< RbVector<double> > {
        
    public:
        MultivariateNormalDistribution(const TypedDagNode< RbVector<double> > *m, const TypedDagNode<MatrixReal>* cov, const TypedDagNode<MatrixReal>* prec, const TypedDagNode<double>* sc);
        virtual                                                         ~MultivariateNormalDistribution(void);                                                //!< Virtual destructor
        
        
        // public member functions
        MultivariateNormalDistribution*                                 clone(void) const;
        
        MatrixReal                                                      computeContrasts(void);
        RevLanguage::RevPtr<RevLanguage::RevVariable>                   executeProcedure(const std::string &n, const std::vector<DagNode*> args, bool &f);  //!< execute the procedure
        void                                                            clampAt(size_t i, double v);
        double                                                          computeLnProbability(void);
        void                                                            redrawValue(void);

    protected:
        // Parameter management functions
        void                                                            swapParameterInternal(const DagNode *oldP, const DagNode *newP);                    //!< Swap a parameter
        
    private:
        
        // members
        const TypedDagNode< RbVector<double> >*                         mean;
        const TypedDagNode< MatrixReal >*                               covariance;
        const TypedDagNode< MatrixReal >*                               precision;
        const TypedDagNode< double >*                                   scale;

        std::vector<bool>                                               observed;
        std::vector<double>                                             observations;
        
    };
    
}

#endif
