#ifndef WishartDistribution_H
#define WishartDistribution_H

#include "TypedDistribution.h"
#include "MatrixReal.h"
#include "TypedDagNode.h"

namespace RevBayesCore {
class DagNode;


/**
 * @brief Wishart distribution class.
 *
 * The Gamma distribution represents a family of distributions
 * defined on the positive real numbers. The Wishart distribution is a multivariate generalization of the Gamma distribution.
 * The Wishart distribution has 2 parameters:
 *
 * @param V A square scaling matrix.
 * @param Df The degrees of freedom
 * Instances of this class can be associated to stochastic variables.
 */


    
    class WishartDistribution : public TypedDistribution<MatrixReal>   {
        
    public:
        
        // inverse Wishart distribution of parameter sigma0 et df degrees of freedom
        WishartDistribution(const TypedDagNode<MatrixReal> *insigma0, const TypedDagNode<long>* indf);
        // specialized version: inverse Wishart of parameter sigma0=kappa*I and df degrees of freedom
        WishartDistribution(const TypedDagNode<long>* indim, const TypedDagNode<double> *inkappa, const TypedDagNode<long>* indf);
        
        virtual                                            ~WishartDistribution(void) {}
        
        // public member functions

        WishartDistribution*                                clone(void) const;                                                          //!< Create an independent clone
        double                                              computeLnProbability(void);
        void                                                redrawValue(void);
        
        long                                                getDF() const {return df->getValue();}
        
        
    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Swap a parameter
        
    private:

        // members
        
        const TypedDagNode<MatrixReal>*                     omega0;     //!<The scaling matrix
        const TypedDagNode<double>*                         kappa;      //!<A value for the diagonal elements of the scaling matrix
        const TypedDagNode<long>*                           df;         //!<The degrees of freedom
        const TypedDagNode<long>*                           dim;        //!<The number of dimensions of the scaling matrix
        
    };
    
}


#endif /* defined(WishartDistribution_H) */
