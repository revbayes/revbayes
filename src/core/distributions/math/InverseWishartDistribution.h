#ifndef InverseWishartDistribution_H
#define	InverseWishartDistribution_H


#include "TypedDistribution.h"
#include "MatrixReal.h"
#include "TypedDagNode.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class RbVector;
    


/**
 * @brief Inverse Wishart distribution class.
 *
 * The Wishart distribution represents a family of distributions
 * defined on the real positive definite matrices. The Inverse Wishart distribution has 2 parameters:
 * @param V a scaling matrix
 * @param df the degrees of freedom.
 * Instances of this class can be associated to stochastic variables.
 *
 */

    class InverseWishartDistribution : public TypedDistribution<MatrixReal>   {
        
    public:
        
        InverseWishartDistribution(const TypedDagNode<MatrixReal> *insigma0, const TypedDagNode<std::int64_t>* indf);
        InverseWishartDistribution(const TypedDagNode<RbVector<double> > *inkappaVector, const TypedDagNode<std::int64_t>* indf);
        InverseWishartDistribution(const TypedDagNode<std::int64_t>* indim, const TypedDagNode<double> *inkappa, const TypedDagNode<std::int64_t>* indf);

        virtual                                            ~InverseWishartDistribution(void) {}
        
        // public member functions

        InverseWishartDistribution*                         clone(void) const;                                                          //!< Create an independent clone
        double                                              computeLnProbability(void);
        void                                                redrawValue(void);
        
        int                                                 getDF() const {return (int)df->getValue();}
        
        const TypedDagNode<MatrixReal>*                     getSigma0() const {return sigma0;}
        const TypedDagNode<RbVector<double> >*              getKappaVector() const {return kappaVector;}
        const TypedDagNode<double>*                         getKappa() const {return kappa;}

        
    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Swap a parameter

    private:

        const TypedDagNode<MatrixReal>*                     sigma0;             //!< a scaling matrix
        const TypedDagNode<RbVector<double> >*              kappaVector;        //!< A vector with the values of the diagonal of the scaling matrix
        const TypedDagNode<double>*                         kappa;              //!< A value with the value on the diagonal of the scaling matrix
        const TypedDagNode<std::int64_t>*                            df;                //!< The degrees of freedom
        const TypedDagNode<std::int64_t>*                            dim;               //!< The number of dimensions on the scaling matrix
                
    };
    
}


#endif	/* INVERSEWISHARTDISTRIBUTION_H */

