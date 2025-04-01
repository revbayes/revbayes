#ifndef DecomposedInverseWishartDistribution_H
#define	DecomposedInverseWishartDistribution_H

#include "TypedDistribution.h"
#include "MatrixReal.h"
#include "TypedDagNode.h"



namespace RevBayesCore {
class DagNode;
template <class valueType> class RbVector;
    
    /**
     * @brief Decomposed inverse Wishart distribution class.
     *
     * The Decomposed inverse Wishart distribution represents a family of distributions
     * on the ... The decomposed inverse Wishart distribution has ?? parameters:
     *   ??
     *   ??
     * Instances of this class can be associated to stochastic variables.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Nicolas Lartillot)
     * @since 2014-11-21, version 1.0
     *
     */
    class DecomposedInverseWishartDistribution : public TypedDistribution<MatrixReal>   {
        
    public:
                                                    DecomposedInverseWishartDistribution(const TypedDagNode<MatrixReal> *insigma0, const TypedDagNode<std::int64_t>* indf);
                                                    DecomposedInverseWishartDistribution(const TypedDagNode<RbVector<double> > *inkappaVector, const TypedDagNode<std::int64_t>* indf);
                                                    DecomposedInverseWishartDistribution(const TypedDagNode<std::int64_t>* indim, const TypedDagNode<double> *inkappa, const TypedDagNode<std::int64_t>* indf);
        virtual                                    ~DecomposedInverseWishartDistribution(void) {}
        DecomposedInverseWishartDistribution*       clone(void) const;                                                          //!< Create an independent clone
        double                                      computeLnProbability(void);
        void                                        redrawValue(void);
        std::int64_t                                        getDF(void) const {return df->getValue();}
        
    protected:
        void                                        swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Swap a parameter

    private:
        const TypedDagNode<MatrixReal>*             sigma0;
        const TypedDagNode<RbVector<double> >*      kappaVector;
        const TypedDagNode<double>*                 kappa;
        const TypedDagNode<std::int64_t>*                    df;
        const TypedDagNode<std::int64_t>*                    dim;
    };
    
}

#endif

