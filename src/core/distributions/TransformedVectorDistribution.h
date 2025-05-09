#ifndef TransformedVectorDistribution_H
#define TransformedVectorDistribution_H

#include "TypedDistribution.h"

namespace RevBayesCore {


    /**
     * This class implements the distribution for a transformed random variable.
     * If x ~ base_dist, then this class implements the distribution for f(x).
     *
     * We need to know f(x), f_inverse(x), and f'(x) -- the derivative of f.
     *
     * By specifying these three functions in the constructor, we can use a single class
     * for multiple different transformations.
     *
     * @copyright Copyright 2024-
     * @author The RevBayes Development Core Team (Ben Redelings)
     * @since 2024-07-23, version 1.0
     */
    class TransformedVectorDistribution : public TypedDistribution< RbVector<double> >
    {

    public:
        // Allow returning nothing in case the input value is invalid.
        typedef std::function<std::optional<std::vector<double>>(const std::vector<const DagNode*>&, const std::vector<double>&)> func_t;

        typedef std::function<std::optional<double>(const std::vector<const DagNode*>&, const std::vector<double>&)> jacobian_t;

        typedef std::function<std::optional<double>(const std::vector<const DagNode*>&, double)> scalar_func_t;

        // constructor(s)
        TransformedVectorDistribution(std::unique_ptr<TypedDistribution<RbVector<double>>>& vp, func_t F, func_t FI, jacobian_t LFP, const std::vector<const DagNode*>& p = {});
        TransformedVectorDistribution(std::unique_ptr<TypedDistribution<RbVector<double>>>& vp, scalar_func_t F, scalar_func_t FI, scalar_func_t LFP, const std::vector<const DagNode*>& p = {});
        TransformedVectorDistribution(const TransformedVectorDistribution &d);

        // public member functions
        TransformedVectorDistribution*                      clone(void) const override;                                                             //!< Create an independent clone
        double                                              computeLnProbability(void) override;
        void                                                redrawValue(void) override;
        void                                                getAffected(RbOrderedSet<DagNode *> &affected, const DagNode* affecter) override;
        void                                                setValue(RbVector<double> *v, bool f=false) override;                                   //!< Set the current value, e.g. attach an observation (clamp)

    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP) override;               //!< Swap a parameter

        // keep specialization for derived classes
        void                                                keepSpecialization(const DagNode* affecter) override;
        void                                                restoreSpecialization(const DagNode *restorer) override;
        void                                                touchSpecialization(const DagNode *toucher, bool touchAll) override;

    private:

        // helper methods
        void                                                simulate();

        // the transformation
        func_t                                              f = nullptr;
        func_t                                              f_inverse = nullptr;
        jacobian_t                                          log_f_prime = nullptr;
        std::vector<const DagNode*>                         transform_params;                                                                       //!< Parameters used by the transformation.

        // the base distribution
	std::unique_ptr<TypedDistribution<RbVector<double>>>  base_dist;
    };

}

#endif
