#ifndef TransformedDistribution_H
#define TransformedDistribution_H

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
    class TransformedDistribution : public TypedDistribution< double >
    {

	// Allow returning nothing in case the input value is invalid.
	typedef std::function<std::optional<double>(double)> func_t;

    public:
        // constructor(s)
        TransformedDistribution(const TypedDistribution<double>& vp, func_t F, func_t FI, func_t LFP, const std::vector<const DagNode*>& p = {});
        TransformedDistribution(const TransformedDistribution &d);

        // public member functions
        TransformedDistribution*                            clone(void) const;                                                                                  //!< Create an independent clone
        double                                              computeLnProbability(void);
        void                                                redrawValue(void);

    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Swap a parameter
        
        
    private:
        
        // helper methods
        void                                                simulate();

	// the transformation
	func_t                                              f = nullptr;
	func_t                                              f_inverse = nullptr;
	func_t                                              log_f_prime = nullptr;

        // the base distribution
	std::unique_ptr<TypedDistribution<double>>          base_dist;
    };
    
}

#endif
