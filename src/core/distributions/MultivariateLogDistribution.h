#ifndef MultivariateLogDistribution_H
#define MultivariateLogDistribution_H

#include "TypedDistribution.h"

namespace RevBayesCore {
    
    /**
     * This class takes a multivariate distribution and transforms it by exponentiating each element.
     * The resulting distribution has distribution `dist` in log-space.
     *
     * @copyright Copyright 2024-
     * @author The RevBayes Development Core Team (Ben Redelings)
     * @since 2024-05-20, version 1.0
     */

    class MultivariateLogDistribution : public TypedDistribution< RbVector<double> > {
        
    public:
        // constructor(s)
        MultivariateLogDistribution(const TypedDistribution< RbVector<double> >& vp);
        MultivariateLogDistribution(const MultivariateLogDistribution &d);

        // public member functions
        MultivariateLogDistribution*                        clone(void) const;                                                                      //!< Create an independent clone
        double                                              computeLnProbability(void);
        void                                                redrawValue(void);
        
    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Swap a parameter
        
    private:
        
        // helper methods
        void                                                simulate();
        
        // private members
	std::unique_ptr<TypedDistribution<RbVector<double>>>  dist;
    };
    
}

#endif
