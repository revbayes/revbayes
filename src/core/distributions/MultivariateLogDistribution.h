#ifndef MultivariateLogDistribution_H
#define MultivariateLogDistribution_H

#include "TypedDistribution.h"

namespace RevBayesCore {
    
    
    /**
     * This class implements a generic mixture distribution between several possible values.
     *
     * This mixture can be considered as a multinomial distribution. We specify a vector of probabilities
     * and a vector of values. Then, a value drawn from this distribution takes each value corresponding to
     * its probability.
     * The values are already of the correct mixture type. You may want to apply a mixture allocation move
     * to change between the current value. The values themselves change automatically when the input parameters change.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-11-18, version 1.0
     */
    class MultivariateLogDistribution : public TypedDistribution< RbVector<double> > {
        
    public:
        // constructor(s)
        MultivariateLogDistribution(const TypedDistribution< RbVector<double> >& vp);
        MultivariateLogDistribution(const MultivariateLogDistribution &d);

        // public member functions
        MultivariateLogDistribution*                        clone(void) const;                                                                                  //!< Create an independent clone
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
