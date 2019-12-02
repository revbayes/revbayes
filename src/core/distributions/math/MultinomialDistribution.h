#ifndef MultinomialDistribution_H
#define MultinomialDistribution_H

#include "Simplex.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"

#include <vector>

namespace RevBayesCore {
    
    /**
     * @brief Multinomial distribution class.
     *
     * The Multinomial distribution represents a family of distributions
     * defined on a vector of real number. The Multinomial distribution has 2 parameters:
     *   n .. the number of tries
     *   p .. the vector of category probabilities
     * Instances of this class can be associated to stochastic variables.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Michael Landis)
     * @since 2014-05-06, version 1.0
     *
     */
    class MultinomialDistribution : public TypedDistribution< RbVector<long> > {
        
    public:
        MultinomialDistribution(const TypedDagNode< Simplex > *p, const TypedDagNode<long> *n);
        virtual                                            ~MultinomialDistribution(void);                                                //!< Virtual destructor
        
        // public member functions
        MultinomialDistribution*                            clone(void) const;                                                          //!< Create an independent clone
        double                                              computeLnProbability(void);
        void                                                redrawValue(void);
        
    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Swap a parameter
        
    private:
        
        // members
        const TypedDagNode<Simplex>*                        p;
        const TypedDagNode<long>*                           n;
    };
    
}


#endif
