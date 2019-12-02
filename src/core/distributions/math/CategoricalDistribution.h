#ifndef CategoricalDistribution_H
#define CategoricalDistribution_H

#include "Simplex.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore {
    
    /**
     * \brief Categorical distribution
     *
     * This class implements the categorical distribution. The number
     * of categories in the distribution is determined by the length,
     * k, of the simplex parameter probs.
     *
     * The sample space is considered to be { 1, 2, 3, ..., k }.
     * This is the way the distribution is set up in BUGS, and it
     * is also how it is described on Wikipedia. By using a
     * sample space that is 1-offset instead of 0-offset, we also
     * avoid the need to subtract one from a categorical value
     * when using it as an index, which should be a common
     * scenario.
     *
     * \author Fredrik Ronquist
     */

    class CategoricalDistribution : public TypedDistribution<long> {
        
    public:
        CategoricalDistribution(const TypedDagNode< Simplex >* p);
        virtual                                            ~CategoricalDistribution(void);                                              //!< Virtual destructor
        
        // public member functions
        CategoricalDistribution*                            clone(void) const;                                                          //!< Create an independent clone
        double                                              computeLnProbability(void);                                                 //!< Compute probability of current value
        void                                                redrawValue(void);                                                          //!< Draw new value from prior
        
    protected:

        // parameter management function
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Swap a parameter
        
    private:
        
        // member
        const TypedDagNode< Simplex >*                      probs;                                                                      //!< Simplex with category probabilities
    };
    
}

#endif
