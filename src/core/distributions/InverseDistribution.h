#ifndef InverseDistribution_h
#define InverseDistribution_h

#include "TypedDistribution.h"

namespace RevBayesCore {
    
    /**
     * Provides the inverse of the probability of a distribution.
     *
     * @copyright Copyright 2024-
     * @author Martin R. Smith
     * @since 2024-07-11, version 1.2.5
     *
     */
    template<typename valType>
    class InverseDistribution :  public TypedDistribution< valType > {
     
    public:
        // constructor(s)
        InverseDistribution(const TypedDistribution<valType>& vp);
        InverseDistribution(const InverseDistribution &d);

        // public member functions
        InverseDistribution*                                clone(void) const;                                                       //!< Create an independent clone
        double                                              computeLnProbability(void);
        void                                                redrawValue(void);
        
    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);        //!< Swap a parameter
        
        
    private:
        
        // helper methods
        void                                                simulate();
        
        // private members
        std::unique_ptr<TypedDistribution<valType>>         dist;
    };
}


#endif /* Func_inverse_h */
