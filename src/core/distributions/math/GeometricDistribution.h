#ifndef GeometricDistribution_H
#define GeometricDistribution_H

#include "TypedDistribution.h"
#include "TypedDagNode.h"

namespace RevBayesCore {
    
    /**
     * @brief Geometric distribution class.
     *
     * The geometric distribution represents ...
     * Instances of this class can be associated to stochastic variables.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     *
     */
    class GeometricDistribution : public TypedDistribution<long> {
        
    public:
        GeometricDistribution(const TypedDagNode<double> *p);
        virtual                                            ~GeometricDistribution(void);                                                  //!< Virtual destructor
        
        // public member functions
        GeometricDistribution*                              clone(void) const;                                                          //!< Create an independent clone
        double                                              computeLnProbability(void);
        void                                                redrawValue(void);

    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Swap a parameter
        
    private:
        
        // members
        const TypedDagNode<double>*                          p;
        
    };
    
}

#endif
