#ifndef HalfNormalDistribution_H
#define HalfNormalDistribution_H

#include "ContinuousDistribution.h"

 /**
     * @brief Half-Normal Distribution class.
     *
     * The Half-Normal distribution represents a family of distributions
     * defined on the postive numbers. The Half-Normal Distribution is a special case of the Folded Normal Distribution with a location parameter of 0.
     *  Instances of this class can be associated to stochastic variables.
     * The Half-Normal distribution has 2 parameters:
     * @param Offest This parameter represents the location parameter of the Folded Normal Distribution.
     * @param Sigma This paramer repesents the variance. In the distribution constructor, the parameter is the standard deviation, the square root of variance.
     *
     */



namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;
    
    class HalfNormalDistribution : public ContinuousDistribution {
        
    public:
        HalfNormalDistribution(const TypedDagNode<double> *o, const TypedDagNode<double> *s);
        virtual                                            ~HalfNormalDistribution(void);                                                   //!< Virtual destructor
        
        // public member functions
        double                                              cdf(void) const;                                                                  //!< Cummulative density function
        HalfNormalDistribution*                        clone(void) const;                                                          //!< Create an independent clone
        double                                              computeLnProbability(void);
        double                                              getMax(void) const;
        double                                              getMin(void) const;
        double                                              quantile(double p) const;                                                  //!< Qu
        void                                                redrawValue(void);
        
    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Swap a parameter
        
    private:
        
        // members
        const TypedDagNode<double>*                         offset;
        const TypedDagNode<double>*                         stDev;

    };
    
}

#endif
