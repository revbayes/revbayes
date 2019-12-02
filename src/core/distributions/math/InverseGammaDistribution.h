//
//  InverseGammaDistribution.h
//  rb_mlandis
//
//  Created by Michael Landis on 3/19/13.
//  Copyright (c) 2013 Michael Landis. All rights reserved.
//

#ifndef __rb_mlandis__InverseGammaDistribution__
#define __rb_mlandis__InverseGammaDistribution__


#include "ContinuousDistribution.h"
#include "TypedDagNode.h"

namespace RevBayesCore {
    
    class InverseGammaDistribution : public ContinuousDistribution {
        
    public:
        InverseGammaDistribution(const TypedDagNode<double> *sd, const TypedDagNode<double> *r);
        virtual                                            ~InverseGammaDistribution(void);                                                  //!< Virtual destructor
        
        // public member functions
        double                                              cdf(void) const;                                                                  //!< Cummulative density function
        InverseGammaDistribution*                           clone(void) const;                                                          //!< Create an independent clone
        double                                              computeLnProbability(void);
        double                                              getMax(void) const;
        double                                              getMin(void) const;
        double                                              quantile(double p) const;                                                       //!< Qu
        void                                                redrawValue(void);
        
        
    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Swap a parameter
        
    private:
        
        // members
        const TypedDagNode<double>*                          shape;
        const TypedDagNode<double>*                          rate;
        
    };
    
}


#endif /* defined(__rb_mlandis__InverseGammaDistribution__) */
