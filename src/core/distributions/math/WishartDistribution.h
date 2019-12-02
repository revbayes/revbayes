//
//  WishartDistribution.h
//  revbayes
//
//  Created by Nicolas Lartillot on 2014-03-26.
//  Copyright (c) 2014 revbayes team. All rights reserved.
//

#ifndef __revbayes__WishartDistribution__
#define __revbayes__WishartDistribution__

#include <iostream>



#include "TypedDistribution.h"
#include "MatrixReal.h"
#include "TypedDagNode.h"

namespace RevBayesCore {
    
    class WishartDistribution : public TypedDistribution<MatrixReal>   {
        
    public:
        
        // inverse Wishart distribution of parameter sigma0 et df degrees of freedom
        WishartDistribution(const TypedDagNode<MatrixReal> *insigma0, const TypedDagNode<long>* indf);
        // specialized version: inverse Wishart of parameter sigma0=kappa*I and df degrees of freedom
        WishartDistribution(const TypedDagNode<long>* indim, const TypedDagNode<double> *inkappa, const TypedDagNode<long>* indf);
        
        virtual                                            ~WishartDistribution(void) {}
        
        // public member functions

        WishartDistribution*                                clone(void) const;                                                          //!< Create an independent clone
        double                                              computeLnProbability(void);
        void                                                redrawValue(void);
        
        long                                                getDF() const {return df->getValue();}
        
        
    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Swap a parameter
        
    private:

        // members
        
        const TypedDagNode<MatrixReal>*                     omega0;
        const TypedDagNode<double>*                         kappa;
        const TypedDagNode<long>*                           df;
        const TypedDagNode<long>*                           dim;
        
    };
    
}


#endif /* defined(__revbayes__WishartDistribution__) */
