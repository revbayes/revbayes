/*
 * File:   ExponentialErrorDistribution.h
 * Author: DavidCerny
 *
 * Created on October 8, 2019
 */

#ifndef EXPONENTIALERRORDISTRIBUTION_H
#define EXPONENTIALERRORDISTRIBUTION_H

#include <cstddef>

#include "TypedDistribution.h"
#include "AverageDistanceMatrix.h"
#include "TypedFunction.h"
#include "TypedDagNode.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class RbVector;
    
    class ExponentialErrorDistribution : public TypedDistribution<AverageDistanceMatrix>   {
        
    public:
        
        // Exponential error distribution on an average distance matrix with parameter lambda
        ExponentialErrorDistribution(const TypedDagNode<AverageDistanceMatrix> *adm, const TypedDagNode<double> *l);
        
        // Exponential error distribution on a distance matrix dm = adm.getDistanceMatrix() with parameter lambda
        ExponentialErrorDistribution(const TypedDagNode<DistanceMatrix> *adm, const TypedDagNode<double> *l);

        virtual                                     ~ExponentialErrorDistribution(void) {}
        
        // Public member functions

        ExponentialErrorDistribution*               clone(void) const;                                               //!< creates an independent clone
        double                                      computeLnProbability(void);
        void                                        redrawValue(void);
        
    protected:
        // Parameter management functions
        void                                        swapParameterInternal(const DagNode *oldP, const DagNode *newP); //!< swaps a parameter

    private:
        // Members
        const TypedDagNode<AverageDistanceMatrix>*  avgDistMat;
        const TypedDagNode<DistanceMatrix>*         distMat;
        const TypedDagNode<double>*                 lambda;
                
    };
    
}


#endif    /* EXPONENTIALERRORDISTRIBUTION_H */
