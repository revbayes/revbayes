//
//  WangLandauMcmc.h
//  rb_mlandis
//
//  Created by Michael Landis on 8/12/13.
//  Copyright (c) 2013 Michael Landis. All rights reserved.
//

#ifndef __rb_mlandis__WangLandauMcmc__
#define __rb_mlandis__WangLandauMcmc__

#include "Mcmc.h"
#include <vector>
#include <string>

namespace RevBayesCore
{
    class WangLandauMcmc : public Mcmc
    {
        
    public:
        WangLandauMcmc(const Model& m, const RbVector<Move> &moves, const RbVector<Monitor> &mons, bool ca=true, double ch=1.0, size_t ci=0, size_t ns=1e6, double se=10e-5, double sc=0.25, double ub=0.0, double lb=1000.0, size_t np=20);
        
         virtual void               nextCycle(bool advanceCycle);
        
    protected:

        
    private:
        
        void initializeSteps(void);
        void initializePartitionLowerBound(void);
        void updateProportions(double x, int t);
        void resetProportions(void);
        void updateBias(void);

        std::vector<double> bias;                   // theta
        std::vector<double> proportions;            // nu
        std::vector<double> steps;                  // gamma
        std::vector<double> partitionLowerBound;    // xi
        
        double proportionConstant;

        size_t stepIndex;
        size_t numSteps;
        double stepEpsilon;
        double stepConstant;
        
        size_t partitionIndex;
        size_t numPartitions;
        double upperBoundEnergy;
        double lowerBoundEnergy;
    };
}

#endif /* defined(__rb_mlandis__WangLandauMcmc__) */
