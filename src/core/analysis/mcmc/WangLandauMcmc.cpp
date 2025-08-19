//
//  WangLandauMcmc.cpp
//  rb_mlandis
//
//  Created by Michael Landis on 8/12/13.
//  Copyright (c) 2013 Michael Landis. All rights reserved.
//

#include <cstddef>
#include <cmath>
#include <vector>

#include "RbConstants.h"
#include "RandomNumberFactory.h"
#include "WangLandauMcmc.h"
#include "Mcmc.h"
#include "Model.h"
#include "Move.h"
#include "MoveSchedule.h"
#include "RbVector.h"

namespace RevBayesCore { class Monitor; }

using namespace RevBayesCore;

WangLandauMcmc::WangLandauMcmc(const Model& m, const RbVector<Move> &moves, const RbVector<Monitor> &mons, bool ca, double ch, size_t ci, size_t ns, double se, double sc, double ub, double lb, size_t np) : Mcmc(m,moves,mons), numSteps(ns), stepEpsilon(se), stepConstant(sc), numPartitions(np), upperBoundEnergy(ub), lowerBoundEnergy(lb)
{
    
    // initialize vectors
    initializeSteps();
    initializePartitionLowerBound();
    proportions.resize(numPartitions, 0.0);
    bias.resize(numPartitions, 1.0);
    
}

void WangLandauMcmc::initializeSteps(void)
{
    bool isGeqEpsilon = true;

    // piecewise function for step sizes
    for (size_t n = 0; n < numSteps; n++)
    {
        double v = ( isGeqEpsilon ? pow(stepConstant,-n) : 1.0 / n );
        if (isGeqEpsilon && v < stepEpsilon)
            isGeqEpsilon = false;
        steps.push_back(v);
    }
}

void WangLandauMcmc::initializePartitionLowerBound(void)
{
    
    double db = (upperBoundEnergy - lowerBoundEnergy) / (numPartitions - 1);
    
    double v = upperBoundEnergy;
    for (size_t i = 0; i < (numPartitions - 1); i++)
    {
        v -= db;
        partitionLowerBound.push_back(v);
    }
    partitionLowerBound.push_back(RbConstants::Double::neginf);
}

void WangLandauMcmc::updateProportions(double x, int t)
{
    // which partition is x in?
    partitionIndex = 0;
    for (size_t i = 0; i < numPartitions; i++)
    {
        if (x > partitionLowerBound[i])
        {
            partitionIndex = i;
            break;
        }
    }
    
    // update partitions
    for (size_t i = 0; i < numPartitions; i++)
        proportions[i] = ((t - 1) * proportions[i]) / t;
    
    // further update found partition
    proportions[partitionIndex] += 1.0 / t;
    
    // "flat histogram" resets proportions as needed
    resetProportions();
}

void WangLandauMcmc::resetProportions(void)
{
    double maxValue = 0.0;
    
    for (size_t i = 0; i < numPartitions; i++)
    {
        if (proportions[i] > maxValue)
        {
            maxValue = proportions[i];
        }
    }
    
    double lhs = fabs(maxValue - (1.0/numPartitions));
    double rhs = proportionConstant/numPartitions;
    if (lhs < rhs)
    {
        stepIndex++;
        proportions.resize(numPartitions,0.0);
    }
}

void WangLandauMcmc::updateBias(void)
{
    
    // update bias given current state's partition
    bias[partitionIndex] += steps[stepIndex] * (1.0 - (1.0/numPartitions) );
    
    // renormalize vector
    double sum = 0.0;
    for (size_t i = 0; i < numPartitions; i++)
        sum += bias[i];
    for (size_t i = 0; i < numPartitions; i++)
        bias[i] /= sum;
}


void WangLandauMcmc::nextCycle(bool advanceCycle) {
    
    size_t proposals = size_t( round( schedule->getNumberMovesPerIteration() ) );
    for (size_t i=0; i<proposals; i++)
    {
        // Get the move
        Move& the_move = schedule->nextMove( generation );
        
        the_move.performMcmcStep( 1.0, chain_likelihood_heat, 1.0 );
    }
    
    
    // advance gen cycle if needed (i.e. run()==true, burnin()==false)
    if (advanceCycle)
        generation++;
    
}
