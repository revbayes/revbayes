/**
 * @file
 * This file contains the declaration of DPPTableValueUpdate, which performs the DPP move based on Neal (2000) Algorithm 8
 *
 * @brief Declaration of DPPTableValueUpdate
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date: 2012-05-11 14:54:35 +0200 (Fri, 11 May 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2012-07-20, version 1.0
 *
 * $Id: DPPTableValueUpdate.h $
 */

#ifndef DPPTableValueUpdate_H
#define DPPTableValueUpdate_H

#include <ostream>

#include "DirichletProcessPriorDistribution.h"
#include "AbstractGibbsMove.h"
#include "RbVector.h"
#include "SimpleProposal.h"

namespace RevBayesCore {
    
    template <class valueType>
    class DPPTableValueUpdate : public AbstractGibbsMove {
        
    public:
        DPPTableValueUpdate(StochasticNode< RbVector<valueType> >* v, SimpleProposal<valueType> *p, double w );                                                                      //!< Internal constructor
        
        // Basic utility functions
        DPPTableValueUpdate<valueType>*                         clone(void) const;                                                                  //!< Clone object
        const std::string&                                      getMoveName(void) const;                                                            //!< Get the name of the move for summary printing
        
    protected:
        void													performGibbsMove(void);                                                            //!< Perform move
        void                                                    swapNodeInternal(DagNode *oldN, DagNode *newN);
        
    private:
        double													getLnProbabilityForMove(void);
        int														findTableIDForVal(std::vector<valueType> tvs, valueType val);
        void													dppNormalizeVector(std::vector<double> &v);
        int														findElementNewTable(double u, std::vector<double> lnProb);
        double                                                  safeExponentiation(double x);
        
        StochasticNode< RbVector<valueType> >*                  variable;
        SimpleProposal<valueType>*                              proposal;

    };
}

#include "NormalizeVectorFunction.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbMathVector.h"

#include <cmath>

template <class valueType>
RevBayesCore::DPPTableValueUpdate<valueType>::DPPTableValueUpdate(StochasticNode< RbVector<valueType> > *v, SimpleProposal<valueType>* p, double w) : AbstractGibbsMove( w ),
    variable( v ),
    proposal( p )
{
    addNode( variable );
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
template <class valueType>
RevBayesCore::DPPTableValueUpdate<valueType>* RevBayesCore::DPPTableValueUpdate<valueType>::clone( void ) const
{
    
    return new DPPTableValueUpdate<valueType>( *this );
}


template <class valueType>
const std::string& RevBayesCore::DPPTableValueUpdate<valueType>::getMoveName( void ) const
{
    static std::string name = "DPPTableValueUpdate";
    
    return name;
}


/** Perform the move */
template <class valueType>
void RevBayesCore::DPPTableValueUpdate<valueType>::performGibbsMove( void )
{
    
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    
    DirichletProcessPriorDistribution<valueType>& dist = static_cast<DirichletProcessPriorDistribution<valueType> &>( variable->getDistribution() );
    
    dist.getTableParameters();
    
    dist.createRestaurantVectors();
    int numTables = dist.getNumberOfCategories();
    int numElements = dist.getNumberOfElements();
    std::vector<valueType> tableVals = dist.getTableParameters();
    std::vector<int> allocVec = dist.getElementAllocationVec();
    std::vector<int> numPerTab = dist.getNumElemPerTable();
    std::vector<valueType>& elementVals = variable->getValue();
    
    TypedDistribution<valueType>* g0 = dist.getBaseDistribution();
    
    for (int i=0; i<numTables; i++)
    {
        
        // get old lnL
        double oldLnl = getLnProbabilityForMove();
        
        valueType storedValue = tableVals[i];
        double ln_hastings_ratio = proposal->propose( tableVals[i] );
//        tableVals[i] = newValue;
        
        // Assign new value to elements
        for (int j=0; j<numElements; j++)
        {
            if ( allocVec[j] == i )
            {
                elementVals[j] = tableVals[i];
            }
        }
        
        g0->getValue() = tableVals[i]; // new
        double priorRatio = g0->computeLnProbability();
        g0->getValue() = storedValue; // old
        priorRatio -= g0->computeLnProbability();
        
        variable->touch();
        double newLnl = getLnProbabilityForMove();
        double lnProbRatio = newLnl - oldLnl;
        
        double r = safeExponentiation(priorRatio + lnProbRatio + ln_hastings_ratio);
        double u = rng->uniform01();
        if ( u < r ) //accept
        {
            variable->keep();
        }
        else // reject
        {
            for (int j=0; j<numElements; j++)
            {
                if (allocVec[j] == i)
                {
                    elementVals[j] = storedValue;
                }
            }
            
            tableVals[i] = storedValue;
            variable->touch();
            variable->keep();
        }
    }
    
    dist.createRestaurantVectors();
    
}



template <class valueType>
double RevBayesCore::DPPTableValueUpdate<valueType>::safeExponentiation(double x)
{
    
    if (x < -300.0)
        return 0.0;
    else if (x > 0.0)
        return 1.0;
    else
        return exp(x);
}


template <class valueType>
void RevBayesCore::DPPTableValueUpdate<valueType>::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    variable = static_cast<StochasticNode< RbVector<valueType> > * >( newN );
    
}


template <class valueType>
double RevBayesCore::DPPTableValueUpdate<valueType>::getLnProbabilityForMove(void)
{
    
    RbOrderedSet<DagNode*> affected;
    variable->initiateGetAffectedNodes( affected );
    double lnProb = 0.0;
    for (RbOrderedSet<DagNode*>::iterator it = affected.begin(); it != affected.end(); ++it)
    {
        double lp = (*it)->getLnProbability();
        lnProb += lp;
    }
    
    return lnProb;
}

template <class valueType>
int RevBayesCore::DPPTableValueUpdate<valueType>::findTableIDForVal(std::vector<valueType> tvs, valueType val) {
    
    for (size_t j=0; j<tvs.size(); j++)
    {
        if (tvs[j] == val)
        {
            // this is an existing table
            return int(j);
        }
    }
    
    return -1; // new table
}


template <class valueType>
void RevBayesCore::DPPTableValueUpdate<valueType>::dppNormalizeVector(std::vector<double> &v)
{
    
    size_t n = v.size();
    double lnC = v[0];
    for (size_t i=1; i<n; i++)
    {
        if (v[i] > lnC)
        {
            lnC = v[i];
        }
    }
    
    for (size_t i=0; i<n; i++)
    {
        v[i] -= lnC;
    }
    
    double sum = 0.0;
    for (size_t i=0; i<n; i++)
    {
        if ( v[i] < -300.0 )
        {
            v[i] = 0.0;
        }
        else
        {
            v[i] = exp( v[i] );
        }
        sum += v[i];
        
    }
    
    for (size_t i=0; i<n; i++)
    {
        v[i] /= sum;
    }
    
}

template <class valueType>
int RevBayesCore::DPPTableValueUpdate<valueType>::findElementNewTable(double u, std::vector<double> lnProb)
{
    
    for (size_t j = 0; j < lnProb.size(); j++)
    {
        u -= lnProb[j];
        if (u < 0.0){
            return (int)j;
        }
    }
    
    return -1;
}



#endif

