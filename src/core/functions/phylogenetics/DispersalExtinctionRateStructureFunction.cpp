//
//  DispersalExtinctionRateStructureFunction.cpp
//  revbayes-proj
//
//  Created by Michael Landis on 3/3/15.
//  Copyright (c) 2015 Michael Landis. All rights reserved.
//

#include "DispersalExtinctionRateStructureFunction.h"

#include <cmath>
#include <vector>
#include <map>

#include "Cloner.h"
#include "RbVectorImpl.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

DispersalExtinctionRateStructureFunction::DispersalExtinctionRateStructureFunction(const TypedDagNode<RbVector<RbVector<double> > >* dr, TypedDagNode<RbVector<double> >* er, TypedDagNode<std::int64_t>* rs) : TypedFunction<RbVector<double> >( new RbVector<double>() ),
    dispersalRates( dr ),
    extinctionRates( er ),
    maxRangeSize(rs)
{
    // add the lambda parameter as a parent
    addParameter( dispersalRates );
    addParameter( extinctionRates );
    addParameter( maxRangeSize );
    
    numCharacters = extinctionRates->getValue().size();
    num_states = (int)(pow(2,numCharacters));
    
    (*value).resize(num_states*num_states-num_states, 0.0);
    
    
    makeBits();
    makeTransitions();
    
    update();
}


DispersalExtinctionRateStructureFunction::~DispersalExtinctionRateStructureFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



DispersalExtinctionRateStructureFunction* DispersalExtinctionRateStructureFunction::clone( void ) const
{
    return new DispersalExtinctionRateStructureFunction( *this );
}


void DispersalExtinctionRateStructureFunction::keep( const DagNode *affecter)
{
    //delegate to base class
    TypedFunction< RbVector<double> >::keep( affecter );
    
}

void DispersalExtinctionRateStructureFunction::makeBits(void)
{
    bits = std::vector<std::vector<unsigned> >(num_states, std::vector<unsigned>(numCharacters, 0));
    for (size_t i = 1; i < num_states; i++)
    {
        size_t n = i;
        for (size_t j = 0; j < numCharacters; j++)
        {
            bits[i][j] = n % 2;
            n /= 2;
            if (n == 0)
                break;
        }
    }
    for (size_t i = 0; i < num_states; i++)
    {
        inverseBits[ bits[i] ] = (unsigned)i;
    }
}

void DispersalExtinctionRateStructureFunction::makeTransitions(void)
{
    transitions.resize(num_states);
    lossOrGain.resize(num_states);
    transitionAreas.resize(num_states);
    
    // populate integer-valued transitions between states
    for (size_t i = 1; i < num_states; i++)
    {
        std::vector<unsigned> b = bits[i];
        
        // each row has b.size() events (excluding i==0)
        for (size_t j = 0; j < b.size(); j++)
        {
            std::vector<unsigned> tmp = b;
            
            // change the range cfg at area j
            tmp[j] = (b[j] == 0 ? 1 : 0);
            
            // store integer-valued event
            transitions[i].push_back(inverseBits[tmp]);
            
            // is event a gain or a loss?
            lossOrGain[i].push_back(tmp[j]);
            
            std::vector<unsigned> a;
            // store dispersal event source areas
            if (tmp[j]==1)
            {
                for (size_t k = 0; k < b.size(); k++)
                {
                    if (b[k]==1)
                    {
                        a.push_back((unsigned)k);
                    }
                }
            }
            // extinction events pushes only the lost area
            else
            {
                a.push_back((unsigned)j);
            }
            transitionAreas[i].push_back(a);
        }
    }
    
    
}

void DispersalExtinctionRateStructureFunction::reInitialized( void )
{
    ;//    *value = tau->getValue();
}


void DispersalExtinctionRateStructureFunction::restore( const DagNode *restorer)
{
    //delegate to base class
    TypedFunction< RbVector<double> >::restore( restorer );
}


void DispersalExtinctionRateStructureFunction::touch(const DagNode *toucher)
{
    
    //delegate to base class
    TypedFunction< RbVector<double> >::touch( toucher );
    
}


void DispersalExtinctionRateStructureFunction::update( void )
{
    
    const RbVector<RbVector<double> >& dr = dispersalRates->getValue();
    const RbVector<double>& er            = extinctionRates->getValue();
    
    std::vector<std::vector<unsigned> >::iterator it;
    std::vector<unsigned>::iterator jt;
    
//    std::cout << "\n";
    int offset = 0;
    for (size_t i = 1; i < transitions.size(); i++)
    {
        int n = 0;
        for (size_t j = 0; j < bits[i].size(); j++)
            n += bits[i][j];
        bool maxSize = n >= maxRangeSize->getValue();
        
//        for (size_t j = 0; j < bits[i].size(); j++)
//            std::cout << bits[i][j];
//        std::cout << " : ";
//        bool lower = true;
        for (size_t j = 0; j < transitions[i].size(); j++)
        {
            double v = 0.0;
            // extinction
            if (lossOrGain[i][j] == 0)
            {
                v = er[ transitionAreas[i][j][0] ];
            }
            // dispersal
            else if (lossOrGain[i][j] == 1 && !maxSize)
            {
                for (size_t k = 0; k < transitionAreas[i][j].size(); k++)
                {
                    v += dr[ transitionAreas[i][j][k] ][ j ];
                }
            }

            // store value
            offset = (int)(-i-1);
            if (i > transitions[i][j]) // && lower)
            {
                offset++;
//                offset--;
//                lower = false;
            }            
//            std::cout << "\n";
//            std::cout << num_states*i+transitions[i][j]+offset << " = " << num_states*i << " + " << transitions[i][j] << " + " << offset << " " << v << "\n";
            (*value)[num_states*i + transitions[i][j] + offset] = v;
           
//            std::cout << "  ";
//            for (size_t k = 0; k < bits[i].size(); k++)
//                std::cout << bits[transitions[i][j]][k];
//            std::cout << ":" << v << " ";
        }
//        std::cout << "\n";
    }
    
}



void DispersalExtinctionRateStructureFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == dispersalRates)
    {
        dispersalRates = static_cast<const TypedDagNode<RbVector<RbVector<double> > >* >( newP );
    }
    else if (oldP == extinctionRates)
    {
        extinctionRates = static_cast<const TypedDagNode<RbVector<double> >* >( newP );
    }
    else if (oldP == maxRangeSize)
    {
        maxRangeSize = static_cast<const TypedDagNode<std::int64_t>* >( newP );
    }
}
