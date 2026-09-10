#include "DispersalExtinctionRootStructureFunction.h"

#include <vector>
#include <map>

#include "RbVector.h"
#include "RbVectorImpl.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

DispersalExtinctionRootStructureFunction::DispersalExtinctionRootStructureFunction(TypedDagNode<RbVector<double> >* rf, TypedDagNode<Simplex>* rs) : TypedFunction<Simplex>( new Simplex() ),
    root_frequencies( rf ),
    rangeSize(rs)
{
    // add the lambda parameter as a parent
    addParameter( root_frequencies );
    addParameter( rangeSize );

    num_states = root_frequencies->getValue().size();
    numCharacters = 0;
    size_t n = num_states;
    
    do {
        n /= 2;
        numCharacters++;
    } while(n != 1);

    (*value).resize(num_states, 0.0);
    
    
    makeBits();
    makeIdxByRangeSize();
   
    update();
}


DispersalExtinctionRootStructureFunction::~DispersalExtinctionRootStructureFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



DispersalExtinctionRootStructureFunction* DispersalExtinctionRootStructureFunction::clone( void ) const
{
    return new DispersalExtinctionRootStructureFunction( *this );
}


void DispersalExtinctionRootStructureFunction::makeBits(void)
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
    for (size_t i = 0; i < (size_t)num_states; i++)
    {
        inverseBits[ bits[i] ] = (unsigned)i;
    }
}

void DispersalExtinctionRootStructureFunction::makeIdxByRangeSize(void)
{
    idxByRangeSize.resize(numCharacters + 1);
    for (size_t i = 0; i < num_states; i++)
    {
        int n = 0;
        for (size_t j = 0; j < bits[i].size(); j++)
            n += bits[i][j];
        idxByRangeSize[n].push_back((unsigned)i);
    }
}

void DispersalExtinctionRootStructureFunction::reInitialized( void )
{
    ;//    *value = tau->getValue();
}


void DispersalExtinctionRootStructureFunction::update( void )
{
    const RbVector<double>& rf = root_frequencies->getValue();
    size_t minRangeSize = 1;
    
    // get normalized range sizes, s.t. largest range size value = 1
    const std::vector<double>& rs = rangeSize->getValue();
    double largestRangeSize = rs[0];
    for (size_t i = minRangeSize; i <= numCharacters; i++)
    {
        if (rs[i-1] > largestRangeSize)
        {
            largestRangeSize = rs[i-1];
        }
    }
    std::vector<double> normalizedRangeSize = rs;
    for (size_t i = 0; i < rs.size(); i++)
    {
        normalizedRangeSize[i] = rs[i] / largestRangeSize;
    }
    
    
    std::vector<std::vector<unsigned> >::iterator it;
    std::vector<unsigned>::iterator jt;
    
    *value = RbVector<double>(num_states, 0.0);
    
    double sum = 0.0;
    for (size_t i = minRangeSize; i <= numCharacters; i++)
    {
        for (std::vector<unsigned>::iterator it = idxByRangeSize[i].begin(); it != idxByRangeSize[i].end(); it++)
        {
            (*value)[ *it ] = rf[ *it ] * normalizedRangeSize[i-1];
            sum += rf[ *it ];
        }
    }
    for (size_t i = minRangeSize; i <= numCharacters; i++)
    {
        for (std::vector<unsigned>::iterator it = idxByRangeSize[i].begin(); it != idxByRangeSize[i].end(); it++)
        {
            (*value)[ *it ] = (*value)[ *it ] / sum;
        }
    }
}



void DispersalExtinctionRootStructureFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    if (oldP == root_frequencies)
    {
        root_frequencies = static_cast<const TypedDagNode<RbVector<double> >* >( newP );
    }
    else if (oldP == rangeSize)
    {
        rangeSize = static_cast<const TypedDagNode<Simplex>* >( newP );
    }
}
