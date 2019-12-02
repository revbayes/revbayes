#include "MatrixReal.h"
#include "RateGeneratorSequence.h"
#include "RateGenerator.h"
#include "RateMatrix.h"
#include "RateMatrix_JC.h"
#include "RbVector.h"
#include "RbException.h"
#include "RbMathMatrix.h"
#include "TransitionProbabilityMatrix.h"

#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>

using namespace RevBayesCore;

/** Construct rate matrix with n states */
RateGeneratorSequence::RateGeneratorSequence(size_t ns, size_t nc) : RateGenerator(ns)
{
    
//    numStates            = ns;
    num_characters        = nc;
    needsUpdate          = true;
    
}


/** Destructor */
RateGeneratorSequence::~RateGeneratorSequence(void)
{
    
    
}

size_t RateGeneratorSequence::getNumberOfCharacters( void ) const
{
    return num_characters;
}

std::ostream& RevBayesCore::operator<<(std::ostream& o, const RateGeneratorSequence& x)
{
    std::streamsize previousPrecision = o.precision();
    std::ios_base::fmtflags previousFlags = o.flags();
    
    o.setf(previousFlags);
    o.precision(previousPrecision);
    
    return o;
}

void RateGeneratorSequence::updateMap(void)
{
    if (needsUpdate)
    {
        ; // do nothing ...
        
        needsUpdate = false;
    }
}



