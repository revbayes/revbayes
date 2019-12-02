#include "RateMatrix.h"
#include "RbException.h"
#include "RbMathMatrix.h"
#include "TypedDagNode.h"

#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>

using namespace RevBayesCore;

/** Construct rate matrix with n states */
RateMatrix::RateMatrix(size_t n) : RateGenerator(n)
{

}





/** Destructor */
RateMatrix::~RateMatrix(void)
{
    
}


RateMatrix& RateMatrix::assign(const Assignable &m)
{
    const RateMatrix *rm = dynamic_cast<const RateMatrix*>(&m);
    if ( rm != NULL )
    {
        return operator=(*rm);
    }
    else
    {
        throw RbException("Could not assign rate matrix.");
    }
}
