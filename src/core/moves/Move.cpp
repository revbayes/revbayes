#include "Move.h"

#include <cmath>

#include "RbConstants.h"
#include "RbException.h"

using namespace RevBayesCore;

Move::Move()
{
    
}


Move::~Move()
{
    
}

double Move::getMoveTuningParameter(void) const
{
    return RbConstants::Double::nan;
}

void Move::setMoveTuningParameter(double tp)
{
    throw RbException()<<"setProposalTuningParameter: tuning not implemented for "<<getMoveName();
}

void Move::autoTune(void)
{
    throw RbException()<<"autoTune: auto-tuning not implemented for "<<getMoveName();
}

std::ostream& RevBayesCore::operator<<(std::ostream& o, const Move& x)
{
    o << "Move";
    
    return o;
}
