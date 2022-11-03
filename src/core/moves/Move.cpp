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

std::optional<double> Move::getMoveTuningParameter(void) const
{
    return {};
}

void Move::setMoveTuningParameter(double tp)
{
    throw RbException()<<"setProposalTuningParameter: tuning not implemented for "<<getMoveName();
}

void Move::autoTune(void)
{
    throw RbException()<<"autoTune: auto-tuning not implemented for "<<getMoveName();
}

bool Move::isTunable() const
{
    return getMoveTuningParameter().has_value();
}

std::ostream& RevBayesCore::operator<<(std::ostream& o, const Move& x)
{
    o << "Move";
    
    return o;
}
