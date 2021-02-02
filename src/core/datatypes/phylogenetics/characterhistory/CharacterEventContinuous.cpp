#include "CharacterEventContinuous.h"

#include <sstream> // IWYU pragma: keep

#include "RbException.h"
#include "Cloneable.h"

using namespace RevBayesCore;


CharacterEventContinuous::CharacterEventContinuous(void) : CharacterEvent(),
    value()
{
    
}


CharacterEventContinuous::CharacterEventContinuous(size_t i, double v, double a, size_t t) : CharacterEvent(i,a,t),
    value()
{
    value.push_back( v );
}


CharacterEventContinuous::CharacterEventContinuous(size_t i, const std::vector<double> &v, double a, size_t t) : CharacterEvent(i,a,t),
    value( v )
{

}


CharacterEventContinuous::CharacterEventContinuous(const CharacterEventContinuous& c) : CharacterEvent(c),
    value(c.value)
{
    
}

CharacterEventContinuous::~CharacterEventContinuous(void)
{
    
}

CharacterEventContinuous* CharacterEventContinuous::clone( void ) const
{
    return new CharacterEventContinuous( *this );
}


double CharacterEventContinuous::getState( size_t index ) const
{
    if ( index >= value.size() )
    {
        std::stringstream ss;
        ss << "Cannot access the " << (index+1) << "-th value of this continuous character event because it has only " << value.size() << " events.";

        throw RbException( ss.str() );
    }
    
    return value[index];
}

std::string CharacterEventContinuous::getStateStr(void) const
{
    std::stringstream ss;
    // want to do this based on DiscreteDataType...
    //    ss << (char)(70+state);
    if ( value.size() == 1 )
    {
        ss << value[0];
    }
    else if ( value.size() > 1 )
    {
        ss << "[" << value[0];
        for (size_t i=1; i<value.size(); ++i)
        {
            ss << ", " << value[i];
        }
        ss << "]";
    }
    return ss.str();
}


void CharacterEventContinuous::resize( size_t new_size )
{

    value.resize( new_size );
    
}


void CharacterEventContinuous::setState(double s, size_t index)
{
    if ( index >= value.size() )
    {
        std::stringstream ss;
        ss << "Cannot set the " << (index+1) << "-th value of this continuous character event because it has only " << value.size() << " events.";
        
        throw RbException( ss.str() );
    }
    
    value[index] = s;
}

