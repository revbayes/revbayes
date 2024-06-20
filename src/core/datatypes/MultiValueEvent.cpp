#include "MultiValueEvent.h"

#include <ostream>
#include <string>

#include "RbException.h"

using namespace RevBayesCore;

MultiValueEvent::MultiValueEvent(void)
{
    
}


bool MultiValueEvent::operator==(const MultiValueEvent &mve) const
{
    
    if ( num_events != mve.num_events )
    {
        return false;
    }
    
    if ( names.size() != mve.names.size() || values.size() != mve.values.size() )
    {
        return false;
    }
    
    for (size_t i=0; i<values.size(); ++i)
    {
        if ( names[i] != mve.names[i] )
        {
            return false;
        }
        
        if ( values[i].size() != mve.values[i].size() )
        {
            return false;
        }
        
        for (size_t j=0; j<values[i].size(); ++j)
        {
            if ( values[i][j] != mve.values[i][j] )
            {
                return false;
            }
        }
        
    }
    
    return true;
}


bool MultiValueEvent::operator<(const MultiValueEvent &mve) const
{
    
    if ( num_events < mve.num_events )
    {
        return true;
    }
    
    if ( names.size() < mve.names.size() || values.size() < mve.values.size() )
    {
        return true;
    }
    
    for (size_t i=0; i<values.size(); ++i)
    {
        if ( names[i] < mve.names[i] )
        {
            return true;
        }
        
        if ( values[i].size() < mve.values[i].size() )
        {
            return true;
        }
        
        for (size_t j=0; j<values[i].size(); ++j)
        {
            if ( values[i][j] < mve.values[i][j] )
            {
                return true;
            }
        }
        
    }
    
    return false;
}


bool MultiValueEvent::operator!=(const MultiValueEvent &mve) const
{
    return (operator==(mve) == false);
}



void MultiValueEvent::addValues(const RbVector<double> &v, const std::string &n)
{
    names.push_back(n);
    values.push_back( v );
}



void MultiValueEvent::clear(void)
{
    names.clear();
    values.clear();
}



MultiValueEvent* MultiValueEvent::clone(void) const
{
    return new MultiValueEvent( *this );
}



const std::string& MultiValueEvent::getName(size_t i) const
{
    return names[i];
}



long MultiValueEvent::getNumberOfEvents(void) const
{
    return num_events;
}



size_t MultiValueEvent::getNumberOfValues(void) const
{
    return values.size();
}



RbVector<double>& MultiValueEvent::getValues(size_t i)
{
    return values[i];
}



const RbVector<double>& MultiValueEvent::getValues(size_t i) const
{
    return values[i];
}



RbVector<double>& MultiValueEvent::getValues(const std::string &n)
{
    size_t index = names.size();
    for (size_t i=0; i<names.size(); ++i)
    {
        if ( names[i] == n )
        {
            index = i;
            break;
        }
    }
    
    if ( index == names.size() )
    {
        throw RbException( "Could not find values with name '" + n + "'.");
    }
    
    return values[index];
}



const RbVector<double>& MultiValueEvent::getValues(const std::string &n) const
{
    size_t index = names.size();
    for (size_t i=0; i<names.size(); ++i)
    {
        if ( names[i] == n )
        {
            index = i;
            break;
        }
    }
    
    if ( index == names.size() )
    {
        throw RbException( "Could not find values with name '" + n + "'.");
    }
    
    return values[index];
}


void MultiValueEvent::setNumberOfEvents(long n)
{
    num_events = n;
}



void MultiValueEvent::setValues(const RbVector<double> &v, const std::string &n)
{
    
    size_t index = names.size();
    for (size_t i=0; i<names.size(); ++i)
    {
        if ( names[i] == n )
        {
            index = i;
            break;
        }
    }
    
    if ( index == names.size() )
    {
        throw RbException( "Could not find values with name '" + n + "'.");
    }
    
    values[index] = v;
}



std::ostream& RevBayesCore::operator<<(std::ostream& o, const RevBayesCore::MultiValueEvent& x)
{
    
    for ( size_t j=0; j<x.getNumberOfValues(); ++j )
    {
        const std::string &n = x.getName(j);
        const RbVector<double> &v = x.getValues(j);
        o << n << ":\t[";
        for (size_t i=0; i<v.size(); ++i)
        {
            o << v[i];
            if ( i != (v.size()-1) )
            {
                o << ", ";
            }
        }
        o << "]";
        if ( j != x.getNumberOfValues()-1 )
        {
            o << "\n";
        }
    }
    
    return o;
}
