#include "RbBitSet.h"

#include <string>
#include <sstream> // IWYU pragma: keep

#include "RbException.h"

using namespace RevBayesCore;



RbBitSet::RbBitSet(void) :
    num_set_bits( 0 )
{

}


RbBitSet::RbBitSet(size_t n, bool def) :
    value(n, def),
    num_set_bits( 0 )
{

}


bool RbBitSet::operator[](size_t i) const
{
    // get the internal value
    return value[i];
}


/** Equals comparison */
bool RbBitSet::operator==(const RbBitSet& x) const
{

    return x.value == value;
}

/** Not-Equals comparison */
bool RbBitSet::operator!=(const RbBitSet& x) const
{

    return operator==(x) == false;
}


/** Less than comparison */
bool RbBitSet::operator<(const RbBitSet& x) const
{

    return x.value < value;
}

/** Less than or equals to comparison */
bool RbBitSet::operator<=(const RbBitSet& x) const
{

    return x.value <= value;
}

/** Greater than comparison */
bool RbBitSet::operator>(const RbBitSet& x) const
{

    return x.value > value;
}

/** Greater than or equals to comparison */
bool RbBitSet::operator>=(const RbBitSet& x) const
{

    return !(x.value < value);
}

/** Bitwise and */
RbBitSet RbBitSet::operator&(const RbBitSet& x) const
{
    if (x.value.size() != value.size())
    {
        throw(RbException("Cannot and RbBitSets of unequal size"));
    }
    RbBitSet r(value.size());
    for (size_t i = 0; i < value.size(); i++)
    {
        if ( value[i] && x.value[i] )
        {
            r.set(i);
        }
    }
    return r;
}

/** Bitwise or */
RbBitSet RbBitSet::operator|(const RbBitSet& x) const
{
    if (x.value.size() != value.size())
    {
        throw(RbException("Cannot or RbBitSets of unequal sizes"));
    }
    RbBitSet r(value.size());
    for (size_t i = 0; i < value.size(); i++)
    {
        if ( value[i] || x.value[i] )
        {
            r.set(i);
        }
    }
    return r;
}

/** Bitwise xor */
RbBitSet RbBitSet::operator^(const RbBitSet& x) const
{
    if (x.value.size() != value.size())
    {
        throw(RbException("Cannot xor RbBitSets of unequal size"));
    }
    RbBitSet r(value.size());
    for (size_t i = 0; i < value.size(); i++)
    {
        if ( value[i] != x.value[i] )
        {
            r.set(i);
        }
    }
    return r;
}

/** Unary not */
RbBitSet& RbBitSet::operator~()
{
    for (size_t i = 0; i < value.size(); i++)
    {
        value[i] = !value[i];
    }

    num_set_bits = value.size() - num_set_bits;

    return *this;
}

/** Bitwise and assignment */
RbBitSet& RbBitSet::operator&=(const RbBitSet& x)
{
    if (x.value.size() != value.size())
    {
        throw(RbException("Cannot and RbBitSets of unequal size"));
    }

    *this = *this & x;

    return *this;
}

/** Bitwise or assignment */
RbBitSet& RbBitSet::operator|=(const RbBitSet& x)
{
    if (x.value.size() != value.size())
    {
        throw(RbException("Cannot or RbBitSets of unequal size"));
    }

    *this = *this | x;

    return *this;
}


void RbBitSet::clear(void)
{
    // reset the bitset
    value = std::vector<bool>(value.size(),false);
    num_set_bits = 0;
}

bool RbBitSet::empty(void) const
{
    return value.empty();
}

void RbBitSet::flip(size_t i)
{
    value[i] = ( value[i] == false );
    num_set_bits += ( value[i] ? 1 : -1 );
}

void RbBitSet::flip()
{
    for (size_t i = 0; i < value.size(); i++)
    {
        value[i] = ( value[i] == false );
        num_set_bits += ( value[i] ? 1 : -1 );
    }
}

size_t RbBitSet::getFirstSetBit( void ) const
{
    size_t index = 0;
    while ( index < value.size() && value[index] == false )
    {
        ++index;
    }

    return index;
}

size_t RbBitSet::getHash( void )
{
    return std::hash<std::vector<bool> >{}(value);
}

size_t RbBitSet::getNumberSetBits( void ) const
{
    // get the internal value
    return num_set_bits;
}


bool RbBitSet::isSet(size_t i) const
{
    // get the internal value
    return value[i];
}

void RbBitSet::resize(size_t size)
{
    value.resize(size);
}

void RbBitSet::set(size_t i)
{

    if ( i >= value.size() )
    {
        std::ostringstream ss;
        ss << i;
        throw RbException("Index " + ss.str() +" out of bounds in bitset. This will likely cause unexpected behavior.");
    }

    if (value[i] == false)
    {
        ++num_set_bits;
    }

    // set the internal value
    value[i] = true;
}


size_t RbBitSet::size(void) const
{
    // get the size from the actual bitset
    return value.size();
}


void RbBitSet::unset(size_t i)
{
    if (value[i] == true)
    {
        --num_set_bits;
    }

    // set the internal value
    value[i] = false;
}


std::ostream& RevBayesCore::operator<<(std::ostream& o, const RevBayesCore::RbBitSet& x)
{

    o << "[";
    for (size_t i=0; i<x.size(); ++i)
    {
        o << ( x.isSet(i) ? "1" : "0");
    }
    o << "]";


    return o;
}
