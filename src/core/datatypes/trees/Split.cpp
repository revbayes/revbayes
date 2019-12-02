#include "Split.h"
#include "RbVectorUtilities.h"

#include <algorithm>
#include <iostream>
#include <sstream>


using namespace RevBayesCore;


/**
 * Default constructor required by the revlanguage code.
 */
Split::Split( void ) :
    bitset()
{

}

/**
 * Constructor from vector of taxa.
 */
Split::Split( const std::vector<Taxon> &c, const std::vector<Taxon> &n ) :
    bitset()
{

    std::vector<Taxon> c_unset = c;
    bitset = RbBitSet(n.size(),false);

    for (size_t i=0; i<n.size(); ++i)
    {
      for (std::vector<Taxon>::iterator jt=c_unset.begin(); jt!=c_unset.end(); jt++)
      {
        if (n[i] == (*jt))
        {
          bitset.set(i);
          c_unset.erase(jt);
          break;
        }
      }
    }

    if ( c_unset.size() != 0 )
    {
      throw(RbException("Attempt to create split with taxon not in vector of taxa."));
    }

    // Polarize if needed, for comparison we forbid splits with >50% 1s
    if ( bitset.getNumberSetBits() > size_t(std::floor(bitset.size()/2.0)) )


    {
      bitset = ~bitset;
    }
    // For a split with 50% 1s, we make sure the first bit is a 1 (00110 -> 11001)
    else if ( size_t(bitset.size() % 2) == 0 && (bitset.getNumberSetBits() == size_t(std::floor(bitset.size()/2.0))) )
    {
      if ( bitset[0] == 0 )
      {
        bitset = ~bitset;
      }
    }

}

/**
 * Constructor from bitset representation of vector of taxa.
 */
Split::Split( RbBitSet &c ) :
    bitset()
{

    bitset = c;

    // Polarize if needed, for comparison we forbid splits with >50% 1s
    if ( bitset.getNumberSetBits() > size_t(std::floor(bitset.size()/2.0)) )


    {
      bitset = ~bitset;
    }
    // For a split with 50% 1s, we make sure the first bit is a 1 (00110 -> 11001)
    else if ( size_t(bitset.size() % 2) == 0 && (bitset.getNumberSetBits() == size_t(std::floor(bitset.size()/2.0))) )
    {
      if ( bitset[0] == 0 )
      {
        bitset = ~bitset;
      }
    }

}

Split& Split::operator=(const Split &s)
{
  if (this != &s)
  {
    bitset  = s.bitset;
  }

  return *this;

}
/**
 * Equals operator, we outsource comparison to std::pair<RbBitSet,RbBitSet>.
 * Only if both clades in both subsplits are equal are two Splits are equal.
 */
bool Split::operator==(const Split &s) const
{
  return ( s.getBitset() == bitset );
}


/**
 * Not equals operator that uses the equals operator.
 */
bool Split::operator!=(const Split &s) const
{
  return ( s.getBitset() != bitset );
}

/**
 * Less than operator.
 */
bool Split::operator<(const Split &s) const
{
  return ( s.getBitset() < bitset );
}

/**
 * Less than or equal to operator
 */
bool Split::operator<=(const Split &s) const
{
  return ( s.getBitset() <= bitset );
}

/**
 * Greater than operator.
 */
bool Split::operator>(const Split &s) const
{
  return ( s.getBitset() > bitset );
}

/**
 * Greater than operator or equal to operator
 */
bool Split::operator>=(const Split &s) const
{
  return ( s.getBitset() >= bitset );
}

/**
 * Utility function to access paired bitset
 *
 * \return    The bitset pair.
 */
RbBitSet Split::getBitset( void ) const
{

    return bitset;
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself
 */
Split* Split::clone(void) const
{
    return new Split(*this);
}

/**
 * Is split compatible with this one?
 * \return    true/false
 */
bool Split::isCompatible(const Split &s) const
{
    if ( s.bitset.size() != bitset.size() )
    {
      throw(RbException("Cannot compare subsplits with unequal taxon sizes."));
    }

    // Check if splits are compatible, handled as 4-gamete test
    bool has_00 = false;
    bool has_01 = false;
    bool has_10 = false;
    bool has_11 = false;
    int patterns = 0;

    for (size_t i=0; i<bitset.size(); ++i)
    {
      if ( bitset[i] == 0 )
      {
        if ( s.bitset[i] == 0 )
        {
          if ( !has_00 )
          {
            has_00 = true;
            ++patterns;
          }
        }
        else
        {
          if ( !has_01 )
          {
            has_01 = true;
            ++patterns;
          }
        }
      }
      else
      {
        if ( s.bitset[i] == 0 )
        {
          if ( !has_10 )
          {
            has_10 = true;
            ++patterns;
          }
        }
        else
        {
          if ( !has_11 )
          {
            has_11 = true;
            ++patterns;
          }
        }
      }

      if ( patterns == 4 )
      {
        return false;
      }
    }

    return true;
}

/**
 * Get the number of taxa with 1s.
 *
 * \return       Size of the taxon name vector in both clades combined.
 */
size_t Split::size(void) const
{
      return bitset.getNumberSetBits();
}

/**
 * Write the value of this Split as a string.
 *
 * \return    A single string containing the entire Split.
 */
std::string Split::toString( void ) const
{
    std::string s = "[";

    for (size_t i=0; i<bitset.size(); ++i)
    {
        s += ( bitset.isSet(i) ? "1" : "0");
    }
    s += "]";

    return s;
}


std::ostream& RevBayesCore::operator<<(std::ostream& o, const Split& x) {

    o << x.toString();

    return o;
}
