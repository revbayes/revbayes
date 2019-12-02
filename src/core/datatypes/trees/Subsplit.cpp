#include "Subsplit.h"
#include "RbVectorUtilities.h"

#include <algorithm>
#include <iostream>
#include <sstream>
#include <cmath>


using namespace RevBayesCore;


/**
 * Default constructor required by the revlanguage code.
 */
Subsplit::Subsplit( void ) :
    bitset(),
    is_fake()
{

}

/**
 * Constructor from single taxon (as a clade)
 */
Subsplit::Subsplit( const std::vector<Taxon> &c1, const std::vector<Taxon> &n ) :
    bitset(),
    is_fake( true )
    {
        // TODO: having a single taxon instead of a vector of size 1 may be faster

        // Check that clade contains only single taxon
        if ( c1.size() > 1 )
        {
          throw(RbException("Cannot create fake subsplit from clade of multiple species."));
        }

        bool found = false;

        // bitset representation
        RbBitSet clade_1_bitset = RbBitSet(n.size(),false);
        for (size_t i=0; i<n.size(); ++i)
        {
          if (n[i] == c1[0])
          {
            clade_1_bitset.set(i);
            found = true;
            break;
          }
        }

        if (!found)
        {
          throw(RbException("Attempt to create subsplit with taxon not in vector of taxa."));
        }

        bitset.first = clade_1_bitset;
        bitset.second = clade_1_bitset;
    }


/**
 * Constructor from vectors of taxa.
 */
Subsplit::Subsplit( const std::vector<Taxon> &c1, const std::vector<Taxon> &c2, const std::vector<Taxon> &n ) :
    bitset(),
    is_fake( false )
{
    // bitset representations and check that X and Y are disjoint
    bool disjoint = true;

    std::vector<Taxon> c1_unset = c1;
    std::vector<Taxon> c2_unset = c2;
    RbBitSet clade_1_bitset = RbBitSet(n.size(),false);
    RbBitSet clade_2_bitset = RbBitSet(n.size(),false);

    for (size_t i=0; i<n.size(); ++i)
    {
      for (std::vector<Taxon>::iterator jt=c1_unset.begin(); jt!=c1_unset.end(); jt++)
      {
        if (n[i] == (*jt))
        {
          clade_1_bitset.set(i);
          c1_unset.erase(jt);
          break;
        }
      }
      for (std::vector<Taxon>::iterator jt=c2_unset.begin(); jt!=c2_unset.end(); jt++)
      {
        if (n[i] == (*jt))
        {
          clade_2_bitset.set(i);
          c2_unset.erase(jt);
          break;
        }
      }
      if (clade_1_bitset.isSet(i) && clade_2_bitset.isSet(i))
      {
        disjoint = false;
      }
      if (c1_unset.size() == 0 && c2_unset.size() == 0)
      {
        break;
      }
    }

    if ( c1_unset.size() != 0 || c2_unset.size() != 0 )
    {
      throw(RbException("Attempt to create subsplit with taxon not in vector of taxa."));
    }

    if (!disjoint)
    {
      throw(RbException("Cannot create subsplit from non-disjoint clades."));
    }

    if (clade_1_bitset.getNumberSetBits() == 0 || clade_2_bitset.getNumberSetBits() == 0)
    {
      // std::cout << "Invalid subsplit: " << clade_1_bitset << "|" << clade_2_bitset << std::endl;
      throw(RbException("Cannot create subsplit when one clade is empty."));
    }


    // Order clades
    if ( clade_1_bitset < clade_2_bitset )
    {
      bitset.first = clade_1_bitset;
      bitset.second = clade_2_bitset;
    }
    else
    {
      bitset.first = clade_2_bitset;
      bitset.second = clade_1_bitset;
    }

    // Check if we made a fake subsplit
    if ( clade_1_bitset.getNumberSetBits() == 1 && clade_2_bitset.getNumberSetBits() == 1 && clade_1_bitset.getFirstSetBit() == clade_2_bitset.getFirstSetBit() )
    {
      is_fake = true;
    }
}

/**
 * Constructor from bitsets.
 */
Subsplit::Subsplit( RbBitSet &clade_1_bitset, RbBitSet &clade_2_bitset ) :
    bitset(),
    is_fake( false )
{
    // Check if we are called to make a fake subsplit
    if ( clade_1_bitset.getNumberSetBits() == 1 && clade_2_bitset.getNumberSetBits() == 1 && clade_1_bitset.getFirstSetBit() == clade_2_bitset.getFirstSetBit() )
    {
      is_fake = true;
    }

    // bitset representations and check that X and Y are disjoint
    bool disjoint = true;
    for (size_t i=0; i<clade_1_bitset.size(); ++i)
    {
      if ( clade_1_bitset[i] && clade_2_bitset[i] )
      {
        disjoint = false;
        break;
      }
    }

    if (!disjoint && !is_fake)
    {
      throw(RbException("Cannot create subsplit from non-disjoint clades."));
    }

    if (clade_1_bitset.getNumberSetBits() == 0 || clade_2_bitset.getNumberSetBits() == 0)
    {
      std::cout << "Invalid subsplit: " << clade_1_bitset << "|" << clade_2_bitset << std::endl;
      throw(RbException("Cannot create subsplit when one clade is empty."));
    }

    // Order clades
    if ( clade_1_bitset < clade_2_bitset )
    {
      bitset.first = clade_1_bitset;
      bitset.second = clade_2_bitset;
    }
    else
    {
      bitset.first = clade_2_bitset;
      bitset.second = clade_1_bitset;
    }

}

Subsplit& Subsplit::operator=(const Subsplit &s)
{
  if (this != &s)
  {
    bitset  = s.bitset;
    is_fake = s.is_fake;
  }

  return *this;

}
/**
 * Equals operator, we outsource comparison to std::pair<RbBitSet,RbBitSet>.
 * Only if both clades in both subsplits are equal are two Subsplits are equal.
 */
bool Subsplit::operator==(const Subsplit &s) const
{
  bool tf;

  // If first set bits of both pairs are the same, the Subsplits may be equal (and we have to check)
  // If not, then we know that they must be different
  // Checking this first speeds the process if most comparisons are between non-equal Subsplits
  if ( s.getYBitset().getFirstSetBit() == bitset.first.getFirstSetBit() && s.getZBitset().getFirstSetBit() == bitset.second.getFirstSetBit() )
  {
    tf = (s.getBitset() == bitset);
  }
  else
  {
    tf = false;
  }
  return ( tf );
}


/**
 * Not equals operator that uses the equals operator.
 */
bool Subsplit::operator!=(const Subsplit &s) const
{
  bool tf;

  // If first set bits of both pairs are the same, the Subsplits may be equal (and we have to check)
  // If not, then we know that they must be different
  // Checking this first speeds the process if most comparisons are between non-equal Subsplits
  if ( s.getYBitset().getFirstSetBit() == bitset.first.getFirstSetBit() && s.getZBitset().getFirstSetBit() == bitset.second.getFirstSetBit() )
  {
    tf = (s.getBitset() != bitset);
  }
  else
  {
    tf = true;
  }
  return ( tf );
}

/**
 * Less than operator.
 */
bool Subsplit::operator<(const Subsplit &s) const
{
  return ( s.getBitset() < bitset );
}

/**
 * Less than or equal to operator
 */
bool Subsplit::operator<=(const Subsplit &s) const
{
  return ( s.getBitset() <= bitset );
}

/**
 * Greater than operator.
 */
bool Subsplit::operator>(const Subsplit &s) const
{
  return ( s.getBitset() > bitset );
}

/**
 * Greater than operator or equal to operator
 */
bool Subsplit::operator>=(const Subsplit &s) const
{
  return ( s.getBitset() >= bitset );
}


/**
 * Make subsplit into a single clade object
 *
 * \return    The clade.
 */
Clade Subsplit::asClade( const std::vector<Taxon>& taxa ) const
{

    RbBitSet total_bitset = asCladeBitset();
    Clade total = Clade(total_bitset, taxa);
    return total;
}

/**
 * Make subsplit into a bitset for a single clade object
 *
 * \return    The clade.
 */
RbBitSet Subsplit::asCladeBitset( void ) const
{
    RbBitSet y = bitset.first;
    RbBitSet z = bitset.second;
    RbBitSet total = y | z;

    return total;
}

/**
 * Make subsplit into a bitset for a split, which is like a clade but polarized
 *
 * \return    The clade.
 */
RbBitSet Subsplit::asSplitBitset( void ) const
{
    RbBitSet y = bitset.first;
    RbBitSet z = bitset.second;
    RbBitSet split = y | z;

    // Polarize if needed, for comparison we forbid splits with >50% 1s
    if ( split.getNumberSetBits() > size_t(std::floor(split.size()/2.0)) )
    {
      split = ~split;
    }

    // For a split with 50% 1s, we make sure the first bit is a 1 (00110 -> 11001)
    else if ( size_t(split.size() % 2) == 0 && (split.getNumberSetBits() == size_t(std::floor(split.size()/2.0))) )
    {
      if ( split[0] == 0 )
      {
        split = ~split;
      }
    }

    return split;
}

/*
 * Gives us all the parent-child subsplit pairs that we can get from the edge leading from parent to child subsplit for all rootings of this tree.
 * The edge must be internal.
 *
 * return: vector of parent-child pairs
 */
std::vector<std::pair<Subsplit,Subsplit> > Subsplit::doVirtualRootingNonRootParent(const Subsplit &parent, const Subsplit &child) const
{
    std::vector<std::pair<Subsplit,Subsplit> > all_orientations; // Each subsplit in this vector will be in the order parent-child

    std::pair<Subsplit,Subsplit> orientation; // Use this for each new s,t pair; will be in the order parent-child

    // Define the sets of taxa we need

    // Children clades of child subsplit
    RbBitSet s_y = child.getYBitset();
    RbBitSet s_z = child.getZBitset();

    // Child subsplit
    RbBitSet s = s_y | s_z;

    // Everyone but child subsplit
    RbBitSet not_s = s; ~not_s;

    // Whole clade of parent
    RbBitSet t = parent.asCladeBitset();

    // Sister group of s in t
    RbBitSet t_not_s = t & not_s;

    // All taxa not in the clade defined by parent
    RbBitSet not_t = t; ~not_t;

    // The root split induced by rooting to this edge
    Subsplit root_to_this_edge = Subsplit(s,not_s);

    // The child subsplit induced by rooting to a descendant of s
    Subsplit reversed_child = Subsplit(t_not_s,not_t);

    // std::cout << "In Subsplit::doVirtualRootingNonRootParent" << std::endl;
    // std::cout << "parent = " << parent << std::endl;
    // std::cout << "child = " << child << std::endl;
    // std::cout << "s_y = " << s_y << std::endl;
    // std::cout << "s_z = " << s_z << std::endl;
    // std::cout << "s = " << s << std::endl;
    // std::cout << "t = " << t << std::endl;
    // std::cout << "not_t = " << not_t << std::endl;
    // std::cout << "not_s = " << not_s << std::endl;
    // std::cout << "t_not_s = " << t_not_s << std::endl;

    // Case 1: as-is (virtual root is in direction of current rooting of tree)
    orientation.first = parent;
    orientation.second = child;
    all_orientations.push_back(orientation);

    // Case 2: virtual root in s-and-not-t (aka sister clade to s)
    orientation.first = Subsplit(not_t,s);
    all_orientations.push_back(orientation);

    // Case 3: virtual root on the edge between t and s, track s as descendant
    orientation.first = root_to_this_edge;
    all_orientations.push_back(orientation);

    // Cases 4 and 5 only apply if s is not a tip
    if ( s.getNumberSetBits() > 1 )
    {
      // Case 4: virtual root in s_y
      orientation.first = Subsplit(s_z,not_s);
      orientation.second = reversed_child;
      all_orientations.push_back(orientation);

      // Case 5: virtual root in s_z
      orientation.first = Subsplit(s_y,not_s);
      all_orientations.push_back(orientation);
    }
    else
    {
      Subsplit dummy_subsplit = Subsplit();
      orientation.first = dummy_subsplit;
      orientation.second = dummy_subsplit;

      all_orientations.push_back(orientation);
      all_orientations.push_back(orientation);
    }

    // Case 6: virtual root on the edge between t and s, track not-s as descendant
    orientation.first = root_to_this_edge;
    orientation.second = reversed_child;
    all_orientations.push_back(orientation);

    return all_orientations;
}

/*
 * Gives us all the parent-child subsplit pairs that we can get from the edge leading from parent to child subsplit for all rootings of this tree.
 * The edge must have the root as a parent.
 *
 * return: vector of parent-child pairs
 */
std::vector<std::pair<Subsplit,Subsplit> >  Subsplit::doVirtualRootingRootParent(const Subsplit &sister1, const Subsplit &sister2, const Subsplit &child) const
{
  std::vector<std::pair<Subsplit,Subsplit> > all_orientations; // Each subsplit in this vector will be in the order parent-child

  std::pair<Subsplit,Subsplit> orientation; // Use this for each new s,t pair; will be in the order parent-child

  // Define the sets of taxa we need

  // Children clades of child subsplit
  RbBitSet s_y = child.getYBitset();
  RbBitSet s_z = child.getZBitset();

  // Child subsplit
  RbBitSet s = s_y | s_z;

  // Everyone but child subsplit
  RbBitSet not_s = s; ~not_s;

  // Sister (sibling) 1's clade
  RbBitSet sib1 = sister1.asCladeBitset();

  // Sister (sibling) 2's clade
  RbBitSet sib2 = sister2.asCladeBitset();

  // The root split induced by rooting to this edge
  Subsplit root_to_this_edge = Subsplit(s,not_s);

  // The child subsplit induced by rooting to a descendant of s
  Subsplit reversed_child = Subsplit(sib1,sib2);

  // std::cout << "In Subsplit::getAllParentChildGivenNewRoot" << std::endl;
  // std::cout << "s = " << s << std::endl;
  // std::cout << "s_y = " << s_y << std::endl;
  // std::cout << "s_z = " << s_z << std::endl;
  // std::cout << "sister 1 = " << sib1 << std::endl;
  // std::cout << "sister 2 = " << sib2 << std::endl;

  // Case 1: virtual root is in sister1
  orientation.first = Subsplit(sib2,s);
  orientation.second = child;
  all_orientations.push_back(orientation);

  // Case 2: virtual root is in sister2
  orientation.first = Subsplit(sib1,s);
  all_orientations.push_back(orientation);

  // Case 3: virtual root on the edge subtending s, track s as descendant
  orientation.first = root_to_this_edge;
  all_orientations.push_back(orientation);

  // Cases 4 and 5 only apply if s is not a tip
  if ( s.getNumberSetBits() > 1 )
  {
    // Case 4: virtual root in s_y
    orientation.first = Subsplit(s_z,not_s);
    orientation.second = reversed_child;
    all_orientations.push_back(orientation);

    // Case 5: virtual root in s_z
    orientation.first = Subsplit(s_y,not_s);
    all_orientations.push_back(orientation);
  }
  else
  {
    Subsplit dummy_subsplit = Subsplit();
    orientation.first = dummy_subsplit;
    orientation.second = dummy_subsplit;

    all_orientations.push_back(orientation);
    all_orientations.push_back(orientation);
  }

  // Case 6: virtual root on the edge subtending s, track not-s as descendant
  orientation.first = root_to_this_edge;
  orientation.second = reversed_child;
  all_orientations.push_back(orientation);

  return all_orientations;
}

/**
 * Utility function to access paired bitset
 *
 * \return    The bitset pair.
 */
std::pair<RbBitSet,RbBitSet> Subsplit::getBitset( void ) const
{

    return bitset;
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself
 */
Subsplit* Subsplit::clone(void) const
{
    return new Subsplit(*this);
}

/**
 * Get clade Y.
 *
 * \return    Clade Y.
 */
Clade Subsplit::getY( const std::vector<Taxon>& taxa ) const
{
    Clade y = Clade(bitset.first, taxa);
    return y;
}

RbBitSet Subsplit::getYBitset( void ) const
{

    return bitset.first;
}

/**
 * Get clade Z.
 *
 * \return    Clade Z.
 */
Clade Subsplit::getZ( const std::vector<Taxon>& taxa ) const
{
    if (!is_fake)
    {
      Clade z = Clade(bitset.second, taxa);
      return z;
    }
    else
    {
      throw(RbException("Cannot access subsplit clade Z from a fake subsplit."));
    }
}

RbBitSet Subsplit::getZBitset( void ) const
{

    return bitset.second;
}

/**
 * Is subsplit s a child of this subsplit's clade Y?
 * \return    true/false
 */
bool Subsplit::isChildOfY(const Subsplit &s) const
{
    return s.asCladeBitset() == bitset.first;
}

/**
 * Is subsplit s a child of this subsplit's clade Z?
 * \return    true/false
 */
bool Subsplit::isChildOfZ(const Subsplit &s) const
{
    return s.asCladeBitset() == bitset.second;
}

/**
 * Is subsplit compatible with this one?
 * A subsplit s is compatible with another subsplit t if the clades in s sum to be one of the clades of t
 * This requires that s has all the taxa in t's Y or Z splits, no more and no fewer
 * \return    true/false
 */
bool Subsplit::isCompatible(const Subsplit &s) const
{
    if ( s.bitset.first.size() != bitset.first.size() )
    {
      throw(RbException("Cannot compare subsplits with unequal taxon sizes."));
    }

    // Check if s is equal to t's Y split
    RbBitSet s_total = s.asCladeBitset();
    if ( s_total == bitset.first )
    {
        return true;
    }
    // Check if s is equal to t's Z split
    else if (!is_fake && s_total == bitset.second) // Don't compare to clade Z if there is no clade Z
    {
      return true;
    }
    return false;
}

/**
 * Is subsplit a fake subsplit that contains only a single taxon?
 *
 * \return    true/false
 */
bool Subsplit::isFake() const
{
    return is_fake;
}

/**
 * If this clade (that this subsplit represents) is one side of a root split, make the full root subsplit.
 *
 * \return    true/false
 */
Subsplit Subsplit::rootSplitFromClade() const
{
  RbBitSet c = asCladeBitset();
  RbBitSet c2 = c; ~c2;

  return(Subsplit(c,c2));
}

/**
 * If this clade is one side of a root split, make the full root subsplit.
 *
 * \return    true/false
 */
// Subsplit Subsplit::rootSplitFromClade(RbBitSet &c) const
// {
//   RbBitSet c2 = c; ~c2;
//
//   return(Subsplit(c,c2));
// }

/**
 * Get the number of taxa contained in this Subsplit.
 *
 * \return       Size of the taxon name vector in both clades combined.
 */
size_t Subsplit::size(void) const
{
    if (is_fake)
    {
      return(bitset.first.getNumberSetBits());
    }
    else
    {
      return bitset.first.getNumberSetBits() + bitset.second.getNumberSetBits();
    }
}

/**
 * Are the two splits in a subsplit disjoint (nonoverlapping)?
 *
 * \return    true/false
 */
bool Subsplit::splitsAreDisjoint() const
{
    if ( !is_fake ) // A fake split cannot and need not be disjoint
    {
      for (size_t i=0; i < bitset.first.size(); ++i)
      {
        if (bitset.first[i] && bitset.second[i])
        {
          return false;
        }
      }
    }

    return true;
}


/**
 * Write the value of this Subsplit as a string.
 *
 * \return    A single string containing the entire Subsplit.
 */
std::string Subsplit::toString( void ) const
{
    std::string s;
    if (is_fake)
    {
      std::string y = "[";
      for (size_t i=0; i<bitset.first.size(); ++i)
      {
          y += ( bitset.first.isSet(i) ? "1" : "0");
      }
      y += "]";

      s = y + "|" + y;
    }
    else
    {
      std::string y = "[";
      for (size_t i=0; i<bitset.first.size(); ++i)
      {
          y += ( bitset.first.isSet(i) ? "1" : "0");
      }
      y += "]";

      std::string z = "[";
      for (size_t i=0; i<bitset.second.size(); ++i)
      {
          z += ( bitset.second.isSet(i) ? "1" : "0");
      }
      z += "]";

      s = y + "|" + z;
    }
    return s;
}

/*
 * Gives us the parent-child subsplit that will replace this one on the equivalent edge in a tree rerooted to a specific edge
 *
 * return: parent-child pair
 */
std::pair<Subsplit,Subsplit> Subsplit::virtualRoot(const Subsplit &parent, const Subsplit &child, int case_number) const
{

    std::pair<Subsplit,Subsplit> flipped; // This will be in the order parent-child

    // // Define the sets of taxa we need
    //
    // // Children clades of child subsplit
    // RbBitSet s_y = child.bitset.first;
    // RbBitSet s_z = child.bitset.second;
    //
    // // Child subsplit
    // RbBitSet s = s_y | s_z;
    //
    // // Everyone but child subsplit
    // RbBitSet not_s = ~s;
    //
    // // Whole clade of parent
    // RbBitSet t = parent.asCladeBitset();
    //
    // // Sister group of s in t
    // RbBitSet t_not_s = t & not_s;
    //
    // // All taxa not in the clade defined by parent
    // RbBitSet not_t = ~t;
    //
    // // Re-root to a branch in clade y of subsplit s
    // if (root_to == "s_y")
    // {
    //   flipped.first = Subsplit(s_z,not_s);
    //   flipped.second = Subsplit(t_not_s,not_t);
    // }
    // // Re-root to a branch in clade z of subsplit s
    // else if (root_to == "s_z")
    // {
    //   flipped.first = Subsplit(s_y,not_s);
    //   flipped.second = Subsplit(t_not_s,not_t);
    // }
    // // Re-root to a branch in clade not-s of subsplit t (sister of s in t)
    // else if (root_to == "t_{-s}")
    // {
    //   flipped.first = Subsplit(not_t,s);
    //   flipped.second = child;
    // }
    // // Re-root to this edge, track clade s as child
    // else if (root_to == "e1")
    // {
    //   flipped.first = Subsplit(s,not_s);
    //   flipped.second = child;
    // }
    // // Re-root to this edge, track clade of everyone but s as child
    // else if (root_to == "e2")
    // {
    //   flipped.first = Subsplit(s,not_s);
    //   flipped.second = Subsplit(t_not_s,not_t);
    // }
    // else
    // {
    //   throw(RbException("Invalid re-rooting in reorientParentChild"));
    // }

    return flipped;
}

std::ostream& RevBayesCore::operator<<(std::ostream& o, const Subsplit& x) {

    o << x.toString();

    return o;
}
