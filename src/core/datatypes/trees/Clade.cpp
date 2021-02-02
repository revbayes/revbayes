#include "Clade.h"

#include <algorithm>
#include <iostream>
#include <string>
#include <utility>

#include "RbVectorUtilities.h"
#include "RbException.h"

using std::vector;
using std::set;

using namespace RevBayesCore;


/**
 * Default constructor required by the revlanguage code.
 * We use an empty string and an age of 0.0 for
 * this default object.
 */
Clade::Clade( void ) :
    age( 0.0 ),
    num_missing( 0 ),
    taxa(),
    is_negative_constraint(false),
    is_optional_match(false)
{
    
}


/**
 * Constructor with a single taxon.
 */
Clade::Clade( const Taxon &t, const RbBitSet &b ) :
    age( 0.0 ),
    bitset( b ),
    num_missing( b.size() > 1 ? b.size() - 1 : 0 ),
    taxa(),
    is_negative_constraint(false),
    is_optional_match(false)
{
    
    taxa.push_back( t );
}


/**
 * Default constructor that instantiates the object.
 * Additionally, we sort the vector of taxon names.
 *
 * \param[in]   n    The vector containing the taxon names.
 */
Clade::Clade(const std::vector<Taxon> &n, const RbBitSet &b) :
    age( 0.0 ),
    bitset( b ),
    num_missing( b.size() > n.size() ? int(b.size()) - int(n.size()) : 0 ),
    taxa( n ),
    is_negative_constraint(false),
    is_optional_match(false)
{
    
    VectorUtilities::sort( taxa );
}


/**
 * Default constructor that instantiates the object.
 * Additionally, we sort the vector of taxon names.
 *
 * \param[in]   n    The vector containing the taxon names.
 */
Clade::Clade(const std::set<Taxon> &n, const RbBitSet &b) :
    age( 0.0 ),
    bitset( b ),
    num_missing( b.size() > n.size() ? int(b.size()) - int(n.size()) : 0 ),
    taxa(),
    is_negative_constraint(false),
    is_optional_match(false)
{
    
    for (std::set<Taxon>::const_iterator it=n.begin(); it!=n.end(); ++it)
    {
        taxa.push_back( *it );
    }
    
    VectorUtilities::sort( taxa );
}


/**
 * Default constructor that instantiates the object.
 * Additionally, we sort the vector of taxon names.
 *
 * \param[in]   n    The vector containing the taxon names.
 */
Clade::Clade(const RbBitSet &b, const std::vector<Taxon> &n) :
    age( 0.0 ),
    bitset( b ),
    num_missing( b.size() - b.getNumberSetBits() ),
    is_negative_constraint(false),
    is_optional_match(false)
{

    for (size_t i = 0; i < b.size(); i++)
    {
        if ( b.isSet(i) )
        {
            taxa.push_back(n[i]);
        }
    }
}


/**
 * Overloaded equals operator.
 * Only if we have the exact same taxon names then these two clades are equal.
 */
bool Clade::operator==(const Clade &c) const
{
    
    if ( c.size() != taxa.size() )
    {
        return false;
    }
    
    // Sebastian (10/19/2015): We cannot use the clade age for comparison because
    //                         otherwise we cannot find the same clade in different trees.
//    if ( c.getAge() != age )
//    {
//        return false;
//    }
    
    for (size_t i = 0; i < taxa.size(); ++i)
    {
        if ( taxa[i] != c.getTaxon(i) )
        {
            return false;
        }
    }
    
    return true;
}


/**
 * Not equals operator that uses the equals operator.
 */
bool Clade::operator!=(const Clade &c) const 
{
    return !operator==( c );
}


/**
 * Less than operator so that we can sort the clades.
 */
bool Clade::operator<(const Clade &c) const
{
    
    if ( taxa.size() < c.size() )
    {
        return true;
    }
    else if ( taxa.size() > c.size() )
    {
        return false;
    }
    
    for (size_t i = 0; i < taxa.size(); ++i)
    {
        
        if ( taxa[i] < c.getTaxon(i) )
        {
            return true;
        }
        else if ( taxa[i] > c.getTaxon(i) )
        {
            return false;
        }
        
    }
    
    return false;
}


/**
 * Less than operator so that we can sort the clades.
 */
bool Clade::operator<=(const Clade &c) const
{
    return operator<( c ) || operator==( c );
}


/**
 * Less than operator so that we can sort the clades.
 */
bool Clade::operator>(const Clade &c) const
{
    return operator<( c ) == false && operator==( c ) == false;
}


/**
 * Less than operator so that we can sort the clades.
 */
bool Clade::operator>=(const Clade &c) const
{
    return operator>( c ) == false;
}



/**
 * Get the const-iterator to the first taxon name.
 */
std::vector<Taxon>::const_iterator Clade::begin(void) const
{
    return taxa.begin();
}


/**
 * Get the iterator to the first taxon name.
 */
std::vector<Taxon>::iterator Clade::begin(void)
{
    return taxa.begin();
}


/**
 * Get the const-iterator after the last taxon name.
 */
std::vector<Taxon>::const_iterator Clade::end(void) const
{
    return taxa.end();
}


/**
 * Get the iterator after the last taxon name.
 */
std::vector<Taxon>::iterator Clade::end(void)
{
    return taxa.end();
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself 
 */
Clade* Clade::clone(void) const 
{
    return new Clade(*this);
}


/**
 * Add a taxon to the list.
 *
 * \param[in]    t    The new taxon.
 */
void Clade::addTaxon(const Taxon &t)
{
    taxa.push_back( t );
}


/**
 * Get the clade age.
 *
 * \return    The stored age.
 */
double Clade::getAge( void ) const 
{
    
    return age;
}


/**
  * Get all taxa as a bitset.
  *
  * \return       The bitset.
  */
const RbBitSet& Clade::getBitRepresentation( void ) const
{
    return bitset;
}


/**
 * Get the mrca taxon.
 *
 * \return       The mrca taxon
 *
 */
const std::set<Taxon>& Clade::getMrca(void) const
{
    return mrca;
}


/**
 * Get number of missing taxa.
 *
 * \return       The number of missing taxa.
 */
int Clade::getNumberMissingTaxa( void ) const
{
    return num_missing;
}

/**
 * Get number of missing taxa.
 *
 * \return       The number of missing taxa.
 */
size_t Clade::getNumberOfTaxa( void ) const
{
    return taxa.size();
}


/**
 * Get the vector of optional clade constraints.
 *
 * \return       The optional clade constraints
 */

std::vector<Clade> Clade::getOptionalConstraints(void) const
{
    return optional_constraints;
}


/**
 * Get all taxon names.
 *
 * \return       The vector of taxon names.
 */
std::vector<Taxon>& Clade::getTaxa( void )
{
    return taxa;
}


/**
 * Get all taxon names.
 *
 * \return       The vector of taxon names.
 */
const std::vector<Taxon>& Clade::getTaxa( void ) const
{
    return taxa;
}


/**
 * Get the taxon name at position i.
 *
 * \param[in]    i    The index for the taxon name we are interested in.
 *
 * \return       The name of the taxon.
 */
const Taxon& Clade::getTaxon(size_t i) const
{
    return taxa[i];
}


/**
 * Get the taxon name at position i.
 *
 * \param[in]    i    The index for the taxon name we are interested in.
 *
 * \return       The name of the taxon.
 */
const std::string& Clade::getTaxonName(size_t i) const
{
    return taxa[i].getName();
}

/**
 * Is this clade a negative clade constraint (used with e.g. dnConstrainedTopology).
 *
 * \return       The true/false value of whether the clade is a negative constraint.
 */
bool Clade::isNegativeConstraint(void) const
{
    return is_negative_constraint;
}


/**
 * Is this clade a negative clade constraint (used with e.g. dnConstrainedTopology).
 *
 * \return       The true/false value of whether the clade is a negative constraint.
 */
bool Clade::isOptionalMatch(void) const
{
    return is_optional_match;
}


/**
 * Reset the bitset.
 *
 * \param[in]    map    The map between taxon names and index.
 *
 */
void Clade::resetTaxonBitset(const std::map<std::string, size_t> map)
{
    
    bitset = RbBitSet( map.size() );
    
    for (size_t i = 0; i < taxa.size(); ++i)
    {
        const std::map<std::string, size_t >::const_iterator& it = map.find( taxa[i].getName() );
        if ( it != map.end() )
        {
            bitset.set( it->second );
        }
        else
        {
            throw RbException("Missing taxon with name '" + taxa[i].getName() + "'.");
        }
    }
    
}



/**
 * Set the age of the clade.
 *
 * \param[in]    a  The age of the clade.
 *
 */
void Clade::setAge(double a)
{
    age = a;
}


/**
 * Set the bitset of the clade.
 *
 * \param[in]    b  The bitset representation of this clade.
 *
 */
void Clade::setBitRepresentation( const RbBitSet &b )
{
    bitset = b;
}


/**
 * Set the mrca taxa. Must be taxa already contained in the clade.
 *
 * \param[in]    t      The taxa to be set as the mrca
 *
 */
void Clade::setMrca(const std::set<Taxon>& t)
{
    mrca = t;
}


/**
 * Set the number of missing taxa.
 *
 * \param[in]    n      The number of missing taxa.
 *
 */
void Clade::setNumberMissingTaxa(int n)
{
    num_missing = n;
}


/**
 * Set the taxon age at position i.
 *
 * \param[in]    i    The index for the taxon we are interested in.
 * \param[in]    age  The age of the taxon to set.
 *
 */
void Clade::setTaxonAge(size_t i, double age)
{
    taxa[i].setAge(age);
}

/**
 * Set flag for negative clade constraint.
 *
 * \param[in]    tf   Flag indicating if clade is a negative constraint.
 *
 */
void Clade::setNegativeConstraint(bool tf)
{
    is_negative_constraint = tf;
}

/**
 * Set flag for negative clade constraint.
 *
 * \param[in]    tf   Flag indicating if clade is a negative constraint.
 *
 */
void Clade::setOptionalMatch(bool tf)
{
    is_optional_match = tf;
}



/**
 * Set optional clade constraints, e.g. dnConstrainedTopology must satisfy one of any clades in the set of clades
 *
 * \param[in]   c   Vector of optional clade constraints
 *
 */
void Clade::setOptionalConstraints(std::vector<Clade> c)
{
    optional_constraints = c;
}


/**
 * Get the number of taxa contained in this clade.
 *
 * \return       Size of the taxon name vector.
 */
size_t Clade::size(void) const 
{
    return taxa.size();
}


/**
 * Write the value of this clade as a string.
 *
 * \return    A single string containing the entire clade.
 */
std::string Clade::toString( void ) const
{
    std::string s = "{";
    
    for (size_t i = 0; i < taxa.size(); ++i)
    {
        if ( i > 0 )
        {
            s += ",";
        }
        s += taxa[i].getName();
        if ( std::find(mrca.begin(), mrca.end(), taxa[i]) != mrca.end() )
        {
            s += "*";
        }
    }
    s += "}";
    
    return s;
}

namespace RevBayesCore
{

std::ostream& operator<<(std::ostream& o, const Clade& x) {
   
    o << x.toString();
   
    return o;
}

void set_ages_for_constraint_top(Clade& clade, const vector<Taxon>& taxa)
{
    // set the ages of each of the taxa in the constraint
    for (size_t j = 0; j < clade.size(); ++j)
    {
        bool found = false;
        for (auto& taxon: taxa)
        {
            if ( taxon.getName() == clade.getTaxonName(j) )
            {
                clade.setTaxonAge(j, taxon.getAge());
                found = true;
                break;
            }
        }
        if (not found)
            throw RbException("set_ages_for_constraint: can't find taxon " + clade.getTaxonName(j) + " in full taxon set!");
    }
}

void set_ages_for_constraint(Clade& clade, const vector<Taxon>& taxa)
{
    // set the ages of each of the taxa in the constraint
    set_ages_for_constraint_top( clade, taxa );

    // set ages for optional constraints
    std::vector<Clade> optional_constraints = clade.getOptionalConstraints();
    for (auto& optional_constraint: optional_constraints)
        set_ages_for_constraint_top( optional_constraint, taxa );

    clade.setOptionalConstraints( optional_constraints );
}

bool clade_nested_within(const Clade& clade1, const Clade& clade2)
{
    set<Taxon> taxa1;
    for(auto& taxon: clade1.getTaxa())
        taxa1.insert(taxon);

    set<Taxon> taxa2;
    for(auto& taxon: clade2.getTaxa())
        taxa2.insert(taxon);

    return std::includes(clade2.begin(), clade2.end(), clade1.begin(), clade1.end());
}

bool clades_overlap(const Clade& clade1, const Clade& clade2)
{
    set<Taxon> taxa1;
    for(auto& taxon: clade1.getTaxa())
        taxa1.insert(taxon);

    set<Taxon> taxa2;
    for(auto& taxon: clade2.getTaxa())
        taxa2.insert(taxon);

    auto i = taxa1.begin();
    auto j = taxa2.begin();
    while (i != taxa1.end() && j != taxa2.end())
    {
      if (*i == *j)
        return true;
      else if (*i < *j)
        ++i;
      else
        ++j;
    }
    return false;
}

bool clades_conflict(const Clade& clade1, const Clade& clade2)
{
    return clades_overlap(clade1, clade2) and
        (not clade_nested_within(clade1, clade2)) and
        (not clade_nested_within(clade2, clade2));
}

}
