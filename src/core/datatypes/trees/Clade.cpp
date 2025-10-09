#include "Clade.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <utility>

#include "StringUtilities.h"
#include "RbVectorUtilities.h"
#include "RbException.h"

using namespace RevBayesCore;

using std::set;
using std::string;
using std::vector;

/**
 * Default constructor required by the revlanguage code.
 * We use an empty string and an age of 0.0 for
 * this default object.
 */
Clade::Clade( void )
{
    
}


/**
 * Constructor with a single taxon.
 */
Clade::Clade( const Taxon &t, const RbBitSet &b ) :
    bitset( b ),
    num_missing( b.size() > 1 ? b.size() - 1 : 0 )
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
    bitset( b ),
    num_missing( b.size() > n.size() ? int(b.size()) - int(n.size()) : 0 ),
    taxa( n )
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
    bitset( b ),
    clade_name( "" ),
    num_missing( b.size() > n.size() ? int(b.size()) - int(n.size()) : 0 )
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
    bitset( b ),
    num_missing( b.size() - b.count() )
{

    for (size_t i = 0; i < b.size(); i++)
    {
        if ( b.test(i) )
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
    
    // Sebastian (20210519): We currently do not compare the clade names anymore.
//    if ( c.clade_name != clade_name )
//    {
//        return false;
//    }
    
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
    
    // Sebastian (20210519): We don't compare the clade name anymore but only it's other members
//    if ( c.clade_name != clade_name )
//    {
//        return c.clade_name < clade_name;
//    }
    
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
 * Does the provided clade conflicts with this clade? A conflict is if the clades overlap but are non-nested.
 *
 * \param[in]    c    The clade.
 *
 * \return True/False whether the clades are nested.
 */
bool Clade::conflicts(const Clade& c) const
{
    return overlaps(c) && ( isNestedWithin(c) == false ) && ( c.isNestedWithin(*this) == false );
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
 * Get the clade name.
 *
 * \return       The name of the clade
 *
 */
const std::string& Clade::getCladeName(void) const
{
    return clade_name;
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

const std::vector<Clade>& Clade::getOptionalConstraints(void) const
{
    assert(optional_constraints.has_value());
    return *optional_constraints;
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
 * Is this clade an optional constraint (used with e.g. dnConstrainedTopology).
 * An "optional" constraint means that at least one of the referenced clades must be true.
 *
 * \return       The true/false value of whether the clade is an optional constraint.
 */
bool Clade::hasOptionalConstraints(void) const
{
    return optional_constraints.has_value();
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
 * Is the provided clade nested within me? It is only nested if all it's taxa are nested within me.
 *
 * \param[in]    c    Theclade.
 *
 * \return       True/False, if there is an overlap.
 */
bool Clade::isNestedWithin(const Clade& c) const
{
    size_t N = taxa.size();

    // do a quick test if the other clade has more taxa
    // in that case it could never be nested.
    if ( N < c.size() )
    {
        return false;
    }

    bool nested = true;
    for ( size_t i=0; i<c.size(); ++i)
    {
        const Taxon& t = c.getTaxon(i);
        bool found = false;
        for ( size_t j=0; j<N; ++j)
        {
            if ( taxa[j] == t )
            {
                found = true;
                break;
            }
        }
        
        if ( found == false )
        {
            nested = false;
            break;
        }
    }

    return nested;
}


/**
 * Does the provided clade overlaps with myself? An overlap is if we share at least one taxon.
 *
 * \param[in]    c    The clade.
 *
 * \return       True/False, if there is an overlap.
 */
bool Clade::overlaps(const Clade& c) const
{
    size_t N = taxa.size();
    for ( size_t i=0; i<c.size(); ++i)
    {
        const Taxon& t = c.getTaxon(i);
        for ( size_t j=0; j<N; ++j)
        {
            if ( taxa[j] == t )
            {
                return true;
            }
        }
    }

    return false;
}

/**
 * Compute the set of overlapping taxa between clade c and myself.
 *
 * \param[in]    c    The clade.
 *
 * \return       The set of overlapping taxa.
 *
 * We use this method to give a more helpful error message
 *  when we claim that clades conflict.
 */
set<Taxon> Clade::intersection(const Clade& c) const
{
    set<Taxon> both;
    size_t N = taxa.size();
    for ( size_t i=0; i<c.size(); ++i)
    {
        const Taxon& t = c.getTaxon(i);
        for ( size_t j=0; j<N; ++j)
        {
            if ( taxa[j] == t )
            {
                both.insert(t);
            }
        }
    }

    return both;
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
            throw RbException() << "Missing taxon with name '" << taxa[i].getName() << "'.";
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
 * Set the clade name.
 *
 * \param[in]    n      The new name
 *
 */
void Clade::setCladeName(const std::string& n)
{
    clade_name = n;
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
 * Set optional clade constraints, e.g. dnConstrainedTopology must satisfy one of any clades in the set of clades
 *
 * \param[in]   c   Vector of optional clade constraints
 *
 */
void Clade::setOptionalConstraints(const std::vector<Clade>& c)
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
    std::string s;

    if ( hasOptionalConstraints() == true )
    {
        vector<string> cstrings;
        for(auto& c: getOptionalConstraints())
        {
            cstrings.push_back(c.toString());
        }


        s = "(" + StringUtilities::join(cstrings, " OR ") + ")";
    }
    else
    {
        vector<string> tstrings;
        for (size_t i = 0; i < taxa.size(); ++i)
        {
            string t = taxa[i].getName();
            if ( std::find(mrca.begin(), mrca.end(), taxa[i]) != mrca.end() )
            {
                t += "*";
            }
            tstrings.push_back(t);
        }

        s  = "{" + StringUtilities::join(tstrings,", ") + "}";
    }

    if (isNegativeConstraint())
        s = "NOT " + s;

    return s;
}

namespace RevBayesCore
{

std::ostream& operator<<(std::ostream& o, const Clade& x) {
   
    o << x.toString();
   
    return o;
}

void Clade::setAges(const std::vector<Taxon>& taxa)
{
    size_t N = size();
    // set the ages of each of the taxa in the constraint
    for (size_t j = 0; j < N; ++j)
    {
        bool found = false;
        for (size_t i=0; i<taxa.size(); ++i)
        {
            const Taxon& t = taxa[i];
            if ( t.getName() == getTaxonName(j) )
            {
                setTaxonAge(j, t.getAge());
                found = true;
                break;
            }
        }
        if ( found == false )
        {
            throw RbException() << "set_ages_for_constraint: can't find taxon " << getTaxonName(j) << " in full taxon set!";
        }
    }
}

// std::sort needs a "strict weak ordering" that returns true
//    only if c1 is ordered with respect to c2, and c1 < c2.
//
// cladeWithin only provides a partial ordering, so we can't use std::sort with that.
// We would need to use a topological sort, which is more complex.
// Real numbers plus NaNs also provide only a partial ordering, so we can't just pretend that
//   an NaN is a number.
//
// Its not really clear how to handle negative constraints or optional constraints.
// They must all be equivalent, we decide to put them last.

// We really should separate the idea of constraints (CladeConstraint,
// NotCladeConstraint, OptionalCladeConstraint) from the idea of clades.

bool cladeBefore(const Clade& c1, const Clade& c2)
{
    if (&c1 == &c2)
        return false;

    // Negative constraints and Optional constraints are ordered after all regular constraints.
    else if (not (c1.isNegativeConstraint() or c1.hasOptionalConstraints()) and (c2.isNegativeConstraint() or c2.hasOptionalConstraints()))
        return true;
    else if (c1.isNegativeConstraint() or c1.hasOptionalConstraints() or c2.isNegativeConstraint() or c2.hasOptionalConstraints())
        return false;

    else
    {
        // Clades without ages are equivalent.
        // They have to come first or last, we pick first.
        bool c1_has_age = not std::isnan(c1.getAge());
        bool c2_has_age = not std::isnan(c2.getAge());

        if (not c1_has_age and     c2_has_age) return true;
        if (    c1_has_age and not c2_has_age) return false;

        return c1.getAge() < c2.getAge();
    }
}

bool cladeSmaller(const Clade& c1, const Clade& c2)
{
    if (&c1 == &c2)
        return false;

    // Negative constraints and Optional constraints are ordered after all regular constraints.
    else if (not (c1.isNegativeConstraint() or c1.hasOptionalConstraints()) and (c2.isNegativeConstraint() or c2.hasOptionalConstraints()))
        return true;
    else if (c1.isNegativeConstraint() or c1.hasOptionalConstraints() or c2.hasOptionalConstraints() or c2.hasOptionalConstraints())
        return false;

    else
        return c1.getTaxa().size() < c2.getTaxa().size();
}


bool cladeWithin(const Clade& c1, const Clade& c2)
{
    if (&c1 == &c2)
        return false;

    // Negative constraints and Optional constraints are ordered after all regular constraints.
    else if (not (c1.isNegativeConstraint() or c1.hasOptionalConstraints()) and (c2.isNegativeConstraint() or c2.hasOptionalConstraints()))
        return true;
    else if (c1.isNegativeConstraint() or c1.hasOptionalConstraints() or c2.isNegativeConstraint() or c2.hasOptionalConstraints())
        return false;

    // c2 is nested with c1
    else if (c1.isNestedWithin(c2))
        return false;
    // c1 is nested with c2, but not the reverse
    else if (c2.isNestedWithin(c1))
        return true;
    else if (c1.overlaps(c2))
    {
        auto both = c1.intersection(c2);
        RbException e;
        e<<"Cannot simulate tree: conflicting monophyletic clade constraints:\n"
         <<"  Clade1 = "<<c1<<"\n"
         <<"  Clade2 = "<<c2<<"\n"
         <<"  Overlap = ";
        for(auto& taxon: both)
            e<<taxon.getName()<<" ";
        throw e;
    }
    else
    {
        // unordered!
        return false;
    }
}

}

