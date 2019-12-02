#include "Taxon.h"

using namespace RevBayesCore;


/**
 * Default constructor.
 */
Taxon::Taxon( void ) :
    age_range(  ),
    name( "" ),
    species_name( "" )
{
    
}


/**
 * Constructor simply initiating the object and its members.
 *
 * \param[in]    n     The name of the taxon.
 */
Taxon::Taxon(const std::string &n) :
    age_range(  ),
    name( n ),
    species_name( n )
{
    
}


/**
 * Equals operator.
 * We check the species name and the individuals name.
 */
bool Taxon::operator==(const RevBayesCore::Taxon &t) const
{
    
    if ( species_name != t.species_name )
    {
        return false;
    }
    
    if ( name != t.name)
    {
        return false;
    }
    
//    if ( age_range != t.age_range)
//    {
//        return false;
//    }

    return true;
}


/**
 * Not equals operator. We simply invert the result of the equals operation.
 */
bool Taxon::operator!=(const RevBayesCore::Taxon &t) const
{
    
    return !operator==(t);
}


/**
 * Less-than operator.
 * We check first the species name and then the indidivuals name.
 */
bool Taxon::operator<(const RevBayesCore::Taxon &t) const
{
    
    if ( species_name < t.species_name)
    {
        return true;
    }
    else if ( species_name > t.species_name )
    {
        return false;
    }
    
    if ( name < t.name )
    {
        return true;
    }
    else if ( name > t.name )
    {
        return false;
    }
    
    // by default return false.
    return false;
}



/**
 * Less-than or equals operator.
 */
bool Taxon::operator<=(const RevBayesCore::Taxon &t) const
{
    
    return operator<(t) || operator==(t);
}


/**
 * Greater-than operator.
 * We check first the species name and then the indidivuals name.
 */
bool Taxon::operator>(const RevBayesCore::Taxon &t) const
{
    
    return operator<=(t) == false;
}



/**
 * Greater-than or equals operator.
 */
bool Taxon::operator>=(const RevBayesCore::Taxon &t) const
{
    
    return operator>(t) || operator==(t);
}


/**
 * Get the age for this taxon.
 *
 * \return    The age.
 */
double Taxon::getAge( void ) const
{
    return age_range.getMin();
}


/**
 * Get the date info for this taxon.
 *
 * \return    The date.
 */
const TimeInterval& Taxon::getAgeRange( void ) const
{
    return age_range;
}


/**
 * Get the name info for this taxon.
 *
 * \return    The name.
 */
const std::string& Taxon::getName( void ) const
{
    return name;
}


/**
 * Get the species name for this taxon.
 *
 * \return    The species name.
 */
const std::string& Taxon::getSpeciesName( void ) const
{
    return species_name;
}


/**
 * Set the age for this taxon.
 *
 * \param[in]    a     The age.
 */
void Taxon::setAge(double a)
{
    age_range.setMin(a);
    age_range.setMax(a);
}


/**
 * Set the date info for this taxon.
 *
 * \param[in]    d     The date.
 */
void Taxon::setAgeRange( const TimeInterval &d )
{
    age_range = d;
}


/**
 * Set the name info for this taxon.
 *
 * \param[in]    n     The name.
 */
void Taxon::setName( const std::string &n )
{
    name = n;
}


/**
 * Set the species name for this taxon.
 *
 * \param[in]    n     The species name.
 */
void Taxon::setSpeciesName( const std::string &sn )
{
    species_name = sn;
}


std::ostream& RevBayesCore::operator<<(std::ostream& o, const Taxon& x)
{
    // Sebastian: We will not write out the species name or date anymore
    // These info need to be queried specifically
//    o << x.getName() << ":" << x.getspecies_name() << ":" << x.getDate();
    o << x.getName();
    return o;
}
