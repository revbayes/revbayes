#include "Taxon.h"

#include <string>

using namespace RevBayesCore;


/**
 * Default constructor.
 */
Taxon::Taxon( void ) :
    name( "" ),
    species_name( "" ),
    extinct(false)
{
    
}


/**
 * Constructor simply initiating the object and its members.
 *
 * \param[in]    n     The name of the taxon.
 */
Taxon::Taxon(const std::string &n) :
    name( n ),
    species_name( n ),
    extinct(false)
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
 * Add an occurrence age for this taxon, with uncertainty.
 *
 * \param[in]    d     The age range.
 */
void Taxon::addOccurrence( const TimeInterval &d )
{
    occurrences[d]++;

    if ( d.getMax() > age_range.getMax() )
    {
        age_range.setMax( d.getMax() );
    }
    if ( d.getMin() < age_range.getMin() )
    {
        age_range.setMin( d.getMin() );
    }
}


/**
 * Get the extinct status flag.
 *
 * \return    The flag
 */
bool Taxon::isExtinct(void) const
{
    return extinct;
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
 * Get the occurrences for this taxon.
 *
 * \return    A map of occurrence age ranges to counts.
 */
const std::map<TimeInterval, size_t>& Taxon::getOccurrences( void ) const
{
    return occurrences;
}


/**
 * Get the age range for this taxon across all occurrences.
 *
 * \return    The age range.
 */
const TimeInterval& Taxon::getAgeRange( void ) const
{
    return age_range;
}


/**
 * Get the max age for this taxon across all occurrences.
 *
 * \return    The age.
 */
double Taxon::getMaxAge( void ) const
{
    return age_range.getMax();
}


/**
 * Get the min age for this taxon across all occurrences.
 *
 * \return    The age.
 */
double Taxon::getMinAge( void ) const
{
    return age_range.getMin();
}

void Taxon::setMinAge( double a )
{
    age_range.setMin( a );
}


void Taxon::setMaxAge( double a )
{
    age_range.setMax( a );
}



/**
* Get a JSON-formatted string description of this object.
*
* \return    The JSON-formatted string.
*/
const std::string Taxon::getJsonRespresentation(void) const {

    std::string jsonStr = "";
    jsonStr += "{\"Taxon\": {";
    jsonStr += "\"name\": \"" + name + "\", ";
    jsonStr += "\"speciesName\": \"" + species_name + "\", ";
    jsonStr += "\"TimeInterval\": {\"minAge\": " + std::to_string(age_range.getMin()) + ", ";
    jsonStr += "\"maxAge\": " + std::to_string(age_range.getMax()) + "}}";
    return jsonStr;
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
    
    if ( age_range.getMax() < a )
    {
        age_range.setMax(a);
        age_range.setMin(a);
    }
    else
    {
        age_range.setMin(a);
        // Sebastian (20210519): We should not automatically set the max age
        // because otherwise we always get a range of size 0.
//        age_range.setMax(a);
    }
}


/**
 * Set the age range for this taxon.
 *
 * \param[in]    d     The age range.
 */
void Taxon::setAgeRange( const TimeInterval &d )
{
    age_range = d;
}


/**
 * Set the extinct status flag.
 *
 * \param[in]    a     The age.
 */
void Taxon::setExtinct(bool e)
{
    extinct = e;
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
 * \param[in]    sn     The species name.
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
