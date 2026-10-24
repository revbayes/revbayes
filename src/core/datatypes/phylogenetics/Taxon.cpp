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

    // The first occurrence establishes an absent constraint. An explicit range, including
    // [0,0], still contributes to the interval covering all occurrences.
    age_range = age_range ? coveringInterval(*age_range, d) : d;
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
 * Get the occurrences for this taxon.
 *
 * \return    A map of occurrence age ranges to counts.
 */
const std::map<TimeInterval, size_t>& Taxon::getOccurrences( void ) const
{
    return occurrences;
}


/**
 * Get the explicitly established age range, if any, for this taxon.
 *
 * \return    The age range.
 */
const std::optional<TimeInterval>& Taxon::getAgeRange( void ) const
{
    return age_range;
}


// Update a copy before installing it, so validation does not first materialize an absent range.
void Taxon::setMinAge( double a )
{
    // COMPAT(age-range): A single-endpoint update starts from the old default [0,0].
    // Revisit whether callers must establish a complete interval first.
    auto range = age_range.value_or(TimeInterval(0.0, 0.0));
    range.setMin(a);
    age_range = range;
}


// Update a copy before installing it, so validation does not first materialize an absent range.
void Taxon::setMaxAge( double a )
{
    // COMPAT(age-range): A single-endpoint update starts from the old default [0,0].
    // Revisit whether callers must establish a complete interval first.
    auto range = age_range.value_or(TimeInterval(0.0, 0.0));
    range.setMax(a);
    age_range = range;
}



/**
* Get a JSON-formatted string description of this object.
*
* \return    The JSON-formatted string.
*/
const std::string Taxon::getJsonRespresentation(void) const {

    // COMPAT(age-range): Preserve the serialized endpoints; absence is not exported yet.
    const auto range = age_range.value_or(TimeInterval(0.0, 0.0));
    std::string jsonStr = "";
    jsonStr += "{\"Taxon\": {";
    jsonStr += "\"name\": \"" + name + "\", ";
    jsonStr += "\"speciesName\": \"" + species_name + "\", ";
    jsonStr += "\"TimeInterval\": {\"minAge\": " + std::to_string(range.getMin()) + ", ";
    jsonStr += "\"maxAge\": " + std::to_string(range.getMax()) + "}}";
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
