#include "ContinuousTaxonData.h"

#include <cstdio>
#include <ostream>

#include "RbException.h"
#include "RbMathLogic.h"
#include "Cloneable.h"
#include "RbConstants.h" // IWYU pragma: keep

using namespace RevBayesCore;


///**
// * Default constructor.
// * Does nothing except instanciating the object.
// */
//ContinuousTaxonData::ContinuousTaxonData(void) : AbstractTaxonData( Taxon("") ),
//    sequence()
//{
//    
//}


/**
 * Constructor with taxon name.
 * Does nothing except instanciating the object.
 */
ContinuousTaxonData::ContinuousTaxonData(const Taxon &t) : AbstractTaxonData( t ),
    sequence()
{
    
}


/**
 * Subscript const operator for convenience access.
 *
 * \param[in]    i    The position of the character.
 *
 * \return            A non-const reference to the character
 */
double& ContinuousTaxonData::operator[](size_t i)
{
    
    if (i >= sequence.size())
    {
        throw RbException("Index out of bounds");
    }
    
    return sequence[i];
}


/**
 * Subscript const operator for convenience access.
 *
 * \param[in]    i    The position of the character.
 *
 * \return            A const reference to the character
 */
const double& ContinuousTaxonData::operator[](size_t i) const
{
    
    if (i >= sequence.size())
    {
        throw RbException("Index out of bounds");
    }
    
    return sequence[i];
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the object.
 */
ContinuousTaxonData* ContinuousTaxonData::clone( void ) const
{
    
    return new ContinuousTaxonData(*this);
}


/**
 * Add another character data object to this character data object.
 *
 * \param[in]    obsd    The CharacterData object that should be added.
 */
void ContinuousTaxonData::concatenate(const AbstractTaxonData &obsd)
{
    
    const ContinuousTaxonData* rhs = dynamic_cast<const ContinuousTaxonData* >( &obsd );
    if ( rhs == NULL )
    {
        throw RbException("Adding wrong character data type into TaxonData!!!");
    }
    
    
    concatenate( *rhs );
}


/**
 * Add another character data object to this character data object.
 *
 * \param[in]    obsd    The CharacterData object that should be added.
 */
void ContinuousTaxonData::concatenate(const ContinuousTaxonData &obsd)
{
    
    sequence.insert( sequence.end(), obsd.sequence.begin(), obsd.sequence.end() );
    
}


/**
 * Push back a new character.
 * 
 * \param[in]    newChar    The new character.
 */
void ContinuousTaxonData::addCharacter(const double &newChar)
{
    
    sequence.push_back( newChar );
    isResolved.push_back( true );
}


/**
 * Push back a new character.
 * 
 * \param[in]    newChar    The new character.
 */
void ContinuousTaxonData::addCharacter(const double &newChar, const bool tf)
{
    
    sequence.push_back( newChar );
    isResolved.push_back( tf );
}


/**
 * Push back a new character.
 * 
 * \param[in]    index    The position character.
 */
double& ContinuousTaxonData::getCharacter(size_t index)
{
    
    if (index >= sequence.size())
    {
        throw RbException("Index out of bounds");
    }
    
    return sequence[index];
}


/**
 * Get-operator for convenience access.
 *
 * \param[in]    index    The position of the character.
 *
 * \return            A const reference to the character
 */
const double& ContinuousTaxonData::getCharacter(size_t index) const
{
    
    if (index >= sequence.size())
    {
        throw RbException("Index out of bounds");
    }
    
    return sequence[index];
}


/**
* Get a JSON-formatted string description of this object
*
* \return            A JSON-formatted string
*/
std::string ContinuousTaxonData::getJsonRepresentation(void) const {

    std::string jsonStr = "";
    
    jsonStr += "{\"ContinuousTaxonData\": {";
    jsonStr += taxon.getJsonRespresentation();
    jsonStr += ", \"charData\": [";
    for (int i=0; i<sequence.size(); i++)
        {
        char tempCStr[20];
        snprintf(tempCStr, sizeof(tempCStr), "%.10e", sequence[i]);
        std::string tempStr = tempCStr;
        jsonStr += tempStr;
        if (i + 1 < sequence.size())
            jsonStr += ",";
        }
    jsonStr += "]";
    jsonStr += "}}";
    
    return jsonStr;
}

/**
 * Get the number of character stored in this object
 *
 * \return            The number of characters
 */
size_t ContinuousTaxonData::getNumberOfCharacters(void) const 
{
    
    return sequence.size();
}


/**
 * Computes the percentage of the sequences that is missing.
 *
 * \return            Percentage of missing characters.
 */
double ContinuousTaxonData::getPercentageMissing( void ) const
{
    double numMissing = 0.0;
    for (size_t i = 0; i < sequence.size(); ++i)
    {
        if ( RevBayesCore::RbMath::isNan(sequence[i]) )
        {
            ++numMissing;
        }
    }
    
    return numMissing / sequence.size();
}

std::string ContinuousTaxonData::getStateLabels(void)
{

    return "";
}

std::string ContinuousTaxonData::getStringRepresentation(size_t idx) const
{

    if ( RevBayesCore::RbMath::isNan(sequence[idx]) )
    {
        return "-";
    }
    
    char tempCStr[20];
    snprintf(tempCStr, sizeof(tempCStr), "%1.2lf", sequence[idx]);
    std::string tempStr = tempCStr;
    return tempStr;
}


bool ContinuousTaxonData::isCharacterResolved(size_t idx) const
{

    if (idx >= isResolved.size())
    {
        throw RbException("Index out of bounds");
    }
    
    return isResolved[idx];
}


/**
 * Determines whether the sequences completely missing.
 *
 * \return            True (missing) or false (observed).
 */
bool ContinuousTaxonData::isSequenceMissing( void ) const
{
    
    for (size_t i = 0; i < sequence.size(); ++i)
    {
        if ( RevBayesCore::RbMath::isNan(sequence[i]) == false )
        {
            return false;
        }
    }
    
    return true;
}


/**
 * Remove characters.
 *
 */
void ContinuousTaxonData::removeCharacters(const std::set<size_t> &idx)
{
    
    size_t alreadyRemoved = 0;
    for (std::set<size_t>::const_iterator it = idx.begin(); it != idx.end(); ++it)
    {
        size_t index = *it + alreadyRemoved;
        sequence.erase(sequence.begin() + index);
        ++alreadyRemoved;
    }
}


/**
 * Determines whether the sequences completely missing.
 *
 * \return            True (missing) or false (observed).
 */
void ContinuousTaxonData::setAllCharactersMissing( void )
{
    
    for (size_t i = 0; i < sequence.size(); ++i)
    {
        sequence[i] = RbConstants::Double::nan;
        isResolved[i] = false;
    }
    
}


std::ostream& RevBayesCore::operator<<(std::ostream& o, const RevBayesCore::ContinuousTaxonData& x) {
    
    o << x.getTaxonName() << ":" << std::endl;
    for (size_t i = 0; i < x.getNumberOfCharacters(); ++i) 
    {
        if ( i > 0 )
        {
            o << ", ";
        }
        o << x[i];
    }
    o << std::endl;
    
    return o;
}

