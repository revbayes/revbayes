#ifndef DiscreteTaxonData_H
#define DiscreteTaxonData_H

#include "AbstractDiscreteTaxonData.h"
#include "DiscreteCharacterState.h"
#include "RbOptions.h"


#include <string>
#include <vector>

namespace RevBayesCore {

    template<class charType>
    class DiscreteTaxonData : public AbstractDiscreteTaxonData {

    public:
                                                        DiscreteTaxonData(const Taxon &t);                                  //!< Set type spec of container from type of elements

        charType&                                       operator[](size_t i);                                               //!< Index op allowing change
        const charType&                                 operator[](size_t i) const;                                         //!< Const index op

        // implemented methods of the Cloneable interface
        DiscreteTaxonData<charType>*                    clone(void) const;

        // TaxonData functions
        void                                            addCharacter(const CharacterState &newChar );                       //!< Push back a new character
        void                                            addCharacter(const DiscreteCharacterState &newChar );               //!< Push back a new character
        void                                            addCharacter(const charType &newChar );                             //!< Push back a new character
        void                                            addCharacter(const CharacterState &newChar, bool tf);               //!< Push back a new character
        void                                            addCharacter(const DiscreteCharacterState &newChar, bool tf);       //!< Push back a new character
        void                                            addCharacter(const charType &newChar, bool tf);                     //!< Push back a new character
        DiscreteTaxonData<NaturalNumbersState>*         combineCharacters(const AbstractDiscreteTaxonData &d) const;        //!< Concatenate sequences
        DiscreteTaxonData<NaturalNumbersState>*         combineCharacters(const DiscreteTaxonData &d) const;                //!< Concatenate sequences
        void                                            concatenate(const AbstractTaxonData &d);                            //!< Concatenate sequences
        void                                            concatenate(const AbstractDiscreteTaxonData &d);                    //!< Concatenate sequences
        void                                            concatenate(const DiscreteTaxonData &d);                            //!< Concatenate sequences
        const charType&                                 getCharacter(size_t index) const;                                   //!< Get the character at position index
        charType&                                       getCharacter(size_t index);                                         //!< Get the character at position index (non-const to return non-const character)
        std::string                                     getJsonRepresentation(void) const;
        size_t                                          getNumberOfCharacters(void) const;                                  //!< How many characters
        double                                          getPercentageMissing(void) const;                                   //!< Returns the percentage of missing data for this sequence
        std::string                                     getStringRepresentation(size_t idx) const;
        std::string                                     getStateLabels(void);                                               //!< Get the possible state labels
        bool                                            isCharacterResolved(size_t idx) const;                              //!< Returns whether the character is fully resolved (e.g., "A" or "1.32") or not (e.g., "AC" or "?")
        bool                                            isSequenceMissing(void) const;                                      //!< Returns whether the contains only missing data or has some actual observations
        void                                            removeCharacters(const std::set<size_t> &i);                        //!< Remove all the characters with a given index
        void                                            setAllCharactersMissing(void);                                      //!< Set all characters as missing
        
    private:

        std::vector<charType>                           sequence;
        RbBitSet                                        is_resolved;

    };

    // Global functions using the class
    template<class charType>
    std::ostream&                                       operator<<(std::ostream& o, const DiscreteTaxonData<charType>& x);  //!< Overloaded output operator

}


#include "CharacterState.h"
#include "RbException.h"
#include "NaturalNumbersState.h"


/**
 * Constructor with taxon name.
 * Does nothing except instanciating the object.
 */
template<class charType>
RevBayesCore::DiscreteTaxonData<charType>::DiscreteTaxonData(const Taxon &t) : AbstractDiscreteTaxonData( t ),
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
template<class charType>
charType& RevBayesCore::DiscreteTaxonData<charType>::operator[](size_t i)
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
template<class charType>
const charType& RevBayesCore::DiscreteTaxonData<charType>::operator[](size_t i) const
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
template<class charType>
RevBayesCore::DiscreteTaxonData<charType>* RevBayesCore::DiscreteTaxonData<charType>::clone( void ) const
{

    return new DiscreteTaxonData<charType>(*this);
}


/**
 * Add another character data object to this character data object.
 *
 * \param[in]    obsd    The CharacterData object that should be added.
 */
template<class charType>
RevBayesCore::DiscreteTaxonData<RevBayesCore::NaturalNumbersState>* RevBayesCore::DiscreteTaxonData<charType>::combineCharacters(const AbstractDiscreteTaxonData &obsd) const
{
    
    const DiscreteTaxonData<charType>* rhs = dynamic_cast<const DiscreteTaxonData<charType>* >( &obsd );
    if ( rhs == NULL )
    {
        throw RbException("Combining wrong character data type into TaxonData!!!");
    }
    
    
    return combineCharacters( *rhs );
}


/**
 * Add another character data object to this character data object.
 *
 * \param[in]    obsd    The CharacterData object that should be added.
 */
template<class charType>
RevBayesCore::DiscreteTaxonData<RevBayesCore::NaturalNumbersState>* RevBayesCore::DiscreteTaxonData<charType>::combineCharacters(const DiscreteTaxonData<charType> &obsd) const
{
    
    charType first_char             = getCharacter( 0 );
    charType second_char            = obsd.getCharacter( 0 );
    
    size_t num_states_first         = first_char.getNumberOfStates();
    size_t num_states_second        = second_char.getNumberOfStates();
    
    NaturalNumbersState new_state   = NaturalNumbersState(num_states_first * num_states_second);
    
    int k=0;
    bool first_state_set = false;
    for ( size_t i=0; i<num_states_first; ++i)
    {
        k = int(i*num_states_second);
        if ( first_char.isStateSet( i ) == true || first_char.isMissingState() == true || first_char.isGapState() == true )
        {
            for ( size_t j=0; j<num_states_second; ++j)
            {
                if ( second_char.isStateSet( j ) == true || second_char.isMissingState() == true || second_char.isGapState() == true )
                {
                    if ( first_state_set == false )
                    {
                        new_state.setStateByIndex(k);
                        first_state_set = true;
                    }
                    else
                    {
                        new_state.addState(k);
                    }
                }
                ++k;
            }
        }
    }

    DiscreteTaxonData<NaturalNumbersState> *new_taxon_data = new DiscreteTaxonData<NaturalNumbersState>( taxon );
    new_taxon_data->addCharacter( new_state );
    
    return new_taxon_data;
}



/**
 * Add another character data object to this character data object.
 *
 * \param[in]    obsd    The CharacterData object that should be added.
 */
template<class charType>
void RevBayesCore::DiscreteTaxonData<charType>::concatenate(const AbstractTaxonData &obsd)
{

    const DiscreteTaxonData<charType>* rhs = dynamic_cast<const DiscreteTaxonData<charType>* >( &obsd );
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
template<class charType>
void RevBayesCore::DiscreteTaxonData<charType>::concatenate(const AbstractDiscreteTaxonData &obsd)
{

    const DiscreteTaxonData<charType>* rhs = dynamic_cast<const DiscreteTaxonData<charType>* >( &obsd );
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
template<class charType>
void RevBayesCore::DiscreteTaxonData<charType>::concatenate(const DiscreteTaxonData<charType> &obsd)
{

    sequence.insert( sequence.end(), obsd.sequence.begin(), obsd.sequence.end() );

}



/**
 * Push back a new character.
 *
 * \param[in]    newChar    The new character.
 */
template<class charType>
void RevBayesCore::DiscreteTaxonData<charType>::addCharacter( const CharacterState &newChar )
{

#   ifdef ASSERTIONS_ALL
    if ( dynamic_cast<const charType&>( newChar ) == NULL )
    {
        throw RbException("Inserting wrong character type into TaxonData!!!");
    }
#   endif

    addCharacter( static_cast<const charType &>(newChar) );
}


/**
 * Push back a new character.
 *
 * \param[in]    newChar    The new character.
 */
template<class charType>
void RevBayesCore::DiscreteTaxonData<charType>::addCharacter( const DiscreteCharacterState &newChar )
{

#   ifdef ASSERTIONS_ALL
    if ( dynamic_cast<const charType&>( newChar ) == NULL )
    {
        throw RbException("Inserting wrong character type into TaxonData!!!");
    }
#   endif

    addCharacter( static_cast<const charType &>(newChar) );

}


/**
 * Push back a new character.
 *
 * \param[in]    newChar    The new character.
 */
template<class charType>
void RevBayesCore::DiscreteTaxonData<charType>::addCharacter( const charType &newChar )
{

    sequence.push_back( newChar );
    is_resolved.push_back(true);
}


/**
 * Push back a new character.
 *
 * \param[in]    newChar    The new character.
 */
template<class charType>
void RevBayesCore::DiscreteTaxonData<charType>::addCharacter( const CharacterState &newChar, bool tf )
{

#   ifdef ASSERTIONS_ALL
    if ( dynamic_cast<const charType&>( newChar ) == NULL )
    {
        throw RbException("Inserting wrong character type into TaxonData!!!");
    }
#   endif

    addCharacter( static_cast<const charType &>(newChar), tf );
}


/**
 * Push back a new character.
 *
 * \param[in]    newChar    The new character.
 */
template<class charType>
void RevBayesCore::DiscreteTaxonData<charType>::addCharacter( const DiscreteCharacterState &newChar, bool tf )
{

#   ifdef ASSERTIONS_ALL
    if ( dynamic_cast<const charType&>( newChar ) == NULL )
    {
        throw RbException("Inserting wrong character type into TaxonData!!!");
    }
#   endif

    addCharacter( static_cast<const charType &>(newChar), tf );

}


/**
 * Push back a new character.
 *
 * \param[in]    newChar    The new character.
 */
template<class charType>
void RevBayesCore::DiscreteTaxonData<charType>::addCharacter( const charType &newChar, bool tf )
{

    sequence.push_back( newChar );
    is_resolved.push_back(tf);
}



/**
 * Push back a new character.
 *
 * \param[in]    index    The position character.
 */
template<class charType>
charType& RevBayesCore::DiscreteTaxonData<charType>::getCharacter(size_t index)
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
 * \param[in]    index  The position of the character.
 *
 * \return              A non-const reference to the character
 */
template<class charType>
const charType& RevBayesCore::DiscreteTaxonData<charType>::getCharacter(size_t index) const
{

    if (index >= sequence.size())
    {
        throw RbException("Index out of bounds");
    }

    return sequence[index];
}


/**
 * Get the number of character stored in this object
 *
 * \return            The number of characters
 */
template<class charType>
size_t RevBayesCore::DiscreteTaxonData<charType>::getNumberOfCharacters(void) const
{

    return sequence.size();
}


/**
* Obtain a JSON-formatted string describing the object.
*
* \return            JSON-formatted string.
*/
template<class charType>
std::string RevBayesCore::DiscreteTaxonData<charType>::getJsonRepresentation(void) const {

    std::string jsonStr = "";
    
    jsonStr += "{\"DiscreteTaxonData\": ";
    jsonStr += taxon.getJsonRespresentation();
    jsonStr += ", \"charData\": [";
    for (int i=0; i<sequence.size(); i++)
        {
        jsonStr += "\"" + sequence[i].getStringValue() + "\"";
        if (i + 1 < sequence.size())
            jsonStr += ",";
        }
    jsonStr += "]";
    
    jsonStr += "}}";
    
    return jsonStr;
}


/**
 * Computes the percentage of the sequences that is missing.
 *
 * \return            Percentage of missing characters.
 */
template<class charType>
double RevBayesCore::DiscreteTaxonData<charType>::getPercentageMissing( void ) const
{
    double numMissing = 0.0;
    for (size_t i = 0; i < sequence.size(); ++i)
    {
        if ( sequence[i].isMissingState() == true || sequence[i].isGapState() == true )
        {
            ++numMissing;
        }
    }

    return numMissing / sequence.size();
}


template<class charType>
std::string RevBayesCore::DiscreteTaxonData<charType>::getStateLabels(void)
{

    if (sequence.size() == 0)
    {
        return "";
    }

    return sequence[0].getStateLabels();
}



template<class charType>
std::string RevBayesCore::DiscreteTaxonData<charType>::getStringRepresentation(size_t idx) const
{

    return sequence[idx].getStringValue();
}


/**
 * Determines whether the character at idx is "resolved," meaning there is a fully-resolved character state at that position.
 *
 * \return            True (resolved) or false (unresolved).
 */
template<class charType>
bool RevBayesCore::DiscreteTaxonData<charType>::isCharacterResolved(size_t idx) const
{

    if (idx >= is_resolved.size())
    {
        throw RbException("Index out of bounds");
    }

    return is_resolved[idx];
}


/**
 * Determines whether the sequences completely missing.
 *
 * \return            True (missing) or false (observed).
 */
template<class charType>
bool RevBayesCore::DiscreteTaxonData<charType>::isSequenceMissing( void ) const
{

    for (size_t i = 0; i < sequence.size(); ++i)
    {
        if ( sequence[i].isMissingState() == false && sequence[i].isGapState() == false )
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
template<class charType>
void RevBayesCore::DiscreteTaxonData<charType>::removeCharacters(const std::set<size_t> &idx)
{

//    size_t alreadyRemoved = 0;
//    for (std::set<size_t>::const_iterator it = idx.begin(); it != idx.end(); ++it)
//    {
//        size_t index = *it - alreadyRemoved;
//        sequence.erase(sequence.begin() + index);
//        ++alreadyRemoved;
//    }

    std::vector<charType> included;
    for (size_t i = 0; i < sequence.size(); ++i)
    {
        if ( idx.find(i) == idx.end() )
        {
            included.push_back( sequence[i] );
        }
    }

    sequence = included;

}


/**
 * Determines whether the sequences completely missing.
 *
 */
template<class charType>
void RevBayesCore::DiscreteTaxonData<charType>::setAllCharactersMissing( void )
{

    for (size_t i = 0; i < sequence.size(); ++i)
    {
        sequence[i].setMissingState( true );
    }

}

template<class charType>
std::ostream& RevBayesCore::operator<<(std::ostream& o, const DiscreteTaxonData<charType>& x)
{

    o << x.getTaxonName() << ":" << std::endl;
    for (size_t i = 0; i < x.getNumberOfCharacters(); ++i)
    {
        o << x[i];
    }
    o << std::endl;

    return o;
}


#endif
