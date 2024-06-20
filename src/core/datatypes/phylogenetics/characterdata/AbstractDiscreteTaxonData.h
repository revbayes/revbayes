#ifndef AbstractDiscreteTaxonData_H
#define AbstractDiscreteTaxonData_H

#include <cstddef>
#include <iosfwd>

#include "AbstractTaxonData.h"
#include "DiscreteCharacterState.h"
#include "Taxon.h"

namespace RevBayesCore {
    
    class NaturalNumbersState;
    template<typename charType>
    class DiscreteTaxonData;
    class CharacterState;
    
    /**
     * Abstract class for all taxon objects.
     *
     * This abstract class provides the base class for all taxon data objects.
     * A taxon data object contains a vector of character objects and defines additional
     * convenience functions to access the data.
     * Note that this class is a pure interface and thus contains only pure virtual functions!
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2013-04-15, version 1.0
     */
    class AbstractDiscreteTaxonData : public AbstractTaxonData {
        
    public:
        virtual                                                    ~AbstractDiscreteTaxonData(void) {}
        
        virtual DiscreteCharacterState&                             operator[](size_t i) = 0;                                           //!< Index op allowing change
        virtual const DiscreteCharacterState&                       operator[](size_t i) const = 0;                                     //!< Const index op
        
        virtual AbstractDiscreteTaxonData*                          clone(void) const = 0;
        
        // AbstractTaxonData functions
        virtual void                                                addCharacter(const CharacterState &newChar ) = 0;                   //!< Push back a new character
        virtual void                                                addCharacter(const DiscreteCharacterState &newChar ) = 0;           //!< Push back a new character
        virtual void                                                addCharacter(const CharacterState &newChar, bool tf) = 0;           //!< Push back a new character
        virtual void                                                addCharacter(const DiscreteCharacterState &newChar, bool tf) = 0;   //!< Push back a new character
        virtual DiscreteTaxonData<NaturalNumbersState>*             combineCharacters(const AbstractDiscreteTaxonData &d) const = 0;    //!< Concatenate sequences
        virtual void                                                concatenate(const AbstractTaxonData &d) = 0;                        //!< Concatenate sequences
        virtual void                                                concatenate(const AbstractDiscreteTaxonData &d) = 0;                //!< Concatenate sequences
        virtual const DiscreteCharacterState&                       getCharacter(size_t index) const = 0;                               //!< Get the character at position index
        virtual DiscreteCharacterState&                             getCharacter(size_t index) = 0;                                     //!< Get the character at position index (non-const to return non-const character)
        virtual size_t                                              getNumberOfCharacters(void) const = 0;                              //!< How many characters
        virtual double                                              getPercentageMissing(void) const = 0;                               //!< Returns the percentage of missing data for this sequence
        virtual std::string                                         getStateLabels(void) = 0;                                           //!< Get the possible state labels
        virtual bool                                                isCharacterResolved(size_t idx) const = 0;                          //!< Returns whether the character is fully resolved (e.g., "A" or "1.32") or not (e.g., "AC" or "?")
        
    protected:
        AbstractDiscreteTaxonData(const Taxon &t);                                                                                 //!< Default constructor
        
    };
    
    // Global functions using the class
    std::ostream&                       operator<<(std::ostream& o, const AbstractDiscreteTaxonData& x);          //!< Overloaded output operator
    
}

#endif
