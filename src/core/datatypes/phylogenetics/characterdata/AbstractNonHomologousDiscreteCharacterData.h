#ifndef AbstractNonHomologousDiscreteCharacterData_H
#define AbstractNonHomologousDiscreteCharacterData_H

#include <cstddef>
#include <vector>
#include <iosfwd>

#include "NonHomologousCharacterData.h"

namespace RevBayesCore {
class AbstractDiscreteTaxonData;
class AbstractTaxonData;
class DiscreteCharacterState;
    
    /**
     * Abstract class for all DISCRETE character data objects.
     *
     * The DISCRETE character data class is derived from the interface character data
     * and simply specifies that all derived classes must contain discrete characters.
     * This simplifies return values because it is known that all are derived from discrete characters
     * and there are several function that only work on discrete characters and not on
     * continuous characters.
     *
     * Note that this class is a pure interface and thus contains only pure virtual functions!
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-02-16, version 1.0
     */
    class AbstractNonHomologousDiscreteCharacterData : public NonHomologousCharacterData {
        
        
    public:
        
        virtual                                                    ~AbstractNonHomologousDiscreteCharacterData(void) {}
        
        // Overloaded operators
        virtual const AbstractTaxonData&                            operator[](size_t i) const = 0;                                             //!< Subscript operator (const)
        
        // methods of the Cloneable interface
        virtual AbstractNonHomologousDiscreteCharacterData*         clone(void) const = 0;
        
        // CharacterData functions
        virtual const DiscreteCharacterState&                       getCharacter(size_t tn, size_t cn) const = 0;                               //!< Return a reference to a character element in the character matrix
        virtual std::string                                         getDataType(void) const = 0;                                                //!< Return the data type of this character data matrix
        virtual std::vector<size_t>                                 getNumberOfCharacters(void) const = 0;                                                              //!< Number of characters
        virtual size_t                                              getNumberOfCharacters(size_t i) const = 0;                                                          //!< Number of characters
        virtual AbstractDiscreteTaxonData&                          getTaxonData(size_t tn) = 0;                                                //!< Return a reference to a sequence in the character matrix
        virtual const AbstractDiscreteTaxonData&                    getTaxonData(size_t tn) const = 0;                                          //!< Return a reference to a sequence in the character matrix
        virtual AbstractDiscreteTaxonData&                          getTaxonData(const std::string &tn) = 0;                                    //!< Return a reference to a sequence in the character matrix
        virtual const AbstractDiscreteTaxonData&                    getTaxonData(const std::string &tn) const = 0;                              //!< Return a reference to a sequence in the character matrix
        
    protected:
        AbstractNonHomologousDiscreteCharacterData() {}                                          //!< Constructor requires character type
        
    };
    
    // Global functions using the class
    std::ostream&                                   operator<<(std::ostream& o, const AbstractNonHomologousDiscreteCharacterData& x);        //!< Overloaded output operator
    
}

#endif
