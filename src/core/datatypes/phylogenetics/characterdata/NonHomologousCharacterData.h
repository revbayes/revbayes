#ifndef NonHomologousCharacterData_H
#define NonHomologousCharacterData_H

#include <cstddef>
#include <vector>
#include <iosfwd>

#include "AbstractCharacterData.h"

namespace RevBayesCore {
class AbstractTaxonData;
    
    /**
     * Abstract class for all character data objects.
     *
     * This abstract class provides the base class for all character data objects.
     * A character data object contains a vector of TaxonData objects and defines additional
     * convenience functions to access the data stored in these TaxonData objects.
     * Note that this class is a pure interface and thus contains only pure virtual functions!
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2013-04-15, version 1.0
     */
    class NonHomologousCharacterData : public AbstractCharacterData {
        
    public:
        virtual                                ~NonHomologousCharacterData(void) {}
        
        // Overloaded operators
        virtual const AbstractTaxonData&        operator[](size_t i) const = 0;                                                                     //!< Subscript operator (const)

        // methods of the Cloneable interface
        virtual NonHomologousCharacterData*     clone(void) const = 0;
        
        // CharacterData functions
        virtual std::vector<size_t>             getNumberOfCharacters(void) const = 0;                                                              //!< Number of characters
        virtual size_t                          getNumberOfCharacters(size_t i) const = 0;                                                          //!< Number of characters
        bool                                    isHomologyEstablished(void) const { return  false; }                                                //!< Returns whether the homology of the characters has been established
        virtual std::string                     getDataType(void) const = 0;                                                    //!< Return the data type of this character data matrix
        
    protected:
                                                NonHomologousCharacterData(void) {}                                                                 //!< Constructor requires character type
        
    };
    
    // Global functions using the class
    std::ostream&                               operator<<(std::ostream& o, const NonHomologousCharacterData& x);                                   //!< Overloaded output operator
    
}

#endif
