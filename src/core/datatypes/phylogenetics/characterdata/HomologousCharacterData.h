#ifndef HomologousCharacterData_H
#define HomologousCharacterData_H

#include <cstddef>
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
    class HomologousCharacterData : public AbstractCharacterData {
        
    public:
        virtual                                ~HomologousCharacterData(void) {}
        
        // Overloaded operators
        virtual const AbstractTaxonData&        operator[](size_t i) const = 0;                                             //!< Subscript operator (const)
        
        // methods of the Cloneable interface
        virtual HomologousCharacterData*        clone(void) const = 0;
        
        // methods of the Serializable interface
        void                                    writeToFile(const path &dir, const std::string &fn) const = 0;

        // pure virtual CharacterData functions
        virtual void                            concatenate(const AbstractCharacterData &d, std::string type = "") = 0;         //!< Concatenate two sequences
        virtual void                            concatenate(const HomologousCharacterData &d, std::string type = "") = 0;       //!< Concatenate two sequences
        virtual void                            excludeAllCharacters(void) = 0;                                                 //!< Exclude all characters
        virtual void                            excludeCharacter(size_t i) = 0;                                                 //!< Exclude character
        virtual std::string                     getDataType(void) const = 0;                                                    //!< Return the data type of this character data matrix
        virtual size_t                          getNumberOfCharacters(void) const = 0;                                          //!< Number of characters
        virtual size_t                          getNumberOfIncludedCharacters(void) const = 0;                                  //!< Number of characters
        virtual void                            includeCharacter(size_t i) = 0;                                                 //!< Include character
        virtual bool                            isCharacterExcluded(size_t i) const = 0;                                        //!< Is the character excluded
        virtual bool                            isCharacterResolved(size_t txIdx, size_t chIdx) const = 0;                      //!< Returns whether the character is fully resolved (e.g., "A" or "1.32") or not (e.g., "AC" or "?")
        virtual bool                            isCharacterResolved(const std::string &tn, size_t chIdx) const = 0;             //!< Returns whether the character is fully resolved (e.g., "A" or "1.32") or not (e.g., "AC" or "?")
        virtual void                            removeExcludedCharacters(void) = 0;                                              //!< Remove all the excluded characters
        virtual void                            restoreCharacter(size_t i) = 0;                                                 //!< Restore character
        
        // CharacterData functions
        bool                                    isHomologyEstablished(void) const { return true; }                              //!< Returns whether the homology of the characters has been established
        
    protected:
                                                HomologousCharacterData(void) {}                                                //!< Constructor requires character type
        
    };
    
    // Global functions using the class
    std::ostream&                               operator<<(std::ostream& o, const HomologousCharacterData& x);                    //!< Overloaded output operator
    
}

#endif
