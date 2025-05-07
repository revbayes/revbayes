#ifndef ContinuousTaxonData_H
#define ContinuousTaxonData_H

#include <cstddef>
#include <vector>
#include <iosfwd>
#include <set>

#include "AbstractTaxonData.h"
#include "Taxon.h"

namespace RevBayesCore {
    
    class ContinuousTaxonData : public AbstractTaxonData {
        
    public:
//                                              ContinuousTaxonData(void);                                          //!< Set type spec of container from type of elements
                                                ContinuousTaxonData(const Taxon &t);                                //!< Set type spec of container from type of elements
        
        double&                                 operator[](size_t i);                                               //!< Index op allowing change
        const double&                           operator[](size_t i) const;                                         //!< Const index op
        
        // implemented methods of the Cloneable interface
        ContinuousTaxonData*                    clone(void) const;
        
        // TaxonData functions
        void                                    addCharacter(const double &newChar);                                //!< Push back a new character
        void                                    addCharacter(const double &newChar, const bool tf);                 //!< Push back a new character
//      void                                    combineCharacters(const AbstractTaxonData &d);                      //!< Concatenate sequences
//      void                                    combineCharacters(const ContinuousTaxonData &d);                    //!< Concatenate sequences
        void                                    concatenate(const AbstractTaxonData &d);                            //!< Concatenate sequences
        void                                    concatenate(const ContinuousTaxonData &d);                          //!< Concatenate sequences
        const double&                           getCharacter(size_t index) const;                                   //!< Get the character at position index
        double&                                 getCharacter(size_t index);                                         //!< Get the character at position index (non-const to return non-const character)
        std::string                             getJsonRepresentation(void) const;
        size_t                                  getNumberOfCharacters(void) const;                                  //!< How many characters
        double                                  getPercentageMissing(void) const;                                   //!< Returns the percentage of missing data for this sequence
        std::string                             getStringRepresentation(size_t idx) const;
        std::string                             getStateLabels(void);                                               //!< Get the possible state labels
        bool                                    isCharacterResolved(size_t idx) const;                              //!< Returns whether the character is fully resolved (e.g., "A" or "1.32") or not (e.g., "AC" or "?")
        bool                                    isSequenceMissing(void) const;                                      //!< Returns whether the contains only missing data or has some actual observations
        void                                    removeCharacters(const std::set<size_t> &i);                        //!< Remove all the characters with a given index
        void                                    setAllCharactersMissing(void);                                      //!< Set all characters as missing
        
    private:
        
        std::vector<double>                     sequence;
        std::vector<bool>                       isResolved;
        
    };
    
    // Global functions using the class
    std::ostream&                               operator<<(std::ostream& o, const ContinuousTaxonData& x);          //!< Overloaded output operator
    
}

#endif
