#ifndef AbstractTaxonData_H
#define AbstractTaxonData_H

#include <cstddef>
#include <set>
#include <iosfwd>

#include "Cloneable.h"
#include "Taxon.h"

namespace RevBayesCore {
    
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
    class AbstractTaxonData : public Cloneable {
    
    public:
        virtual                                ~AbstractTaxonData(void) {}
        
        // methods of the Cloneable interface
        virtual AbstractTaxonData*              clone(void) const = 0;
        
        bool                                    operator==(const AbstractTaxonData &rm) const { return this == &rm; }
        bool                                    operator!=(const AbstractTaxonData &rm) const { return !operator==(rm); }
        bool                                    operator<(const AbstractTaxonData &rm) const { return this < &rm; }
        bool                                    operator<=(const AbstractTaxonData &rm) const { return operator<(rm) || operator==(rm); }

        
        // AbstractTaxonData functions
        virtual void                            concatenate(const AbstractTaxonData &d) = 0;                        //!< Concatenate sequences
        virtual std::string                     getJsonRepresentation(void) const = 0;                              //!< Return a JSON represenation of the taxon data
        virtual size_t                          getNumberOfCharacters(void) const = 0;                              //!< How many characters
        virtual double                          getPercentageMissing(void) const = 0;                               //!< Returns the percentage of missing data for this sequence
        const Taxon&                            getTaxon(void) const;                                               //!< Return the name of the character vector
        const std::string&                      getTaxonName(void) const;                                           //!< Return the name of the character vector
        virtual bool                            isCharacterResolved(size_t idx) const = 0;                          //!< Returns whether the character is fully resolved (e.g., "A" or "1.32") or not (e.g., "AC" or "?")
        virtual bool                            isSequenceMissing(void) const = 0;                                  //!< Returns whether the contains only missing data or has some actual observations
        virtual void                            removeCharacters(const std::set<size_t> &i) = 0;                    //!< Remove all the characters with a given index
        virtual void                            setAllCharactersMissing(void) = 0;                                  //!< Set all characters as missing
        void                                    setTaxon(const Taxon &tn);                                          //!< Set the taxon name
        virtual std::string                     getStringRepresentation(size_t idx) const = 0;
        virtual std::string                     getStateLabels(void) = 0;                                           //!< Get the possible state labels
    
    protected:
                                                AbstractTaxonData(const Taxon &t);                                  //!< Default constructor

        Taxon                                   taxon;                                                              //!< Name of the taxon for this vector of characters

    };

}

#endif
