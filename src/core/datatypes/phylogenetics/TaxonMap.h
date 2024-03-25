#ifndef TaxonMap_H
#define TaxonMap_H

#include <cstddef>
#include <vector>
#include <map>
#include <iosfwd>

#include "Taxon.h"

namespace RevBayesCore {
    
    class Tree;
    
    
    /**
     * @brief TaxonMap object representing an individual.
     *
     * This class represents a TaxonMap.
     * A TaxonMap simply consists of a TaxonMap name and a date when this TaxonMap has been observed.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-03-18, version 1.0
     */
    class TaxonMap {
        
    public:
        
        TaxonMap(void);                                //!< Default constructor
        TaxonMap(const Tree &t);                                                            //!< Constructor from tree object
        virtual                             ~TaxonMap() {}
        
//        bool                                operator==(const TaxonMap &t) const;           //!< Equals operators
//        bool                                operator!=(const TaxonMap &t) const;           //!< Not-quals operators
//        bool                                operator<(const TaxonMap &t) const;            //!< Less-than operators
//        bool                                operator<=(const TaxonMap &t) const;           //!< Less-than operators
//        bool                                operator>(const TaxonMap &t) const;            //!< Less-than operators
//        bool                                operator>=(const TaxonMap &t) const;           //!< Less-than operators
        
        // public methods
        void                                addTaxon(const Taxon &t);                       //!< Get the age for this TaxonMap.
        bool                                hasTaxon(const Taxon &t) const;                 //!< Get the i-th Taxon
        size_t                              getTaxonIndex(const Taxon &t) const;            //!< Get the i-th Taxon

        const Taxon&                        getTaxon(size_t i) const;                       //!< Get the i-th Taxon

        int                                 size() const;                                   //!< How many taxa in the map?
        
        std::string                         print() const;

    private:
        
        std::vector<Taxon>                  taxa;
        std::map<Taxon, size_t>             taxa_map;
        
    };
    
    // Global functions using the class
    std::ostream&                       operator<<(std::ostream& o, const TaxonMap& x);                                         //!< Overloaded output operator
}

#endif
