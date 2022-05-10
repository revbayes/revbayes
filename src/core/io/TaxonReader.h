#ifndef TaxonReader_H
#define TaxonReader_H

#include "DelimitedDataReader.h"
#include "Taxon.h"

namespace RevBayesCore {
    
    
    
    /**
     * Reader for taxon names mapped to other information including species name and sampling dates.
     *
     * This reader is a simple file reader of a delimited file, e.g., by tab-stops.
     * In the first column should be the taxon name and in the second column the date of sampling.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-03-18, version 1.0
     *
     */
    class TaxonReader : public DelimitedDataReader {
        
    public:
        
        TaxonReader(const std::string &fn, const std::string& d="");                        //!< Constructor
        
        const std::vector<Taxon>&   getTaxa(void) const;                        //!< Get the taxa.

        
    protected:
        
        std::vector<Taxon>          taxa;
    };
}

#endif
