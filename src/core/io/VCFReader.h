#ifndef VCFReader_H
#define VCFReader_H

#include "DelimitedDataReader.h"
#include "AbstractHomologousDiscreteCharacterData.h"
#include "HomologousDiscreteCharacterData.h"
#include "BinaryState.h"
#include "DnaState.h"

#include <string>
#include <vector>

namespace RevBayesCore {
    
    
    /**
     * Reader for VCF files.
     *
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2015-03-03, version 1.0
     *
     */
    class VCFReader : public DelimitedDataReader {

    public:
        
        enum PLOIDY { HAPLOID, DIPLOID, POLYPLOID };
        enum UNKOWN_TREATMENT { MISSING, REFERENCE, ALTERNATIVE };
        
        VCFReader(const std::string &fn, PLOIDY p, UNKOWN_TREATMENT u, bool read_data=true);
        
        HomologousDiscreteCharacterData<DnaState>*              readDNAMatrix( void );
        HomologousDiscreteCharacterData<BinaryState>*           readBinaryMatrix( void );
        void                                                    convertToCountsFile( const std::string& fn, const RbVector<Taxon>& taxa_list, const std::string& type );

    protected:
        
        std::vector<size_t>                                     extractStateIndices(std::string alleles, const std::string& type);

        
        std::string                                             filename;
        PLOIDY                                                  ploidy;
        UNKOWN_TREATMENT                                        unkown_treatment;

    };
    
}

#endif
