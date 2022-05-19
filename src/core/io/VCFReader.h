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
        
        VCFReader(const std::string &fn, PLOIDY p, UNKOWN_TREATMENT u);
        
        HomologousDiscreteCharacterData<DnaState>*              readDNAMatrix( void );
        HomologousDiscreteCharacterData<BinaryState>*           readBinaryMatrix( void );

    protected:
        
        std::string                                             filename;
        PLOIDY                                                  ploidy;
        UNKOWN_TREATMENT                                        unkown_treatment;

    };
    
}

#endif
