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
    class VCFReader : public DelimitedDataReader, public Cloneable {
        
    public:
        
        enum PLOIDY { HAPLOID, DIPLOID, POLYPLOID };
        enum UNKOWN_TREATMENT { MISSING, REFERENCE, ALTERNATIVE };
        
        VCFReader(const std::string &fn, PLOIDY p, UNKOWN_TREATMENT u, bool read_data=true);
        
        VCFReader*                                              clone(void) const;
        HomologousDiscreteCharacterData<DnaState>*              readDNAMatrix( void );
        HomologousDiscreteCharacterData<BinaryState>*           readBinaryMatrix( void );
        void                                                    convertToCountsFile( const std::string& fn, const RbVector<Taxon>& taxa_list, const std::string& type, long thinning, long skip );
        RbVector<long>                                          convertToSFS( const RbVector<Taxon>& taxa_list );
        void                                                    computeMonomorphicVariableStatistics( const std::string& fn, const RbVector<Taxon>& taxa_list);

    protected:
        void                                                    mapSpeciesNames(const RbVector<Taxon>& taxa_list, std::vector<std::string>& species_names, std::map<std::string, size_t>& species_names_to_index, std::vector< std::vector<size_t> >& indices_of_taxa_per_species );
        std::vector<size_t>                                     extractStateIndices(std::string alleles, const std::string& type);
        
        std::string                                             filename;
        PLOIDY                                                  ploidy;
        UNKOWN_TREATMENT                                        unkown_treatment;
        
        bool                                                    statistics_ready;
        std::vector< std::vector<size_t> >                      mono_in_A_var_in_B;
        std::vector< std::vector<size_t> >                      mono_in_both_equal;
        std::vector< std::vector<size_t> >                      mono_in_both_diff;
        std::vector< std::vector<size_t> >                      var_in_both;

    };
    
}

#endif
