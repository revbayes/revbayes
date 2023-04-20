#ifndef CountFileToNaturalNumbersConverter_H
#define CountFileToNaturalNumbersConverter_H

#include <vector>
#include <string>
#include "RbFileManager.h"

namespace RevBayesCore {
    
    
    /**
     * Converts a count file to a natural numbers type file
     * the natural numbers are PoMo states
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Rui Borges)
     * @since 2021-07-17, version 1.0
     *
     */
    class CountFileToNaturalNumbersConverter {
        
    public:
        CountFileToNaturalNumbersConverter();

        void               cfconverter( const path &fi, const size_t n_individuals, const path &fo );

        //const size_t     getNumberOfAlleles( void );
        //const size_t     getNumberOfIndividuals( void );

        const size_t          getState(std::string counts, size_t n_alleles, size_t n_individuals, std::vector<std::string>& vector_edges, std::vector<int>& matrix_edges );
        const size_t          sample_weight(size_t M, size_t m, size_t N);
        const size_t          sample_edge(size_t allele_index, std::vector<int>& vector);
        const size_t          get_index(std::vector<std::string>& vector, std::string element);

    protected:



        
    };
    
}
#endif
