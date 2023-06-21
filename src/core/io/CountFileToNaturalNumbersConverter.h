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

        void                    cfconverter( const path &fi, const size_t n_individuals, const path &fo );

        //size_t                  getNumberOfAlleles( void );
        //size_t                  getNumberOfIndividuals( void );

        size_t                  getState(const std::string& counts, size_t n_alleles, size_t n_individuals, const std::vector<std::string>& vector_edges, const std::vector<int>& matrix_edges );
        size_t                  sampleWeight(size_t M, size_t m, size_t N);
        size_t                  sampleEdge(size_t allele_index, const std::vector<int>& vector);
        size_t                  getIndex(const std::vector<std::string>& vector, const std::string& element);

    protected:


        
    };
    
}
#endif
