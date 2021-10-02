#ifndef MatrixReader_H
#define MatrixReader_H

#include "DelimitedDataReader.h"
#include "MatrixReal.h"

namespace RevBayesCore {
	
	
	/**
	 * Reader for distance matrices from Phylip format files.
	 * The first line contains the number of tips.
	 * The following lines contain pairwise distances.
	 *
	 *
	 *
	 * @copyright Copyright 2009-
	 * @author The RevBayes Development Core Team (Bastien Boussau)
	 * @since 2015-03-03, version 1.0
	 *
	 */
	class MatrixReader : public DelimitedDataReader {
		
    public:
        MatrixReader(const std::string &fn, std::string d="", size_t ls=0);
        
        const MatrixReal&   							    getMatrix(void);

        
        
    protected:
        
        // protected member only accessible for derived classes
        MatrixReal                                          matrix;
		
	};
	
}

#endif
