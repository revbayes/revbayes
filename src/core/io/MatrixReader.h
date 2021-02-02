#ifndef MatrixReader_H
#define MatrixReader_H

#include <stddef.h>
#include <iosfwd>

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
	class MatrixReader {
		
    public:
        MatrixReader(const std::string &fn, char d='\t', size_t ls=0);
        
        void                                                readData( size_t ls);
        const std::string&                                  getFilename(void);
        const MatrixReal&   							    getMatrix(void);

        
        
    protected:
        
        // protected member only accessible for derived classes
        std::string                                         filename;
        char                                                delimiter;
		MatrixReal                                          matrix;
		
	};
	
}

#endif
