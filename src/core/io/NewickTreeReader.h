#ifndef NewickTreeReader_H
#define NewickTreeReader_H

#include <string>
#include <vector>
#include "RbFileManager.h"

namespace RevBayesCore {
    
    class Tree;
    
    /**
     * Newick tree reader.
     *
     * The newick tree reader provides convenience methods for reading trees in Newick format.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-01-29, version 1.0
     *
     */
    class NewickTreeReader {
        
    public:
        NewickTreeReader();                                                                          //!< Default constructor.
        
        std::vector<Tree*>*                 readBranchLengthTrees(const path &fn);                   //!< Read a set of trees with branch lengths in newick format from a file.
                
    };
    
}

#endif
