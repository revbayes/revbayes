/**
 * @file
 * This file contains the declaration of NewickConvert, which is our simple class to
 * convert trees written in Newick format into our own tree format.
 *
 * @brief Declaration of NclReader
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date: 2012-07-12 11:08:07 +0200 (Thu, 12 Jul 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2009-08-16, version 1.0
 * @extends DAGNode
 *
 * $Id: NclReader.h 1674 2012-07-12 09:08:07Z hoehna $
 */

#ifndef NewickConverter_H
#define NewickConverter_H

#include "RbFileManager.h"

#include <vector>
#include <iosfwd>

namespace RevBayesCore {

    class BranchHistory;
    class CharacterHistoryDiscrete;
    class Tree;
    class TopologyNode;

    class NewickConverter {

    public:
        NewickConverter();
        virtual                    ~NewickConverter();
    
        Tree*                       convertFromNewick(const std::string &n);
        CharacterHistoryDiscrete*   convertSimmapFromNewick(const std::string &n );


    private:
        TopologyNode*           createNode(const std::string &n, std::vector<TopologyNode*> &nodes, std::vector<double> &brlens);
        TopologyNode*           createSimmapNode(const std::string &n, std::vector<TopologyNode*> &nodes, std::vector<double> &brlens, std::vector<BranchHistory*> &histories);
    };

}

#endif
