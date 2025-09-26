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


#include <vector>
#include <iosfwd>
#include <optional>

namespace RevBayesCore {

    class Tree;
    class TopologyNode;

    // Why keep this class at all?
    // Hypothetically, in the future we could us it to store options to the Newick conversion process.
    class NewickConverter
    {

    public:
        NewickConverter();
        virtual                 ~NewickConverter();
        
        
        Tree*                   convertFromNewick(const std::string &n);
    };

    std::optional<std::pair<TopologyNode*, int>> parseTree(const std::string& input, int start_pos);
    std::optional<std::pair<TopologyNode*, int>> parseSubTree(const std::string& input, int start_pos);
    std::optional<int> checkChar(const std::string& input, int start_pos, char c);
    std::optional<std::pair<char, int>> parseChar(const std::string& input, int start_pos);
    std::optional<std::pair<char, int>> parseQuotedChar(const std::string& input, int start_pos);
    std::optional<std::pair<std::string, int>> parseQuotedName(const std::string& input, int start_pos);
    std::optional<std::pair<std::string, int>> parseUnquotedName(const std::string& input, int start_pos);
    std::optional<std::pair<TopologyNode*, int>> parseInternal(const std::string& input, int start_pos);
    std::optional<std::pair<TopologyNode*, int>> parseLeaf(const std::string& input, int start_pos);
    std::optional<std::pair<TopologyNode*, int>> parseBranch(const std::string& input, int start_pos);
    std::optional<std::pair<std::vector<TopologyNode*>, int>> parseBranchSet(const std::string& input, int start_pos);
    std::optional<std::pair<std::string, int>> parseName(const std::string& input, int start_pos);
    std::optional<std::pair<char, int>> parseUnquotedChar(const std::string& input, int start_pos);
}

#endif
