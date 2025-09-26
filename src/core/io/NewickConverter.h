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

    // A ParseResult<T> is either
    // (i) Success -> a value and a new offset (non-empty optional).
    // (i) Failure -> no information (empty optional)
    template <typename T>
    using ParseResult = std::optional<std::pair<T,int>>;

    // A parser is a function that takes a string and an offset and returns a ParseResult<T>
    template <typename T>
    using Parser = ParseResult<T>(const std::string&,int);

    // Why keep this class at all?
    // Hypothetically, in the future we could us it to store options to the Newick conversion process.
    class NewickConverter
    {

    public:
        NewickConverter();
        virtual                 ~NewickConverter();
        
        
        Tree*                   convertFromNewick(const std::string &n);
    };

    ParseResult<TopologyNode*> parseTree(const std::string& input, int start_pos);
    ParseResult<TopologyNode*> parseSubTree(const std::string& input, int start_pos);
    ParseResult<char> checkChar(const std::string& input, int start_pos, char c);
    ParseResult<char> parseChar(const std::string& input, int start_pos);
    ParseResult<char> parseQuotedChar(const std::string& input, int start_pos);
    ParseResult<std::string> parseQuotedName(const std::string& input, int start_pos);
    ParseResult<std::string> parseUnquotedName(const std::string& input, int start_pos);
    ParseResult<TopologyNode*> parseInternal(const std::string& input, int start_pos);
    ParseResult<TopologyNode*> parseLeaf(const std::string& input, int start_pos);
    ParseResult<TopologyNode*> parseBranch(const std::string& input, int start_pos);
    ParseResult<std::vector<TopologyNode*>> parseBranchSet(const std::string& input, int start_pos);
    ParseResult<std::string> parseName(const std::string& input, int start_pos);
    ParseResult<char> parseUnquotedChar(const std::string& input, int start_pos);
}

#endif
