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
#include <optional>
#include <cassert>
#include <string>
#include <variant>

namespace RevBayesCore {

    class BranchHistory;
    class CharacterHistoryDiscrete;
    class Tree;
    class TopologyNode;

    template <typename T>
    struct ParseSuccess
    {
        T value;       // holds parsed thing
        int next_pos;  // index in the string where parsing to continue if success = true
        ParseSuccess(const T& t, int i):value(t), next_pos(i) {}
    };

    struct ParseFail
    {
        int err_pos;            // where it failed
        bool soft;              // soft means we didn't find something, hard means we found it but its malformed.
        ParseFail(int i, bool s=true): err_pos(i), soft(s) {}
    };


    // A ParseResult<T> is either
    // (i) Success -> a value and a new offset (non-empty optional).
    // (i) Failure -> no information (empty optional)
    // A parser is a function that takes a string and an offset and returns a ParseResult<T>
    //using Parser = ParseResult<T>(const std::string&,int);
    template <typename T>
    class ParseResult: public std::variant<ParseSuccess<T>, ParseFail>
    {
    public:

        using std::variant<ParseSuccess<T>, ParseFail>::variant;

        bool success() const
        {
            return std::holds_alternative<ParseSuccess<T>>(*this);
        }

        const ParseFail& as_failure() const
        {
            assert(not success());
            return std::get<ParseFail>(*this);
        }

        bool failure() const
        {
            return not success();
        }
        
        bool hard_failure() const
        {
            if (success()) return false;
            return not as_failure().soft;
        }
        
        bool soft_failure() const
        {
            if (success()) return false;
            return as_failure().soft;
        }
        
        ParseFail as_hard_failure() const
        {
            assert(not success());
            auto f = as_failure();
            f.soft = false;
            return f;
        }

        operator bool() const
        {
            return success();
        }

        T& value()
        {
            assert(success());
            return std::get<ParseSuccess<T>>(*this).value;
        }

        const T& value() const
        {
            assert(success());
            return std::get<ParseSuccess<T>>(*this).value;
        }

        int& next_pos()
        {
            assert(success());
            return std::get<ParseSuccess<T>>(*this).next_pos;
        }

        int next_pos() const
        {
            assert(success());
            return std::get<ParseSuccess<T>>(*this).next_pos;
        }

        int& err_pos()
        {
            assert(not success());
            return std::get<ParseFail>(*this).err_pos;
        }

        int err_pos() const
        {
            assert(not success());
            return std::get<ParseFail>(*this).err_pos;
        }
    };

    // Why keep this class at all?
    // Hypothetically, in the future we could us it to store options to the Newick conversion process.
    class NewickConverter
    {

    public:
        NewickConverter();
        virtual                    ~NewickConverter();
        
        Tree*                       convertFromNewick(const std::string &n);
        CharacterHistoryDiscrete*   convertSimmapFromNewick(const std::string &n );


    private:
        TopologyNode*           createSimmapNode(const std::string &n, std::vector<TopologyNode*> &nodes, std::vector<double> &brlens, std::vector<BranchHistory*> &histories);
    };

    void postProcessNewick(Tree* tree);
    std::vector<Tree*> readNewicks(const std::string& input);

    ParseResult<std::vector<Tree*>> parseTrees(const std::string& input, int start_pos);
    ParseResult<Tree*> parseTree(const std::string& input, int start_pos);
    ParseResult<TopologyNode*> parseSubTree(const std::string& input, int start_pos);
    ParseResult<char> checkChar(const std::string& input, int start_pos, char c);
    ParseResult<char> parseChar(const std::string& input, int start_pos);
    ParseResult<char> parseQuotedChar(const std::string& input, int start_pos);
    ParseResult<std::pair<std::optional<double>,std::vector<std::string>>> parseBranch(const std::string& input, int start_pos);
    ParseResult<std::string> parseQuotedName(const std::string& input, int start_pos);
    ParseResult<std::string> parseUnquotedName(const std::string& input, int start_pos);
    ParseResult<std::vector<TopologyNode*>> parseDescendants(const std::string& input, int start_pos);
    ParseResult<std::string> parseName(const std::string& input, int start_pos);
    ParseResult<char> parseUnquotedChar(const std::string& input, int start_pos);
}

#endif
