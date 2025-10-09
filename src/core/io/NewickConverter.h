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
#include <cassert>
namespace RevBayesCore {

    class Tree;
    class TopologyNode;

    // A ParseResult<T> is either
    // (i) Success -> a value and a new offset (non-empty optional).
    // (i) Failure -> no information (empty optional)
    // A parser is a function that takes a string and an offset and returns a ParseResult<T>
    //using Parser = ParseResult<T>(const std::string&,int);
    template <typename T>
    class ParseResult {
        bool success;      // true if parsing worked
        T value_;           // holds parsed thing if success = true
        int next_pos_;      // index in the string where parsing to continue if success = true
        std::string err_message_; // informative error if success = false
        int err_pos_;       // where it failed
    public:
        // define constructors
        ParseResult(bool b, const T& t, int i):success(true), value_(t), next_pos_(i){

        }
        ParseResult(const std::string& s, int i):success(false), err_message_(s), err_pos_(i){

        }
        ~ParseResult(){

        }
        operator bool() const{
            return success;
        }

        T& value(){
            assert(success);
            return value_;
        }
        const T& value() const{
            assert(success);
            return value_;
        }
        int& next_pos(){
            assert(success);
            return next_pos_;
        }
        int next_pos() const{
            assert(success);
            return next_pos_;
        }

        std::string& err_message(){
            assert(!success);
            return err_message_;
        }
        const std::string& err_message() const{
            assert(!success);
            return err_message_;
        }
        int& err_pos(){
            assert(!success);
            return err_pos_;
        }
        int err_pos() const{
            assert(!success);
            return err_pos_;
        }
    };

    template <typename T>
    ParseResult<T> ParseSuccess(const T& t, int i){    
        return ParseResult<T>(true, t, i);
    }
    template <typename T>
    ParseResult<T> ParseFail(const std::string& s, int i){    
        return ParseResult<T>(s, i);
    }
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
