#include "NewickConverter.h"

#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <vector>
#include <string>

#include "RbException.h"
#include "TopologyNode.h"
#include "Tree.h"

using std::optional;
using std::pair;
using std::string;
using std::vector;

namespace RevBayesCore
{

NewickConverter::NewickConverter()
{

}


NewickConverter::~NewickConverter()
{

}

Tree* NewickConverter::convertFromNewick(std::string const &newick)
{

    // create and allocate the tree object
    Tree *t = new Tree();

    // construct the tree starting from the root
    auto result = parseTree(newick, 0);
    if (not result){
        throw RbException()<<result.err_message()<< " at location " << result.err_pos();
    }
    if (result.next_pos() != newick.size())
        throw RbException()<<"Junk at end of newick string: '"<<newick.substr(result.next_pos())<<"'";

    // set up the tree
    auto root = result.value();

    t->setRoot( root, true );

    // try to set node indices from attributes
    t->tryReadIndicesFromParameters(true);

    // trees with 2-degree root nodes should not be rerooted
    t->setRooted( root->getNumberOfChildren() == 2 );
    
    // make this tree first a branch length tree
    // hence, tell the root to use branch lengths and not ages (with recursive call)
    root->setUseAges(false, true);
    
    // return the tree, the caller is responsible for destruction
    return t;
}

/* We are currently using this Newick grammar:
  
    Tree → Subtree ";"
    Subtree → Leaf | Internal
    Leaf → Name
    Internal → "(" BranchSet ")" Name
    BranchSet → Branch | Branch "," BranchSet
    Branch → Subtree Length
    Name -> QuotedName | UnquotedName | empty
    Length → ":" number | empty

But we should probably switch to this Newick grammar:

    Tree → Subtree ";"
    Subtree → [Descendants] [Name] [Length]
    Descendants = "(" Subtree ("," Subtree)* ")"
    Name -> QuotedName | UnquotedName
    Length -> ":" [Double]
    QuotedName -> "'" QuotedChar* "'"
    UnquotedName -> UnquotedChar*
    UnquotedChar = not punctuation

When we put the spacing in, we should get:

    Tree → SPACE<tree> Subtree ";" SPACE<ignore>
    Subtree → SPACE<node> [Descendants] SPACE<node> [Name] SPACE<node> (Length | empty<no length, no attributes>)
    Descendants = "(" Subtree % "," ")"
    Length -> ":" SPACE<branch> [Double] SPACE<branch>

where SPACE<x> means that we read whitepspace and newick comments, and the comments are attached to x.
So for example SPACE<tree> means that any comments found are considered to refer to the tree.
*/

// succeeds if there is another character in the input, otherwise fails
ParseResult<char> parseChar(const std::string& input, int start_pos)
{
    // if reading beyond end of string return null
    if (start_pos >= input.size())
        return ParseFail("End of input",start_pos);
    else
    {
        char c = input[start_pos];
        return ParseSuccess(c,start_pos+1);
    }
}

string escape(char c)
{
    if (c == '\'')
        return "\\'";
    else if (c == '\n')
        return "\\n";
    else if (c == '\r')
        return "\\r";
    else if (c == '\t')
        return "\\t";
    else
        return string(1,c);
}

// returns a new start_pos if the input contains another character, and it is "c"
ParseResult<char> checkChar(const std::string& input, int start_pos, char c)
{
    // if reading beyond end of string return null
    if (auto maybe_char = parseChar(input, start_pos))
    {
        if (maybe_char.value() == c)
            return maybe_char;
        else
            return ParseFail(std::string("Expected '") + escape(c) + "'", start_pos);
    }
    else
        return ParseFail(std::string("Reached end of input looking for '") + escape(c) + "'", start_pos);
}


ParseResult<string> parseNewickComment(const std::string& input, int start_pos)
{
    // 1. First check for an open bracket
    auto check_open_bracket = checkChar(input, start_pos, '[');
    if (not check_open_bracket)
        return check_open_bracket.as_failure();

    // FIXME: only some comments are supposed to be computer-readable.
    // Comments that start with [& should be saved.
    // Also some others... [%U, [%R, and [&&NHX:
    // Currently we don't try to decide here what to save and what to discard, and just pass back all comments.

    // 2. Then read characters until an end bracket
    std::string newick_comment;
    while(auto check_char = parseChar(input, start_pos))
    {
        start_pos = check_char.next_pos();

        if (check_char.value() != ']')
            newick_comment += check_char.value();
        else
            return ParseSuccess(newick_comment, start_pos);
    }

    // 3. If we get here, then ending bracket is missing!
    return ParseFail("Newick comment is missing final ']'", start_pos);
}

// Which values should count as whitespace?

// whiteChar -> [ \t\n\r]    
ParseResult<char> parseWhiteChar(const std::string& input, int start_pos)
{
    // 1. Try to read a character
    if (auto check_char = parseChar(input, start_pos))
    {
        // 2. If it is whitespace, return the character
        if (strchr(" \t\n\r", check_char.value()))
            return check_char;
        else
            return ParseFail("Not whitespace", start_pos);
    }
    else
        return check_char;
}

// OneWhitespace -> NewickComment OR WhiteChar    
ParseResult<optional<string>> parseOneWhitespace(const std::string& input, int start_pos)
{
    // 1. First try to read a newick comment
    if (auto check_comment = parseNewickComment(input, start_pos))
    {
        return ParseSuccess<optional<string>>({check_comment.value()}, check_comment.next_pos());
    }
    // 2. Then try to read a whitespace character
    else if (auto check_char = parseWhiteChar(input, start_pos))
    {
        // If we found a whitespace char, return success with no comment
        return ParseSuccess<optional<string>>({}, check_comment.next_pos());
    }
    // 3. If no comment and no whitespace character, return null
    else
        return ParseFail("Expected white space ", start_pos);
}

// Whitespace -> OneWhitespace*
ParseResult<vector<string>> parseWhitespace(const std::string& input, int start_pos)
{
    vector<string> newick_comments;
    while(auto check = parseOneWhitespace(input, start_pos))
    {
        start_pos = check.next_pos();
        if (check.value())
            newick_comments.push_back(*check.value());
    }
    return ParseSuccess(newick_comments, start_pos);
}


// Tree -> Subtree ";"
ParseResult<TopologyNode*> parseTree(const std::string& input, int start_pos)
{
    // 1. Get the Subtree
    auto check_subtree = parseSubTree(input,start_pos);

    // Return error message if we failed to find a subtree.
    if (not check_subtree) return check_subtree;  

    // 2. Check the semicolon
    auto check_semi = checkChar(input, check_subtree.next_pos(), ';');

    // Return error message if we failed to find a semicolon.
    if (not check_semi) return check_semi.as_failure();

    return ParseSuccess(check_subtree.value(), check_semi.next_pos());
}

// subtree -> [Descendants] [Name]
ParseResult<TopologyNode*> parseSubTree(const std::string& input, int start_pos)
{
    auto node = new TopologyNode;

    // Parse Descendants and add children to node
    if (auto maybe_children = parseDescendants(input, start_pos)) {
        start_pos = maybe_children.next_pos();
        for(auto& child: maybe_children.value())
        {
            node->addChild(child);
            child->setParent(node);
        }
    }
    else if (maybe_children.hard_failure())
        return maybe_children.as_failure();
    
    // Read Name
    if (auto maybe_name = parseName(input, start_pos)) {
        start_pos = maybe_name.next_pos();
        node -> setName(maybe_name.value());
    }
    else if (maybe_name.hard_failure())
        return maybe_name.as_failure();

    //add name and children
    return ParseSuccess(node, start_pos); 
}

ParseResult<double> parseDouble(const std::string& input, int start_pos)
{
    char* endptr;
    double value = std::strtod(input.c_str() + start_pos, &endptr);

    errno = 0;
    start_pos = endptr - input.c_str();
    if (errno == 0)
        return ParseSuccess(value, start_pos);
    else
        return ParseFail("Failure reading floating point number", start_pos);
}

// Length -> ":" number | empty
ParseResult<optional<double>> parseLength(const std::string& input, int start_pos)
{
    auto maybe_colon = checkChar(input, start_pos, ':');
    if (not maybe_colon)
        return ParseSuccess<optional<double>>({}, start_pos);

    auto maybe_length = parseDouble(input, maybe_colon.next_pos());
    if (not maybe_length)
        return ParseSuccess<optional<double>>({}, maybe_colon.next_pos());

    return ParseSuccess<optional<double>>(maybe_length.value(), maybe_length.next_pos());
}


// Branch -> Subtree Length
ParseResult<TopologyNode*> parseBranch(const std::string& input, int start_pos)
{
    auto maybe_subtree = parseSubTree(input, start_pos);
    if (not maybe_subtree) return maybe_subtree.as_failure();

    auto maybe_length = parseLength(input, maybe_subtree.next_pos());
    if (not maybe_length) return maybe_length.as_failure();

    auto node = maybe_subtree.value();
    auto length = maybe_length.value();
    if (length)
        node->setBranchLength(*length);

    return ParseSuccess(node, maybe_length.next_pos());
}

// Descendants -> "(" > Branch ("," > Branch)* > ")"
ParseResult<std::vector<TopologyNode*>> parseDescendants(const std::string& input, int start_pos)
{
    // Read left parenthesis
    if (auto maybe_lparen = checkChar(input, start_pos, '('))
        start_pos = maybe_lparen.next_pos();
    else
        return maybe_lparen.as_failure();

    vector<TopologyNode*> branches;
    auto maybe_branch = parseBranch(input, start_pos);
    if (not maybe_branch)
        return maybe_branch.as_hard_failure();
    branches.push_back(maybe_branch.value());
    start_pos = maybe_branch.next_pos();

    while(auto maybe_comma = checkChar(input, start_pos, ','))
    {
        start_pos = maybe_comma.next_pos();
        auto maybe_branch2 = parseBranch(input, start_pos);
        if (not maybe_branch2)
            return maybe_branch2.as_hard_failure();
        branches.push_back(maybe_branch2.value());
        start_pos = maybe_branch2.next_pos();
    }


    // Read right parenthesis
    if (auto maybe_rparen = checkChar(input, start_pos, ')')) {
        start_pos = maybe_rparen.next_pos();
    }
    else
        return maybe_rparen.as_hard_failure();

    return ParseSuccess(branches, start_pos);
}

// this function is matching (not a quote) or (two quotes)
// QuotedChar -> [^'] | ''
ParseResult<char> parseQuotedChar(const std::string& input, int start_pos){
    assert(start_pos>=0);

    if (auto maybe_char1 = parseChar(input, start_pos)){
        if (maybe_char1.value() == '\'') {
            // two quotes returns quote
            if (auto maybe_char2 = checkChar(input, maybe_char1.next_pos(), '\''))
                // new_start_pos2 will ALWAYS be start_pos + 2
                return ParseSuccess<char>('\'', maybe_char2.next_pos());
            // one quote fails
            else
                return ParseFail("Single quote not followed by another single quote", start_pos);
        }
        else
            //new_start_pos will ALWAYS be start_pos+1
            return ParseSuccess<char>(maybe_char1.value(), maybe_char1.next_pos());
    }
    else
        return maybe_char1;
}

//no quotes 
ParseResult<char> parseUnquotedChar(const std::string& input, int start_pos)
{
    assert(start_pos>=0);

    if (auto maybe_char = parseChar(input, start_pos))
    {
        char c = maybe_char.value();
        // c is _
        if (c == '_')
            return ParseSuccess(' ', maybe_char.next_pos());
        // c is a legal character that is not _
        else if (!strchr("()[]':;, ", c)) 
            return ParseSuccess(c, maybe_char.next_pos());
        // c is an illegal character (punctuation)
        else 
            return ParseFail("Illegal character", start_pos);
    }
    else
        return maybe_char;
}

// QuotedName -> ' + QuotedChar* + '
// quotedName -> char('\'') > many(quotedChar) > char('\'')
ParseResult<std::string> parseQuotedName(const std::string& input, int start_pos)
{
    auto maybe_quote1 = checkChar(input, start_pos, '\'');
    if (not maybe_quote1)
        return maybe_quote1.as_failure();
    else
        start_pos = maybe_quote1.next_pos();

    std::string name; 
    while (auto check = parseQuotedChar(input, start_pos))
    {
        name += check.value();
        start_pos = check.next_pos();
    }

    auto maybe_quote2 = checkChar(input, start_pos, '\'');
    if (not maybe_quote2)
        return maybe_quote2.as_hard_failure();
    else
        start_pos = maybe_quote2.next_pos();

    return ParseSuccess(name, start_pos);
}

// UnquotedName -> UnquotedChar+
// unquotedName -> some(unquotedChar)
ParseResult<std::string> parseUnquotedName(const std::string& input, int start_pos){
    std::string name;
    if (auto check = parseUnquotedChar(input, start_pos))
    {
        name += check.value();
        start_pos = check.next_pos();
    }
    //if not return null
    else
        return check.as_failure();

    //continue check
    while (auto check = parseUnquotedChar(input, start_pos)){
        name += check.value();
        start_pos = check.next_pos();
    }

    return ParseSuccess(name, start_pos);
}

// Name -> QuotedName | UnquotedName | empty
ParseResult<std::string> parseName(const std::string& input, int start_pos){
    if (auto maybe_quoted = parseQuotedName(input, start_pos); maybe_quoted.success() or maybe_quoted.hard_failure())
        return maybe_quoted;
    else 
        return parseUnquotedName(input, start_pos);
}

//AdmixtureTree* NewickConverter::getAdmixtureTreeFromNewick(std::string const &n)
//{
//
//    // create and allocate the tree object
//    AdmixtureTree *t = new AdmixtureTree();
//
//    std::vector<TopologyNode*> nodes;
//    std::vector<double> brlens;
//
//    // construct the tree starting from the root
//    //TopologyNode *root = createNode( n, nodes, brlens );
//
//    // convert to AdmixtureNode*
//    std::vector<AdmixtureNode*> adm_nodes;
//    for (size_t i = 0; i < nodes.size(); i++)
//        adm_nodes.push_back(static_cast<AdmixtureNode*>(nodes[i]));
//
//    // set up the tree
//    t->setRoot(adm_nodes[adm_nodes.size()-1]);
//
//    // set the branch lengths
//    //for (size_t i = 0; i < nodes.size(); ++i) {
//    t->setAgesFromBrlens(brlens);
//    ;//t->setBranchLength(nodes[i]->getIndex(), brlens[i]);
//    // }
//
//    // return the tree, the caller is responsible for destruction
//    return t;
//}
}
