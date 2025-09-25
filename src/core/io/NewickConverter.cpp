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

// A parser is a function that takes a string and an offset, and returns either
// (i) failure (empty optional)
// (ii) a value and a new offset (non-empty optional).
template <typename T>
using Parser = std::optional<std::pair<T,int>>(const std::string&,int);

Tree* NewickConverter::convertFromNewick(std::string const &n)
{

    // create and allocate the tree object
    Tree *t = new Tree();

    std::vector<TopologyNode*> nodes;
    std::vector<double> brlens;

    // create a string-stream and throw the string into it
    std::stringstream ss (std::stringstream::in | std::stringstream::out);
    ss << n;

    // ignore white spaces
    std::string trimmed = "";
    char c;
    while ( ss.good() )
    {
        // check for EOF
        int c_int = ss.get();
        if (c_int != EOF)
        {
            c = char( c_int );
            if ( c != ' ')
            {
                trimmed += c;
            }
        }
    }

    // construct the tree starting from the root
    TopologyNode *root = createNode( trimmed, nodes, brlens );

    // set up the tree
    t->setRoot( root, true );

    // try to set node indices from attributes
    t->tryReadIndicesFromParameters(true);

    // set the branch lengths
    for (size_t i = 0; i < nodes.size(); ++i)
    {
        t->getNode( nodes[i]->getIndex() ).setBranchLength( brlens[i] );
    }

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

Comments:

    - We could simplify by doing `Subtree -> [ "(" BranchSet ")" ] Name`, where [..] means optional.
    - Newick comments are considered whitespace.  They can go anywhere whitespace can go.
    - If we explicitly consider whitespace, then

          Leaf -> Name

      looks like

          Leaf -> Whitespace Name Whitespace
      
      This is a bit uglier to read.
          
    - ???

*/

// succeeds if there is another character in the input, otherwise fails
std::optional<std::pair<char, int>> parseChar(const std::string& input, int start_pos){
    // if reading beyond end of string return null
    if (start_pos >= input.size()){
        return {};
    }
    char c = input[start_pos];
    return optional<pair<char, int>>(pair<char,int>(c,start_pos+1));
}

// returns a new start_pos if the input contains another character, and it is "c"
std::optional<int> checkChar(const std::string& input, int start_pos, char c){
    // if reading beyond end of string return null
    if (auto check=parseChar(input, start_pos)){
        auto [c2, new_start_pos] = *check;
        if (c2 == c){
            return new_start_pos;
        }
        else{
            return {};
        }
    }
    else{
        return {};
    }
}


optional<pair<string,int>> parseNewickComment(const std::string& input, int start_pos)
{
    // 1. First check for an open bracket
    auto check_open_bracket = checkChar(input, start_pos, '[');
    if (not check_open_bracket)
        return {};

    // FIXME: only some comments are supposed to be computer-readable.
    // Comments that start with [& should be saved.
    // Also some others... [%U, [%R, and [&&NHX:
    // Currently we don't try to decide here what to save and what to discard, and just pass back all comments.

    // 2. Then read characters until an end bracket
    std::string newick_comment;
    while(auto check_char = parseChar(input, start_pos))
    {
        auto& [c, new_start_pos] = *check_char;
        start_pos = new_start_pos;
        if (c != ']')
            newick_comment += 'c';
        else
            return {{newick_comment, start_pos}};
    }

    // 3. If we get here, then ending bracket is missing!
    return {};
}

optional<pair<char,int>> parseWhiteChar(const std::string& input, int start_pos)
{
    // 1. Try to read a character
    auto check_char = parseChar(input, start_pos);
    if (not check_char) return {};
    auto& [c, new_start_pos] = *check_char;

    // 2. If it is whitespace, return the character
    if (c == ' ' or c == '\t' or c == '\n' or c == '\r')
        return {{c, new_start_pos}};
    // 3. If it is not whitespace, return null
    else
        return {};

    // QUESTION: Which characters should count as whitespace?
}
    
// OneWhitespace -> NewickComment OR WhiteChar    
optional<pair<optional<string>,int>> parseOneWhitespace(const std::string& input, int start_pos)
{
    // 1. First try to read a newick comment
    if (auto check_comment = parseNewickComment(input, start_pos))
    {
        auto& [comment, new_start_pos] = *check_comment;
        return {{{comment}, new_start_pos}};
    }
    // 2. Then try to read a whitespace character
    else if (auto check_char = parseWhiteChar(input, start_pos))
    {
        auto& [c, new_start_pos] = *check_char;
        // If we found a whitespace char, return success with no comment
        return {{{},new_start_pos}};
    }
    // 3. If no comment and no whitespace character, return null
    else
        return {};
}

// Whitespace -> OneWhitespace*
optional<pair<vector<string>,int>> parseWhitespace(const std::string& input, int start_pos)
{
    vector<string> newick_comments;
    while(auto check = parseOneWhitespace(input, start_pos))
    {
        auto& [maybe_comment,new_start_pos] = *check;
        start_pos = new_start_pos;
        if (maybe_comment)
            newick_comments.push_back(*maybe_comment);
    }
    return {{newick_comments, start_pos}};
}


// Tree -> Subtree ";"
std::optional<std::pair<TopologyNode*, int>> parseTree(const std::string& input, int start_pos)
{
    // 1. Get the Subtree
    auto check_subtree = parseSubTree(input,start_pos);
    if (not check_subtree)
        return {};
    auto& [subtree, new_start_pos] = *check_subtree;

    // 2. Check the semicolon
    auto check_semi = checkChar(input, new_start_pos, ';');
    if (not check_semi)
        return {};

    int new_start_pos2 = *check_semi;

    return {{subtree, new_start_pos2}};
}

// subtree -> internal OR leaf
std::optional<std::pair<TopologyNode*, int>> parseSubTree(const std::string& input, int start_pos){

    if (auto check = parseInternal(input, start_pos))
        return check;
    
    else {
        return parseLeaf(input, start_pos);
    }
}

// Internal -> '(' BranchSet ')' Name
std::optional<std::pair<TopologyNode*, int>> parseInternal(const std::string& input, int start_pos){
    auto node = new TopologyNode;
    // Check if we have left parenthesis
    if (auto check = checkChar(input, start_pos, '('))
        start_pos = check.value();

    else
        return {};
    //add children to node: TODO
    if (auto check = parseBranchSet(input, start_pos)){
        auto [children, new_start_pos] = check.value();
        start_pos = new_start_pos;
    }

    else
        return {};

    // Check for right parenthesis
    if (auto check = checkChar(input, start_pos, ')'))
        start_pos = check.value();

    else
        return {};
    
    if (auto check = parseName(input, start_pos)){
        auto[name, new_start_pos] = check.value();
        start_pos = new_start_pos;
        node -> setName(name);
    }

    else
        return {}; 

    // construct new topology node with children `children` and name `name`
    
    //add name and children
    return {{node, start_pos}}; 
}

// Length -> ":" number | empty
std::optional<std::pair<optional<double>, int>> parseLength(const std::string& input, int start_pos)
{
    // This parser should always succeed, but might return an empty value.
    return {};
}


// Branch -> Subtree Length
std::optional<std::pair<TopologyNode*, int>> parseBranch(const std::string& input, int start_pos)
{
    return {};
}

// BranchSet -> Branch (, Branch)*
std::optional<std::pair<std::vector<TopologyNode*>, int>> parseBranchSet(const std::string& input, int start_pos)
{
    std::vector<TopologyNode*> branches;
    return {};
}

// Leaf -> name or empty
std::optional<std::pair<TopologyNode*, int>> parseLeaf(const std::string& input, int start_pos)
{
    //adding new node
    auto node = new TopologyNode;
    if (auto check = parseName(input, start_pos)){
        auto [name, new_start_pos] = check.value();
        start_pos = new_start_pos;
        node->setName(name);
    }
    //add name and children
    return {{node, start_pos}}; 
}

// this function is matching (not a quote) or (two quotes)
// 
std::optional<std::pair<char, int>> parseQuotedChar(const std::string& input, int start_pos){
    assert(start_pos>=0);
    if (auto check = parseChar(input, start_pos)){
        auto [c, new_start_pos] = check.value();
        if (c == '\'') {
            // two quotes returns quote
            if (auto new_start_pos2 = checkChar(input, new_start_pos, '\'')){
                // new_start_pos2 will ALWAYS be start_pos + 2
                return optional<pair<char,int>>(pair<char,int>(c,*new_start_pos2));
            }
            // one quote fails
            else{
                return {};
            }
        }
        else {
            //new_start_pos will ALWAYS be start_pos+1
            return optional<pair<char,int>>(pair<char,int>(c,new_start_pos));
        }
    }

    else{
        return {};
    }
}

//no quotes 
std::optional<std::pair<char, int>> parseUnquotedChar(const std::string& input, int start_pos){
    assert(start_pos>=0);
    if (auto check = parseChar(input, start_pos)){
        auto [c, new_start_pos] = check.value();
        //check if c is _
        if (c == '_'){
            return optional<pair<char,int>>(pair<char,int>(' ',new_start_pos));
        }
        //is c an illegal character (punctuation)
        if (!strchr("()[]':;, ", c)) {
            return optional<pair<char,int>>(pair<char,int>(c,new_start_pos));
        }
        else {
            return {};
        }
    }

    else{
        return {};
    }
}

std::optional<std::pair<std::string, int>> parseQuotedName(const std::string& input, int start_pos){
    if (!checkChar(input, start_pos, '\'')){
        return {};
    }
    else{
        start_pos++;
    }
    std::string name; 
    while (auto check = parseQuotedChar(input, start_pos)){
        auto [c, new_start_pos] = check.value();
        name += c;
        start_pos = new_start_pos;
    }
    if (!checkChar(input, start_pos, '\'')){
        return {};
    }
    else{
        start_pos++;
    }
    return optional<pair<std::string, int>>(pair<std::string,int>(name, start_pos));
}

std::optional<std::pair<std::string, int>> parseUnquotedName(const std::string& input, int start_pos){
    std::string name; 
    //checks for unquoted character
    if (auto check = parseUnquotedChar(input, start_pos)){
        auto [c, new_start_pos] = check.value();
        name += c;
        start_pos = new_start_pos;
    }
    //if non return null
    else{
        return {};
    }
    //continue check
    while (auto check = parseUnquotedChar(input, start_pos)){
        auto [c, new_start_pos] = check.value();
        name += c;
        start_pos = new_start_pos;
    }
    return optional<pair<std::string, int>>(pair<std::string,int>(name, start_pos));
}

// Name -> QuotedName | UnquotedName | empty
std::optional<std::pair<std::string, int>> parseName(const std::string& input, int start_pos){
    if (auto check = parseQuotedName(input, start_pos))
        return check;
    
    else {
        return parseUnquotedName(input, start_pos);
    }

}
// This routine has 4 copies of attribute parsing from comments -- 2 for node attributes, and 2 for branch attributes.
// Probably "index" should only be handled in node attributes.  And perhaps only allowed there too, since its a magic attribute.

TopologyNode* NewickConverter::createNode(const std::string &n, std::vector<TopologyNode*> &nodes, std::vector<double> &brlens) {
 // create a string-stream and throw the string into it
 std::stringstream ss (std::stringstream::in | std::stringstream::out);
   ss << n;
   char c = ' ';
   ss.get(c);
   // the initial character has to be '('
   //   fixme: actually, the string 'a;' is valid newick.
    if ( c != '('){
         throw RbException() << "Error while converting Newick tree. We expected an opening parenthesis, but didn't get one. Problematic string: " << n;
    }
    TopologyNode *node = new TopologyNode();
    while ( ss.good() && ss.peek() != ')' )
    {

        TopologyNode *childNode;
        if (ss.peek() == '(' )
        {
            // we received an internal node
            int depth = 0;
            std::string child = "";
            do
            {
                ss.get(c);
                child += c;
                if ( c == '(' )
                {
                    depth++;
                }
                else if ( c == ')' )
                {
                    depth--;
                }

            } while ( ss.good() && depth > 0 );

            // construct the child node
            childNode = createNode( child, nodes, brlens );
        }
        else
        {
            // construct the node
            childNode = new TopologyNode();
        }

        // set the parent child relationship
        node->addChild( childNode );
        childNode->setParent( node );

        // read the optional label
        std::string lbl = "";
        while ( ss.good() && (c = char( ss.peek() ) ) != ':' && c != '[' && c != ';' && c != ',' && c != ')')
        {
            lbl += char( ss.get() );
        }
        childNode->setName( lbl );

        // read the optional node parameters
        if ( ss.peek() == '[' )
        {

            do
            {

                ss.ignore();
                // ignore the '&' before parameter name
                if ( ss.peek() == '&')
                {
                    ss.ignore();
                }

                // read the parameter name
                std::string paramName = "";
                while ( ss.good() && (c = char( ss.peek() ) ) != '=' && c != ',' && c != ']')
                {
                    paramName += char( ss.get() );
                }

                // ignore the equal sign between parameter name and value
                if ( ss.peek() == '=')
                {
                    ss.ignore();
                }

                // read the parameter name
                std::string paramValue = "";
                while ( ss.good() && (c = char( ss.peek() ) ) != ']' && c != ',' && c != ':')
                {
                    paramValue += char( ss.get() );
                }

                childNode->addNodeParameter_(paramName, paramValue);

            } while ( (c = char( ss.peek() ) ) == ',' );

            // ignore the final ']'
            if ( (c = char( ss.peek( ) ) ) == ']' )
            {
                ss.ignore();
            }

        }

        // read the optional branch length
        if ( ss.peek() == ':' )
        {
            ss.ignore();
            std::string time = "";
            while ( ss.good() && (c = char( ss.peek( ) ) ) != ';' && c != ',' && c != ')' && c != '[' )
            {
                time += char( ss.get() );
            }

            std::istringstream stm;
            stm.str(time);
            double d;
            stm >> d;
            nodes.push_back( childNode );
            brlens.push_back( d );
        }
        else
        {
            nodes.push_back( childNode );
            brlens.push_back( 0.0 );
        }

        // read the optional branch parameters
        if ( char( ss.peek() ) == '[' )
        {

            do
            {

                ss.ignore();

                // ignore the '&' before parameter name
                if ( char( ss.peek() ) == '&')
                {
                    ss.ignore();
                }

                // read the parameter name
                std::string paramName = "";
                while ( ss.good() && (c = char( ss.peek() )) != '=' && c != ',')
                {
                    paramName += char( ss.get() );
                }

                // ignore the equal sign between parameter name and value
                if ( char( ss.peek() ) == '=')
                {
                    ss.ignore();
                }

                // read the parameter name
                std::string paramValue = "";
                while ( ss.good() && (c = char( ss.peek() )) != ']' && c != ',' && c != ':')
                {
                    paramValue += char( ss.get() );
                }

                childNode->addBranchParameter(paramName, paramValue);

            } while ( (c = char( ss.peek() )) != ']' );

            ss.ignore();

        }


        // skip comma
        if ( char( ss.peek() ) == ',' )
        {
            ss.ignore();
        }

        // skip comma
        if ( char( ss.peek() ) == ';' )
        {
            // Avoid infinite loop.
            throw RbException()<<"Not enough closing parentheses!";
        }
    }

    // remove closing parenthesis
    ss.ignore();

    // read the optional label, checking for EOF = '\377'
    std::string lbl = "";
    while ( ss.good() && (c = char( ss.peek() )) != ':' && c != ';' && c != ',' && c != '[' && c != '\377')
    {
        lbl += char( ss.get() );
    }
    node->setName( lbl );

    // read the optional node parameters
    if ( char( ss.peek() ) == '[' )
    {

        do
        {

            ss.ignore();

            // ignore the '&' before parameter name
            if ( char( ss.peek() ) == '&')
            {
                ss.ignore();
            }

            // read the parameter name
            std::string paramName = "";
            while ( ss.good() && (c = char( ss.peek() )) != '=' && c != ',')
            {
                paramName += char( ss.get() );
            }

            // ignore the equal sign between parameter name and value
            if ( char( ss.peek() ) == '=')
            {
                ss.ignore();
            }

            // read the parameter name
            std::string paramValue = "";
            while ( ss.good() && (c = char( ss.peek() )) != ']' && c != ',' && c != ':')
            {
                paramValue += char( ss.get() );
            }

            node->addNodeParameter_(paramName, paramValue);

        } while ( (c = char( ss.peek() )) == ',' );

        // ignore the final ']'
        if ( (c = char( ss.peek() )) == ']' )
        {
            ss.ignore();
        }

    }

    // read the optional branch length
    if ( char( ss.peek() ) == ':' )
    {
        ss.ignore();
        std::string time = "";
        while ( ss.good() && (c = char( ss.peek() )) != ';' && c != ',' && c != '[' )
        {
            time += char( ss.get() );
        }

        std::istringstream stm;
        stm.str(time);
        double d;
        stm >> d;
        nodes.push_back( node );
        brlens.push_back( d );
    }
    else
    {
        nodes.push_back( node );
        brlens.push_back( 0.0 );
    }



    // read the optional branch parameters
    if ( char( ss.peek() ) == '[' )
    {

        do
        {

            ss.ignore();

            // ignore the '&' before parameter name
            if ( char( ss.peek() ) == '&')
            {
                ss.ignore();
            }

            // read the parameter name
            std::string paramName = "";
            while ( ss.good() && (c = char( ss.peek() )) != '=' && c != ',')
            {
                paramName += char( ss.get() );
            }

            // ignore the equal sign between parameter name and value
            if ( char( ss.peek() ) == '=')
            {
                ss.ignore();
            }

            // read the parameter name
            std::string paramValue = "";
            while ( ss.good() && (c = char( ss.peek() )) != ']' && c != ',' && c != ':')
            {
                paramValue += char( ss.get() );
            }

            node->addBranchParameter(paramName, paramValue);

        } while ( (c = char( ss.peek() )) != ']' );

        ss.ignore();

    }


    return node;
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
