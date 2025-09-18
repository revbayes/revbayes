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
    - ???

*/

// Tree -> Subtree ";"
std::optional<std::pair<TopologyNode*, int>> parseTree(const std::string& input, int start_pos)
{
    // 1. Get the Subtree
    auto check_subtree = parseSubTree(input,start_pos);
    if (not check_subtree)
        return {};

    auto& [subtree, new_start_pos] = *check_subtree;

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

// Internal -> '(' BranchSet ')' Name
std::optional<std::pair<TopologyNode*, int>> parseInternal(const std::string& input, int start_pos){
    // Check if we have left parenthesis
    if (auto check = checkChar(input, start_pos, '('))
        start_pos = check.value();

    else
        return {};

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
    }

    else
        return {}; 

    // construct new topology node with children `children` and name `name`
    auto node = new TopologyNode;
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


// Leaf -> Name
std::optional<std::pair<TopologyNode*, int>> parseLeaf(const std::string& input, int start_pos)
{
    return {};
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
    return {};
}

// Name -> QuotedName | UnquotedName | empty
std::optional<std::pair<std::string, int>> parseName(const std::string& input, int start_pos){
    //handle newick escaping hear, read newick minus parsing syntax?
    //* skip any whitespace
    // quoted name or non-quoted name, make two new functions for this
    // parse quoted name: if failed, instead of {}, call parse-unquoted name
    // return call parse unquoted
    
    // check if starting position is out of bounds
    // i think static cast is important to have .size be an int
    if (start_pos < 0 || start_pos >= static_cast<int>(input.size())) {
        return std::nullopt;
    }

    // skip whitespace
    int pos = start_pos;
    while (pos < static_cast<int>(input.size())) {
        char c = input[pos];
        if (c != ' ' && c != '\t' && c != '\r' && c != '\n') {
            break;
        }
        ++pos;
    }
    if (pos >= static_cast<int>(input.size())) {
        return std::nullopt;
    }


    //unquoted case, is character important character
    auto isdelim = [](char c) {
        switch (c) {
            case '(': case ')': case ',': case ':': case ';':
            case '[': case ']':
            case ' ':
                return true;
            default:
                return false;
        }
    };

    std::string out;
    while (pos < static_cast<int>(input.size())) {
        char c = input[pos];
        if (isdelim(c)) break;  // stop at a delimiter
        out.push_back(c);
        ++pos;
    }

    if (out.empty()) {
        return std::nullopt;
    }

    return std::make_pair(out, pos);

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
