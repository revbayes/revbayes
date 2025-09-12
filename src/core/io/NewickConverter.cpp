#include "NewickConverter.h"

#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <vector>
#include <string>

#include "RbException.h"
#include "TopologyNode.h"
#include "Tree.h"

using namespace RevBayesCore;

NewickConverter::NewickConverter()
{

}


NewickConverter::~NewickConverter()
{

}



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

// subtree -> internal OR leaf
std::optional<std::pair<TopologyNode*, int>> NewickConverter::parseSubTree(const std::string input, int start_pos){

    if (auto check = parseInternal(input, start_pos))
        return check;
    
    else {
        return parseLeaf(input, start_pos);
    }
}

//
std::optional<int> NewickConverter::parseChar(const std::string input, int start_pos, char c){
    // if reading beyond end of string return null
    if (start_pos >= n.size()){
        return {};
    }
    else if (input[start_pos] != c) {
        return {};
    }
    else return start_pos+1;
}
// Internal -> '(' BranchSet ')' Name
std::optional<std::pair<TopologyNode*, int>> NewickConverter::parseInternal(const std::string& input, int start_pos){
    // Check if we have left parenthesis
    if (auto check = parseChar(input, start_pos, '('))
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
    if (auto check = parseChar(input, start_pos, ')'))
        start_pos = check.value();

    else
        return {};
    
    if (auto check = parseName(input, start_pos)){
        auto[name, new_start_pos] = check.value();
        start_pos = new_start_pos;
    }

    else
        return {}; 

    // construct new topology node with children children and name name
}

std::optional<std::pair<std::string, int>> NewickConverter::parseName(const std::string& input, int start_pos){
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

    // is name quoted?
    if (input[pos] == '\'') {          
        ++pos; // index past quote               
        std::string out;
        // need to verify this loop is working as intended
        // idea is if we seee a single quote it may be an internal quote not the terminal one
        while (pos < static_cast<int>(input.size())) {
            char c = input[pos];
            if (c == '\'') {
                if (pos + 1 < static_cast<int>(input.size()) && input[pos + 1] == '\'') {
                    // just learned about push_back may not work
                    out.push_back('\'');
                    pos += 2;
                    continue;
                } else {
                    // Closing quote found
                    ++pos;
                    return std::make_pair(out, pos);
                }
            }
            out.push_back(c);
            ++pos;
        }

        return std::nullopt; // no closing quote found
    }
    //unquoted case, is character important character
    auto isdelim = [](char c) {
        switch (c) {
            case '(': case ')': case ',': case ':': case ';':
            case '[': case ']':
            case ' '
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

TopologyNode* NewickConverter::createNode(const std::string &n, int& start_pos, std::vector<TopologyNode*> &nodes, std::vector<double> &brlens) {
    
    // if reading beyond end of string return null
    if (start_pos >= n.size()){
        return nullptr;                                                                                                                                                                        
    }
    char c = n[start_pos];

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
