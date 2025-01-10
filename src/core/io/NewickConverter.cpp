#include "NewickConverter.h"

#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <vector>
#include <string>

#include "BranchHistoryDiscrete.h"
#include "CharacterHistoryDiscrete.h"
#include "RbException.h"
#include "StringUtilities.h"
#include "TopologyNode.h"
#include "Tree.h"
#include "TreeUtilities.h"

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


CharacterHistoryDiscrete* NewickConverter::convertSimmapFromNewick(const std::string& n)
{
    bool reindex = true;

    // create and allocate the tree object
    Tree *t = new Tree();

    std::vector<TopologyNode*> nodes;
    std::vector<double> brlens;
    std::vector<BranchHistory*> histories;

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
    TopologyNode *root = createSimmapNode( trimmed, nodes, brlens, histories );
    nodes.push_back( root );
    brlens.push_back( 0.0 );
    BranchHistoryDiscrete* root_history = new BranchHistoryDiscrete(1,0,0);
    histories.push_back( root_history );

    // set up the tree
    t->setRoot( root, reindex );
    
    std::vector<BranchHistory*> sorted_histories = std::vector<BranchHistory*>(histories.size(), NULL);

    // set the branch lengths
    for (size_t i = 0; i < nodes.size(); ++i)
    {
        t->getNode( nodes[i]->getIndex() ).setBranchLength( brlens[i] );
        sorted_histories[ nodes[i]->getIndex() ] = histories[i];
    }

    // make all internal nodes bifurcating
    // this is important for fossil trees which have sampled ancestors
//    t->makeInternalNodesBifurcating( reindex, true );

    // trees with 2-degree root nodes should not be rerooted
    t->setRooted( root->getNumberOfChildren() == 2 );
    
    // make this tree first a branch length tree
    // hence, tell the root to use branch lengths and not ages (with recursive call)
    root->setUseAges(false, true);
    root->setUseAges(true, true);

    double max_depth = root->getMaxDepth();

    // recursive creation of the tree
    TreeUtilities::setAgesRecursively( *root, max_depth );
    
    // update the ages of the character change events
    for (size_t i = 0; i < nodes.size(); ++i)
    {
        double this_node_age = nodes[i]->getAge();
        BranchHistory* this_history = histories[i];
        const std::multiset<CharacterEvent*,CharacterEventCompare>& old_history = this_history->getHistory();
        std::multiset<CharacterEvent*,CharacterEventCompare> new_history;
        std::multiset<CharacterEvent*,CharacterEventCompare>::const_iterator it = old_history.begin();
        for (; it != old_history.end(); ++it)
        {
            CharacterEvent* this_event = (*it)->clone();
            this_event->setAge( this_event->getAge() + this_node_age );
            new_history.insert( this_event );
        }
        this_history->setHistory( new_history );
    }

//    // if this tree is ultrametric, then we should use ages and transform so!
//    if ( t->isUltrametric() == true )
//    {
//        root->setUseAges(true, true);
//    }
    
    CharacterHistoryDiscrete* new_char_hist = new CharacterHistoryDiscrete();
    new_char_hist->setTree( t );
    new_char_hist->setHistories( sorted_histories );

    // return the tree, the caller is responsible for destruction
    return new_char_hist;
}


// This routine has 4 copies of attribute parsing from comments -- 2 for node attributes, and 2 for branch attributes.
// Probably "index" should only be handled in node attributes.  And perhaps only allowed there too, since its a magic attribute.

TopologyNode* NewickConverter::createNode(const std::string &n, std::vector<TopologyNode*> &nodes, std::vector<double> &brlens)
{

    // create a string-stream and throw the string into it
    std::stringstream ss (std::stringstream::in | std::stringstream::out);
    ss << n;

    char c = ' ';
    ss.get(c);

    // the initial character has to be '('
    //   fixme: actually, the string 'a;' is valid newick.
    if ( c != '(')
    {
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



TopologyNode* NewickConverter::createSimmapNode(const std::string &n, std::vector<TopologyNode*> &nodes, std::vector<double> &brlens, std::vector<BranchHistory*> &histories)
{

    // create a string-stream and throw the string into it
    std::stringstream ss (std::stringstream::in | std::stringstream::out);
    ss << n;

    char c = ' ';
    ss.get(c);

    // the initial character has to be '('
    if ( c != '(')
    {
        throw RbException("Error while converting Newick tree. We expected an opening parenthesis, but didn't get one. Problematic string: " + n);
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
            childNode = createSimmapNode( child, nodes, brlens, histories );
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
        while ( ss.good() && (c = char( ss.peek() ) ) != ':' && c != '[' && c != '{' && c != ';' && c != ',' && c != ')')
        {
            lbl += char( ss.get() );
        }
        childNode->setName( lbl );

        // read the branch length
        if ( ss.peek() == ':' )
        {
            ss.ignore();
            
            // the total branch time
            double time = 0.0;
            BranchHistoryDiscrete* history = new BranchHistoryDiscrete(1,0,0);
            
            // read the node character history
            if ( ss.peek() == '{' )
            {
                std::vector<size_t> states;
                std::vector<double> times;
                do
                {
                    ss.ignore();

                    // read the state
                    std::string state_str = "";
                    while ( ss.good() && (c = char( ss.peek() ) ) != ',' && c != '}')
                    {
                        state_str += char( ss.get() );
                    }
                    size_t state = StringUtilities::asIntegerNumber( state_str );
                    states.push_back( state );

                    // ignore the equal sign between parameter name and value
                    if ( ss.peek() == ',')
                    {
                        ss.ignore();
                    }

                    // read the parameter name
                    std::string duration_str = "";
                    while ( ss.good() && (c = char( ss.peek() ) ) != ';' && c != '}' && c != ':' )
                    {
                        duration_str += char( ss.get() );
                    }
                    double duration = atof( duration_str.c_str() );
                    time += duration;
                    times.push_back( time );

                } while ( (c = char( ss.peek() ) ) == ';' || c == ':' );

                // ignore the final '}'
                if ( (c = char( ss.peek( ) ) ) == '}' )
                {
                    ss.ignore();
                }
                
                for ( size_t i=1; i<states.size(); ++i )
                {
                    CharacterEventDiscrete* this_hist_event = new CharacterEventDiscrete(0, states[i-1], times[i-1]);
                    history->addEvent( this_hist_event );
                }
                
                // add the parent and child histories
                std::vector<CharacterEvent*> child_characters;
                CharacterEventDiscrete* this_child_char = new CharacterEventDiscrete(0, states[0], 1.0);
                child_characters.push_back( this_child_char );
                history->setChildCharacters( child_characters );

                std::vector<CharacterEvent*> parent_characters;
                CharacterEventDiscrete* this_parent_char = new CharacterEventDiscrete(0, states[states.size()-1], 0.0);
                parent_characters.push_back( this_parent_char );
                history->setParentCharacters( parent_characters );

            }
            
            nodes.push_back( childNode );
            brlens.push_back( time );
            histories.push_back( history );
        }
        else
        {
            nodes.push_back( childNode );
            brlens.push_back( 0.0 );
        }


        // skip comma
        if ( char( ss.peek() ) == ',' )
        {
            ss.ignore();
        }

        // skip semi-colon
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
    while ( ss.good() && (c = char( ss.peek() )) != ':' && c != ';' && c != ',' && c != '{' && c != '\377')
    {
        lbl += char( ss.get() );
    }
    node->setName( lbl );

//    // read the optional branch length
//    if ( char( ss.peek() ) == ':' )
//    {
//        
//        ss.ignore();
//            
//        // the total branch time
//        std::string time = "";
//            
//        // read the node character history
//        if ( ss.peek() == '{' )
//        {
//            std::vector<size_t> states;
//            std::vector<double> times;
//            do
//            {
//                ss.ignore();
//                
//                // read the state
//                std::string state_str = "";
//                while ( ss.good() && (c = char( ss.peek() ) ) != ',' && c != '}')
//                {
//                    state_str += char( ss.get() );
//                }
//                size_t state = StringUtilities::asIntegerNumber( state_str );
//                states.push_back( state );
//
//                // ignore the equal sign between parameter name and value
//                if ( ss.peek() == ',')
//                {
//                    ss.ignore();
//                }
//
//                // read the parameter name
//                std::string duration_str = "";
//                while ( ss.good() && (c = char( ss.peek() ) ) != ';' && c != '}' )
//                {
//                    duration_str += char( ss.get() );
//                }
//                double duration = atof( duration_str.c_str() );
//                times.push_back( duration );
//                time += duration;
//
//            } while ( (c = char( ss.peek() ) ) == ';' );
//
//            // ignore the final '}'
//            if ( (c = char( ss.peek( ) ) ) == '}' )
//            {
//                ss.ignore();
//            }
//        }
//    }
//    else
//    {
//        nodes.push_back( node );
//        brlens.push_back( 0.0 );
//    }


    return node;
}


