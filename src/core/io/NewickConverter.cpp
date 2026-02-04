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
    // construct the tree starting from the root
    auto result = parseTree(newick, 0);
    if (not result){
        throw RbException()<<result.err_message()<< " at location " << result.err_pos();
    }
    if (result.next_pos() != newick.size())
        throw RbException()<<"Junk at end of newick string: '"<<newick.substr(result.next_pos())<<"'";

    // set up the tree
    auto tree = result.value();

    // try to set node indices from attributes
    tree->tryReadIndicesFromParameters(true);

    // trees with 2-degree root nodes should not be rerooted
    tree->setRooted( tree->getRoot().getNumberOfChildren() == 2 );
    
    // make this tree first a branch length tree
    // hence, tell the root to use branch lengths and not ages (with recursive call)
    tree->getRoot().setUseAges(false, true);
    
    // return the tree, the caller is responsible for destruction
    return tree;
}


/* We are currently using this Newick grammar:
  
    Tree → Subtree ";"
    Subtree → [Descendants] [Name] [Branch]
    Descendants = "(" Subtree ("," Subtree)* ")"
    Name -> QuotedName | UnquotedName
    Branch -> ":" [Double]
    QuotedName -> "'" QuotedChar* "'"
    UnquotedName -> UnquotedChar*
    UnquotedChar = not punctuation

When we put the spacing in, we get:

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
    else
        start_pos = check_open_bracket.next_pos();

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
        return ParseSuccess<optional<string>>({}, check_char.next_pos());
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

std::string trim(const std::string& s) {
    // Define all characters considered "whitespace"
    const std::string WHITESPACE = " \n\r\t\f\v";

    // Find the first character that is NOT whitespace
    size_t start = s.find_first_not_of(WHITESPACE);
    
    // If the string is all whitespace, return an empty string
    if (start == std::string::npos) {
        return "";
    }

    // Find the last character that is NOT whitespace
    size_t end = s.find_last_not_of(WHITESPACE);

    // Return the substring between start and end
    return s.substr(start, end - start + 1);
}

/// Split by commas that are not inside { and }, then trim whitespace
/// and add non-empty chunks to the output list.
std::vector<std::string> split_comment(const std::string& comment)
{
// "[&B=2,C=4]" -> {"B=2", "C=4" }
// "[&A=1]" -> { "A=1"}
// "[  A=1  ]" -> { "A=1"}
// "[  A=1 , R  ]" -> { "A=1", "R"}
// "[  A=1 ,]" -> { "A=1"}
// "[ , ]" -> { }
// "[,]" -> { }

    // The comment does NOT include the square brackets
    assert(not comment.ends_with(']'));

    // If the comment doesn't start with & or % then ignore it.
    if (not (comment.starts_with('&') or comment.starts_with('%'))) return {};

    // Note that this function does NOT handle NHX comments.
    vector<std::string> chunks;
    int start = 1;
    int depth = 0;

    // Handle chunks that are not the last chunk
    for(int i=start;i<comment.size();i++)
    {
        char c = comment[i];

        if (c == ',' and depth == 0)
        {
            int length = i-start;
            auto chunk = trim( comment.substr(start,length) );
            if (not chunk.empty())
                chunks.push_back( chunk );
            start = i+1;
        }

        // This is from BEAST.  Thanks for making our lives more complicated.
        else if (c == '{')
            depth++;

        else if (c == '}')
            depth--;
    }

    // Handle the last chunk
    auto last_chunk = trim( comment.substr(start) );
    if (not last_chunk.empty())
        chunks.push_back( last_chunk );
    
    return chunks;
}
    
// Tree -> SPACE<tree> Subtree ";"
ParseResult<Tree*> parseTree(const std::string& input, int start_pos)
{
    vector<string> tree_comments;

    // 1. Skip whitespace and record tree_comments
    if (auto maybe_whitespace = parseWhitespace(input, start_pos))
    {
        tree_comments = std::move(maybe_whitespace.value());
        start_pos = maybe_whitespace.next_pos();
    }
    else if (maybe_whitespace.hard_failure())
        return maybe_whitespace.as_hard_failure();
    
    // 2. Get the Subtree
    TopologyNode* root = nullptr;
    if (auto check_subtree = parseSubTree(input,start_pos))
    {
        start_pos = check_subtree.next_pos();
        root = check_subtree.value();
    }
    else
        return check_subtree.as_failure();

    // 3. Check the semicolon
    if (auto check_semi = checkChar(input, start_pos, ';'))
        start_pos = check_semi.next_pos();
    else
        return check_semi.as_failure();;

    // 4. Create the tree
    Tree* tree = new Tree();
    
    tree->setRoot( root, true );

    for(auto& comment: tree_comments)
        for(auto& chunk: split_comment(comment))
            tree->addParameter_(chunk);

    return ParseSuccess(tree, start_pos);
}

// subtree -> SPACE [Descendants] SPACE [Name] SPACE [Branch]
ParseResult<TopologyNode*> parseSubTree(const std::string& input, int start_pos)
{
    auto node = new TopologyNode;
    vector<string> node_comments;

    // Skip whitespace and record node_comments
    if (auto maybe_whitespace = parseWhitespace(input, start_pos))
    {
        node_comments = std::move(maybe_whitespace.value());
        start_pos = maybe_whitespace.next_pos();
    }
    else if (maybe_whitespace.hard_failure())
        return maybe_whitespace.as_hard_failure();
    
    // Parse [Descendants]
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
    
    // Skip whitespace and record node_comments
    if (auto maybe_whitespace = parseWhitespace(input, start_pos))
    {
        auto& new_comments = maybe_whitespace.value();
        std::move(new_comments.begin(), new_comments.end(), std::back_inserter(node_comments));
        start_pos = maybe_whitespace.next_pos();
    }
    else if (maybe_whitespace.hard_failure())
        return maybe_whitespace.as_hard_failure();
    
    // Parse [Name]
    if (auto maybe_name = parseName(input, start_pos)) {
        start_pos = maybe_name.next_pos();
        node -> setName(maybe_name.value());
    }
    else if (maybe_name.hard_failure())
        return maybe_name.as_failure();

    // Skip whitespace and record node_comments
    if (auto maybe_whitespace = parseWhitespace(input, start_pos))
    {
        auto& new_comments = maybe_whitespace.value();
        std::move(new_comments.begin(), new_comments.end(), std::back_inserter(node_comments));
        start_pos = maybe_whitespace.next_pos();
    }
    else if (maybe_whitespace.hard_failure())
        return maybe_whitespace.as_hard_failure();
    
    // Parse [Branch]
    if (auto maybe_branch = parseBranch(input, start_pos))
    {
        start_pos = maybe_branch.next_pos();
        auto& [length, branch_comments] = maybe_branch.value();
        if (length)
            node->setBranchLength(*length);

        // Add branch comments
        for(auto& comment: branch_comments)
            for(auto& chunk: split_comment(comment))
                node->addBranchParameter_(chunk);
    }
    else if (maybe_branch.hard_failure())
        return maybe_branch.as_failure();

    // Add node comments
    for(auto& comment: node_comments)
        for(auto& chunk: split_comment(comment))
            node->addNodeParameter_(chunk);
    
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

// Length -> ":" SPACE [number] SPACE | empty
ParseResult<pair<optional<double>,vector<string>>> parseBranch(const std::string& input, int start_pos)
{
    optional<double> length;
    vector<string> comments;

    // Check for colon
    if (auto maybe_colon = checkChar(input, start_pos, ':'))
        start_pos = maybe_colon.next_pos();
    else
        return ParseSuccess<pair<optional<double>,vector<string>>>({length, comments}, start_pos);

    // Skip whitespace and record comments
    if (auto maybe_whitespace = parseWhitespace(input, start_pos))
    {
        comments = std::move(maybe_whitespace.value());
        start_pos = maybe_whitespace.next_pos();
    }
    else if (maybe_whitespace.hard_failure())
        return maybe_whitespace.as_hard_failure();
    
    // Parse Length
    if (auto maybe_length = parseDouble(input, start_pos))
    {
        start_pos = maybe_length.next_pos();
        length = maybe_length.value();
    }
    else if (maybe_length.hard_failure())
        return maybe_length.as_hard_failure();
    
    // Skip whitespace and record comments
    if (auto maybe_whitespace = parseWhitespace(input, start_pos))
    {
        auto& new_comments = maybe_whitespace.value();
        std::move(new_comments.begin(), new_comments.end(), std::back_inserter(comments));
        start_pos = maybe_whitespace.next_pos();
    }
    else if (maybe_whitespace.hard_failure())
        return maybe_whitespace.as_hard_failure();
    
    return ParseSuccess<pair<optional<double>, vector<string>>>({length, comments}, start_pos);
}

// Descendants -> "(" > Branch ("," > Branch)* > ")"
ParseResult<std::vector<TopologyNode*>> parseDescendants(const std::string& input, int start_pos)
{
    // Read left parenthesis
    if (auto maybe_lparen = checkChar(input, start_pos, '('))
        start_pos = maybe_lparen.next_pos();
    else
        return maybe_lparen.as_failure();

    vector<TopologyNode*> subtrees;
    auto maybe_subtree = parseSubTree(input, start_pos);
    if (not maybe_subtree)
        return maybe_subtree.as_hard_failure();
    subtrees.push_back(maybe_subtree.value());
    start_pos = maybe_subtree.next_pos();

    while(auto maybe_comma = checkChar(input, start_pos, ','))
    {
        start_pos = maybe_comma.next_pos();
        auto maybe_subtree2 = parseSubTree(input, start_pos);
        if (not maybe_subtree2)
            return maybe_subtree2.as_hard_failure();
        subtrees.push_back(maybe_subtree2.value());
        start_pos = maybe_subtree2.next_pos();
    }


    // Read right parenthesis
    if (auto maybe_rparen = checkChar(input, start_pos, ')')) {
        start_pos = maybe_rparen.next_pos();
    }
    else
        return maybe_rparen.as_hard_failure();

    return ParseSuccess(subtrees, start_pos);
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
        if (c == '_' and false)
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
    
    CharacterHistoryDiscrete* new_char_hist = new CharacterHistoryDiscrete();
    new_char_hist->setTree( t );
    new_char_hist->setHistories( sorted_histories );
    
    size_t max_obs_state = new_char_hist->getMaxObservedState();
    new_char_hist->setNumberOfStates( max_obs_state );

    // return the tree, the caller is responsible for destruction
    return new_char_hist;
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


    return node;
}

}
