#include <stddef.h>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include "DelimitedCharacterDataWriter.h"
#include "HomologousDiscreteCharacterData.h"
#include "NclReader.h"
#include "TreeDiscreteCharacterData.h"
#include "AbstractHomologousDiscreteCharacterData.h"
#include "Cloneable.h"
#include "RbException.h"
#include "RbFileManager.h"
#include "TopologyNode.h"
#include "Tree.h"
#include "StandardState.h" // IWYU pragma: keep

namespace RevBayesCore { class AbstractCharacterData; }


using namespace RevBayesCore;

TreeDiscreteCharacterData::TreeDiscreteCharacterData() :
        character_data(NULL)
{
    
}


TreeDiscreteCharacterData* TreeDiscreteCharacterData::clone( void ) const
{
    return new TreeDiscreteCharacterData( *this );
}



AbstractHomologousDiscreteCharacterData& TreeDiscreteCharacterData::getCharacterData( void )
{
    
    return *character_data;
}



const AbstractHomologousDiscreteCharacterData& TreeDiscreteCharacterData::getCharacterData( void ) const
{

    return *character_data;
}


std::vector<double> TreeDiscreteCharacterData::getTimeInStates( void )
{
    return time_in_states;
}


bool TreeDiscreteCharacterData::hasCharacterData( void ) const
{

    return character_data != NULL;
}


/**
 * Initialize this object from a file
 *
 * \param[in]   idx    The site at which we want to know if it is constant?
 */
void TreeDiscreteCharacterData::initFromFile(const std::string &dir, const std::string &fn)
{
    RbFileManager fm = RbFileManager(dir, fn + ".nex");
    
    // get an instance of the NCL reader
    NclReader reader = NclReader();
    
    std::string myFileType = "nexus";
    std::string dType = character_data->getDataType();
    
    std::string suffix = "|" + dType;
    suffix += "|noninterleaved";
    myFileType += suffix;
    
    // read the content of the file now
    std::vector<AbstractCharacterData*> m_i = reader.readMatrices( fm.getFullFileName(), myFileType );
    
    if ( m_i.size() < 1 )
    {
        const std::set<std::string> &warnings = reader.getWarnings();
        for (std::set<std::string>::iterator it = warnings.begin(); it != warnings.end(); ++it)
        {
            std::cerr << "NCL-Warning:\t\t" << *it << std::endl;
        }
        throw RbException("Could not read character data matrix from file \"" + fm.getFullFileName() + "\".");
    }
    
    HomologousDiscreteCharacterData<StandardState> *coreM = static_cast<HomologousDiscreteCharacterData<StandardState> *>( m_i[0] );
    
    character_data = coreM;
    
}


/**
 * Initialize this object from a file
 *
 * \param[in]   s    The value as string
 */
void TreeDiscreteCharacterData::initFromString(const std::string &s)
{
    throw RbException("Cannot initialize a tree with a discrete character data matrix from a string.");
}

void TreeDiscreteCharacterData::setCharacterData( AbstractHomologousDiscreteCharacterData *d )
{
	if ( character_data != NULL )
	{
		delete character_data;
	}
    character_data = d;
}


void TreeDiscreteCharacterData::setTimeInStates( std::vector<double> t )
{
    time_in_states = t;
}


void TreeDiscreteCharacterData::writeToFile(const std::string &dir, const std::string &fn) const
{
    // do not write a file if the tree is invalid
    if (this->getNumberOfTips() > 1)
    {
        RbFileManager fm = RbFileManager(dir, fn + ".newick");
        fm.createDirectoryForFile();
        
        // open the stream to the file
        std::fstream o;
        o.open( fm.getFullFileName().c_str(), std::fstream::out);

        // write the value of the node
        o << getNewickRepresentation();
        o << std::endl;

        // close the stream
        o.close();
       
        // many SSE models use NaturalNumber states, which are incompatible
        // with the NEXUS format, so write the tips states to a separate
        // tab-delimited file
        fm = RbFileManager(dir, fn + ".tsv");
        RevBayesCore::DelimitedCharacterDataWriter writer; 
        writer.writeData(fm.getFullFileName(), *character_data);

        // write the character history's time spent in each state
        fm = RbFileManager(dir, fn + "_time_in_states.tsv");
        o.open( fm.getFullFileName().c_str(), std::fstream::out);
        for (size_t i = 0; i < time_in_states.size(); i++)
        {
            o << time_in_states[i]; 
            if (i < (time_in_states.size() - 1))
            {
                o << "\t";
            }
        }
        o << std::endl;
        o.close();
    }
}


void TreeDiscreteCharacterData::setTree(const Tree &t)
{
    
    nodes.clear();
    delete root;
    root = NULL;
    
    num_tips    = t.getNumberOfTips();
    num_nodes   = t.getNumberOfNodes();
    rooted      = t.isRooted();
    
    TopologyNode* newRoot = t.getRoot().clone();
    
    // set the root. This will also set the nodes vector
    // do not reorder node indices when copying (WP)
    setRoot(newRoot, false);
}
