#ifndef JointConditionalJointConditionalAncestralStateMonitor_H
#define JointConditionalJointConditionalAncestralStateMonitor_H

#include "AbstractHomologousDiscreteCharacterData.h"
#include "StateDependentSpeciationExtinctionProcess.h"
#include "GeneralizedLineageHeterogeneousBirthDeathSamplingProcess.h"
#include "VariableMonitor.h"
#include "Tree.h"
#include "TypedDagNode.h"
#include "StochasticNode.h"
#include "DnaState.h"

#include <fstream>
#include <string>
#include <vector>
#include <typeinfo>

namespace RevBayesCore {
    
    /**
     * @brief Declaration and implementation of the JointConditionalAncestralStateMonitor class.
     *
     * @file
	 * Declaration and implementation of the JointConditionalAncestralStateMonitor class which
	 * monitors the ancestral states of each internal node and
	 * prints their value into a file.
     *
     * @copyright Copyright 2014-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @version 1.0
     *
     */
    template<class characterType>
    class JointConditionalAncestralStateMonitor : public VariableMonitor {
        
    public:
        // Constructors and Destructors
		JointConditionalAncestralStateMonitor(TypedDagNode<Tree> *t, StochasticNode<AbstractHomologousDiscreteCharacterData>* ch, unsigned long g, const std::string &fname, const std::string &del, bool wt, bool wss);                                  //!< Constructor
		JointConditionalAncestralStateMonitor(StochasticNode<Tree>* ch, unsigned long g, const std::string &fname, const std::string &del, bool wt, bool wss);
        virtual ~JointConditionalAncestralStateMonitor(void);
        
        JointConditionalAncestralStateMonitor*          clone(void) const;                                                  //!< Clone the object
        
        // functions you may want to overwrite
        virtual void                                    monitorVariables(unsigned long gen);                                //!< Monitor at generation gen
        virtual void                                    printFileHeader(void);                                              //!< Print header
        virtual void                                    swapNode(DagNode *oldN, DagNode *newN);

    private:
        
        // members
		TypedDagNode<Tree>*                             tree;
		StochasticNode<AbstractHomologousDiscreteCharacterData>*            ctmc;
        StochasticNode<Tree>*                                               cdbdp;
        bool                                            withTips;
        bool                                            withStartStates;
    };
    
}


#include "DagNode.h"
#include "Model.h"
#include "Monitor.h"
#include "RbFileManager.h"
#include "PhyloCTMCSiteHomogeneous.h"
#include "AbstractPhyloCTMCSiteHomogeneous.h"

using namespace RevBayesCore;


/* Constructor for state dependent birth death process */
template<class characterType>
JointConditionalAncestralStateMonitor<characterType>::JointConditionalAncestralStateMonitor(StochasticNode<Tree>* ch, unsigned long g, const std::string &fname, const std::string &del, bool wt, bool wss) :
    VariableMonitor(ch, g, fname, del, false, false, false),
    cdbdp( ch ),
    withTips( wt ),
    withStartStates( wss )
{
    ctmc = NULL;
    
    // the cdbdp is both the tree and character evolution model
    addVariable( cdbdp );
    tree = cdbdp;
}


/* Constructor for CTMC */
template<class characterType>
JointConditionalAncestralStateMonitor<characterType>::JointConditionalAncestralStateMonitor(TypedDagNode<Tree> *t, StochasticNode<AbstractHomologousDiscreteCharacterData>* ch, unsigned long g, const std::string &fname, const std::string &del, bool wt, bool wss) :
    VariableMonitor(ch, g, fname, del, false, false, false),
    tree( t ),
    ctmc( ch ),
    withTips( wt ),
    withStartStates( wss )
{
    cdbdp = NULL;
    
    // the cdbdp is both the tree and character evolution model
    addVariable( tree );
    addVariable( ctmc );
    
}


/**
 * Destructor.
 */
template<class characterType>
JointConditionalAncestralStateMonitor<characterType>::~JointConditionalAncestralStateMonitor()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself
 */
template<class characterType>
JointConditionalAncestralStateMonitor<characterType>* JointConditionalAncestralStateMonitor<characterType>::clone(void) const
{
    
    return new JointConditionalAncestralStateMonitor<characterType>(*this);
}



/**
 * Monitor value at given generation.
 *
 * \param[in]   gen    The current generation.
 */
template<class characterType>
void JointConditionalAncestralStateMonitor<characterType>::monitorVariables(unsigned long gen)
{
    auto& separator = to<SeparatorFormat>(format)->separator;
    
    size_t num_sites;
    size_t num_nodes;
        
        
    // get the distribution for the character
    AbstractPhyloCTMCSiteHomogeneous<characterType>          *dist_ctmc    = NULL;
    StateDependentSpeciationExtinctionProcess                *dist_bd      = NULL;
    GeneralizedLineageHeterogeneousBirthDeathSamplingProcess *dist_glhbdsp = NULL;
    if ( ctmc != NULL )
    {
        dist_ctmc = static_cast<AbstractPhyloCTMCSiteHomogeneous<characterType>* >( &ctmc->getDistribution() );
        num_sites = dist_ctmc->getValue().getNumberOfCharacters();
        num_nodes = tree->getValue().getNumberOfNodes();
    }
    else
    {
        dist_bd = dynamic_cast<StateDependentSpeciationExtinctionProcess*>( &cdbdp->getDistribution() ); // this does!
        if ( dist_bd == NULL )
        {
        	dist_glhbdsp = dynamic_cast<GeneralizedLineageHeterogeneousBirthDeathSamplingProcess*>( &cdbdp->getDistribution() );
        }
        num_sites = 1;
        num_nodes = tree->getValue().getNumberOfNodes();
    }
        
    std::vector<std::vector<characterType> > startStates(num_nodes, std::vector<characterType>(num_sites));
    std::vector<std::vector<characterType> > endStates(num_nodes, std::vector<characterType>(num_sites));
    
    // draw ancestral states
    if ( ctmc != NULL )
    {
        dist_ctmc->drawJointConditionalAncestralStates(startStates, endStates);
    }
    else if ( dist_bd != NULL )
    {
        std::vector<size_t> startStatesIndexes(num_nodes);
        std::vector<size_t> endStatesIndexes(num_nodes);
        dist_bd->drawJointConditionalAncestralStates(startStatesIndexes, endStatesIndexes);
        
        // Let us check first the type of the data and the one we expect.
        characterType tmp = characterType();
        if ( dist_bd->getCharacterData().getTaxonData(0)[0].getDataType() != tmp.getDataType() )
        {
            throw RbException() << "The character type in the ancestral state monitor does not match. \" The data has type " << dist_bd->getCharacterData().getTaxonData(0)[0].getDataType() << "\" but the monitor expected \"" << tmp.getDataType() << "\".";
        }
            
        // now give as an object that we can clone.
        // this is necessary because otherwise we would not know the state labels or size of the character
        characterType *tmp_char = dynamic_cast< characterType* >( dist_bd->getCharacterData().getTaxonData(0)[0].clone() );
            
        for (size_t i = 0; i < startStatesIndexes.size(); i++)
        {
            startStates[i][0] = characterType( *tmp_char );
            startStates[i][0].setStateByIndex(startStatesIndexes[i]);
            startStates[i][0].setMissingState(false);
            endStates[i][0]   = characterType( *tmp_char );
            endStates[i][0].setStateByIndex(endStatesIndexes[i]);
            endStates[i][0].setMissingState(false);
        }
            
        delete tmp_char;
            
    }
    else if ( dist_glhbdsp != NULL )
    {
        std::vector<size_t> startStatesIndexes(num_nodes);
        std::vector<size_t> endStatesIndexes(num_nodes);

    	dist_glhbdsp->drawJointConditionalAncestralStates(startStatesIndexes, endStatesIndexes);

        // Let us check first the type of the data and the one we expect.
        characterType tmp = characterType();
        if ( dist_glhbdsp->getCharacterData().getTaxonData(0)[0].getDataType() != tmp.getDataType() )
        {
            throw RbException() << "The character type in the ancestral state monitor does not match. \" The data has type " << dist_glhbdsp->getCharacterData().getTaxonData(0)[0].getDataType() << "\" but the monitor expected \"" << tmp.getDataType() << "\".";
        }

        // now give as an object that we can clone.
        // this is necessary because otherwise we would not know the state labels or size of the character
        characterType *tmp_char = dynamic_cast< characterType* >( dist_glhbdsp->getCharacterData().getTaxonData(0)[0].clone() );

        for (size_t i = 0; i < startStatesIndexes.size(); i++)
        {
            startStates[i][0] = characterType( *tmp_char );
            startStates[i][0].setStateByIndex(startStatesIndexes[i]);
            startStates[i][0].setMissingState(false);
            endStates[i][0]   = characterType( *tmp_char );
            endStates[i][0].setStateByIndex(endStatesIndexes[i]);
            endStates[i][0].setMissingState(false);
        }

        delete tmp_char;

    }
    else
    {
    	throw RbException("No valid distribution found.");
    }
        
    // print ancestral states
    const std::vector<TopologyNode*>& nds = tree->getValue().getNodes();
    for (int i = 0; i < nds.size(); i++)
    {
        size_t node_index = nds[i]->getIndex();
        
        // start states
        if (withStartStates)
        {
            // add a separator before every new element
            out_stream << separator;

            // print out ancestral states....
            for (int j = 0; j < startStates[node_index].size(); j++)
            {
                if (j != 0)
                    out_stream << ",";
                out_stream << startStates[node_index][j].getStringValue();
            }
        }
            
        // end states
        if ( withTips || !nds[i]->isTip() )
        {
            out_stream << separator;
            for (int j = 0; j < endStates[node_index].size(); j++)
            {
                if (j != 0)
                    out_stream << ",";
                out_stream << endStates[node_index][j].getStringValue();
            }
        }
    }
    
}


/**
 * Print header for monitored values
 */
template<class characterType>
void JointConditionalAncestralStateMonitor<characterType>::printFileHeader()
{
    auto& separator = to<SeparatorFormat>(format)->separator;
    
	std::vector<TopologyNode*> nodes = tree->getValue().getNodes();
	
	// iterate through all tree nodes and make header with node index
	for (int i = 0; i < tree->getValue().getNumberOfNodes(); i++)
    {
		TopologyNode* nd = nodes[i];
        size_t node_index = nd->getIndex();
		if (withStartStates)
        {
            out_stream << separator;
            out_stream << "start_" << node_index+1;
        }
        
        // end states
        if ( withTips || !nd->isTip() )
        {
			out_stream << separator;
			out_stream << "end_" << node_index+1;
		}
    }
    
}


template<class characterType>
void JointConditionalAncestralStateMonitor<characterType>::swapNode(DagNode *oldN, DagNode* newN)
{
	
	if ( oldN == tree )
	{
		tree = static_cast< TypedDagNode<Tree> *>( newN );
	}
	else if ( oldN == ctmc )
	{
		ctmc = static_cast< StochasticNode<AbstractHomologousDiscreteCharacterData> *>( newN );
	}
    else if ( oldN == cdbdp )
    {
        cdbdp = static_cast< StochasticNode<Tree> *>( newN );
        tree = static_cast< StochasticNode<Tree> *>( newN );
    }
    
    Monitor::swapNode( oldN, newN );
	
}

#endif 
