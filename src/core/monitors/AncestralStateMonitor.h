#ifndef AncestralStateMonitor_H
#define AncestralStateMonitor_H

#include "VariableMonitor.h"
#include "Tree.h"
#include "TypedDagNode.h"
#include "StochasticNode.h"
#include "DnaState.h"

#include <fstream>
#include <string>
#include <vector>
#include <typeinfo>

#define FAILED 0
#define PASSED 1
#define NOT_CHECKED 2

namespace RevBayesCore {
    
    /**
     * @brief Declaration and implementation of the AncestralStateMonitor class.
     *
     * @file 
	 * Declaration and implementation of the AncestralStateMonitor class which
	 * monitors the ancestral states of each internal node and
	 * prints their value into a file.
     *
     * @copyright Copyright 2014-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @version 1.0
     *
     */
    template<class characterType> 
    class AncestralStateMonitor : public VariableMonitor {
	
    public:
        // Constructors and Destructors
		AncestralStateMonitor(TypedDagNode<Tree> *t, DagNode* &ch, unsigned long g, const std::string &fname, const std::string &del);                                  //!< Constructor
		
        AncestralStateMonitor(const AncestralStateMonitor &m);
        virtual ~AncestralStateMonitor(void);
        
        AncestralStateMonitor*              clone(void) const;                                                  //!< Clone the object
        
        // Monitor functions
        void                                monitorVariables(unsigned long gen);                                //!< Monitor at generation gen
        void                                printFileHeader(void);                                              //!< Print header
        
        // getters and setters
		void								swapNode(DagNode *oldN, DagNode *newN);
		
    private:
        
        // members
		TypedDagNode<Tree>*                 tree;
		DagNode*                            ctmc;
		bool                                stochastic_nodes_only;
    };
    
}


#include "DagNode.h"
#include "Model.h"
#include "Monitor.h"
#include "RbFileManager.h"
#include "PhyloCTMCSiteHomogeneous.h"
#include "AbstractPhyloCTMCSiteHomogeneous.h"

using namespace RevBayesCore;


/* Constructor */
template<class characterType>
AncestralStateMonitor<characterType>::AncestralStateMonitor(TypedDagNode<Tree> *t, DagNode* &ch, unsigned long g, const std::string &fname, const std::string &del) : VariableMonitor(std::vector<DagNode*>(), g, fname, del, false, false, false),
    tree( t ),
    ctmc( ch ),
    stochastic_nodes_only( false )
{
	
	nodes.push_back( tree );
	nodes.push_back( ctmc );
    
    // tell the nodes that we have a reference to it (avoids deletion)
    tree->incrementReferenceCount();
    ctmc->incrementReferenceCount();
    
    // tell the ctmc to use marginal likelihoods
    StochasticNode<PhyloCTMCSiteHomogeneous<characterType> > *char_stoch = (StochasticNode<PhyloCTMCSiteHomogeneous<characterType> >*) ctmc;
    PhyloCTMCSiteHomogeneous<characterType> *dist = (PhyloCTMCSiteHomogeneous<characterType>*) &char_stoch->getDistribution();
    dist->setUseMarginalLikelihoods( true );

}


/**
 * Copy constructor.
 */
template<class characterType>
AncestralStateMonitor<characterType>::AncestralStateMonitor( const AncestralStateMonitor &m) : VariableMonitor( m ),
    tree( m.tree ),
    ctmc( m.ctmc ),
    stochastic_nodes_only( m.stochastic_nodes_only )
{
    
}


/**
 * Destructor.
 */
template<class characterType>
AncestralStateMonitor<characterType>::~AncestralStateMonitor()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself 
 */
template<class characterType>
AncestralStateMonitor<characterType>* AncestralStateMonitor<characterType>::clone(void) const 
{
    
    return new AncestralStateMonitor<characterType>(*this);
}


/** 
 * Monitor value at given generation.
 *
 * \param[in]   gen    The current generation.
 */
template<class characterType>
void AncestralStateMonitor<characterType>::monitorVariables(unsigned long gen)
{
    auto& separator = to<SeparatorFormat>(format)->separator;
    
    // convert 'ctmc' from a DagNode to a StochasticNode so that we can call ctmc->getDistribution()
    StochasticNode<PhyloCTMCSiteHomogeneous<characterType> > *char_stoch = (StochasticNode<PhyloCTMCSiteHomogeneous<characterType> >*) ctmc;
		
    // now  get the TypedDistribution and cast it as an AbstractSiteHomogeneousMixtureCharEvoModel distribution
    PhyloCTMCSiteHomogeneous<characterType> *dist = (PhyloCTMCSiteHomogeneous<characterType>*) &char_stoch->getDistribution();
		
    // update the marginal node likelihoods
    dist->updateMarginalNodeLikelihoods();
        
    std::vector<TopologyNode*> nodes = tree->getValue().getNodes();
		
    // loop through all tree nodes
    for (int i = 0; i < tree->getValue().getNumberOfNodes(); i++)
    {
			
        // we need to print values for all internal nodes
        // we assume that tip nodes always precede internal nodes
        TopologyNode* the_node = nodes[i];
			
        if ( the_node->isTip() == false )
        {
				
            // add a separator before every new element
            out_stream << separator;
				
            // for each node print
            // site1,site2,site3
				
            std::vector<characterType> ancestralStates = dist->drawAncestralStatesForNode( *the_node );
                
            bool has_cladogenic_states = false; // dist->hasCladogenicStates();
            std::vector<characterType> cladogenicStates;
            if ( has_cladogenic_states == true )
            {
                ; // cladogenicStates = dist->drawCladogenicStatesForNode( *the_node );
            }
				
            // print out ancestral states....
            for (int j = 0; j < ancestralStates.size(); j++)
            {
                out_stream << ancestralStates[j].getStringValue();
                if (j != ancestralStates.size()-1 && typeid(ancestralStates[j]) != typeid(DnaState) )
                {
                    out_stream << ",";
                }
            }
            
        }
        
    }
    
}


/** 
 * Print header for monitored values 
 */
template<class characterType>
void AncestralStateMonitor<characterType>::printFileHeader()
{
    auto& separator = to<SeparatorFormat>(format)->separator;
	
	std::vector<TopologyNode*> nodes = tree->getValue().getNodes();
	
	// iterate through all tree nodes and make header with node index
	for (int i = 0; i < tree->getValue().getNumberOfNodes(); i++)
    {
		TopologyNode* the_node = nodes[i];
		if ( the_node->isTip() == false )
        {
			// add a separator before every new element
			out_stream << separator;
			
			// print the node index + 1 to be consistent with Rev language one-based indexes

			out_stream << the_node->getIndex() + 1;
		}
    }
    
}



template<class characterType>
void AncestralStateMonitor<characterType>::swapNode(DagNode *oldN, DagNode* newN)
{
	
	if ( oldN == tree ) 
	{
		tree = static_cast< TypedDagNode<Tree> *>( newN );
	}
	else if ( oldN == ctmc )
	{
		ctmc = newN;
	}

    VariableMonitor::swapNode( oldN, newN );
	
}




#endif

