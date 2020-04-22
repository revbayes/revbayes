//
//  BirthDeathStochasticCharacterMappingMonitor.hpp
//  revbayes-proj
//
//  Created by Michael Landis on 4/5/20.
//  Copyright © 2020 Michael Landis. All rights reserved.
//

#ifndef BirthDeathStochasticCharacterMappingMonitor_h
#define BirthDeathStochasticCharacterMappingMonitor_h

#include "AbstractHomologousDiscreteCharacterData.h"
#include "StateDependentSpeciationExtinctionProcess.h"
#include "VariableMonitor.h"
#include "TypedDagNode.h"
#include "StochasticNode.h"

#include <fstream>
#include <string>
#include <vector>
#include <typeinfo>


namespace RevBayesCore {
    
    /**
     * @brief Declaration and implementation of the BirthDeathStochasticCharacterMappingMonitor class.
     *
     * @file
     * Declaration and implementation of the BirthDeathStochasticCharacterMappingMonitor class which
     * monitors samples of character histories drawn from the state-dependent birth death process
     * and PhyloCTMC and prints their value into a file.
     *
     */
    template<class characterType>
    class BirthDeathStochasticCharacterMappingMonitor : public VariableMonitor {
        
    public:
        
        // Constructors and Destructors
        // BirthDeathStochasticCharacterMappingMonitor(StochasticNode<Tree>* ch, unsigned long g, const std::string &fname, bool is, bool sd, const std::string &del);
        // BirthDeathStochasticCharacterMappingMonitor(TypedDagNode<Tree> *t, StochasticNode<AbstractHomologousDiscreteCharacterData>* ch, unsigned long g, const std::string &fname, bool is, const std::string &del);
        BirthDeathStochasticCharacterMappingMonitor(StochasticNode<AbstractHomologousDiscreteCharacterData>* ch, StochasticNode<Tree>* t, unsigned long g, const std::string &fname, bool is, bool sd, const std::string &del, size_t idx);
        BirthDeathStochasticCharacterMappingMonitor(const BirthDeathStochasticCharacterMappingMonitor &m);
        virtual ~BirthDeathStochasticCharacterMappingMonitor(void);
        
        BirthDeathStochasticCharacterMappingMonitor*              clone(void) const;                                                  //!< Clone the object
        
        // Monitor functions
        void                                            monitorVariables(unsigned long gen);                                 //!< Monitor at generation gen
        void                                            printFileHeader(void);                                              //!< Print header
        
        // getters and setters
        void                                            swapNode(DagNode *oldN, DagNode *newN);
        
    private:
        
        // members
        
//        StochasticNode<Tree>*                           cdbdp;                                                              //!< The character dependent birth death process we are monitoring
        StochasticNode<Tree>*                           tree;
        StochasticNode<AbstractHomologousDiscreteCharacterData>*            ctmc;
        bool                                            include_simmaps;                                                    //!< Should we print out SIMMAP/phytools compatible character histories?
        bool                                            use_simmap_default;                                                 //!< Should we generate stochastic mappings for extinct/unsampled subclades?
        size_t                                          index;
        
    };
    
}


#include "AbstractPhyloCTMCSiteHomogeneous.h"
#include "ConstantRateBirthDeathProcess.h"
#include "DagNode.h"
#include "Model.h"
#include "Monitor.h"
#include "RbFileManager.h"

using namespace RevBayesCore;

/* Constructor for CTMC */
template<class characterType>
BirthDeathStochasticCharacterMappingMonitor<characterType>::BirthDeathStochasticCharacterMappingMonitor(StochasticNode<AbstractHomologousDiscreteCharacterData>* ch, StochasticNode<Tree>* t, unsigned long g, const std::string &fname, bool is, bool sd, const std::string &del, size_t idx) :
VariableMonitor(ch, g, fname, del, false, false, false),
ctmc( ch ),
tree( t ),
include_simmaps( is ),
use_simmap_default( sd ),
index(idx)
{
//    cdbdp = NULL;
    
    AbstractPhyloCTMCSiteHomogeneous<characterType> *ctmc_dist = NULL;
    AbstractBirthDeathProcess *tree_dist = NULL;
    ctmc_dist = static_cast<AbstractPhyloCTMCSiteHomogeneous<characterType>* >( &ctmc->getDistribution() );
    tree_dist = static_cast<AbstractBirthDeathProcess* >( &t->getDistribution() );
    
    addVariable( tree );
    addVariable( ctmc );
}


/**
 * Copy constructor.
 */
template<class characterType>
BirthDeathStochasticCharacterMappingMonitor<characterType>::BirthDeathStochasticCharacterMappingMonitor( const BirthDeathStochasticCharacterMappingMonitor &m) : VariableMonitor( m ),
tree( m.tree ),
//cdbdp( m.cdbdp ),
ctmc( m.ctmc ),
include_simmaps( m.include_simmaps ),
use_simmap_default( m.use_simmap_default ),
index( m.index )
{
    
}


/**
 * Destructor.
 */
template<class characterType>
BirthDeathStochasticCharacterMappingMonitor<characterType>::~BirthDeathStochasticCharacterMappingMonitor()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself
 */
template<class characterType>
BirthDeathStochasticCharacterMappingMonitor<characterType>* BirthDeathStochasticCharacterMappingMonitor<characterType>::clone(void) const
{
    
    return new BirthDeathStochasticCharacterMappingMonitor(*this);
}


/**
 * Monitor value at given generation.
 *
 * \param[in]   gen    The current generation.
 */
template<class characterType>
void BirthDeathStochasticCharacterMappingMonitor<characterType>::monitorVariables(unsigned long gen)
{
    
    
    // create empty vector for simmap strings
    std::vector<std::string> character_histories;
    
    Tree tree_restore = Tree(tree->getValue());
    AbstractHomologousDiscreteCharacterData* data_restore = ctmc->getValue().clone();
    AbstractHomologousDiscreteCharacterData* data_complete = ctmc->getValue().clone();

    // get the distribution for the character
    AbstractPhyloCTMCSiteHomogeneous<characterType> *ctmc_dist = NULL;
    
    
    if ( ctmc != NULL && tree != NULL )
    {
        ctmc_dist = static_cast<AbstractPhyloCTMCSiteHomogeneous<characterType>* >( &ctmc->getDistribution() );
        ConstantRateBirthDeathProcess *bdp_dist = static_cast<ConstantRateBirthDeathProcess*>( &tree->getDistribution() );
        double lnProb_before = ctmc_dist->computeLnProbability();
        
        // draw complete tree
        bdp_dist->simulateHiddenClades();
        
        // add missing taxa
        std::vector<TopologyNode*> nodes = bdp_dist->getValue().getNodes();
 
        for (size_t i = 0; i < nodes.size(); i++) {
            TopologyNode* n = nodes[i];
            if (n->isTip()) {
                std::string taxon_str = n->getName();
//                std::cout << taxon_str << "\n";
                data_complete->addMissingTaxon( taxon_str );
            }
        }

        // get number of nodes
        size_t num_nodes = bdp_dist->getValue().getNumberOfNodes();
        
        // allocate vector for character histories
        character_histories = std::vector<std::string>( num_nodes );
        
        // assign complete matrix to CTMC
        ctmc_dist->setValue( data_complete, true );
 
        // need to update tree value within ctmc
        double lnProb_after = ctmc_dist->computeLnProbability();
//        std::cout << lnProb_before << " " << lnProb_after << "\n";
        
        // draw stochastic mappings for CTMC
        ctmc_dist->drawStochasticCharacterMap( character_histories, index, use_simmap_default );
        
    }

    // print to monitor file
    std::vector<TopologyNode*> nds = tree->getValue().getNodes();
    
    
//    for (int i = 0; i < nds.size(); i++)
//    {
//        
//        size_t node_index = nds[i]->getIndex();
//        
//        // add a separator before every new element
//        out_stream << separator;
//        
//        // print out this branch's character history in the format
//        // used by SIMMAP and phytools
//        out_stream << character_histories[ node_index ];
//        
//    }
    
    
    if ( include_simmaps == true )
    {
        // print out the SIMMAP/phytools compatible newick string as the last column of the log file
        out_stream << separator;
        Tree t = Tree(tree->getValue());
        out_stream << t.getNewickRepresentation();
        out_stream << separator;
        t.clearNodeParameters();
        t.addNodeParameter( "character_history", character_histories, false );
        out_stream << t.getSimmapNewickRepresentation();
    }
    
    
    // restore distribution values
    if ( ctmc != NULL && tree != NULL )
    {
//        AbstractBirthDeathProcess *bdp_dist = static_cast<AbstractBirthDeathProcess*>( &tree->getDistribution() );
        ConstantRateBirthDeathProcess *bdp_dist = static_cast<ConstantRateBirthDeathProcess*>( &tree->getDistribution() );
//        bdp_dist->getValue().initFromString( tree_restore.getNewickRepresentation() );
        
//        bdp_dist->getValue().setRoot( &bdp_dist->getValue().getRoot(), true );
//        TopologyNode* root = &bdp_dist->getValue().getRoot();
//        root->setAge( root->getAge() );
        bdp_dist->getValue().setRoot( &tree->getValue().getRoot(), true );
        bdp_dist->getValue().getNewickRepresentation();
        bdp_dist->getValue().clearTaxonBitSetMap();
        bdp_dist->getValue().getTaxonBitSetMap();
        
        bdp_dist->getValue().pruneTaxaWithoutSampledDescendants();
//        std::cout << bdp_dist->getValue() << "\n";
        // update nodes/taxa?
        
//        Tree tt = bdp_dist->getValue();
//        std::cout << tt << "\n";
//        //        }
//        tree->getValue().setRoot( &tree->getValue().getRoot(), true );
//        tree->getValue().getNewickRepresentation();
//        tree->getValue().clearTaxonBitSetMap();
//        tree->getValue().getTaxonBitSetMap();
        
        ctmc->setValue( data_restore, true );
    }

}


/**
 * Print header for monitored values
 */
template<class characterType>
void BirthDeathStochasticCharacterMappingMonitor<characterType>::printFileHeader()
{
    std::vector<TopologyNode*> nodes = tree->getValue().getNodes();
    
    /*
    // iterate through all tree nodes and make header with node index
    for (int i = 0; i < tree->getValue().getNumberOfNodes(); i++)
    {
        TopologyNode* nd = nodes[i];
        size_t node_index = nd->getIndex();
        
        out_stream << separator;
        out_stream << node_index + 1;
        
    }
    */
    
    out_stream << separator;
    out_stream << "tree";
    if ( include_simmaps == true )
    {
        out_stream << separator;
        out_stream << "simmap";
    }
    
}


template<class characterType>
void BirthDeathStochasticCharacterMappingMonitor<characterType>::swapNode(DagNode *oldN, DagNode* newN)
{
    
    if ( oldN == tree )
    {
        tree = static_cast< StochasticNode<Tree> *>( newN );
    }
    else if ( oldN == ctmc )
    {
        ctmc = static_cast< StochasticNode<AbstractHomologousDiscreteCharacterData> *>( newN );
    }
//    else if ( oldN == cdbdp )
//    {
//        cdbdp = static_cast< StochasticNode<Tree> *>( newN );
//        tree = static_cast< StochasticNode<Tree> *>( newN );
//    }
    
    VariableMonitor::swapNode( oldN, newN );
    
}


#endif /* BirthDeathStochasticCharacterMappingMonitor_hpp */
