//
//  CharacterHistorySummaryMonitor.h
//  revbayes-proj
//
//  Created by Michael Landis on 12/3/18.
//  Copyright © 2018 Michael Landis. All rights reserved.
//

#ifndef CharacterHistorySummaryMonitor_h
#define CharacterHistorySummaryMonitor_h

#include "AbstractCharacterData.h"
#include "TreeHistoryCtmc.h"
#include "BiogeographicCladoEvent.h"
#include "GeneralTreeHistoryCtmc.h"
#include "DagNode.h"
#include "Model.h"
#include "Monitor.h"
#include "BranchHistory.h"
#include "RbException.h"
#include "RbFileManager.h"
#include "StochasticNode.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"

#include <fstream>
#include <cstdint>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace RevBayesCore {
    
    template<class charType>
    class CharacterHistorySummaryMonitor : public Monitor {
        
    public:
        // Constructors and Destructors
        CharacterHistorySummaryMonitor(StochasticNode<AbstractHomologousDiscreteCharacterData>* s, TypedDagNode<Tree> *t, std::uint64_t g, const path &fname, const std::string &del, bool pp=true, bool l=true, bool pr=true, bool ap=false, bool sm=true, bool sne=false, bool ste=true);
        CharacterHistorySummaryMonitor(const CharacterHistorySummaryMonitor& f);
        
        // basic methods
        CharacterHistorySummaryMonitor*          clone(void) const;                                                  //!< Clone the object
        
        // Monitor functions
        void                                monitor(std::uint64_t gen);                                         //!< Monitor at generation gen
        void                                swapNode(DagNode *oldN, DagNode *newN);
        
        // FileMonitor functions
        void                                closeStream(void);                                                  //!< Close stream after finish writing
        void                                openStream(bool reopen);                                            //!< Open the stream for writing
        void                                printHeader(void);                                                  //!< Print header
        
    private:
        std::string                         buildChainSummaryString(std::uint64_t gen);
        std::string                         buildIterationSummaryTable(std::uint64_t gen);
        std::string                         buildIterationSummaryTable(std::uint64_t gen, TopologyNode* n);
        bool                                isCladogenetic( TopologyNode* n );
        
        // the stream to print
        std::fstream                        outStream;
        
        // parameters
        StochasticNode<AbstractHomologousDiscreteCharacterData>* variable;
        TypedDagNode<Tree>*                 tree;
        std::set<DagNode *>                 nodeVariables;
        path                                filename;
        std::string                         separator;
        bool                                posterior;
        bool                                prior;
        bool                                likelihood;
        bool                                append;
        bool                                showMetadata;
        bool                                showNumEvents;
        bool                                showTreeEvents;
        size_t                              num_states;
        
    };
    
}

/* Constructor */
template<class charType>
RevBayesCore::CharacterHistorySummaryMonitor<charType>::CharacterHistorySummaryMonitor(StochasticNode<AbstractHomologousDiscreteCharacterData>* s, TypedDagNode<Tree>* t, std::uint64_t g, const RevBayesCore::path &fname, const std::string &del, bool pp, bool l, bool pr, bool ap, bool sm, bool sne, bool ste) :
Monitor(g,t),
outStream(),
variable(s),
tree( t ),
filename( fname ),
separator( del ),
posterior( pp ),
prior( pr ),
likelihood( l ),
append(ap),
showMetadata(sm),
showNumEvents(sne),
showTreeEvents(ste),
num_states(0)
{
    nodes.push_back(s);
    //    nodes.push_back(t);
    s->incrementReferenceCount();
    
    num_states = static_cast<const DiscreteCharacterState&>(s->getValue().getCharacter(0,0)).getNumberOfStates();
}

template<class charType>
RevBayesCore::CharacterHistorySummaryMonitor<charType>::CharacterHistorySummaryMonitor(const CharacterHistorySummaryMonitor &m) :
Monitor( m ),
outStream( ),
variable(m.variable),
tree( m.tree ),
nodeVariables( m.nodeVariables ),
showMetadata(m.showMetadata),
showNumEvents(m.showNumEvents),
showTreeEvents(m.showTreeEvents),
num_states(m.num_states)
{
    filename    = m.filename;
    separator   = m.separator;
    prior       = m.prior;
    posterior   = m.posterior;
    likelihood  = m.likelihood;
    append      = m.append;
}


/* Clone the object */
template<class charType>
RevBayesCore::CharacterHistorySummaryMonitor<charType>* RevBayesCore::CharacterHistorySummaryMonitor<charType>::clone(void) const {
    
    return new CharacterHistorySummaryMonitor<charType>(*this);
}

template<class charType>
void RevBayesCore::CharacterHistorySummaryMonitor<charType>::closeStream() {
    outStream.close();
}

template<class charType>
std::string RevBayesCore::CharacterHistorySummaryMonitor<charType>::buildChainSummaryString( std::uint64_t gen ) {
    // print the iteration number first
    
    std::stringstream ss;
    ss << gen;
    
    if ( posterior ) {
        // add a separator before every new element
        ss << separator;
        
        const std::vector<DagNode*> &n = model->getDagNodes();
        double pp = 0.0;
        for (std::vector<DagNode*>::const_iterator it = n.begin(); it != n.end(); ++it) {
            pp += (*it)->getLnProbability();
        }
        ss << pp;
    }
    
    if ( likelihood ) {
        // add a separator before every new element
        ss << separator;
        
        const std::vector<DagNode*> &n = model->getDagNodes();
        double pp = 0.0;
        for (std::vector<DagNode*>::const_iterator it = n.begin(); it != n.end(); ++it) {
            if ( (*it)->isClamped() ) {
                pp += (*it)->getLnProbability();
            }
        }
        ss << pp;
    }
    
    if ( prior ) {
        // add a separator before every new element
        ss << separator;
        
        const std::vector<DagNode*> &n = model->getDagNodes();
        double pp = 0.0;
        for (std::vector<DagNode*>::const_iterator it = n.begin(); it != n.end(); ++it) {
            if ( !(*it)->isClamped() ) {
                pp += (*it)->getLnProbability();
            }
        }
        ss << pp;
    }
    
    return ss.str();
}

template<class charType>
std::string RevBayesCore::CharacterHistorySummaryMonitor<charType>::buildIterationSummaryTable( std::uint64_t gen ) {
    //tree->getValue().getRoot().setNewickNeedsRefreshing(true);
    std::string newick = buildIterationSummaryTable( gen, &tree->getValue().getRoot() );
    return newick;
}

/* Build newick string */
template<class charType>
std::string RevBayesCore::CharacterHistorySummaryMonitor<charType>::buildIterationSummaryTable( std::uint64_t gen, TopologyNode* n ) {
    // create the newick string
    std::stringstream o;
    
    // get character history for tree
    TreeHistoryCtmc<charType>* p = static_cast< TreeHistoryCtmc<charType>* >(&variable->getDistribution());
    
    // get branch history for node
    const BranchHistory& bh = p->getHistory(*n);
    
    // get node index
    std::stringstream node_index_ss;
    node_index_ss << n->getIndex();
    
    // get parent index
    std::stringstream parent_index_ss;
    if ( n->isRoot() ) {
        parent_index_ss << "NA";
    }
    else {
        parent_index_ss << n->getParent().getIndex();
    }
    
    // get children indices
    std::vector<TopologyNode*> children = n->getChildren();
    std::stringstream child1_index_ss;
    std::stringstream child2_index_ss;
    if ( n->isTip() ) {
        child1_index_ss << "NA";
        child2_index_ss << "NA";
    }
    else {
        child1_index_ss << children[0]->getIndex();
        child2_index_ss << children[1]->getIndex();
    }
    
    // get branch lengths
    double branch_end_time = n->getAge();
    double branch_start_time = branch_end_time + n->getBranchLength();
    
    std::stringstream branch_end_time_ss;
    branch_end_time_ss << branch_end_time;
    std::stringstream branch_start_time_ss;
    branch_start_time_ss << branch_start_time;
    
    // get prefix string per row
    std::stringstream ss_prefix;
    ss_prefix << buildChainSummaryString( gen ) << separator;
    ss_prefix << node_index_ss.str() << separator;
    ss_prefix << branch_start_time_ss.str() << separator;
    ss_prefix << branch_end_time_ss.str();
    
    std::stringstream ss_suffix;
    ss_suffix << parent_index_ss.str() << separator;
    ss_suffix << child1_index_ss.str() << separator;
    ss_suffix << child2_index_ss.str();
    
    // get event list
    std::multiset<CharacterEvent*,CharacterEventCompare> h = bh.getHistory();
    
    // what events to expect?
    bool has_cladogenetic = isCladogenetic( n );
    bool has_anagenetic = ( h.size() > 0 );
    
    std::vector<CharacterEvent*> start_states = bh.getParentCharacters();
    std::vector<CharacterEvent*> end_states = bh.getChildCharacters();
    
    std::stringstream start_state_ss;
    std::stringstream end_state_ss;
    
    for (size_t i = 0; i < start_states.size(); i++) {
        start_state_ss << static_cast<CharacterEventDiscrete*>( start_states[i] )->getState();
        end_state_ss << static_cast<CharacterEventDiscrete*>( end_states[i] )->getState();
    }
    
    // no_change along or at end of branch
    if (!has_cladogenetic && !has_anagenetic) {
        o << ss_prefix.str() << separator;
        o << start_state_ss.str() << separator;
        o << end_state_ss.str() << separator;
        o << "NA" << separator; // transition_time
        o << "no_change" << separator; // transition_type
        o << ss_suffix.str() << std::endl;
    }
    else {
        
        // anagenetic events along branch
        if (has_anagenetic) {
            std::vector<CharacterEvent*> next_states = start_states;
            std::multiset<CharacterEvent*, CharacterEventCompare>::reverse_iterator it;
            for (it = h.rbegin(); it != h.rend(); it++)
            {
                std::stringstream curr_state_ss;
                std::stringstream next_state_ss;
                
                CharacterEventDiscrete* evt = static_cast<CharacterEventDiscrete*>( *it );
                next_states[ evt->getSiteIndex() ] = evt;
                for (size_t i = 0; i < start_states.size(); i++) {
                    curr_state_ss << static_cast<CharacterEventDiscrete*>( start_states[i] )->getState();
                    next_state_ss << static_cast<CharacterEventDiscrete*>( next_states[i] )->getState();
                }
                
                o << ss_prefix.str() << separator;
                o << curr_state_ss.str() << separator;
                o << next_state_ss.str() << separator;
                o << evt->getAge() << separator;
                o << "anagenetic" << separator;
                o << ss_suffix.str() << std::endl;
                
                start_states = next_states;
            }
        }
        
        // cladogenetic event at end of branch
        if (has_cladogenetic && !n->isTip()) {
            
            for (size_t i = 0; i < n->getNumberOfChildren(); i++) {
                
                // get child states
                TopologyNode* child = &n->getChild(i);
                const BranchHistory& bh_child = p->getHistory( *child );
                
                const std::vector<CharacterEvent*>& child_states = bh_child.getParentCharacters();
                std::stringstream child_state_ss;
                
                // make child state string
                for (size_t i = 0; i < start_states.size(); i++) {
                    child_state_ss << static_cast<CharacterEventDiscrete*>( child_states[i] )->getState();
                }
                
                // make child cladogenesis row string
                o << ss_prefix.str() << separator;
                o << end_state_ss.str() << separator;
                o << child_state_ss.str() << separator;
                o << branch_end_time_ss.str() << separator;
                o << "cladogenetic" << separator;
                o << ss_suffix.str() << std::endl;
            }
        }
    }
    
    // recurse
    for (size_t i = 0; i < n->getNumberOfChildren(); i++) {
        o << buildIterationSummaryTable( gen, &n->getChild(i) );
    }
    
    return o.str();
}

template<class charType>
bool RevBayesCore::CharacterHistorySummaryMonitor<charType>::isCladogenetic( TopologyNode* n ) {
    
    // get character history for tree
    TreeHistoryCtmc<charType>* p = static_cast< TreeHistoryCtmc<charType>* >(&variable->getDistribution());
    
    bool is_cladogenetic = false;
    if ( !n->isTip() ) {
        const TopologyNode& n_child0 = n->getChild(0);
        const TopologyNode& n_child1 = n->getChild(1);
        const std::vector<CharacterEvent*>& end_states          = p->getHistory(*n).getChildCharacters();
        const std::vector<CharacterEvent*>& start_states_child0 = p->getHistory(n_child0).getParentCharacters();
        const std::vector<CharacterEvent*>& start_states_child1 = p->getHistory(n_child1).getParentCharacters();
        
        for (size_t i = 0; i < end_states.size(); i++) {
            size_t s_node = static_cast<CharacterEventDiscrete*>( end_states[i] )->getState();
            size_t s_child0 = static_cast<CharacterEventDiscrete*>( start_states_child0[i] )->getState();
            size_t s_child1 = static_cast<CharacterEventDiscrete*>( start_states_child1[i] )->getState();
            
            if (s_node != s_child0 || s_node != s_child1) {
                is_cladogenetic = true;
                break;
            }
        }
    }
    return is_cladogenetic;
}

/** Monitor value at generation gen */
template<class charType>
void RevBayesCore::CharacterHistorySummaryMonitor<charType>::monitor(std::uint64_t gen) {
    
    // get the printing frequency
    std::uint64_t samplingFrequency = printgen;
    
    if (gen % samplingFrequency == 0) {

        outStream << buildIterationSummaryTable(gen);
        
    }
}


/** open the file stream for printing */
template<class charType>
void RevBayesCore::CharacterHistorySummaryMonitor<charType>::openStream(bool reopen)
{
    RevBayesCore::createDirectoryForFile(filename);
    
    // open the stream to the file
    if (append)
        outStream.open( filename.string(), std::fstream::out | std::fstream::app);
    else
        outStream.open( filename.string(), std::fstream::out);
}

/** Print header for monitored values */
template<class charType>
void RevBayesCore::CharacterHistorySummaryMonitor<charType>::printHeader()
{
    
    // print one column for the iteration number
    outStream << "iteration";
    
    if ( posterior )
    {
        // add a separator before every new element
        outStream << separator;
        outStream << "posterior";
    }
    
    if ( likelihood ) {
        // add a separator before every new element
        outStream << separator;
        outStream << "likelihood";
    }
    
    if ( prior ) {
        // add a separator before every new element
        outStream << separator;
        outStream << "prior";
    }
    
    outStream << separator << "node_index";
    outStream << separator << "branch_start_time";
    outStream << separator << "branch_end_time";
    outStream << separator << "start_state";
    outStream << separator << "end_state";
    outStream << separator << "transition_time";
    outStream << separator << "transition_type";
    outStream << separator << "parent_index";
    outStream << separator << "child1_index";
    outStream << separator << "child2_index";
    outStream << std::endl;

}

template<class charType>
void RevBayesCore::CharacterHistorySummaryMonitor<charType>::swapNode(DagNode *oldN, DagNode *newN) {
    
    if ( oldN == tree )
    {
        tree = static_cast< TypedDagNode< Tree > *>( newN );
    }
    else if (oldN == variable)
    {
        variable = static_cast<StochasticNode<AbstractHomologousDiscreteCharacterData>* >(newN);
    }
    
    // delegate to base class
    Monitor::swapNode(oldN, newN);
}

#endif /* CharacterHistorySummaryMonitor_h */

