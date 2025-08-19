#ifndef RevBayes_development_branch_StochasticCharacterMappingMonitor_h
#define RevBayes_development_branch_StochasticCharacterMappingMonitor_h

#include "AbstractHomologousDiscreteCharacterData.h"
#include "StateDependentSpeciationExtinctionProcess.h"
#include "GeneralizedLineageHeterogeneousBirthDeathSamplingProcess.h"
#include "VariableMonitor.h"
#include "TypedDagNode.h"
#include "StochasticNode.h"

#include <fstream>
#include <string>
#include <vector>
#include <typeinfo>


namespace RevBayesCore {

    /**
     * @brief Declaration and implementation of the StochasticCharacterMappingMonitor class.
     *
     * @file
     * Declaration and implementation of the StochasticCharacterMappingMonitor class which
     * monitors samples of character histories drawn from the state-dependent birth death process
     * and PhyloCTMC and prints their value into a file.
     *
     */
    template<class characterType>
    class StochasticCharacterMappingMonitor : public VariableMonitor {

    public:

        // Constructors and Destructors
        StochasticCharacterMappingMonitor(StochasticNode<Tree>* ch, std::uint64_t g, const std::string &fname, bool is, bool sd, const std::string &del);
//        StochasticCharacterMappingMonitor(TypedDagNode<Tree> *t, StochasticNode<AbstractHomologousDiscreteCharacterData>* ch, std::uint64_t g, const std::string &fname, bool is, const std::string &del);
        StochasticCharacterMappingMonitor(StochasticNode<AbstractHomologousDiscreteCharacterData>* ch, std::uint64_t g, const std::string &fname, bool is, bool sd, const std::string &del, size_t idx);
        StochasticCharacterMappingMonitor(const StochasticCharacterMappingMonitor &m);
        virtual ~StochasticCharacterMappingMonitor(void);

        StochasticCharacterMappingMonitor*              clone(void) const;                                                  //!< Clone the object

        // Monitor functions
        void                                            monitorVariables(std::uint64_t gen);                                 //!< Monitor at generation gen
        void                                            printFileHeader(void);                                              //!< Print header

        // getters and setters
        void                                            swapNode(DagNode *oldN, DagNode *newN);

    private:

        // members
        TypedDagNode<Tree>*                             tree;
        StochasticNode<Tree>*                           cdbdp;                                                              //!< The character dependent birth death process we are monitoring
        StochasticNode<AbstractHomologousDiscreteCharacterData>*            ctmc;
        bool                                            include_simmaps;                                                    //!< Should we print out SIMMAP/phytools compatible character histories?
        bool                                            use_simmap_default;
        size_t                                          index;

    };

}


#include "AbstractPhyloCTMCSiteHomogeneous.h"
#include "DagNode.h"
#include "Model.h"
#include "Monitor.h"
#include "RbFileManager.h"

using namespace RevBayesCore;


/* Constructor for state dependent birth death process */
template<class characterType>
StochasticCharacterMappingMonitor<characterType>::StochasticCharacterMappingMonitor(StochasticNode<Tree>* ch, std::uint64_t g, const std::string &fname, bool is, bool sd, const std::string &del) : VariableMonitor(ch, g, fname, del, false, false, false),
    cdbdp( ch ),
    include_simmaps( is ),
    use_simmap_default( sd ),
    index(0)
{
    ctmc = NULL;

    // the cdbdp is both the tree and character evolution model
    addVariable( cdbdp );
    tree = static_cast< StochasticNode<Tree> *>( cdbdp );

}

//StochasticCharacterMappingMonitor<characterType>::StochasticCharacterMappingMonitor(TypedDagNode<Tree> *t, StochasticNode<AbstractHomologousDiscreteCharacterData>* ch, std::uint64_t g, const std::string &fname, bool is, const std::string &del) : Monitor(g),
/* Constructor for CTMC */
template<class characterType>
StochasticCharacterMappingMonitor<characterType>::StochasticCharacterMappingMonitor(StochasticNode<AbstractHomologousDiscreteCharacterData>* ch, std::uint64_t g, const std::string &fname, bool is, bool sd, const std::string &del, size_t idx) :
    VariableMonitor(ch, g, fname, del, false, false, false),
    ctmc( ch ),
    include_simmaps( is ),
    use_simmap_default( sd ),
    index(idx)
{
    cdbdp = NULL;

    AbstractPhyloCTMCSiteHomogeneous<characterType> *ctmc_dist = NULL;
    ctmc_dist = static_cast<AbstractPhyloCTMCSiteHomogeneous<characterType>* >( &ctmc->getDistribution() );
    tree = const_cast<TypedDagNode<Tree>* >( ctmc_dist->getTree() );

    addVariable( tree );
    addVariable( ctmc );
}


/**
 * Copy constructor.
 */
template<class characterType>
StochasticCharacterMappingMonitor<characterType>::StochasticCharacterMappingMonitor( const StochasticCharacterMappingMonitor &m) : VariableMonitor( m ),
    tree( m.tree ),
    cdbdp( m.cdbdp ),
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
StochasticCharacterMappingMonitor<characterType>::~StochasticCharacterMappingMonitor()
{

}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself
 */
template<class characterType>
StochasticCharacterMappingMonitor<characterType>* StochasticCharacterMappingMonitor<characterType>::clone(void) const
{

    return new StochasticCharacterMappingMonitor(*this);
}


/**
 * Monitor value at given generation.
 *
 * \param[in]   gen    The current generation.
 */
template<class characterType>
void StochasticCharacterMappingMonitor<characterType>::monitorVariables(std::uint64_t gen)
{
    auto& separator = to<SeparatorFormat>(format)->separator;

    size_t num_nodes;

    // get the distribution for the character
    StateDependentSpeciationExtinctionProcess *sse_process = NULL;
    GeneralizedLineageHeterogeneousBirthDeathSamplingProcess *glhbdsp_process = NULL;
    AbstractPhyloCTMCSiteHomogeneous<characterType> *ctmc_dist = NULL;
    if ( ctmc != NULL )
    {
        ctmc_dist = static_cast<AbstractPhyloCTMCSiteHomogeneous<characterType>* >( &ctmc->getDistribution() );
        num_nodes = tree->getValue().getNumberOfNodes();
    }
    else
    {
        sse_process = dynamic_cast<StateDependentSpeciationExtinctionProcess*>( &cdbdp->getDistribution() );
        if ( sse_process == NULL )
        {
        	glhbdsp_process = dynamic_cast<GeneralizedLineageHeterogeneousBirthDeathSamplingProcess*>( &cdbdp->getDistribution() );
        }
        num_nodes = tree->getValue().getNumberOfNodes();
    }
        
    std::vector<std::string> character_histories( num_nodes );
    
    // draw stochastic character map
    if ( ctmc != NULL )
    {
        ctmc_dist->drawStochasticCharacterMap( character_histories, index, use_simmap_default );
    }
    else if ( sse_process != NULL )
    {
        sse_process->drawStochasticCharacterMap( character_histories );
    }
    else
    {
    	glhbdsp_process->drawStochasticCharacterMap( character_histories );
    }

    // print to monitor file
    const std::vector<TopologyNode*>& nds = tree->getValue().getNodes();
    for (int i = 0; i < nds.size(); i++)
    {

        size_t node_index = nds[i]->getIndex();

        // add a separator before every new element
        out_stream << separator;

        // print out this branch's character history in the format
        // used by SIMMAP and phytools
        out_stream << character_histories[ node_index ];
        
    }

    if ( include_simmaps == true )
    {
        // print out the SIMMAP/phytools compatible newick string as the last column of the log file
        out_stream << separator;
        Tree t = Tree(tree->getValue());
        t.clearNodeParameters();
        t.addNodeParameter( "character_history", character_histories, false );
        out_stream << t.getSimmapNewickRepresentation();
    }

}


/**
 * Print header for monitored values
 */
template<class characterType>
void StochasticCharacterMappingMonitor<characterType>::printFileHeader()
{
    auto& separator = to<SeparatorFormat>(format)->separator;

    std::vector<TopologyNode*> nodes = tree->getValue().getNodes();

    // iterate through all tree nodes and make header with node index
    for (int i = 0; i < tree->getValue().getNumberOfNodes(); i++)
    {
        TopologyNode* nd = nodes[i];
        size_t node_index = nd->getIndex();

        out_stream << separator;
        out_stream << node_index + 1;

    }

    if ( include_simmaps == true )
    {
        out_stream << separator;
        out_stream << "simmap";
    }

}


template<class characterType>
void StochasticCharacterMappingMonitor<characterType>::swapNode(DagNode *oldN, DagNode* newN)
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

    VariableMonitor::swapNode( oldN, newN );

}

#endif
