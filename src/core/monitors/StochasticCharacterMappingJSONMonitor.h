#ifndef RevBayes_development_branch_StochasticCharacterMappingJSONMonitor_h
#define RevBayes_development_branch_StochasticCharacterMappingJSONMonitor_h

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
     * @brief Declaration and implementation of the StochasticCharacterMappingJSONMonitor class.
     *
     * @file
     * Declaration and implementation of the StochasticCharacterMappingJSONMonitor class which
     * monitors samples of character histories drawn from the state-dependent birth death process
     * and PhyloCTMC and prints their value into a file.
     *
     */
    template<class characterType>
    class StochasticCharacterMappingJSONMonitor : public VariableMonitor {

    public:

        // Constructors and Destructors
        StochasticCharacterMappingJSONMonitor(StochasticNode<Tree>* ch, unsigned long g, const std::string &fname, bool is, bool sd, const std::string &del);
        StochasticCharacterMappingJSONMonitor(const StochasticCharacterMappingJSONMonitor &m);
        virtual ~StochasticCharacterMappingJSONMonitor(void);

        StochasticCharacterMappingJSONMonitor*              clone(void) const;                                                  //!< Clone the object

        // Monitor functions
        void                                            monitorVariables(unsigned long gen);                                 //!< Monitor at generation gen
        void                                            printFileHeader(void);                                              //!< Print header

        // getters and setters
        void                                            swapNode(DagNode *oldN, DagNode *newN);

    private:

        // members
        TypedDagNode<Tree>*                             tree;
        StochasticNode<Tree>*                           glhbdsp;                                                              //!< The character dependent birth death process we are monitoring
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
StochasticCharacterMappingJSONMonitor<characterType>::StochasticCharacterMappingJSONMonitor(StochasticNode<Tree>* ch, unsigned long g, const std::string &fname, bool is, bool sd, const std::string &del) : VariableMonitor(ch, g, fname, del, false, false, false),
    glhbdsp( ch ),
    index(0)
{

    // the cdbdp is both the tree and character evolution model
    addVariable( glhbdsp );
    tree = static_cast< StochasticNode<Tree> *>( glhbdsp );

}

/**
 * Copy constructor.
 */
template<class characterType>
StochasticCharacterMappingJSONMonitor<characterType>::StochasticCharacterMappingJSONMonitor( const StochasticCharacterMappingJSONMonitor &m) : VariableMonitor( m ),
    tree( m.tree ),
    glhbdsp( m.glhbdsp ),
    index( m.index )
{

}


/**
 * Destructor.
 */
template<class characterType>
StochasticCharacterMappingJSONMonitor<characterType>::~StochasticCharacterMappingJSONMonitor()
{

}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself
 */
template<class characterType>
StochasticCharacterMappingJSONMonitor<characterType>* StochasticCharacterMappingJSONMonitor<characterType>::clone(void) const
{

    return new StochasticCharacterMappingJSONMonitor(*this);
}


/**
 * Monitor value at given generation.
 *
 * \param[in]   gen    The current generation.
 */
template<class characterType>
void StochasticCharacterMappingJSONMonitor<characterType>::monitorVariables(unsigned long gen)
{

    auto& separator = to<SeparatorFormat>(format)->separator;
    size_t num_nodes;

    // get the distribution for the character
    GeneralizedLineageHeterogeneousBirthDeathSamplingProcess *glhbdsp_process = dynamic_cast<GeneralizedLineageHeterogeneousBirthDeathSamplingProcess*>( &glhbdsp->getDistribution() );

    if ( glhbdsp_process == NULL ) {
    	RbException("Error casting as GLHBDSP object");
    }


    // get the JSON-formatted stochastic map
    std::string json_string;
    glhbdsp_process->drawStochasticCharacterMapJSON(json_string);

    // append
    out_stream << separator;
    out_stream << json_string;

}


/**
 * Print header for monitored values
 */
template<class characterType>
void StochasticCharacterMappingJSONMonitor<characterType>::printFileHeader()
{

    auto& separator = to<SeparatorFormat>(format)->separator;
	out_stream << separator;
	out_stream << "stoch_map";

//	std::vector<TopologyNode*> nodes = tree->getValue().getNodes();
//
//    // iterate through all tree nodes and make header with node index
//    for (int i = 0; i < tree->getValue().getNumberOfNodes(); i++)
//    {
//        TopologyNode* nd = nodes[i];
//        size_t node_index = nd->getIndex();
//
//        out_stream << separator;
//        out_stream << node_index + 1;
//
//    }
//
//    if ( include_simmaps == true )
//    {
//        out_stream << separator;
//        out_stream << "simmap";
//    }

}


template<class characterType>
void StochasticCharacterMappingJSONMonitor<characterType>::swapNode(DagNode *oldN, DagNode* newN)
{

    if ( oldN == glhbdsp )
    {
        glhbdsp = static_cast< StochasticNode<Tree> *>( newN );
        tree = static_cast< StochasticNode<Tree> *>( newN );
    }

    VariableMonitor::swapNode( oldN, newN );

}

#endif
