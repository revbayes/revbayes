#include <cmath>
#include <typeinfo>
#include <algorithm>
#include <cstddef>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "DagNode.h"
#include "HomologousDiscreteCharacterData.h"
#include "PosteriorPredictiveSimulation.h"
#include "StateDependentSpeciationExtinctionProcess.h"
#include "RbFileManager.h"
#include "AbstractHomologousDiscreteCharacterData.h"
#include "Cloneable.h"
#include "DiscreteTaxonData.h"
#include "Model.h"
#include "NaturalNumbersState.h"
#include "Parallelizable.h"
#include "RbException.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "StringUtilities.h"
#include "Taxon.h"
#include "Trace.h"
#include "Tree.h"
#include "TreeDiscreteCharacterData.h"

using namespace RevBayesCore;

PosteriorPredictiveSimulation::PosteriorPredictiveSimulation( const Model &m, const std::string &dir, const RbVector<ModelTrace> &t) : Cloneable(), Parallelizable(),
    model( m ),
    directory( dir ),
    traces( t ),
    condition_on_tips( false )
{
    
}


PosteriorPredictiveSimulation::PosteriorPredictiveSimulation( const Model &m, const std::string &dir, const RbVector<ModelTrace> &t, const std::vector<RevBayesCore::AncestralStateTrace> &ast) : Cloneable(), Parallelizable(),
    model( m ),
    directory( dir ),
    traces( t ),
    ancestral_state_traces( ast ),
    condition_on_tips( true )
{
    
}



PosteriorPredictiveSimulation::~PosteriorPredictiveSimulation()
{
    
}



/** Clone function */
PosteriorPredictiveSimulation* PosteriorPredictiveSimulation::clone() const
{
    
    return new PosteriorPredictiveSimulation(*this);
}


void RevBayesCore::PosteriorPredictiveSimulation::run( int thinning )
{
    
    // some general constant variables
    // this is where we need to implement the posterior predictive simulation
    
    size_t n_samples = traces[0].size();
    size_t n_traces = traces.size();
        
    // build a map for the ancestral state trace labels -> tip indices
    std::map<std::string, size_t> ancestral_state_traces_lookup;
    if (condition_on_tips == true)
    {
        for (size_t z = 0; z < ancestral_state_traces.size(); z++)
        {
            ancestral_state_traces_lookup[ ancestral_state_traces[z].getParameterName() ] = z;
        }
    }

    
    std::vector<DagNode*> nodes = model.getDagNodes();

    size_t sim_pid_start = size_t(floor( (double(pid) / num_processes * n_samples ) ) );
    size_t sim_pid_end   = std::max( int(sim_pid_start), int(floor( (double(pid+1) / num_processes * n_samples ) ) - 1) );
    
    size_t index_sample = sim_pid_start;
    while ( index_sample % thinning > 0 ) ++index_sample;
    size_t current_pp_sim = size_t( floor( index_sample / thinning ) );
    
    
    for ( ; index_sample <= sim_pid_end; ++current_pp_sim, index_sample += thinning)
    {
        
        
        // create a new directory name for this simulation
        path sim_directory_name = directory / ("posterior_predictive_sim_" + std::to_string(current_pp_sim + 1));
        
        // now for the numerical parameters
        for ( size_t j=0; j<n_traces; ++j )
        {
            std::string parameter_name = traces[j].getParameterName();
            
            // iterate over all DAG nodes (variables)
            for ( std::vector<DagNode*>::iterator it = nodes.begin(); it!=nodes.end(); ++it )
            {
                DagNode *the_node = *it;
                
                if ( the_node->getName() == parameter_name )
                {
                    // set the value for the variable with the i-th sample
                    the_node->setValueFromString( traces[j].objectAt( index_sample ) );
                }
            
            }
        
        }
        
        // next we need to simulate the data and store it
        // iterate over all DAG nodes (variables)
        for ( std::vector<DagNode*>::iterator it = nodes.begin(); it!=nodes.end(); ++it )
        {
            DagNode *the_node = *it;
            
            if ( the_node->isClamped() == true )
            {
                the_node->setIgnoreRedraw(false);
                
                // check if the PP simulation must condition on sampled tip states
                auto& ptr = the_node->getDistribution();
                if (condition_on_tips == true && typeid( &ptr ) == typeid(StateDependentSpeciationExtinctionProcess))
                {
                    // set the tip states to the values sampled during this iteration
                    AncestralStateTrace* tip_state_trace;
                    StateDependentSpeciationExtinctionProcess* sse = static_cast<StateDependentSpeciationExtinctionProcess*>( &the_node->getDistribution() );
                    std::vector<std::string> tips = sse->getValue().getTipNames();
                    size_t num_states = static_cast<TreeDiscreteCharacterData*>( &sse->getValue() )->getCharacterData().getNumberOfStates();
                    HomologousDiscreteCharacterData<NaturalNumbersState> *tip_data = new HomologousDiscreteCharacterData<NaturalNumbersState>();

                    // read the ancestral state trace
                    for (size_t i = 0; i < tips.size(); ++i)
                    {
                        size_t tip_index = sse->getValue().getTipIndex(tips[i]);
                        std::string tip_index_anc_str = StringUtilities::toString(tip_index + 1);
                        std::string tip_index_end_str = "end_" + StringUtilities::toString(tip_index + 1);
  
                        if (ancestral_state_traces_lookup.find(tip_index_anc_str) != ancestral_state_traces_lookup.end())
                        {
                            size_t idx = ancestral_state_traces_lookup[tip_index_anc_str];
                            tip_state_trace = &ancestral_state_traces[idx];
                        }
                        else if (ancestral_state_traces_lookup.find(tip_index_end_str) != ancestral_state_traces_lookup.end())
                        {
                            size_t idx = ancestral_state_traces_lookup[tip_index_end_str];
                            tip_state_trace = &ancestral_state_traces[idx];
                        }
                        else
                            throw RbException("Can't find tip_state_trace!");
                        const std::vector<std::string>& tip_state_vector = tip_state_trace->getValues();
                        std::string state_str = tip_state_vector[index_sample];
  
                        // create a taxon data object for each tip
                        DiscreteTaxonData<NaturalNumbersState> this_tip_data = DiscreteTaxonData<NaturalNumbersState>(tips[tip_index]);
                        NaturalNumbersState state = NaturalNumbersState(0, num_states);
                        state.setState(state_str);
                        this_tip_data.addCharacter(state);
                        tip_data->addTaxonData(this_tip_data);
                    }
                   
                    // finally set the tip data to the sampled values
                    static_cast<TreeDiscreteCharacterData*>( &sse->getValue() )->setCharacterData(tip_data);
                }
               
                try 
                {
                    // redraw new values
                    the_node->redraw();
                
                    // we need to store the new simulated data
                    the_node->writeToFile(sim_directory_name);
                }
                catch (RbException &e)
                {
                    
                    std::cerr << "Problem in Posterior Predictive Simulation:" << std::endl;
                    std::cerr << e.getMessage() << std::endl;
                    // skip this simulation
                }
                catch (...)
                {
                    
                    std::cerr << "Problem occurred." << std::endl;
                    // skip this simulation
                }
            }
            
        }
        
    } // end for over all samples
    
}
