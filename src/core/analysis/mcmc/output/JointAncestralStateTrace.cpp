#include "JointAncestralStateTrace.h"

#include <boost/lexical_cast.hpp>
#include <cassert>
#include <numeric>
#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <iterator>
#include <map>
#include <sstream>
#include <string>

#include "ProgressBar.h"
#include "RbBitSet.h"
#include "RbException.h"
#include "RbFileManager.h"
#include "RlUserInterface.h"
#include "StringUtilities.h"
#include "TopologyNode.h"
#include "Tree.h"
#include "boost/algorithm/string/trim.hpp"

using namespace RevBayesCore;

JointAncestralStateTrace::JointAncestralStateTrace(std::vector<AncestralStateTrace> at, TraceTree tt ) :
    ancestral_state_traces(at),
    tree_trace(tt),
    burnin(0)
{
    num_sampled_states = ancestral_state_traces[0].size();
    
    if ( tree_trace.size() > 0 && num_sampled_states != tree_trace.size() )
    {
        throw RbException("The tree trace and the ancestral state trace must contain the same number of samples.");
    }
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */

JointAncestralStateTrace* JointAncestralStateTrace::clone(void) const
{
    
    return new JointAncestralStateTrace(*this);
}

void JointAncestralStateTrace::setBurnin(size_t b)
{
    // get the number of ancestral state samples and the number of tree samples
    if ( b >= ancestral_state_traces[0].size() )
    {
        throw RbException("Burnin size is too large for the ancestral state trace.");
    }
    
    burnin = b;
}

/**
 *
 * Helper function for ancestralStateTree() and cladoAncestralStateTree() that collects and sorts joint ancestral state samples.
 *
 */
void JointAncestralStateTrace::collectJointAncestralStateSamples(int site, Tree &final_summary_tree, const std::vector<TopologyNode*> &summary_nodes, std::vector<std::vector<double> > &pp_end, std::vector<std::vector<double> > &pp_start, std::vector<std::vector<std::string> > &end_states, std::vector<std::vector<std::string> > &start_states, bool clado, ProgressBar &progress, bool verbose)
{
    // loop through all ancestral state samples
    
    // check if sampled tree topology matches the summary tree
    
    // if so collect sampled states
    
    // done looping through samples, now sort them
    
    // save the top 3 in start/end_states and pp vectors
}


/**
 *
 * Helper function for ancestralStateTree() and cladoAncestralStateTree() that traverses the tree from root to tips collecting ancestral state samples.
 *
 */
void JointAncestralStateTrace::recursivelyCollectAncestralStateSamples(size_t node_index, std::string map_parent_state, bool root, bool conditional, int site, Tree &final_summary_tree, const std::vector<TopologyNode*> &summary_nodes, std::vector<std::vector<double> > &pp_end, std::vector<std::vector<double> > &pp_start, std::vector<double> &pp_clade, std::vector<std::vector<std::string> > &end_states, std::vector<std::vector<std::string> > &start_states, bool clado, ProgressBar &progress, size_t &num_finished_nodes, bool verbose)
{
    
    size_t parent_node_index = 0;
    if ( root == false )
    {
        parent_node_index = summary_nodes[node_index]->getParent().getIndex();
    }
    
    size_t sample_clade_index = -1;
    size_t parent_sample_clade_index = -1;
    AncestralStateTrace* ancestral_state_trace_end_state = NULL;
    AncestralStateTrace* ancestral_state_trace_start_1 = NULL;
    AncestralStateTrace* ancestral_state_trace_start_2 = NULL;
    AncestralStateTrace* parent_ancestral_state_trace = NULL;
    bool trace_found_end_state = false;
    bool trace_found_start_1 = false;
    bool trace_found_start_2 = false;
    bool parent_trace_found = false;
    size_t num_samples_end = 0;
    size_t num_samples_start_1 = 0;
    size_t num_samples_start_2 = 0;
    size_t num_samples_clade = 0;
    
    std::vector<TopologyNode*> sn = summary_nodes;
    
    
    // build a lookup map for the ancestral state trace labels -> indices
    std::map<std::string, size_t> ancestral_state_traces_lookup;
    for (size_t z = 0; z < ancestral_state_traces.size(); z++)
    {
        ancestral_state_traces_lookup[ ancestral_state_traces[z].getParameterName() ] = z;
    }
    
    // loop through all the ancestral state samples
    for (size_t j = burnin; j < num_sampled_states; ++j)
    {
        if ( verbose == true )
        {
            progress.update( num_finished_nodes * num_sampled_states + num_sampled_states * (j - burnin) / (num_sampled_states - burnin) );
        }
        
        // if necessary, get the sampled tree from the tree trace
        const Tree &sample_tree = usingTreeTrace() ? tree_trace.objectAt( j ) : final_summary_tree;
        const TopologyNode& sample_root = sample_tree.getRoot();
        
         bool parent_sample_clade_found = true;
        
        if ( usingTreeTrace() == true )
        {
            // check if the clade in the summary tree is also in the sampled tree
            //            sample_clade_index = sample_root.getCladeIndex( summary_nodes[node_index] );
            //            parent_sample_clade_index = sample_root.getCladeIndex( summary_nodes[parent_node_index] );
            try
            {
                sample_clade_index = sample_root.getCladeIndex( summary_nodes[node_index] );
            }
            catch(RbException&)
            {
                continue;
            }
            try
            {
                parent_sample_clade_index = sample_root.getCladeIndex( summary_nodes[parent_node_index] );
            }
            catch(RbException&)
            {
                parent_sample_clade_found = false;
            }
            
            // then we must find the ancestral state traces for this sampled node
            trace_found_end_state = false;
            trace_found_start_1 = false;
            trace_found_start_2 = false;
            parent_trace_found = false;
        }
        else
        {
            sample_clade_index = summary_nodes[node_index]->getIndex();
            parent_sample_clade_index = parent_node_index;
        }
        
        // record the states if the sample tree contains the summary node's clade
        
        num_samples_clade += 1;
        
        size_t sample_clade_index_child_1 = 0;
        size_t sample_clade_index_child_2 = 0;
        
        bool found_child_clade_1 = true;
        bool found_child_clade_2 = true;
        
        if ( !summary_nodes[node_index]->isTip() && clado == true )
        {
            const TopologyNode& sample_node = sample_tree.getNode( sample_clade_index );
            try
            {
                sample_clade_index_child_1 = sample_node.getCladeIndex( &summary_nodes[node_index]->getChild(0) );
            }
            catch(RbException&)
            {
                found_child_clade_1 = false;
            }
            try
            {
                sample_clade_index_child_2 = sample_node.getCladeIndex( &summary_nodes[node_index]->getChild(1) );
            }
            catch(RbException&)
            {
                found_child_clade_2 = false;
            }
        }
        
        std::string sample_clade_index_anc_str = StringUtilities::toString(sample_clade_index + 1);
        std::string sample_clade_index_end_str = "end_" + StringUtilities::toString(sample_clade_index + 1);
        std::string sample_clade_index_child_1_start_str = "start_" + StringUtilities::toString(sample_clade_index_child_1 + 1);
        std::string sample_clade_index_child_2_start_str = "start_" + StringUtilities::toString(sample_clade_index_child_2 + 1);
        
        // find the appropriate end state
        if (ancestral_state_traces_lookup.find(sample_clade_index_anc_str) != ancestral_state_traces_lookup.end())
        {
            size_t idx = ancestral_state_traces_lookup[sample_clade_index_anc_str];
            ancestral_state_trace_end_state = &ancestral_state_traces[idx];
            trace_found_end_state = true;
        }
        else if (ancestral_state_traces_lookup.find(sample_clade_index_end_str) != ancestral_state_traces_lookup.end())
        {
            size_t idx = ancestral_state_traces_lookup[sample_clade_index_end_str];
            ancestral_state_trace_end_state = &ancestral_state_traces[idx];
            trace_found_end_state = true;
        }
        
        // find start state traces if necessary
        if ( clado == true && (!summary_nodes[node_index]->isTip()) )
        {
            
            if (ancestral_state_traces_lookup.find(sample_clade_index_child_1_start_str) != ancestral_state_traces_lookup.end() && found_child_clade_1 == true )
            {
                size_t idx = ancestral_state_traces_lookup[sample_clade_index_child_1_start_str];
                ancestral_state_trace_start_1 = &ancestral_state_traces[idx];
                trace_found_start_1 = true;
            }
            
            if (ancestral_state_traces_lookup.find(sample_clade_index_child_2_start_str) != ancestral_state_traces_lookup.end() && found_child_clade_2 == true )
            {
                size_t idx = ancestral_state_traces_lookup[sample_clade_index_child_2_start_str];
                ancestral_state_trace_start_2 = &ancestral_state_traces[idx];
                trace_found_start_2 = true;
            }
        }
        else
        {
            trace_found_start_1 = true;
            trace_found_start_2 = true;
        }
        
        // if we are conditioning on the parent's state we must get the corresponding sample from the parent
        if ( conditional == true && root == false && parent_trace_found == false && parent_sample_clade_found == true )
        {
            for (size_t k = 0; k < ancestral_state_traces.size(); ++k)
            {
                // if we have an ancestral state trace from an anagenetic-only process
                if (ancestral_state_traces[k].getParameterName() == StringUtilities::toString(parent_sample_clade_index + 1))
                {
                    parent_ancestral_state_trace = &ancestral_state_traces[k];
                    parent_trace_found = true;
                    break;
                }
                // if we have an ancestral state trace from a cladogenetic process
                if (ancestral_state_traces[k].getParameterName() == "end_" + StringUtilities::toString(parent_sample_clade_index + 1))
                {
                    parent_ancestral_state_trace = &ancestral_state_traces[k];
                    parent_trace_found = true;
                    break;
                }
            }
        }
        
        // get the sampled ancestral state for this iteration
        const std::vector<std::string>& ancestral_state_vector_end = ancestral_state_trace_end_state->getValues();
        std::string ancestral_state_end = getSiteState( ancestral_state_vector_end[j], site );
        
        // get the sampled ancestral state from the parent node
        bool count_sample = false;
        std::string sampled_parent_state = "";
        if ( conditional == true && root == false && parent_trace_found == true )
        {
            if (not parent_ancestral_state_trace)
                throw RbException("Can't find parent_ancestral_state_trace");

            const std::vector<std::string>& parent_ancestral_state_vector = parent_ancestral_state_trace->getValues();
            sampled_parent_state = getSiteState( parent_ancestral_state_vector[j], site );
            
            // condition on the parent state?
            if ( sampled_parent_state == map_parent_state )
            {
                count_sample = true;
            }
        }
        
        if ( conditional == false || root == true )
        {
            count_sample = true;
        }
        
        // finally add the sample to our vectors of samples
        if ( count_sample == true )
        {
            if ( clado == true && summary_nodes[node_index]->isTip() == false )
            {
                std::string ancestral_state_start_1 = "";
                std::string ancestral_state_start_2 = "";
                size_t child1 = summary_nodes[node_index]->getChild(0).getIndex();
                size_t child2 = summary_nodes[node_index]->getChild(1).getIndex();
                
                // find & store end states
                bool end_state_found = false;
                size_t k_end = 0;
                for (; k_end < pp_end[node_index].size(); k_end++)
                {
                    if (end_states[node_index][k_end] == ancestral_state_end) {
                        end_state_found = true;
                        break;
                    }
                }
                // update the end pp and states vectors
                if ( end_state_found == false )
                {
                    pp_end[node_index].push_back(1.0);
                    end_states[node_index].push_back( ancestral_state_end );
                    
                }
                else
                {
                    pp_end[node_index][k_end] += 1.0;
                }
                
                // find & store start_state_1
                if (trace_found_start_1 && trace_found_start_2)
                {
                    std::vector<std::string> ancestral_state_trace_start_1_vector = ancestral_state_trace_start_1->getValues();
                    ancestral_state_start_1 = getSiteState( ancestral_state_trace_start_1_vector[j], site );
                    
                    bool start_state_1_found = false;
                    size_t k_start_1 = 0;
                    for (; k_start_1 < pp_start[child1].size(); k_start_1++)
                    {
                        if (start_states[child1][k_start_1] == ancestral_state_start_1) {
                            start_state_1_found = true;
                            break;
                        }
                    }
                    // update the end pp and states vectors
                    if ( start_state_1_found == false )
                    {
                        pp_start[child1].push_back(1.0);
                        start_states[child1].push_back( ancestral_state_start_1 );
                        
                    }
                    else
                    {
                        pp_start[child1][k_start_1] += 1.0;
                    }
                    
                    num_samples_start_1 += 1;
                    //                    }
                    //
                    //                    // find & store start_state_1
                    //                    if (trace_found_start_2)
                    //                    {
                    std::vector<std::string> ancestral_state_trace_start_2_vector = ancestral_state_trace_start_2->getValues();
                    ancestral_state_start_2 = getSiteState( ancestral_state_trace_start_2_vector[j], site );
                    
                    bool start_state_2_found = false;
                    size_t k_start_2 = 0;
                    for (; k_start_2 < pp_start[child2].size(); k_start_2++)
                    {
                        if (start_states[child2][k_start_2] == ancestral_state_start_2) {
                            start_state_2_found = true;
                            break;
                        }
                    }
                    // update the end pp and states vectors
                    if ( start_state_2_found == false )
                    {
                        pp_start[child2].push_back(1.0);
                        start_states[child2].push_back( ancestral_state_start_2 );
                        
                    }
                    else
                    {
                        pp_start[child2][k_start_2] += 1.0;
                    }
                    
                    num_samples_start_2 += 1;
                }
            }
            else
            {
                bool state_found = false;
                size_t k = 0;
                for (; k < pp_end[node_index].size(); k++)
                {
                    if ( end_states[node_index][k] == ancestral_state_end )
                    {
                        state_found = true;
                        break;
                    }
                }
                // update the pp and states vectors
                if ( state_found == false )
                {
                    pp_end[node_index].push_back(1.0);
                    end_states[node_index].push_back(ancestral_state_end);
                }
                else
                {
                    pp_end[node_index][k] += 1.0;
                }
            }
            num_samples_end += 1;
        }
    }
    
    // scale posteriors
    for (size_t i = 0; i < pp_end[node_index].size(); i++)
    {
        pp_end[node_index][i] /= num_samples_end;
    }
    if ( !summary_nodes[node_index]->isTip() ) {
        size_t ch_1_idx = summary_nodes[node_index]->getChild(0).getIndex();
        for (size_t i = 0; i < pp_start[ch_1_idx].size(); i++)
        {
            pp_start[ch_1_idx][i] /= num_samples_start_1;
        }
        size_t ch_2_idx = summary_nodes[node_index]->getChild(1).getIndex();
        for (size_t i = 0; i < pp_start[ch_2_idx].size(); i++)
        {
            pp_start[ch_2_idx][i] /= num_samples_start_2;
        }
    }
    
    pp_clade[node_index] = (double)num_samples_clade / (num_sampled_states - burnin);
    
    std::string map_state = "";
    double max_pp = 0.0;
    if ( conditional == true )
    {
        // find the MAP state
        for (size_t i = 0; i < end_states[node_index].size(); i++)
        {
            if (pp_end[node_index][i] > max_pp)
            {
                map_state = end_states[node_index][i];
                max_pp = pp_end[node_index][i];
            }
        }
    }
    
    // update nodes finished for the progress bar
    num_finished_nodes += 1;
    
    // recurse through tree
    std::vector<int> children_indices = summary_nodes[node_index]->getChildrenIndices();
    for (size_t i = 0; i < children_indices.size(); i++)
    {
        recursivelyCollectAncestralStateSamples(children_indices[i], map_state, false, conditional, site, final_summary_tree, summary_nodes, pp_end, pp_start, pp_clade, end_states, start_states, clado, progress, num_finished_nodes, verbose);
    }
    
}


/*
 * This method calculates ancestral character states for the nodes on the input_summary_tree. When summary_stat = MAP, it annotates
 * the summary tree with the posterior probabilities of the 3 most probable states. The method requires a vector of traces containing
 * sampled ancestral states, and optionally works with a trace containing sampled trees that correspond to the
 * ancestral state samples, enabling ancestral states to be estimated over a distribution of trees. In this case the annotated
 * posterior probability for a given character state is the probability of the node existing times the probability of the
 * node being in the character state (see Pagel et al. 2004).
 *
 * This method summarizes both conditional and marginal ancestral states. When conditional = false the method summarizes marginal
 * ancestral states. conditional = true should only be used if the states were sampled using the JointConditionalAncestralStateMonitor.
 * When summarizing conditional ancestral states the posterior probabilities of each state are conditioned on the parent node being
 * in the MAP state.
 */
Tree* JointAncestralStateTrace::ancestralStateTree(const Tree &input_summary_tree, std::string summary_stat, size_t num_states, int site, bool conditional, bool joint, bool verbose )
{
    
    if ( summary_stat != "MAP" && (conditional == true || joint == true) )
    {
        throw RbException("Invalid reconstruction type: mean ancestral states can only be calculated with a marginal reconstruction.");
    }
    
    std::stringstream ss;
    ss << "Compiling " << summary_stat << " ancestral states from " << num_sampled_states << " samples in the ancestral state trace, using a burnin of " << burnin << " samples.\n";
    RBOUT(ss.str());
    
    RBOUT("Calculating ancestral state posteriors...\n");
    
    // allocate memory for the new summary tree
    Tree* final_summary_tree = new Tree( input_summary_tree );
    
    // 2-d vectors to keep the data (posteriors and ancestral states) of the input_summary_tree nodes: [node][data]
    const std::vector<TopologyNode*> &summary_nodes = final_summary_tree->getNodes();
    std::vector<std::vector<double> > pp_end( summary_nodes.size(), std::vector<double>() );
    std::vector<std::vector<double> > pp_start( summary_nodes.size(), std::vector<double>() );
    std::vector<double> pp_clade( summary_nodes.size(), 0.0 );
    
    std::vector<std::vector<std::string> > states( summary_nodes.size(), std::vector<std::string>() );
    
    size_t num_finished_nodes = 0;
    ProgressBar progress = ProgressBar( summary_nodes.size() * num_sampled_states, 0 );
    if ( verbose == true )
    {
        progress.start();
    }
    
    if ( joint == true )
    {
        // collect and sort the joint samples
        collectJointAncestralStateSamples(site, *final_summary_tree, summary_nodes, pp_end, pp_start, states, states, false, progress, verbose);
    }
    else
    {
        // recurse through summary tree and collect ancestral state samples
        size_t node_index = final_summary_tree->getRoot().getIndex();
        recursivelyCollectAncestralStateSamples(node_index, "", true, conditional, site, *final_summary_tree, summary_nodes, pp_end, pp_start, pp_clade, states, states, false, progress, num_finished_nodes, verbose);
    }
    
    if ( verbose == true )
    {
        progress.finish();
    }
    
    // summarize ancestral states
    if (summary_stat == "MAP")
    {
        // find the K most probable ancestral states for each node and add them to the tree as parameters
        std::vector< std::vector<std::string> > anc_state = std::vector< std::vector<std::string> >(num_states, std::vector<std::string>());
        std::vector< std::vector<double> > anc_state_pp = std::vector< std::vector<double> >(num_states, std::vector<double>());
        std::vector<double> anc_state_other_pp;
        
        std::vector<double> posteriors;
        
        for (int node_index = 0; node_index < summary_nodes.size(); ++node_index)
        {
            

            double state_pp = 0.0;
            double other_pp = 0.0;
            double total_node_pp = 0.0;
            
            std::string state = "";
            
            std::multimap<double,std::string, std::greater<double> > state_map;
            
            // loop through all states for this node
            for (int j = 0; j < pp_end[node_index].size(); j++)
            {
                total_node_pp += pp_end[node_index][j];
                state_map.insert( std::pair<double, std::string>(pp_end[node_index][j], states[node_index][j]) );
            }
            
            posteriors.push_back(total_node_pp);
            other_pp = total_node_pp;
//            size_t current_state_count = 0;
            std::multimap<double,std::string>::iterator it=state_map.begin();
            for (size_t j=0; j<num_states; ++j)
            {
                
                if ( it != state_map.end() )
                {
                    state_pp = it->first;
                    other_pp -= state_pp;
                    anc_state[j].push_back( it->second );
                    anc_state_pp[j].push_back(state_pp);
                    ++it;
                }
                else
                {
                    anc_state[j].push_back( "NA" );
                    anc_state_pp[j].push_back(0.0);
                }
            }
            
            if (other_pp > 0.0001)
            {
                anc_state_other_pp.push_back(other_pp);
            }
            else
            {
                anc_state_other_pp.push_back(0.0);
            }
            
        }
        
        final_summary_tree->clearNodeParameters();
        final_summary_tree->addNodeParameter("posterior", posteriors, false);
        for (size_t i=0; i<num_states; ++i)
        {
            final_summary_tree->addNodeParameter("anc_state_"+StringUtilities::to_string(i+1), anc_state[i], false);
            final_summary_tree->addNodeParameter("anc_state_"+StringUtilities::to_string(i+1)+"_pp", anc_state_pp[i], false);
        }
        final_summary_tree->addNodeParameter("anc_state_other_pp", anc_state_other_pp, false);
        
    }
    else
    {
        // calculate the mean and 95% CI and add to the tree as annotation
        double hpd = 0.95;
        std::vector<double> means( summary_nodes.size(), 0.0 );
        std::vector<double> uppers( summary_nodes.size(), 0.0 );
        std::vector<double> lowers( summary_nodes.size(), 0.0 );
        std::vector<double> posteriors( summary_nodes.size(), 0.0 );
        
        for (int i = 0; i < summary_nodes.size(); i++)
        {
            
            double node_pp = 0.0;
            std::vector<double> state_samples;
            
            // loop through all states for this node and collect samples
            for (int j = 0; j < pp_end[i].size(); j++)
            {
                node_pp += pp_end[i][j];
                state_samples.push_back( boost::lexical_cast<double>( states[i][j] ) );
            }
            posteriors[i] = node_pp;
            
            // calculate mean value
            double samples_sum = std::accumulate(state_samples.begin(), state_samples.end(), 0.0);
            double mean = samples_sum / state_samples.size();
            means[i] = mean;
            
            // sort the samples by frequency and calculate interval
            std::sort(state_samples.begin(), state_samples.end());
            size_t interval_start = ( (1.0 - hpd) / 2.0 ) * state_samples.size();
            size_t interval_end = ( 1.0 - (1.0 - hpd) / 2.0 ) * state_samples.size();
            interval_end = (interval_end >= state_samples.size() ? state_samples.size()-1 : interval_end);
            
            double lower = state_samples[interval_start];
            double upper = state_samples[interval_end];
            
            lowers[i] = lower;
            uppers[i] = upper;
            
        }
        
        final_summary_tree->clearNodeParameters();
        final_summary_tree->addNodeParameter("posterior", posteriors, true);
        final_summary_tree->addNodeParameter("mean", means, false);
        final_summary_tree->addNodeParameter("lower_95%_CI", lowers, true);
        final_summary_tree->addNodeParameter("upper_95%_CI", uppers, true);
        
    }
    
    return final_summary_tree;
}


/*
 * This method calculates the MAP ancestral character states for the nodes on the input_tree. This method
 * is identical to the ancestralStateTree function except that is calculates the MAP states resulting from
 * a cladogenetic event, so for each node the MAP state includes the end state and the starting states for
 * the two daughter lineages.
 */
Tree* JointAncestralStateTrace::cladoAncestralStateTree2(const Tree &input_summary_tree, std::string summary_stat, size_t num_states, int site, bool conditional, bool joint, bool verbose )
{
    
    if ( summary_stat != "MAP" && (conditional == true || joint == true) )
    {
        throw RbException("Invalid reconstruction type: mean ancestral states can only be calculated with a marginal reconstruction.");
    }
    
    std::stringstream ss;
    ss << "Compiling " << summary_stat << " ancestral states from " << num_sampled_states << " samples in the ancestral state trace, using a burnin of " << burnin << " samples.\n";
    RBOUT(ss.str());
    
    RBOUT("Calculating ancestral state posteriors...\n");
    
    // allocate memory for the new tree
    Tree* final_summary_tree = new Tree(input_summary_tree);
    
    // 2-d vectors to keep the data (posteriors and states) of the inputTree nodes: [node][data]
    const std::vector<TopologyNode*> &summary_nodes = final_summary_tree->getNodes();
    std::vector<std::vector<double> > pp_end( summary_nodes.size(), std::vector<double>() );
    std::vector<std::vector<double> > pp_start( summary_nodes.size(), std::vector<double>() );
    std::vector<double> pp_clade( summary_nodes.size(), 0.0 );
    std::vector<std::vector<std::string> > end_states( summary_nodes.size(), std::vector<std::string>() );
    std::vector<std::vector<std::string> > start_states( summary_nodes.size(), std::vector<std::string>() );
    
    size_t num_finished_nodes = 0;
    ProgressBar progress = ProgressBar( summary_nodes.size() * num_sampled_states, 0 );
    if ( verbose == true )
    {
        progress.start();
    }
    
    if ( joint == true )
    {
        // collect and sort the joint samples
        collectJointAncestralStateSamples(site, *final_summary_tree, summary_nodes, pp_end, pp_start, end_states, start_states, false, progress, verbose);
    }
    else
    {
        // recurse through summary tree and collect ancestral state samples
        size_t node_index = final_summary_tree->getRoot().getIndex();
        recursivelyCollectAncestralStateSamples(node_index, "", true, conditional, site, *final_summary_tree, summary_nodes, pp_end, pp_start, pp_clade, end_states, start_states, true, progress, num_finished_nodes, verbose);
    }
    
    if ( verbose == true )
    {
        progress.finish();
    }
    
    // summarize ancestral states
    if (summary_stat == "MAP")
    {
        // find the K most probable ancestral states for each node and add them to the tree as annotations
        
        // find the K most probable ancestral states for each node and add them to the tree as parameters
        std::vector< std::vector<std::string> > start_state = std::vector< std::vector<std::string> >(num_states, std::vector<std::string>());
        std::vector< std::vector<double> > start_state_pp = std::vector< std::vector<double> >(num_states, std::vector<double>());
        std::vector<double> start_state_other_pp;
        std::vector< std::vector<std::string> > end_state = std::vector< std::vector<std::string> >(num_states, std::vector<std::string>());
        std::vector< std::vector<double> > end_state_pp = std::vector< std::vector<double> >(num_states, std::vector<double>());
        std::vector<double> end_state_other_pp;
        
        std::vector<double> posteriors;
        
        for (int node_index = 0; node_index < summary_nodes.size(); ++node_index)
        {
            
            double start_state_pp_j = 0.0;
            double start_other_pp = 0.0;
            double start_total_node_pp = 0.0;
            double end_state_pp_j = 0.0;
            double end_other_pp = 0.0;
            double end_total_node_pp = 0.0;
            
//            std::string state = "";
            
            std::multimap<double,std::string, std::greater<double> > start_state_map;
            std::multimap<double,std::string, std::greater<double> > end_state_map;
            
            // loop through all start states for this node
            for (int j = 0; j < pp_start[node_index].size(); j++)
            {
                start_total_node_pp += pp_start[node_index][j];
                start_state_map.insert( std::pair<double, std::string>(pp_start[node_index][j], start_states[node_index][j]) );
            }
            
            // loop through all end states for this node
            for (int j = 0; j < pp_end[node_index].size(); j++)
            {
                end_total_node_pp += pp_end[node_index][j];
                end_state_map.insert( std::pair<double, std::string>(pp_end[node_index][j], end_states[node_index][j]) );
            }
            
            posteriors.push_back(pp_clade[node_index]);
            end_other_pp = end_total_node_pp;
            start_other_pp = start_total_node_pp;
            
            // process start state
            std::multimap<double,std::string>::iterator start_it = start_state_map.begin();
            for (size_t j=0; j<num_states; ++j)
            {
                
                if ( start_it != start_state_map.end() )
                {
                    start_state_pp_j = start_it->first;
                    start_other_pp -= start_state_pp_j;
                    start_state[j].push_back( start_it->second );
                    start_state_pp[j].push_back(start_state_pp_j);
                    ++start_it;
                }
                else
                {
                    start_state[j].push_back( "NA" );
                    start_state_pp[j].push_back(0.0);
                }
            }
            
            if (start_other_pp > 0.0001)
            {
                start_state_other_pp.push_back(start_other_pp);
            }
            else
            {
                start_state_other_pp.push_back(0.0);
            }
            
            // process end state
            std::multimap<double,std::string>::iterator end_it = end_state_map.begin();
            for (size_t j=0; j<num_states; ++j)
            {
                
                if ( end_it != end_state_map.end() )
                {
                    end_state_pp_j = end_it->first;
                    end_other_pp -= end_state_pp_j;
                    end_state[j].push_back( end_it->second );
                    end_state_pp[j].push_back(end_state_pp_j);
                    ++end_it;
                }
                else
                {
                    end_state[j].push_back( "NA" );
                    end_state_pp[j].push_back(0.0);
                }
            }
            
            if (end_other_pp > 0.0001)
            {
                end_state_other_pp.push_back(end_other_pp);
            }
            else
            {
                end_state_other_pp.push_back(0.0);
            }
            
        }
        
        // ordered this way to match original newick annotations for
        // original cladoAncestralStateTree as closely as possible
        final_summary_tree->clearNodeParameters();
        final_summary_tree->addNodeParameter("posterior", posteriors, false);
        for (size_t i=0; i<num_states; ++i)
        {
            final_summary_tree->addNodeParameter("end_state_"+StringUtilities::to_string(i+1), end_state[i], false);
        }
        for (size_t i=0; i<num_states; ++i)
        {
            final_summary_tree->addNodeParameter("end_state_"+StringUtilities::to_string(i+1)+"_pp", end_state_pp[i], false);
            
        }
        final_summary_tree->addNodeParameter("end_state_other_pp", end_state_other_pp, false);
        for (size_t i=0; i<num_states; ++i)
        {
            final_summary_tree->addNodeParameter("start_state_"+StringUtilities::to_string(i+1), start_state[i], false);
        }
        for (size_t i=0; i<num_states; ++i)
        {
            final_summary_tree->addNodeParameter("start_state_"+StringUtilities::to_string(i+1)+"_pp", start_state_pp[i], false);
        }
        final_summary_tree->addNodeParameter("start_state_other_pp", start_state_other_pp, false);
        
        
    }
    else
    {
        // calculate the mean and 95% CI and add to the tree as annotation
        double hpd = 0.95;
        std::vector<double> start_means( summary_nodes.size(), 0.0 );
        std::vector<double> start_uppers( summary_nodes.size(), 0.0 );
        std::vector<double> start_lowers( summary_nodes.size(), 0.0 );
        std::vector<double> end_means( summary_nodes.size(), 0.0 );
        std::vector<double> end_uppers( summary_nodes.size(), 0.0 );
        std::vector<double> end_lowers( summary_nodes.size(), 0.0 );
        std::vector<double> posteriors( summary_nodes.size(), 0.0 );
        
        for (int i = 0; i < summary_nodes.size(); i++)
        {
            
            double node_pp = 0.0;
            std::vector<double> state_samples_end;
            std::vector<double> state_samples_start;
            
            // loop through all states for this node and collect samples
            for (int j = 0; j < pp_end[i].size(); j++)
            {
                node_pp += pp_end[i][j];
                state_samples_end.push_back( boost::lexical_cast<double>( end_states[i][j] ) );
                state_samples_start.push_back( boost::lexical_cast<double>( start_states[i][j] ) );
            }
            posteriors[i] = node_pp;
            
            // calculate mean value for end states
            double samples_sum = std::accumulate(state_samples_end.begin(), state_samples_end.end(), 0.0);
            double mean = samples_sum / state_samples_end.size();
            end_means[i] = mean;
            
            // calculate mean value for start states
            samples_sum = std::accumulate(state_samples_start.begin(), state_samples_start.end(), 0.0);
            mean = samples_sum / state_samples_start.size();
            start_means[i] = mean;
            
            // sort the samples by frequency and calculate interval for the end states
            std::sort(state_samples_end.begin(), state_samples_end.end());
            size_t interval_start = ( (1.0 - hpd) / 2.0 ) * state_samples_end.size();
            size_t interval_end = ( 1.0 - (1.0 - hpd) / 2.0 ) * state_samples_end.size();
            interval_end = (interval_end >= state_samples_end.size() ? state_samples_end.size()-1 : interval_end);
            
            double lower = state_samples_end[interval_start];
            double upper = state_samples_end[interval_end];
            
            end_lowers[i] = lower;
            end_uppers[i] = upper;
            
            // sort the samples by frequency and calculate interval for the start states
            std::sort(state_samples_start.begin(), state_samples_start.end());
            interval_start = ( (1.0 - hpd) / 2.0 ) * state_samples_start.size();
            interval_end = ( 1.0 - (1.0 - hpd) / 2.0 ) * state_samples_start.size();
            interval_end = (interval_end >= state_samples_start.size() ? state_samples_start.size()-1 : interval_end);
            
            lower = state_samples_start[interval_start];
            upper = state_samples_start[interval_end];
            
            start_lowers[i] = lower;
            start_uppers[i] = upper;
        }
        
        final_summary_tree->clearNodeParameters();
        final_summary_tree->addNodeParameter("posterior", posteriors, true);
        final_summary_tree->addNodeParameter("end_mean", end_means, false);
        final_summary_tree->addNodeParameter("end_lower_95%_CI", end_lowers, true);
        final_summary_tree->addNodeParameter("end_upper_95%_CI", end_uppers, true);
        final_summary_tree->addNodeParameter("start_mean", start_means, true);
        final_summary_tree->addNodeParameter("start_lower_95%_CI", start_lowers, true);
        final_summary_tree->addNodeParameter("start_upper_95%_CI", start_uppers, true);
    }
    
    return final_summary_tree;
}


/*
 * This method calculates the MAP ancestral character states for the nodes on the input_tree. This method
 * is identical to the ancestralStateTree function except that is calculates the MAP states resulting from
 * a cladogenetic event, so for each node the MAP state includes the end state and the starting states for
 * the two daughter lineages.
 */
Tree* JointAncestralStateTrace::cladoAncestralStateTree(const Tree &input_summary_tree, std::string summary_stat, size_t num_states, int site, bool conditional, bool joint, bool verbose )
{
    
    if ( summary_stat != "MAP" && (conditional == true || joint == true) )
    {
        throw RbException("Invalid reconstruction type: mean ancestral states can only be calculated with a marginal reconstruction.");
    }
    
    std::stringstream ss;
    ss << "Compiling " << summary_stat << " ancestral states from " << num_sampled_states << " samples in the ancestral state trace, using a burnin of " << burnin << " samples.\n";
    RBOUT(ss.str());
    
    RBOUT("Calculating ancestral state posteriors...\n");
    
    // allocate memory for the new tree
    Tree* final_summary_tree = new Tree(input_summary_tree);
    
    // 2-d vectors to keep the data (posteriors and states) of the inputTree nodes: [node][data]
    const std::vector<TopologyNode*> &summary_nodes = final_summary_tree->getNodes();
    std::vector<std::vector<double> > pp_end( summary_nodes.size(), std::vector<double>() );
    std::vector<std::vector<double> > pp_start( summary_nodes.size(), std::vector<double>() );
    std::vector<double> pp_clade( summary_nodes.size(), 0.0 );
    std::vector<std::vector<std::string> > end_states( summary_nodes.size(), std::vector<std::string>() );
    std::vector<std::vector<std::string> > start_states( summary_nodes.size(), std::vector<std::string>() );
    
    size_t num_finished_nodes = 0;
    ProgressBar progress = ProgressBar( summary_nodes.size() * num_sampled_states, 0 );
    if ( verbose == true )
    {
        progress.start();
    }
    
    if ( joint == true )
    {
        // collect and sort the joint samples
        collectJointAncestralStateSamples(site, *final_summary_tree, summary_nodes, pp_end, pp_start, end_states, start_states, false, progress, verbose);
    }
    else
    {
        // recurse through summary tree and collect ancestral state samples
        size_t node_index = final_summary_tree->getRoot().getIndex();
        recursivelyCollectAncestralStateSamples(node_index, "", true, conditional, site, *final_summary_tree, summary_nodes, pp_end, pp_start, pp_clade, end_states, start_states, true, progress, num_finished_nodes, verbose);
    }
    
    if ( verbose == true )
    {
        progress.finish();
    }
    
    // summarize ancestral states
    if (summary_stat == "MAP")
    {
        // find the 3 most probable ancestral states for each node and add them to the tree as annotations
        
        std::vector<std::string> end_state_1( summary_nodes.size(), std::string("NA") );
        std::vector<std::string> end_state_2( summary_nodes.size(), std::string("NA") );
        std::vector<std::string> end_state_3( summary_nodes.size(), std::string("NA") );
        
        std::vector<double> end_state_1_pp( summary_nodes.size(), 0.0 );
        std::vector<double> end_state_2_pp( summary_nodes.size(), 0.0 );
        std::vector<double> end_state_3_pp( summary_nodes.size(), 0.0 );
        std::vector<double> end_state_other_pp( summary_nodes.size(), 0.0 );
        
        std::vector<std::string> start_state_1( summary_nodes.size(), std::string("NA") );
        std::vector<std::string> start_state_2( summary_nodes.size(), std::string("NA") );
        std::vector<std::string> start_state_3( summary_nodes.size(), std::string("NA") );
        
        std::vector<double> start_state_1_pp( summary_nodes.size(), 0.0 );
        std::vector<double> start_state_2_pp( summary_nodes.size(), 0.0 );
        std::vector<double> start_state_3_pp( summary_nodes.size(), 0.0 );
        std::vector<double> start_state_other_pp( summary_nodes.size(), 0.0 );
        
        std::vector<double> posteriors( summary_nodes.size(), 0.0 );
        
        for (int i = 0; i < summary_nodes.size(); i++)
        {

            
            std::string end_state1 = "";
            std::string end_state2 = "";
            std::string end_state3 = "";
            
            std::string start_child1_state1 = "";
            std::string start_child1_state2 = "";
            std::string start_child1_state3 = "";
            
            
            // end states
            std::vector<double> best_end_pp_node;
            std::vector<std::string> best_end_states_node;
            computeMarginalCladogeneticStateProbs(pp_end[i], end_states[i], best_end_pp_node, best_end_states_node);
            
            end_state_1[ i ] = best_end_states_node[0];
            end_state_2[ i ] = best_end_states_node[1];
            end_state_3[ i ] = best_end_states_node[2];
            end_state_1_pp[ i ] = best_end_pp_node[0];
            end_state_2_pp[ i ] = best_end_pp_node[1];
            end_state_3_pp[ i ] = best_end_pp_node[2];
            
            
            // start state, child 1
            std::vector<double> best_start_pp;
            std::vector<std::string> best_start_states;
            computeMarginalCladogeneticStateProbs(pp_start[i], start_states[i], best_start_pp, best_start_states);
            
            start_state_1[ i ] = best_start_states[0];
            start_state_2[ i ] = best_start_states[1];
            start_state_3[ i ] = best_start_states[2];
            start_state_1_pp[ i ] = best_start_pp[0];
            start_state_2_pp[ i ] = best_start_pp[1];
            start_state_3_pp[ i ] = best_start_pp[2];
            
            // sum best state probabilities
            double total_node_pp = best_end_pp_node[0] + best_end_pp_node[1] + best_end_pp_node[2];
            double end_other_pp = 1.0 - total_node_pp;
            double start_child1_other_pp = 1.0 - best_start_pp[0] - best_start_pp[1] - best_start_pp[2];
            posteriors[i] = pp_clade[i];
            
            
            end_state_other_pp[i] = 0.0;
            start_state_other_pp[i] = 0.0;
            if (end_other_pp > 0.0001)
            {
                end_state_other_pp[i] = end_other_pp;
            }
            if (start_child1_other_pp)
            {
                start_state_other_pp[i] = start_child1_other_pp;
                
            }
            
            if (i == final_summary_tree->getRoot().getIndex())
            {
                start_state_1[i] = end_state_1[i];
                start_state_2[i] = end_state_2[i];
                start_state_3[i] = end_state_3[i];
                
                start_state_1_pp[i] = end_state_1_pp[i];
                start_state_2_pp[i] = end_state_2_pp[i];
                start_state_3_pp[i] = end_state_3_pp[i];
                start_state_other_pp[i] = end_state_other_pp[i];
            }
        }
        
        final_summary_tree->clearNodeParameters();
        final_summary_tree->addNodeParameter("posterior", posteriors, false);
        
        final_summary_tree->addNodeParameter("end_state_1", end_state_1, false);
        final_summary_tree->addNodeParameter("end_state_2", end_state_2, false);
        final_summary_tree->addNodeParameter("end_state_3", end_state_3, false);
        final_summary_tree->addNodeParameter("end_state_1_pp", end_state_1_pp, false);
        final_summary_tree->addNodeParameter("end_state_2_pp", end_state_2_pp, false);
        final_summary_tree->addNodeParameter("end_state_3_pp", end_state_3_pp, false);
        final_summary_tree->addNodeParameter("end_state_other_pp", end_state_other_pp, false);
        
        final_summary_tree->addNodeParameter("start_state_1", start_state_1, false);
        final_summary_tree->addNodeParameter("start_state_2", start_state_2, false);
        final_summary_tree->addNodeParameter("start_state_3", start_state_3, false);
        final_summary_tree->addNodeParameter("start_state_1_pp", start_state_1_pp, false);
        final_summary_tree->addNodeParameter("start_state_2_pp", start_state_2_pp, false);
        final_summary_tree->addNodeParameter("start_state_3_pp", start_state_3_pp, false);
        final_summary_tree->addNodeParameter("start_state_other_pp", start_state_other_pp, false);
        
    }
    else
    {
        // calculate the mean and 95% CI and add to the tree as annotation
        double hpd = 0.95;
        std::vector<double> start_means( summary_nodes.size(), 0.0 );
        std::vector<double> start_uppers( summary_nodes.size(), 0.0 );
        std::vector<double> start_lowers( summary_nodes.size(), 0.0 );
        std::vector<double> end_means( summary_nodes.size(), 0.0 );
        std::vector<double> end_uppers( summary_nodes.size(), 0.0 );
        std::vector<double> end_lowers( summary_nodes.size(), 0.0 );
        std::vector<double> posteriors( summary_nodes.size(), 0.0 );
        
        for (int i = 0; i < summary_nodes.size(); i++)
        {
            
            double node_pp = 0.0;
            std::vector<double> state_samples_end;
            std::vector<double> state_samples_start;
            
            // loop through all states for this node and collect samples
            for (int j = 0; j < pp_end[i].size(); j++)
            {
                node_pp += pp_end[i][j];
                state_samples_end.push_back( boost::lexical_cast<double>( end_states[i][j] ) );
                state_samples_start.push_back( boost::lexical_cast<double>( start_states[i][j] ) );
            }
            posteriors[i] = node_pp;
            
            // calculate mean value for end states
            double samples_sum = std::accumulate(state_samples_end.begin(), state_samples_end.end(), 0.0);
            double mean = samples_sum / state_samples_end.size();
            end_means[i] = mean;
            
            // calculate mean value for start states
            samples_sum = std::accumulate(state_samples_start.begin(), state_samples_start.end(), 0.0);
            mean = samples_sum / state_samples_start.size();
            start_means[i] = mean;
            
            // sort the samples by frequency and calculate interval for the end states
            std::sort(state_samples_end.begin(), state_samples_end.end());
            size_t interval_start = ( (1.0 - hpd) / 2.0 ) * state_samples_end.size();
            size_t interval_end = ( 1.0 - (1.0 - hpd) / 2.0 ) * state_samples_end.size();
            interval_end = (interval_end >= state_samples_end.size() ? state_samples_end.size()-1 : interval_end);
            
            double lower = state_samples_end[interval_start];
            double upper = state_samples_end[interval_end];
            
            end_lowers[i] = lower;
            end_uppers[i] = upper;
            
            // sort the samples by frequency and calculate interval for the start states
            std::sort(state_samples_start.begin(), state_samples_start.end());
            interval_start = ( (1.0 - hpd) / 2.0 ) * state_samples_start.size();
            interval_end = ( 1.0 - (1.0 - hpd) / 2.0 ) * state_samples_start.size();
            interval_end = (interval_end >= state_samples_start.size() ? state_samples_start.size()-1 : interval_end);
            
            lower = state_samples_start[interval_start];
            upper = state_samples_start[interval_end];
            
            start_lowers[i] = lower;
            start_uppers[i] = upper;
        }
        
        final_summary_tree->clearNodeParameters();
        final_summary_tree->addNodeParameter("posterior", posteriors, true);
        final_summary_tree->addNodeParameter("end_mean", end_means, false);
        final_summary_tree->addNodeParameter("end_lower_95%_CI", end_lowers, true);
        final_summary_tree->addNodeParameter("end_upper_95%_CI", end_uppers, true);
        final_summary_tree->addNodeParameter("start_mean", start_means, true);
        final_summary_tree->addNodeParameter("start_lower_95%_CI", start_lowers, true);
        final_summary_tree->addNodeParameter("start_upper_95%_CI", start_uppers, true);
    }
    
    return final_summary_tree;
}

void JointAncestralStateTrace::computeMarginalCladogeneticStateProbs(std::vector<double> pp, std::vector<std::string> states, std::vector<double>& best_pp, std::vector<std::string>& best_states)
{
    best_pp = std::vector<double>(3, 0.0);
    best_states = std::vector<std::string>(3, "NA");
    
    // get probs per state
    std::map<std::string, double> state_sorted_map;
    for (size_t j = 0; j < states.size(); j++) {
        std::map<std::string, double>::iterator it;
        it = state_sorted_map.find( states[j] );
        if (it == state_sorted_map.end()) {
            state_sorted_map[ states[j] ] = 0.0;
        }
        state_sorted_map[ states[j] ] = state_sorted_map[ states[j] ] + pp[j];
    }
    
    // get reverse probs per state
    std::multimap<double, std::string> prob_sorted_map;
    for (std::map<std::string, double>::iterator it = state_sorted_map.begin(); it != state_sorted_map.end(); it++)
    {
        prob_sorted_map.insert( std::pair<double, std::string>( it->second, it->first) );
    }
    
    // populate most three most probable states
    size_t k = 0;
    for (std::multimap<double, std::string>::reverse_iterator rit = prob_sorted_map.rbegin(); rit != prob_sorted_map.rend(); rit++)
    {
        //        std::cout << rit->first << " " << rit->second << "\n";
        best_pp[k] = rit->first;
        best_states[k] = rit->second;
        k++;
        if (k >= 3) break;
    }
    
    // done!
    return;
}

void JointAncestralStateTrace::computeMarginalCladogeneticStateProbs2(std::vector<double> pp, std::vector<std::string> states, std::vector<double>& best_pp, std::vector<std::string>& best_states, size_t num_states)
{
    best_pp = std::vector<double>(num_states, 0.0);
    best_states = std::vector<std::string>(num_states, "NA");
    
    // get probs per state
    std::map<std::string, double> state_sorted_map;
    for (size_t j = 0; j < states.size(); j++) {
        std::map<std::string, double>::iterator it;
        it = state_sorted_map.find( states[j] );
        if (it == state_sorted_map.end()) {
            state_sorted_map[ states[j] ] = 0.0;
        }
        state_sorted_map[ states[j] ] = state_sorted_map[ states[j] ] + pp[j];
    }
    
    // get reverse probs per state
    std::multimap<double, std::string> prob_sorted_map;
    for (std::map<std::string, double>::iterator it = state_sorted_map.begin(); it != state_sorted_map.end(); it++)
    {
        prob_sorted_map.insert( std::pair<double, std::string>( it->second, it->first) );
    }
    
    // populate most three most probable states
    size_t k = 0;
    for (std::multimap<double, std::string>::reverse_iterator rit = prob_sorted_map.rbegin(); rit != prob_sorted_map.rend(); rit++)
    {
        //        std::cout << rit->first << " " << rit->second << "\n";
        best_pp[k] = rit->first;
        best_states[k] = rit->second;
        k++;
        if (k >= num_states) break;
    }
    
    // done!
    return;
}


/**
 *
 * Helper function for characterMapTree() that traverses the tree from root to tips collecting stochastic character map samples.
 *
 */
void JointAncestralStateTrace::recursivelyCollectCharacterMapSamples(size_t node_index, size_t map_parent_state, bool root, bool conditional, Tree &final_summary_tree, const std::vector<TopologyNode*> &summary_nodes, std::vector<std::string> &map_character_history, std::vector<std::string> &map_character_history_posteriors, std::vector<std::string> &map_character_history_shift_prob, ProgressBar &progress, size_t &num_finished_nodes, int NUM_TIME_SLICES, bool verbose, bool bwd_time)
{
    
    double dt = final_summary_tree.getRoot().getMaxDepth() / double(NUM_TIME_SLICES);
    
    size_t sample_clade_index;
    bool trace_found = false;
    bool parent_trace_found = false;
    AncestralStateTrace* ancestralstate_trace = NULL;
    AncestralStateTrace* parent_trace = NULL;
    
    std::vector< std::vector< std::pair<size_t, double> > > branch_maps_all = std::vector< std::vector< std::pair<size_t, double> > >();
    std::vector< std::vector< std::pair<size_t, double> > > branch_maps_conditional = std::vector< std::vector< std::pair<size_t, double> > >();
    
    // loop through all the sampled character histories for this branch
    for (size_t j = burnin; j < num_sampled_states; ++j)
    {
        if ( verbose == true )
        {
            progress.update( num_finished_nodes * num_sampled_states + num_sampled_states * (j - burnin) / (num_sampled_states - burnin) );
        }
        
        // if necessary, get the sampled tree from the tree trace
        const Tree &sample_tree = (usingTreeTrace()) ? tree_trace.objectAt( j ) : final_summary_tree;
        const TopologyNode& sample_root = sample_tree.getRoot();
        
        try
        {
            if ( usingTreeTrace() == true )
            {
                // check if the clade in the summary tree is also in the sampled tree
                sample_clade_index = sample_root.getCladeIndex( summary_nodes[node_index] );
                
                // and we must also find the trace for this node index
                trace_found = false;
            }
            else
            {
                sample_clade_index = summary_nodes[node_index]->getIndex();
            }
            
        }
        catch(RbException&)
        {
            continue;
        }
        
        bool use_sample = true;
        
        // check if we must condition on the parent's end state
        if ( conditional == true && root == false )
        {
            size_t sample_parent_index;
            if ( usingTreeTrace() == true )
            {
                sample_parent_index = sample_tree.getNode( sample_clade_index ).getParent().getIndex();
                parent_trace_found = false;
            }
            else
            {
                sample_parent_index = summary_nodes[sample_clade_index]->getParent().getIndex();
            }
            
            if ( parent_trace_found == false )
            {
                for (size_t k = 0; k < ancestral_state_traces.size(); ++k)
                {
                    if (ancestral_state_traces[k].getParameterName() == StringUtilities::toString(sample_parent_index + 1))
                    {
                        parent_trace = &ancestral_state_traces[k];
                        parent_trace_found = true;
                        break;
                    }
                }
            }
            assert(parent_trace);
            
            // get the sampled character history for the parent for this iteration
            const std::vector<std::string>& parent_vector = parent_trace->getValues();
            std::string character_history = parent_vector[j];
            
            // parse sampled SIMMAP string
            std::vector< std::pair<size_t, double> > parent_branch_map = parseSIMMAPForNode(character_history, bwd_time);
            
            // finally check against the map state of the parent
            size_t parent_end_state = parent_branch_map[ parent_branch_map.size() - 1 ].first;
            if ( parent_end_state != map_parent_state )
            {
                use_sample = false;
            }
        }
        
        // if necessary find the AncestralStateTrace for the sampled node
        if ( trace_found == false )
        {
            for (size_t k = 0; k < ancestral_state_traces.size(); ++k)
            {
                if (ancestral_state_traces[k].getParameterName() == StringUtilities::toString(sample_clade_index + 1))
                {
                    ancestralstate_trace = &ancestral_state_traces[k];
                    trace_found = true;
                    break;
                }
            }
        }
        assert(ancestralstate_trace);
        
        // get the sampled character history for this iteration
        const std::vector<std::string>& ancestralstate_vector = ancestralstate_trace->getValues();
        std::string character_history = ancestralstate_vector[j];
        
        // parse sampled SIMMAP string
        std::vector< std::pair<size_t, double> > this_branch_map = parseSIMMAPForNode(character_history, bwd_time);
        
        if ( use_sample == true )
        {
            branch_maps_conditional.push_back(this_branch_map);
        }
        branch_maps_all.push_back(this_branch_map);
    }
    
    // now summarize maximum a posteriori states in each time slice along the branch
    
    double branch_len = summary_nodes[node_index]->getBranchLength();
    
    // strings to hold the SIMMAP strings for the summarized character history
    std::string branch_map_history = "}";
    std::string branch_map_history_posteriors = "}";
    std::string branch_map_history_shift_prob = "}";
    
    // loop through each time slice along this branch
    double current_time = 0.0;
    size_t current_dt = 1;
    bool finished_branch = false;
    size_t map_state = 0;
    size_t old_map_state = 0;
//    std::vector<double> old_state_pp = std::vector<double>();
    while ( finished_branch == false )
    {
        std::vector< std::vector< std::pair<size_t, double> > > branch_maps;
        if ( current_dt == 1 && conditional == true && root == false )
        {
            branch_maps = branch_maps_conditional;
        }
        else
        {
            branch_maps = branch_maps_all;
        }
        
        current_time = current_dt * dt;
        if ( current_time >= branch_len )
        {
            current_time = branch_len;
            finished_branch = true;
        }
        
        std::vector<size_t> states = std::vector<size_t>();
        std::vector<double> posteriors = std::vector<double>();
        
        double this_shift_prob = 0.0;
        
        // loop through all sampled character histories for this branch
        // and compile each state sampled in this time slice
        for (size_t k = 0; k < branch_maps.size(); k++)
        {
            double sample_time = 0.0;
            
            bool use_sample = true;
            if ( conditional == true && current_dt > 1 && root == false )
            {
                use_sample = false;
            }
            
            // find this time slice in this sampled character history
            for (size_t l = 0; l < branch_maps[k].size(); l++)
            {
                sample_time += branch_maps[k][l].second;
                
                // if we are conditioning on the previous time slice's MAP state
                // check whether or not we should use this sample by checking
                // the state in the last time slice
                // for current_dt == 1 we have already excluded branch_maps whose
                // parent didn't end with the parent MAP end state
                if ( conditional == true && current_dt > 1 )
                {
                    size_t window_size = 1;
                    if ( current_time - (window_size * dt) <= sample_time || l == branch_maps[k].size() - 1)
                    {
                        if ( branch_maps[k][l].first == old_map_state )
                        {
                            use_sample = true;
                        }
                    }
                }
                
                if ( use_sample == true && sample_time < branch_len )
                {
                    if ( (current_time >= sample_time && ((current_dt-1) * dt) < sample_time) )
                    {
                        if ( current_time != dt || true )
                        {
                            if ( (sample_time+1E-5) < branch_len )
                            {
                                this_shift_prob += 1.0;
                            }
                        }
                    }
                }
                
                if ( use_sample == true && ( current_time <= sample_time || l == branch_maps[k].size() - 1 ) )
                {
                    bool state_found = false;
                    for (size_t m = 0; m < states.size(); m++)
                    {
                        if (states[m] == branch_maps[k][l].first)
                        {
                            state_found = true;
                            posteriors[m] += 1.0;
                            break;
                        }
                    }
                    if (state_found == false)
                    {
                        states.push_back(branch_maps[k][l].first);
                        posteriors.push_back(1.0);
                    }
                    break;
                }
            }
        }
        
        
        // scale posteriors
        double sum_posteriors = 0.0;
        for (size_t i = 0; i < posteriors.size(); i++)
        {
            sum_posteriors += posteriors[i];
        }
        for (size_t i = 0; i < posteriors.size(); i++)
        {
            posteriors[i] /= sum_posteriors;
        }
        this_shift_prob /= double(branch_maps.size());
        
        // find the MAP state for this time slice
        double map_state_pp = 0.0;
//        double shift_prob = 0.0;
        for (size_t k = 0; k < states.size(); k++)
        {
            if (posteriors[k] > map_state_pp)
            {
                map_state = states[k];
                map_state_pp = posteriors[k];
            }
            
//            if ( old_state_pp.size() == states.size() )
//            {
//                shift_prob += old_state_pp[k] * (1-posteriors[k]);
//            }
        }
        old_map_state = map_state;
//        old_state_pp = posteriors;
        
        // now add this time slice to the SIMMAP strings
        
        double time_slice_length = current_time - ((current_dt - 1) * dt);
        
        branch_map_history = "," + StringUtilities::toString(time_slice_length) + branch_map_history;
        branch_map_history_posteriors = "," + StringUtilities::toString(time_slice_length) + branch_map_history_posteriors;
        branch_map_history_shift_prob = "," + StringUtilities::toString(time_slice_length) + branch_map_history_shift_prob;
        
        branch_map_history = StringUtilities::toString(map_state) + branch_map_history;
        branch_map_history_posteriors = StringUtilities::toString(int(map_state_pp * 100)) + branch_map_history_posteriors;
        branch_map_history_shift_prob = StringUtilities::toString(int(this_shift_prob * 100)) + branch_map_history_shift_prob;
        
        if (finished_branch == false )
        {
            branch_map_history = ":" + branch_map_history;
            branch_map_history_posteriors = ":" + branch_map_history_posteriors;
            branch_map_history_shift_prob = ":" + branch_map_history_shift_prob;
        }
        
        current_dt++;
    }
    
    // finally finish the SIMMAP strings for this node
    branch_map_history = "{" + branch_map_history;
    branch_map_history_posteriors = "{" + branch_map_history_posteriors;
    branch_map_history_shift_prob = "{" + branch_map_history_shift_prob;
    map_character_history[node_index] = branch_map_history;
    map_character_history_posteriors[node_index] = branch_map_history_posteriors;
    map_character_history_shift_prob[node_index] = branch_map_history_shift_prob;
    
    // update nodes finished for the progress bar
    num_finished_nodes += 1;
    
    // recurse through tree
    std::vector<int> children_indices = summary_nodes[node_index]->getChildrenIndices();
    for (size_t i = 0; i < children_indices.size(); i++)
    {
        recursivelyCollectCharacterMapSamples(children_indices[i], map_state, false, conditional, final_summary_tree, summary_nodes, map_character_history, map_character_history_posteriors, map_character_history_shift_prob, progress, num_finished_nodes, NUM_TIME_SLICES, verbose);
    }
}


/*
 * This method summarizes marginal, joint and conditional character histories over a give summary tree from a distribution of stochastic character map samples.
 */
Tree* JointAncestralStateTrace::characterMapTree(const Tree &input_summary_tree, int NUM_TIME_SLICES, bool conditional, bool joint, bool verbose)
{
    if ( NUM_TIME_SLICES < 1 )
    {
        throw RbException("The number of time slices must be at least 1.");
    }
    
    std::stringstream ss;
    ss << "Calculating maximum a posteriori character map from " << num_sampled_states << " samples in the character map trace, using a burnin of " << burnin << " samples.\n";
    RBOUT(ss.str());
    
    RBOUT("Compiling character maps...\n");
    
    // allocate memory for the new summary tree
    Tree* final_summary_tree = new Tree( input_summary_tree );
    const std::vector<TopologyNode*> &summary_nodes = final_summary_tree->getNodes();
    
    std::vector<std::string> map_character_history( summary_nodes.size(), "" );
    std::vector<std::string> map_character_history_posteriors( summary_nodes.size(), "" );
    std::vector<std::string> map_character_history_shift_prob( summary_nodes.size(), "" );
    
    bool process_active = true;
    ProgressBar progress = ProgressBar( summary_nodes.size() * num_sampled_states, 0 );
    if ( verbose == true && process_active == true )
    {
        progress.start();
    }
    
    // recurse through summary tree and collect ancestral state samples
    size_t node_index = final_summary_tree->getRoot().getIndex();
    size_t num_finished_nodes = 0;
    recursivelyCollectCharacterMapSamples(node_index, 0, true, conditional, *final_summary_tree, summary_nodes, map_character_history, map_character_history_posteriors, map_character_history_shift_prob, progress, num_finished_nodes, NUM_TIME_SLICES, verbose);
    
    if ( verbose == true && process_active == true )
    {
        progress.finish();
    }
    
    // add all the SIMMAP strings to the tree
    final_summary_tree->clearNodeParameters();
    final_summary_tree->addNodeParameter("map_character_history", map_character_history, false);
    final_summary_tree->addNodeParameter("map_character_history_posteriors", map_character_history_posteriors, false);
    final_summary_tree->addNodeParameter("map_character_history_shift_prob", map_character_history_shift_prob, false);
    
    return final_summary_tree;
    
}


/*
 * Helper function that parses a SIMMAP character history for a single branch.
 * These strings represent character histories for a single branch, either
 *
 * (i) in backward time (frontend: use_simmap_default = "True")
 *     e.g., {2,0.1:1,1.9},
 *     which is read {state_2,time_in_state_2:state_1,time_in_state_1},
 *     where left means young (tips) and right means old (root)
 *
 * or
 *
 * (ii) in forward time (frontend: use_simmap_default = "False")
 *      e.g., {1,1.9:2,0.1},
 *      which is read {state_1,time_in_state_1:state_2,time_in_state_2},
 *      where left means old (root) and right means young (tips)
 *
 *
 * We loop through the string from right to left (in the case of i) or
 * from left to right (in the case of ii), and store events in forward time
 * (root to tip).
 *
 * Returns vector of events: [<state_1, time_in_state_1>, <state_2, time_in_state_2>]
 */
std::vector< std::pair<size_t, double> > JointAncestralStateTrace::parseSIMMAPForNode(std::string character_history,
                                                                                      bool bwd_time) {
    
    boost::trim(character_history);
    
    // Now parse the sampled SIMMAP string:
    bool parsed_first_element = false;
    std::vector< std::pair<size_t, double> > this_branch_map = std::vector< std::pair<size_t, double> >();
    std::pair<size_t, double> this_event = std::pair<size_t, double>();
    std::string state = "";
    std::string time = "";
    size_t str_size = character_history.size();
    size_t pos = 1;
    size_t pos_stop_condition = 0;
    
    if (bwd_time) {
        pos = str_size - 1;
        // pos_stop_condition is 0
    }
    else {
        // pos is 1
        pos_stop_condition = str_size;
    }
    
    while (true) {
        if ( ( (bwd_time && pos == (str_size - 1)) || (!bwd_time && pos == 1) ) &&
              std::string(1, character_history[0]).compare("{") != 0 &&
              std::string(1, character_history[pos]).compare("}") != 0
            ) {
            throw RbException("Error while summarizing character maps: trace does not contain valid SIMMAP string.");
        }

        else if (std::string(1, character_history[pos]).compare(",") == 0 ) {
            parsed_first_element = true;
            
            // if backward in time (right-to-left), when it gets to comma ","
            // it must have read through a time
            if (bwd_time) {
                this_event.second = std::atof( time.c_str() );
            }
            
            // if forward in time (left-to-right), when it gets to comma ","
            // it must have read through a state integer
            else {
                this_event.first = std::atoi( state.c_str() );
            }
        }

        else if ( std::string(1, character_history[pos]).compare(":") == 0 ||
                 pos == pos_stop_condition ) {

            if (bwd_time) {
                this_event.first = std::atoi( state.c_str() );
            }
            
            // forward time
            else {
                this_event.second = std::atof( time.c_str() );
            }
            
            // youngest to oldest if bwd_time
            // (oldest to youngest if forward time!!!)
            this_branch_map.push_back( this_event );
            
            // finished whole stochastic mapping string for this node!
            if (pos == pos_stop_condition) {
                break;
            }
            
            else {
                state = "";
                time = "";
                parsed_first_element = false;
            }
        }

        else if ( parsed_first_element == false ) {
            if (bwd_time) {
                time = std::string(1, character_history[pos]) + time;
            }
            
            // forward time
            else {
                state = state + std::string(1, character_history[pos]);
            }
        }

        // parsed_first_element == true
        else {
            if (bwd_time) {
                state = std::string(1, character_history[pos]) + state;
            }
            
            // forward time
            else {
                time = time + std::string(1, character_history[pos]);
            }
        }
        
        if (bwd_time) { pos--; }
        else { pos++; }
    }
    
    return this_branch_map;
}


/**
 *
 * Summarizes sampled character histories for a vector of stochastic character map traces for each node of a given summary tree.
 *
 */
void JointAncestralStateTrace::summarizeCharacterMaps(Tree input_tree, const path& filename, bool verbose, std::string separator, bool bwd_time)
{
    std::vector<TopologyNode*> summary_nodes;
    bool condition_on_tree = false;
    
    if ( input_tree.getNumberOfTips() != 0 )
    {
        summary_nodes = input_tree.getNodes();
        condition_on_tree = true;
    }
    
    if ( condition_on_tree == false && usingTreeTrace() == false )
    {
        throw RbException("You must enter a tree trace, a tree, or both to summarize stochastic character maps.");
    }
    
    // get the number of nodes
    size_t num_nodes = 0;
    if (condition_on_tree == true)
    {
        num_nodes = input_tree.getNumberOfNodes();
    }
    else
    {
        num_nodes = tree_trace.objectAt(0).getNumberOfNodes();
    }
    
    std::stringstream ss;
    ss << "Compiling stochastic character maps from " << num_sampled_states << " samples in the character map trace, using a burnin of " << burnin << " samples.\n";
    RBOUT(ss.str());
    
    ProgressBar progress = ProgressBar( num_nodes * num_sampled_states, 0 );
    if ( verbose == true )
    {
        progress.start();
    }
    
    // open stream for file output
    createDirectoryForFile( filename );
    std::ofstream out( filename.string() );
    
    // write column headers
    out << "iteration" << separator;
    out << "node_index" << separator;
    out << "branch_start_time" << separator;
    out << "branch_end_time" << separator;
    out << "start_state" << separator;
    out << "end_state" << separator;
    out << "transition_time" << separator;
    out << "transition_type" << separator;
    out << "parent_index" << separator;
    out << "child1_index" << separator;
    out << "child2_index" << std::endl;
    
    // get the iteration trace
    AncestralStateTrace* iteration_trace = NULL;
    bool trace_found = false;
    for (size_t k = 0; k < ancestral_state_traces.size(); k++)
    {
        if (ancestral_state_traces[k].getParameterName() == "Iteration")
        {
            iteration_trace = &ancestral_state_traces[k];
            trace_found = true;
            break;
        }
    }
    
    if ( trace_found == false )
    {
        throw RbException("The ancestral state trace is missing the 'Iteration' column.");
    }
    
    const std::vector<std::string>& iteration_vector = iteration_trace->getValues();
    
    // loop through all nodes
    for (size_t i = 0; i < num_nodes; ++i)
    {
        size_t sample_clade_index = i;
        trace_found = false;
        AncestralStateTrace* character_map_trace = NULL;
        std::vector<std::string> ancestralstate_vector;
        
        // loop through all the stochastic character map samples
        for (size_t j = burnin; j < num_sampled_states; ++j)
        {
            
            if ( verbose == true )
            {
                progress.update( i * num_sampled_states + num_sampled_states * (j - burnin) / (num_sampled_states - burnin) );
            }
            
            // if necessary, get the sampled tree from the tree trace
            const Tree &sample_tree = (usingTreeTrace()) ? tree_trace.objectAt( j ) : input_tree;
            
            try
            {
                if ( condition_on_tree == true )
                {
                    if ( usingTreeTrace() == true )
                    {
                        // check if the clade in the summary tree is also in the sampled tree
                        const TopologyNode& sample_root = sample_tree.getRoot();
                        sample_clade_index = sample_root.getCladeIndex( summary_nodes[i] );
                        
                        // and we must also find the trace for this node index
                        trace_found = false;
                    }
                    else
                    {
                        sample_clade_index = summary_nodes[i]->getIndex();
                    }
                }
            }
            catch(RbException&)
            {
                continue;
            }
            
            // if necessary find the AncestralStateTrace for the sampled node
            if ( trace_found == false )
            {
                for (size_t k = 0; k < ancestral_state_traces.size(); k++)
                {
                    if (ancestral_state_traces[k].getParameterName() == StringUtilities::toString(sample_clade_index + 1))
                    {
                        character_map_trace = &ancestral_state_traces[k];
                        ancestralstate_vector = character_map_trace->getValues();
                        trace_found = true;
                        break;
                    }
                }
            }
            
            // get the iteration
            std::string iteration = iteration_vector[j];
            
            // get the sampled character history for this node for this iteration
            std::string character_history = ancestralstate_vector[j];
            
            // parse sampled SIMMAP string
            std::vector< std::pair<size_t, double> > this_branch_map = parseSIMMAPForNode(character_history, bwd_time);
            
            double start_time = sample_tree.getNode( sample_clade_index ).getAge() + sample_tree.getNode( sample_clade_index ).getBranchLength();
            double end_time = sample_tree.getNode( sample_clade_index ).getAge();
            
            double current_time = start_time;
            size_t current_state;
            size_t end_state;
            if ( this_branch_map.size() > 0 ) {
                if (bwd_time) {
                    current_state = this_branch_map[0].first;
                    end_state = this_branch_map[ this_branch_map.size() - 1 ].first;
                }
                // forward time
                else {
                    end_state = this_branch_map[0].first;
                    current_state = this_branch_map[ this_branch_map.size() - 1 ].first;
                }
            }
            else
            {
                throw RbException()<<"There were no sampled character histories for node "<<(sample_clade_index + 1)<<" in the summary tree.";
            }
            
            // write output if there was no change along the branch
            if (this_branch_map.size() == 1)
            {
                // write node index
                out << iteration << separator;
                
                // write node index
                if (condition_on_tree && usingTreeTrace())
                {
                    out << summary_nodes[i]->getIndex() + 1 << separator;
                }
                else
                {
                    out << sample_clade_index + 1 << separator;
                }
                
                // write branch start/end times
                out << start_time << separator;
                out << end_time << separator;
                
                // write start state
                out << current_state << separator;
                
                // write end state
                out << end_state << separator;
                
                // write transition time
                out << "NA" << separator;
                
                out << "no_change" << separator;
                
                // write parent/child indices
                std::string parent_index = "NA";
                std::string child1_index = "NA";
                std::string child2_index = "NA";
                if (condition_on_tree && usingTreeTrace())
                {
                    if (summary_nodes[i]->isRoot() == false)
                    {
                        parent_index = StringUtilities::toString( summary_nodes[i]->getParent().getIndex() + 1 );
                    }
                    if (summary_nodes[i]->isTip() == false)
                    {
                        child1_index = StringUtilities::toString( summary_nodes[i]->getChild( 0 ).getIndex() + 1 );
                        child2_index = StringUtilities::toString( summary_nodes[i]->getChild( 1 ).getIndex() + 1 );
                    }
                }
                else
                {
                    if (sample_tree.getNode( sample_clade_index ).isRoot() == false)
                    {
                        parent_index = StringUtilities::toString( sample_tree.getNode( sample_clade_index ).getParent().getIndex() + 1 );
                    }
                    if (sample_tree.getNode( sample_clade_index ).isTip() == false)
                    {
                        child1_index = StringUtilities::toString( sample_tree.getNode( sample_clade_index ).getChild( 0 ).getIndex() + 1 );
                        child2_index = StringUtilities::toString( sample_tree.getNode( sample_clade_index ).getChild( 1 ).getIndex() + 1 );
                    }
                }
                out << parent_index << separator;
                out << child1_index << separator << child2_index << std::endl;
            }
            
            // write output for each anagenetic transition along the branch in forward time (towards tips)
            for (size_t k = 1; k < this_branch_map.size(); k++)
            {
                if (this_branch_map[k].first != current_state)
                {
                    // write node index
                    out << iteration << separator;
                    
                    // write node index
                    if (condition_on_tree && usingTreeTrace())
                    {
                        out << summary_nodes[i]->getIndex() + 1 << separator;
                    }
                    else
                    {
                        out << sample_clade_index + 1 << separator;
                    }
                    
                    // write branch start/end times
                    out << start_time << separator;
                    out << end_time << separator;
                    
                    // write start state
                    out << current_state << separator;
                    
                    // write end state
                    out << this_branch_map[k].first << separator;
                    
                    // write transition time
                    out << current_time - this_branch_map[k - 1].second << separator;
                    
                    out << "anagenetic" << separator;
                    
                    // write parent/child indices
                    std::string parent_index = "NA";
                    std::string child1_index = "NA";
                    std::string child2_index = "NA";
                    if (condition_on_tree && usingTreeTrace())
                    {
                        if (summary_nodes[i]->isRoot() == false)
                        {
                            parent_index = StringUtilities::toString( summary_nodes[i]->getParent().getIndex() + 1 );
                        }
                        if (summary_nodes[i]->isTip() == false)
                        {
                            child1_index = StringUtilities::toString( summary_nodes[i]->getChild( 0 ).getIndex() + 1 );
                            child2_index = StringUtilities::toString( summary_nodes[i]->getChild( 1 ).getIndex() + 1 );
                        }
                    }
                    else
                    {
                        if (sample_tree.getNode( sample_clade_index ).isRoot() == false)
                        {
                            parent_index = StringUtilities::toString( sample_tree.getNode( sample_clade_index ).getParent().getIndex() + 1 );
                        }
                        if (sample_tree.getNode( sample_clade_index ).isTip() == false)
                        {
                            child1_index = StringUtilities::toString( sample_tree.getNode( sample_clade_index ).getChild( 0 ).getIndex() + 1 );
                            child2_index = StringUtilities::toString( sample_tree.getNode( sample_clade_index ).getChild( 1 ).getIndex() + 1 );
                        }
                    }
                    out << parent_index << separator;
                    out << child1_index << separator << child2_index << std::endl;
                }
                current_state = this_branch_map[k].first;
                current_time = current_time - this_branch_map[k - 1].second;
            }
            
            // now check this node's children's start states to see if there were any cladogenetic transitions
            //
            // will call parseSIMMAPForNode() on children now
            std::vector<int> children_indices = sample_tree.getNode( sample_clade_index ).getChildrenIndices();
            
            for (int k = 0; k < children_indices.size(); k++)
            {
                size_t child_index = children_indices[k];
                
                for (size_t l = 0; l < ancestral_state_traces.size(); l++)
                {
                    if (ancestral_state_traces[l].getParameterName() == StringUtilities::toString(child_index + 1))
                    {
                        character_map_trace = &ancestral_state_traces[l];
                        break;
                    }
                }
                if (not character_map_trace)
                    throw RbException("Couldn't find character_map_trace!");
                
                // get the sampled character history for the child for this iteration
                const std::vector<std::string>& ancestralstate_vector_child = character_map_trace->getValues();
                std::string character_history_child = ancestralstate_vector_child[j];
                
                // parse sampled SIMMAP string
                std::vector< std::pair<size_t, double> > child_branch_map = parseSIMMAPForNode(character_history_child, bwd_time);
                
                // get child's start state
                size_t child_start_state;
                if (bwd_time) {
                    child_start_state = child_branch_map[0].first;
                }
                // forward time
                // child_branch_map should always be reversed here, but sometimes the
                // stoch map logger fails to reverse the input string of the current function!!!
                else {
                    child_start_state = child_branch_map[(child_branch_map.size()-1)].first;
                }
                
                if (end_state != child_start_state)
                {
                    // write node index
                    out << iteration << separator;
                    
                    // write node index
                    if (condition_on_tree && usingTreeTrace())
                    {
                        out << summary_nodes[i]->getIndex() + 1 << separator;
                    }
                    else
                    {
                        out << sample_clade_index + 1 << separator;
                    }
                    
                    // write branch start/end times
                    out << start_time << separator;
                    out << end_time << separator;
                    
                    // write start state
                    out << end_state << separator;
                    
                    // write end state
                    out << child_start_state << separator;
                    
                    // write transition time
                    out << end_time << separator;
                    
                    out << "cladogenetic" << separator;
                    
                    // write parent/child indices
                    std::string parent_index = "NA";
                    std::string child1_index = "NA";
                    std::string child2_index = "NA";
                    if (condition_on_tree && usingTreeTrace())
                    {
                        if (summary_nodes[i]->isRoot() == false)
                        {
                            parent_index = StringUtilities::toString( summary_nodes[i]->getParent().getIndex() + 1 );
                        }
                        if (summary_nodes[i]->isTip() == false)
                        {
                            child1_index = StringUtilities::toString( summary_nodes[i]->getChild( 0 ).getIndex() + 1 );
                            child2_index = StringUtilities::toString( summary_nodes[i]->getChild( 1 ).getIndex() + 1 );
                        }
                    }
                    else
                    {
                        if (sample_tree.getNode( sample_clade_index ).isRoot() == false)
                        {
                            parent_index = StringUtilities::toString( sample_tree.getNode( sample_clade_index ).getParent().getIndex() + 1 );
                        }
                        if (sample_tree.getNode( sample_clade_index ).isTip() == false)
                        {
                            child1_index = StringUtilities::toString( sample_tree.getNode( sample_clade_index ).getChild( 0 ).getIndex() + 1 );
                            child2_index = StringUtilities::toString( sample_tree.getNode( sample_clade_index ).getChild( 1 ).getIndex() + 1 );
                        }
                    }
                    out << parent_index << separator;
                    out << child1_index << separator << child2_index << std::endl;
                }
            }
        }
    }
    
    out.close();
    
    if ( verbose == true )
    {
        progress.finish();
    }
    
}


/*
 * Split a string of sampled states for multiple sites (e.g. "5,12,3") and return the sample for a single site.
 */
std::string JointAncestralStateTrace::getSiteState( const std::string &site_sample, size_t site )
{
    std::vector<std::string> states;
    std::istringstream ss( site_sample );
    std::string state;
    
    while(std::getline(ss, state, ','))
    {
        states.push_back(state);
    }

    if ( states.size() == 0 )
    {
        throw RbException()<<"Problem while getting state for site '"<<site<<"'.";
    }

    if (site >= states.size())
    {
        site = 0;
    }
    return states[site];
}
