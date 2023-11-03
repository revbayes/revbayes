#ifndef JointAncestralStateTrace_H
#define JointAncestralStateTrace_H

#include <cstddef>
#include <iosfwd>
#include <utility>
#include <vector>

#include "TraceTree.h"
#include "Cloneable.h"
#include "Trace.h"

namespace RevBayesCore {
class ProgressBar;
class TopologyNode;
class Tree;

    class JointAncestralStateTrace : public Cloneable {
        
    public:
        
        JointAncestralStateTrace(std::vector<AncestralStateTrace> at, TraceTree tt);
        
        JointAncestralStateTrace*                   clone(void) const;

        Tree*                                       ancestralStateTree(const Tree &inputTree, std::string summary_stat, size_t num_states, int site, bool conditional, bool joint, bool verbose);
        Tree*                                       characterMapTree(const Tree &input_summary_tree, int NUM_TIME_SLICES = 500, bool conditional = false, bool joint = false, bool verbose = true);
        Tree*                                       cladoAncestralStateTree(const Tree &inputTree, std::string summary_stat, int site, bool conditional, bool joint, bool verbose);
        void                                        summarizeCharacterMaps(Tree inputTree, const path& filename, bool verbose, std::string separator);

        size_t                                      getBurnin() const { return burnin; }
        void                                        setBurnin(size_t b);

    private:

        void                                        collectJointAncestralStateSamples(int site, Tree &final_summary_tree, const std::vector<TopologyNode*> &summary_nodes, std::vector<std::vector<double> > &pp_end, std::vector<std::vector<double> > &pp_start, std::vector<std::vector<std::string> > &end_states, std::vector<std::vector<std::string> > &start_states, bool clado, ProgressBar &progress, bool verbose);
        std::string                                 getSiteState( const std::string &site_sample, size_t site );
        void                                        recursivelyCollectAncestralStateSamples(size_t node_index, std::string map_parent_state, bool root, bool conditional, int site, Tree &final_summary_tree, const std::vector<TopologyNode*> &summary_nodes, std::vector<std::vector<double> > &pp_end, std::vector<std::vector<double> > &pp_start, std::vector<double> &pp_clade, std::vector<std::vector<std::string> > &end_states, std::vector<std::vector<std::string> > &start_states, bool clado, ProgressBar &progress, size_t &num_finished_nodes, bool verbose);
        void                                        recursivelyCollectCharacterMapSamples(size_t node_index, size_t map_parent_state, bool root, bool conditional, Tree &final_summary_tree, const std::vector<TopologyNode*> &summary_nodes, std::vector<std::string> &map_character_history, std::vector<std::string> &map_character_history_posteriors, std::vector<std::string> &map_character_history_shift_prob, ProgressBar &progress, size_t &num_finished_nodes, int NUM_TIME_SLICES, bool verbose);
        void                                        computeMarginalCladogeneticStateProbs(std::vector<double> pp, std::vector<std::string> states, std::vector<double>& best_pp, std::vector<std::string>& best_states);

        std::vector< std::pair<size_t, double> >    parseSIMMAPForNode(std::string character_history);

        bool                                        usingTreeTrace(void) const { return tree_trace.size() > 0; };

        std::vector<AncestralStateTrace>            ancestral_state_traces;
        size_t                                      num_sampled_states;
        TraceTree                                   tree_trace;

        size_t                                      burnin;

    };
    

} //end namespace RevBayesCore


#endif
