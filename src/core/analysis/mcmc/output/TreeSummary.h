#ifndef TreeSummary_H
#define TreeSummary_H

#include "Clade.h"
#include "Trace.h"
#include "Tree.h"

namespace RevBayesCore {

    class MatrixReal;
    class TraceTree;

    class TreeSummary {

        /*
         * This struct represents a value/count pair that is sorted by count
         */
        template <class T>
        struct Sample : public std::pair<T, std::int64_t>
        {
            Sample(T t, std::int64_t l) : std::pair<T, std::int64_t>(t,l) {}

            inline bool operator<(const Sample<T>& rhs) const
            {
                if (this->second == rhs.second)
                    return this->first < rhs.first;
                else
                    return this->second < rhs.second;
            }
        };

        /*
         * This struct represents a tree bipartition (split) that can be rooted or unrooted
         */
        struct Split : public std::pair<RbBitSet, std::set<Taxon> >
        {
            Split( RbBitSet b, std::set<Taxon> m, bool r) : std::pair<RbBitSet, std::set<Taxon> >( !r && b[0] ? ~b : b, m) {}

            inline bool operator()(const Sample<Split>& s)
            {
                return (*this) == s.first;
            }
        };

    public:

        /*
         * This struct determines which annotations are reported in the summary tree
         */
        struct AnnotationReport
        {
            bool clade_probs = true;
            bool conditional_clade_ages = false;
            bool conditional_clade_probs = false;
            bool conditional_tree_ages = false;
            bool MAP_parameters = false;
            bool node_ages = true;
            bool mean_node_ages = true;
            double node_ages_HPD = 0.95;
            bool sampled_ancestor_probs = true;
            bool force_positive_branch_lengths = false;
            bool use_outgroup = false;  // Is this unused?
        };


        /*
         * Declaration of the TreeTrace class
         */
        TreeSummary( TraceTree* t, bool c = true );
        TreeSummary( std::vector<TraceTree* > t, bool c = true );
        virtual ~TreeSummary(){}

        TreeSummary*                                   clone(void) const;
        void                                           annotateTree(Tree &inputTree, AnnotationReport report, bool verbose);
        double                                         cladeProbability(const Clade &c, bool verbose);
        double                                         computeEntropy(double credible_interval_size, bool verbose);
        std::vector<double>                            computePairwiseRFDistance(double credible_interval_size, bool verbose);
        std::vector<double>                            computeTreeLengths(void);
        std::vector< std::pair<Tree, std::int64_t> >   getCredibleSetOfTrees(double credible_interval_size, bool verbose=true);
        std::int64_t                                   getTopologyCount(const Tree &t, bool verbose);
        double                                         getTopologyFrequency(const Tree &t, bool verbose);
        std::vector<Clade>                             getUniqueClades(double min_clade_probability=0.05, bool non_trivial_only=true, bool verbose=true);
        std::vector<Tree>                              getUniqueTrees(double credible_interval_size=0.95, bool verbose=true);
        bool                                           isClock(void) const;
        bool                                           isCoveredInInterval(const std::string &v, double credible_interval_size, bool verbose);
        bool                                           isCoveredInInterval(const Tree &t, double credible_interval_size, bool verbose);
        bool                                           isDirty(void) const;
        double                                         jointCladeProbability(const RbVector<Clade> &c, bool verbose);
        double                                         maxdiff(bool verbose);
        Tree*                                          mapTree(AnnotationReport report, bool verbose);
        Tree*                                          mccTree(AnnotationReport report, bool verbose);
        Tree*                                          mrTree(AnnotationReport report, double cutoff, bool verbose);
        void                                           printCladeSummary(std::ostream& o, double min_clade_probability=0.05, bool verbose=true);
        void                                           printTreeSummary(std::ostream& o, double credible_interval_size=0.95, bool verbose=true);
        std::int64_t                                   sampleSize(bool post = false) const;
        void                                           setOutgroup(const Clade &c);

    protected:

        Split                                          collectTreeSample(const TopologyNode&, RbBitSet&, std::string, std::map<Split, std::int64_t>&);
        void                                           enforceNonnegativeBranchLengths(TopologyNode& tree) const;
        TopologyNode*                                  findParentNode(TopologyNode&, const Split &, std::vector<TopologyNode*>&, RbBitSet& ) const;
        void                                           mapContinuous(Tree &inputTree, const std::string &n, size_t paramIndex, double hpd, bool np, bool verbose ) const;
        void                                           mapDiscrete(Tree &inputTree, const std::string &n, size_t paramIndex, size_t num, bool np, bool verbose ) const;
        void                                           mapParameters(Tree &inputTree, bool verbose) const;
        std::int64_t                                   splitCount(const Split &n) const;
        double                                         splitFrequency(const Split &n) const;
        void                                           summarize(bool verbose);

        std::vector<TraceTree* >                       traces;

        bool                                           clock;
        bool                                           rooted;

        bool                                           computed = false;
        std::map<Split, std::int64_t>                  clade_counts;
        std::set<Sample<Split> >                       clade_samples;
        std::map<Taxon, std::int64_t >                 sampled_ancestor_counts;
        std::map<std::string, std::int64_t>            tree_counts;
        std::set<Sample<std::string> >                 tree_samples;

        std::map<Split, std::vector<double> >                           clade_ages;
        std::map<Split, std::map<Split, std::vector<double> > >         conditional_clade_ages;
        std::map<std::string, std::map<Split, std::vector<double> > >   tree_clade_ages;

        boost::optional<Clade>                         outgroup;
    };

}


#endif
