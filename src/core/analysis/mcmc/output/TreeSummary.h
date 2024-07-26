#ifndef TreeSummary_H
#define TreeSummary_H

#include "Clade.h"
#include "Trace.h"
#include "Tree.h"

namespace RevBayesCore {

    class MatrixReal;
    class TraceTree;

    class TreeSummary {

    public:
        
        /*
         * This struct represents a value/count pair that is sorted by count
         */
        template <class T>
        struct Sample : public std::pair<T, long>
        {
            Sample(T t, long l) : std::pair<T, long>(t,l) {}

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

        TreeSummary*                               clone(void) const;
        void                                       annotateTree(Tree &inputTree, AnnotationReport report, bool verbose );
        double                                     cladeProbability(const Clade &c, bool verbose);
        MatrixReal                                 computeConnectivity( double credible_interval_size, const std::string& m, bool verbose );
        double                                     computeEntropy( double credible_interval_size, int num_taxa, bool verbose );
        std::vector<double>                        computePairwiseRFDistance( double credible_interval_size, bool verbose );
        std::vector<double>                        computeTreeLengths(void);
        const std::set<Sample<std::string> >&      getTreeSamples(void) const;
        std::vector<Clade>                         getUniqueClades(double ci=0.95, bool non_trivial_only=true, bool verbose=true);
        std::vector<Tree>                          getUniqueTrees(double ci=0.95, bool verbose=true);
        long                                       getTopologyCount(const Tree &t, bool verbose);
        double                                     getTopologyFrequency(const Tree &t, bool verbose);
        bool                                       isClock(void) const;
        bool                                       isCoveredInInterval(const std::string &v, double size, bool verbose);
        bool                                       isCoveredInInterval(const Tree &t, double size, bool verbose);
        bool                                       isDirty(void) const;
        double                                     jointCladeProbability(const RbVector<Clade> &c, bool verbose);
        double                                     maxdiff(bool verbose);
        Tree*                                      mapTree(AnnotationReport report, bool verbose);
        Tree*                                      mccTree(AnnotationReport report, bool verbose);
        Tree*                                      mrTree(AnnotationReport report, double cutoff, bool verbose);
        void                                       printTreeSummary(std::ostream& o, double ci=0.95, bool verbose=true);
        void                                       printCladeSummary(std::ostream& o, double minP=0.05, bool verbose=true);
        long                                       sampleSize(bool post = false) const;
        void                                       setOutgroup(const Clade &c);

    protected:

        Split                                      collectTreeSample(const TopologyNode&, RbBitSet&, std::string, std::map<Split, long>&);
        void                                       enforceNonnegativeBranchLengths(TopologyNode& tree) const;
        TopologyNode*                              findParentNode(TopologyNode&, const Split &, std::vector<TopologyNode*>&, RbBitSet& ) const;
        double                                     jointSplitFrequency(const std::vector<Split>& s) const;
        void                                       mapContinuous(Tree &inputTree, const std::string &n, size_t paramIndex, double hpd, bool np, bool verbose ) const;
        void                                       mapDiscrete(Tree &inputTree, const std::string &n, size_t paramIndex, size_t num, bool np, bool verbose ) const;
        void                                       mapParameters(Tree &inputTree, bool verbose) const;
        long                                       splitCount(const Split &n) const;
        double                                     splitFrequency(const Split &n) const;
        void                                       summarize(bool verbose);

        std::vector<TraceTree* >                   traces;

        bool                                       clock;
        bool                                       rooted;

        bool                                       computed = false;
        std::map<Split, long>                      clade_counts;
        std::set<Sample<Split> >                   clade_samples;
        std::map<Taxon, long >                     sampled_ancestor_counts;
        std::map<std::string, long>                tree_counts;
        std::set<Sample<std::string> >             tree_samples;

        std::map<Split, std::vector<double> >                           clade_ages;
        std::map<Split, std::map<Split, std::vector<double> > >         conditional_clade_ages;
        std::map<std::string, std::map<Split, std::vector<double> > >   tree_clade_ages;

        boost::optional<Clade>                     outgroup;
    };

}


#endif
