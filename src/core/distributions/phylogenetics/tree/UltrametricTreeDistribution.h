#ifndef UltrametricTreeDistribution_H
#define UltrametricTreeDistribution_H

#include "TraceTree.h"
#include "Tree.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore {
    
    class Clade;
    
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
    
//    /*
//     * This struct represents a tree bipartition (split) that can be rooted or unrooted
//     */
//    struct Split : public std::pair<RbBitSet, std::set<Taxon> >
//    {
//        Split( RbBitSet b, std::set<Taxon> m) : std::pair<RbBitSet, std::set<Taxon> >( b[0] ? ~b : b, m) {}
//
//        inline bool operator()(const Sample<Split>& s)
//        {
//            return (*this) == s.first;
//        }
//    };
    
    /*
     * This struct represents a tree bipartition (split) that can be rooted or unrooted
     */
    struct Split : public RbBitSet
    {
        Split( RbBitSet b ) : RbBitSet( b ) {}
        
        inline bool operator()(const Sample<Split>& s)
        {
            return (*this) == s.first;
        }
    };
    
    class UltrametricTreeDistribution : public TypedDistribution<Tree>, public MemberObject< RbVector<double> > {
        
    public:
        
        enum MEAN { ARITHMETIC, HARMONIC };
        
        UltrametricTreeDistribution(TypedDistribution<Tree>* tp, TypedDistribution<double>* rp, TypedDagNode<double> *ra, TypedDagNode<double> *rbf, const TraceTree &tt, Trace<double>* d = NULL, MEAN m = ARITHMETIC);   //!< Constructor
        UltrametricTreeDistribution(const UltrametricTreeDistribution &d);                                              //!< Copy-constructor
        virtual                                            ~UltrametricTreeDistribution(void);                          //!< Virtual destructor
        
        
        UltrametricTreeDistribution&                        operator=(const UltrametricTreeDistribution &d);            //!< Assignment operator

        // public member functions
        UltrametricTreeDistribution*                        clone(void) const;                                          //!< Create an independent clone
        double                                              computeLnProbability(void);                                 //!< Compute ln prob of current value
        void                                                executeMethod(const std::string &n, const std::vector<const DagNode*> &args, RbVector<double> &rv) const;     //!< Map the member methods to internal function calls
        void                                                redrawValue(void);                                          //!< Draw a new random value from distribution
        
        
    protected:
        // Parameter management functions
        virtual void                                        getAffected(RbOrderedSet<DagNode *>& affected, const DagNode* affecter);                                      //!< get affected nodes
        virtual void                                        keepSpecialization(const DagNode* affecter);
        virtual void                                        restoreSpecialization(const DagNode *restorer);
        virtual void                                        touchSpecialization(const DagNode *toucher, bool touchAll);
        
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Swap a parameter
        
    private:
        
        // helper functions
        void                                                attachTimes(Tree *psi, std::vector<TopologyNode *> &tips, size_t index, const std::vector<double> &times, double T);
        void                                                buildRandomBinaryHistory(std::vector<TopologyNode *> &tips);
        Split                                               collectTreeSample(const TopologyNode& n, RbBitSet& in, std::map<Split, double>& bl);
        Split                                               collectSplits(const TopologyNode& n, RbBitSet& in, std::vector<Split>& s) const;
        double                                              computeBranchRateLnProbability(const Tree &x, const std::string &newick, const std::vector<Split> &s, size_t index) const;
        void                                                computeBranchRates(const Tree &x, const std::string &newick, const std::vector<Split> &s, size_t index, std::vector<double> &rates) const;
        void                                                prepareTreeSamples(const std::vector<Tree>& trees);
        void                                                simulateTree(void);
        
        // members
        TypedDistribution<Tree>*                            tree_prior;
        TypedDistribution<double>*                          rate_prior;
        const TypedDagNode<double>*                         root_age;
        const TypedDagNode<double>*                         root_branch_fraction;
        std::vector<Tree>                                   trees;
        Trace<double>*                                      sample_prior_density;
        
        std::vector<std::string>                            trees_newick;
        std::map<std::string, std::vector<size_t> >         topology_indices;
        std::vector<std::map<Split, double> >               tree_branch_lengths;
        std::string                                         outgroup;
        
        size_t                                              num_samples;
        size_t                                              sample_block_start;
        size_t                                              sample_block_end;
        size_t                                              sample_block_size;

#ifdef RB_MPI
        std::vector<size_t>                                 pid_per_sample;
#endif

        std::vector<double>                                 ln_probs;
        size_t                                              num_taxa;
//        std::vector<Taxon>                                  taxa;
        MEAN                                                mean_method;
    };
    
}

#endif

