#ifndef PopGenInfiniteSites_H
#define PopGenInfiniteSites_H

#include "AbstractHomologousDiscreteCharacterData.h"
#include "BinaryState.h"
#include "Tree.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"


namespace RevBayesCore {
    
    /**
     * @brief Homogeneous distribution of character state evolution along a tree class (PhyloCTMC).
     *
     *
     *
     */
    class PopGenInfiniteSites : public TypedDistribution< AbstractHomologousDiscreteCharacterData > {
        
    public:
        // Note, we need the size of the alignment in the constructor to correctly simulate an initial state
        PopGenInfiniteSites(const TypedDagNode<Tree> *tree, const TypedDagNode< RbVector<Tree> > *trees, const TypedDagNode<double> *rate, const TypedDagNode< RbVector<double> > *rates, size_t ns );                                                                                //!< Copy constructor
        virtual                                                            ~PopGenInfiniteSites(void);                                                     //!< Virtual destructor
        
        // public member functions
        // pure virtual
        virtual PopGenInfiniteSites*                                        clone(void) const;                                                                      //!< Create an independent clone
        
        // non-virtual
        double                                                              computeLnProbability(void);
//        void                                                                executeMethod(const std::string &n, const std::vector<const DagNode*> &args, RbVector<double> &rv) const;     //!< Map the member methods to internal function calls
//        void                                                                executeMethod(const std::string &n, const std::vector<const DagNode*> &args, MatrixReal &rv) const;     //!< Map the member methods to internal function calls
        
        void                                                                redrawValue(void);

        
    protected:
        
        // Parameter management functions.
        void                                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                             //!< Swap a parameter
        
        
        // virtual methods that may be overwritten, but then the derived class should call this methods
//        void                                                                keepSpecialization(DagNode* affecter);
//        void                                                                restoreSpecialization(DagNode *restorer);
//        void                                                                touchSpecialization(DagNode *toucher, bool touchAll);
        
        // virtual methods that you may want to overwrite
        void                                                                compress(void);
        
        
        // members
//        double                                                              ln_prob;
//        double                                                              stored_ln_prob;
//        size_t                                                              num_nodes;
        size_t                                                              num_sites;
        
        // the likelihoods
        
        
        // the data
//        std::vector<size_t>                                                 pattern_counts;
        std::vector<bool>                                                   site_invariant;
        bool                                                                compressed;
        std::vector<size_t>                                                 site_pattern;    // an array that keeps track of which pattern is used for each site
        
        // flags for likelihood recomputation
        bool                                                                touched;
        
        // members
        const TypedDagNode< double >*                                       homogeneous_mutation_rates;
        const TypedDagNode< RbVector< double > >*                           heterogeneous_mutation_rates;
        const TypedDagNode< Tree >*                                         homogeneous_tau;
        const TypedDagNode< RbVector< Tree > >*                             heterogeneous_tau;
        
        
        
    };
    
}


#endif /* PopGenInfinitesSites_hpp */
