#ifndef PhyloBrownianCharacterHistoryProcess_H
#define PhyloBrownianCharacterHistoryProcess_H

#include "AbstractPhyloContinuousCharacterHistoryProcess.h"
#include "TreeChangeEventListener.h"

namespace RevBayesCore {
    
    /**
     * @brief Homogeneous distribution of character state evolution along a tree class (PhyloCTMC).
     *
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2015-01-23, version 1.0
     */
    class PhyloBrownianCharacterHistoryProcess : public AbstractPhyloContinuousCharacterHistoryProcess, public TreeChangeEventListener {
        
    public:
        // Note, we need the size of the alignment in the constructor to correctly simulate an initial state
        PhyloBrownianCharacterHistoryProcess(const TypedDagNode<Tree> *t, size_t nSites);
        virtual                                                            ~PhyloBrownianCharacterHistoryProcess(void);                                                         //!< Virtual destructor
        
        // public member functions
        // pure virtual
        virtual PhyloBrownianCharacterHistoryProcess*                       clone(void) const;                                                                      //!< Create an independent clone
        
        // non-virtual
        void                                                                fireTreeChangeEvent(const TopologyNode &n, const unsigned& m=0);                                             //!< The tree has changed and we want to know which part.
        double                                                              computeLnProbability(void);
        
    protected:
        
        // virtual methods that may be overwritten, but then the derived class should call this methods
        void                                                                computeNodeLnProbability( const TopologyNode &node, size_t node_index );
        virtual void                                                        keepSpecialization(DagNode* affecter);
        void                                                                flagNodeDirty(const TopologyNode& n);
//        void                                                                recursiveComputeLnProbability( const TopologyNode &node, size_t node_index );
        void                                                                resetValue( void );
        virtual void                                                        restoreSpecialization(DagNode *restorer);
        void                                                                simulateInternalNodeStates( void );
        virtual void                                                        simulateRecursively(const TopologyNode& node, std::vector< ContinuousTaxonData > &t);
        virtual void                                                        simulateRecursivelyInternal(const TopologyNode& node, std::vector< ContinuousTaxonData > &t);
        virtual void                                                        touchSpecialization(DagNode *toucher, bool touchAll);

        // Parameter management functions.
        virtual void                                                        swapParameterInternal(const DagNode *oldP, const DagNode *newP);                         //!< Swap a parameter

        // the likelihoods
        std::vector<std::vector<std::vector<double> > >                     partial_likelihoods;
        std::vector<size_t>                                                 active_likelihood;
        
        // convenience variables available for derived classes too
        std::vector<bool>                                                   changed_nodes;
        std::vector<bool>                                                   dirty_nodes;

    private:
                
    };
    
}


#endif
