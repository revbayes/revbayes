#ifndef PhyloBrownianRegressionProcess_H
#define PhyloBrownianRegressionProcess_H

#include "AbstractPhyloBrownianProcess.h"
#include "TreeChangeEventListener.h"

namespace RevBayesCore {
    
    /**
     * @brief A Brownian motion with regression (predictor variables)
     *
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2015-01-23, version 1.0
     */
    class PhyloBrownianRegressionProcess : public AbstractPhyloBrownianProcess, public TreeChangeEventListener {
        
    public:
        // Note, we need the size of the alignment in the constructor to correctly simulate an initial state
        PhyloBrownianRegressionProcess(const TypedDagNode<Tree> *t, const ContinuousCharacterData& pred );
        virtual                                                            ~PhyloBrownianRegressionProcess(void);                                                              //!< Virtual destructor
        
        // public member functions
        // pure virtual
        virtual PhyloBrownianRegressionProcess*                             clone(void) const;                                                                      //!< Create an independent clone
        
        // non-virtual
        void                                                                fireTreeChangeEvent(const TopologyNode &n, const unsigned& m=0);                                             //!< The tree has changed and we want to know which part.
        double                                                              computeLnProbability(void);
        void                                                                setMeanPredictor(const TypedDagNode< double >* r);
        void                                                                setMeanPredictor(const TypedDagNode< RbVector< double > >* r);
        void                                                                setSlope(const TypedDagNode< double >* s);
        void                                                                setSlope(const TypedDagNode< RbVector< double > >* s);
        
    protected:
        
        // virtual methods that may be overwritten, but then the derived class should call this methods
        virtual void                                                        keepSpecialization(const DagNode* affecter);
        void                                                                recursiveComputeLnProbability( const TopologyNode &node, size_t node_index );
        void                                                                recursivelyFlagNodeDirty(const TopologyNode& n);
        void                                                                resetValue( void );
        virtual void                                                        restoreSpecialization(const DagNode *restorer);
        std::vector<double>                                                 simulateRootCharacters(size_t n);
        double                                                              sumRootLikelihood(void);
        virtual void                                                        touchSpecialization(const DagNode *toucher, bool touchAll);

        // Parameter management functions.
        virtual void                                                        swapParameterInternal(const DagNode *oldP, const DagNode *newP);                         //!< Swap a parameter

        // the likelihoods
        std::vector<std::vector<std::vector<double> > >                     partial_likelihoods;
        std::vector<std::vector<std::vector<double> > >                     contrasts;
        std::vector<std::vector<double> >                                   contrast_uncertainty;
        std::vector<std::vector<std::vector<double> > >                     contrast_uncertainty_per_site;
        std::vector<size_t>                                                 active_likelihood;
        
        // convenience variables available for derived classes too
        std::vector<bool>                                                   changed_nodes;
        std::vector<bool>                                                   dirty_nodes;

        bool                                                                use_missing_data;
        
        
    private:
        
        double                                                              getSlope(size_t i) const;
        double                                                              getMeanPredictor(size_t i) const;
        double                                                              getPredictor(size_t n, size_t i) const;
        
        const TypedDagNode< double >*                                       single_mean_predictor;
        const TypedDagNode< RbVector< double > >*                           multiple_mean_predictor;
        const TypedDagNode< double >*                                       single_slope;
        const TypedDagNode< RbVector< double > >*                           multiple_slope;
        
        ContinuousCharacterData                                             predictors;
        size_t                                                              num_predictors;
        
    };
    
}


#endif
