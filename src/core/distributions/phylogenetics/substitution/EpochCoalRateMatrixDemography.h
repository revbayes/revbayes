#ifndef EpochCoalRateMatrixDemography_H
#define EpochCoalRateMatrixDemography_H

#include "RateMatrix_BinaryMutationCoalescent.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"


namespace RevBayesCore {
    
    /**
     * @brief Homogeneous distribution of character state evolution along a tree class (PhyloCTMC).
     *
     *
     *
     */
    class EpochCoalRateMatrixDemography : public TypedDistribution< RbVector<double> >, public MemberObject< RbVector<double> > {
        
    public:
        enum CODING { ALL, NO_MONOMORPHIC, NO_SINGLETONS };
        
        // Note, we need the size of the alignment in the constructor to correctly simulate an initial state
        EpochCoalRateMatrixDemography(const TypedDagNode< RbVector<double> > *ne,
                            const TypedDagNode< RbVector<double> > *et,
                            const TypedDagNode< RbVector<double> > *mu,
                            const TypedDagNode< Simplex >* asfs,
                            long vps,
                            long n_sites,
                            long n_ind,
                            bool f=false,
                            CODING cod=ALL);
        virtual                                            ~EpochCoalRateMatrixDemography(void);                                                     //!< Virtual destructor
        
        // public member functions

        // pure virtual
        virtual EpochCoalRateMatrixDemography*             clone(void) const;                                                                      //!< Create an independent clone
        
        // non-virtual
        double                                             computeLnProbability(void);
        
        void                                               redrawValue(void);

        
    protected:
        
        // Parameter management functions.
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                             //!< Swap a parameter
        void                                                executeMethod(const std::string &n, const std::vector<const DagNode*> &args, RbVector<double> &rv) const;
        
        
    private:
        
        void                                                initialize( void );

        // parameters (estimated)
        const TypedDagNode< Simplex >*                      ancestral_SFS;
        const TypedDagNode< RbVector< double > >*           ne;
        const TypedDagNode< RbVector< double > >*           epoch_times;
        const TypedDagNode< RbVector< double > >*           mu;
        
        std::vector<RateMatrix_BinaryMutationCoalescent>    rate_matrices;

        // members (fixed)
        long                                                num_states;
        long                                                num_sites;
        long                                                num_individuals;
        bool                                                folded;
        CODING                                              coding;

        
    };
    
}


#endif /* PopGenInfinitesSites_hpp */
