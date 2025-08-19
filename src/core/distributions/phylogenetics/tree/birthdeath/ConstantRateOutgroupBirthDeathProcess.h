#ifndef ConstantRateOutgroupBirthDeathProcess_h
#define ConstantRateOutgroupBirthDeathProcess_h

#include <cstddef>
#include <iosfwd>
#include <vector>

#include "AbstractBirthDeathProcess.h"
#include "TreeChangeEventListener.h"
#include "Clade.h"
#include "RbBitSet.h"
#include "RbException.h"
#include "Taxon.h"
#include "TopologyNode.h"

namespace RevBayesCore {
class DagNode;
class Tree;
template <class valueType> class TypedDagNode;
    
    class ConstantRateOutgroupBirthDeathProcess : public AbstractBirthDeathProcess, TreeChangeEventListener {

    public:
        ConstantRateOutgroupBirthDeathProcess(const TypedDagNode<double> *ra, const TypedDagNode<double> *s, const TypedDagNode<double> *e,
                                              const TypedDagNode<double> *ri, const TypedDagNode<double> *ro,
                                              const std::string &ssi, const std::string &sso,
                                              const std::string &cdt, const std::vector<Taxon> &tn, const std::vector<Taxon> &tno);
        
        // public member functions
        ConstantRateOutgroupBirthDeathProcess*              clone(void) const;
        virtual void                                        fireTreeChangeEvent(const TopologyNode &n, const unsigned& m=0);                                                 //!< The tree has changed and we want
        virtual void                                        setValue(Tree *v, bool f=false );
        
    protected:
        // Parameter management functions
        bool                                                isLnProbabilityNonZero(void);
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                //!< Swap a parameter
        virtual void                                        keepSpecialization(const DagNode* affecter);
        virtual void                                        restoreSpecialization(const DagNode *restorer);
        virtual void                                        touchSpecialization(const DagNode *toucher, bool touchAll);
        void                                                recursivelyFlagNodesDirty(const TopologyNode& n);
        RbBitSet                                            recursivelyUpdateClades(const TopologyNode& node);

        
    private:
        TopologyNode*                                       findMrcaNode( Clade c ) const;
        void                                                initializeBitSets(void);
        bool                                                matchesConstraints(void);
        void                                                recursivelyCollectCladeAges( TopologyNode* node, std::vector<double>& ages, TopologyNode* barrier_node=NULL ) const;
//        void                                                recursivelyCollectIngroupAges( TopologyNode* node, std::vector<double>& ages ) const;
        // helper functions
        double                                              computeLnProbabilityTimes(void) const;                                          //!< Compute the log-transformed probability of the current value.
        double                                              computeLnProbabilityDiversifiedSampling(std::vector<double> ages, double presentTime, double samplingProb, size_t num_taxa) const;
        double                                              lnProbNumTaxa(size_t n, double start, double end, bool MRCA) const              { throw RbException("Cannot compute P(nTaxa)."); }
        double                                              simulateDivergenceTime(double origin, double present) const;                    //!< Simulate a speciation event.
        void                                                simulateOutgroupTree( void );
        double                                              pSurvival(double start, double end) const;                                      //!< Compute the probability of
        double                                              pSurvival(double start, double end, double sampling_prob) const;                                      //!< Compute the probability of survival of the process (without incomplete taxon sampling).
        double												pOne(double t, double b, double d, double r) const;
        double												pZero(double t, double b, double d, double r) const;
        double                                              rateIntegral(double t_low, double t_high) const;
        
        // members
        const TypedDagNode<double>*                         lambda;                                                                         //!< The speciation rate.
        const TypedDagNode<double>*                         mu;                                                                             //!< The extinction rate.
        const TypedDagNode<double>*                         rho_ingroup;                                                                    //!< The sampling probability of extant taxa within the ingroup.
        const TypedDagNode<double>*                         rho_outgroup;                                                                   //!< The sampling probability of extant taxa within the outgroup.
        std::string                                         sampling_strategy_ingroup;
        std::string                                         sampling_strategy_outgroup;
        std::vector<Taxon>                                  taxa_outgroup;
        std::vector<Taxon>                                  taxa_ingroup;
        
        Clade                                               clade_ingroup;
        std::vector<RbBitSet>                               active_clades;
        std::vector<RbBitSet>                               stored_clades;
        std::vector<bool>                                   dirty_nodes;
    };
    
    
}

#endif /* ConstantRateOutgroupBirthDeathProcess_h */
