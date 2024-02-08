#ifndef GeneralizedLineageHeterogeneousBirthDeathSamplingProcess_H
#define GeneralizedLineageHeterogeneousBirthDeathSamplingProcess_H

#include "AbstractHomologousDiscreteCharacterData.h"
#include "TreeDiscreteCharacterData.h"
#include "CladogeneticProbabilityMatrix.h"
#include "RateGenerator.h"
#include "MatrixReal.h"
#include "Simplex.h"
#include "Taxon.h"
#include "Tree.h"
#include "TreeDiscreteCharacterData.h"
#include "TreeChangeEventListener.h"
#include "TypedDagNode.h"
#include "RbVector.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlString.h"
#include "Loader.h"

#include <vector>

namespace RevBayesCore {

    /**
     * @file
     * This is the generalized lineage-heterogenous (/state-dependent) birth-death-sampling process.
     * An interface for tensorphylo.
     *
     * Michael R. May and Xavier Meyer 2020/3/08
     *
     */
    class GeneralizedLineageHeterogeneousBirthDeathSamplingProcess : public TypedDistribution<Tree>, public TreeChangeEventListener, public MemberObject< RbVector<double> > {
        
    public:
        GeneralizedLineageHeterogeneousBirthDeathSamplingProcess(const std::vector<Taxon>&     taxa_,
																 const TypedDagNode<double>*   age_,
																 const std::string&            condition_type_,
																 const TypedDagNode<Simplex >* root_frequency_,
																 const size_t                  num_states_,
																 bool                          use_origin_,
																 bool                          zero_indexed_,
																 size_t                        n_proc_,
																 double                        abs_tol_,
																 double                        rel_tol_);

        // pure virtual member functions
        virtual GeneralizedLineageHeterogeneousBirthDeathSamplingProcess* clone(void) const;
        virtual                                                           ~GeneralizedLineageHeterogeneousBirthDeathSamplingProcess(void);                                                              //!< Virtual destructor
        
        double                                                            computeLnProbability(void);
        void                                                              fireTreeChangeEvent(const TopologyNode &n, const unsigned& m=0);                                                 //!< The tree has changed and we want to know which part.
        const AbstractHomologousDiscreteCharacterData&                    getCharacterData() const;
        virtual void                                                      redrawValue(void);
        virtual void                                                      setValue(Tree *v, bool f=false);                                                                    //!< Set the current value, e.g. attach an observation (clamp)
        
        void                                                              drawStochasticCharacterMap(std::vector<std::string>& character_histories, bool use_simmap_default=true);
        void                                                              drawStochasticCharacterMap(std::vector<std::string>& character_histories, std::vector<double>& branch_lambda, std::vector<double>& branch_mu, std::vector<double>& branch_phi, std::vector<double>& branch_delta, std::vector<long>& num_events);
        void                                                              drawJointConditionalAncestralStates(std::vector<size_t>& startStates, std::vector<size_t>& endStates);

        void                                                              dumpModel(std::string file_name);

        // parameter setters
        void                                                              setLambda(const TypedDagNode< RbVector<double> >* param);
        void                                                              setLambda(const TypedDagNode< RbVector< RbVector<double> > >* param, const TypedDagNode< RbVector<double> >* times);
        void                                                              setMu(const TypedDagNode< RbVector<double> >* param);
        void                                                              setMu(const TypedDagNode< RbVector< RbVector<double> > >* param, const TypedDagNode< RbVector<double> >* times);
        void                                                              setPhi(const TypedDagNode< RbVector<double> >* param);
        void                                                              setPhi(const TypedDagNode< RbVector< RbVector<double> > >* param, const TypedDagNode< RbVector<double> >* times);
        void                                                              setDelta(const TypedDagNode< RbVector<double> >* param);
        void                                                              setDelta(const TypedDagNode< RbVector< RbVector<double> > >* param, const TypedDagNode< RbVector<double> >* times);
        void                                                              setUpsilon(const TypedDagNode< RbVector< RbVector<double> > >* param, const TypedDagNode< RbVector<double> >* times);
        void                                                              setGamma(const TypedDagNode< RbVector< RbVector<double> > >* param, const TypedDagNode< RbVector<double> >* times);
        void                                                              setRho(const TypedDagNode< double >* param);
        void                                                              setRho(const TypedDagNode< RbVector< RbVector<double> > >* param, const TypedDagNode< RbVector<double> >* times);
        void                                                              setXi(const TypedDagNode< RbVector< RbVector<double> > >* param, const TypedDagNode< RbVector<double> >* times);
        void                                                              setEta(const TypedDagNode< double >* param);
        void                                                              setEta(const TypedDagNode< RateGenerator >* param);
        void                                                              setEta(const TypedDagNode< RbVector< RateGenerator > >* param, const TypedDagNode< RbVector<double> >* times);
        void                                                              setOmega(const TypedDagNode< CladogeneticProbabilityMatrix >* param);
        void                                                              setOmega(const TypedDagNode< RbVector< CladogeneticProbabilityMatrix > >* param, const TypedDagNode< RbVector<double> >* times);
        void                                                              setZeta(const TypedDagNode< RbVector< MatrixReal > >* param);

    protected:

        // likelihoods
        double current_ln_prob = 0.0;
        double old_ln_prob = 0.0;
        bool   probability_dirty = true;
        bool   tp_can_reset = true;

        // tensorphylo interface
        size_t                                    n_proc;
        TensorPhylo::DistributionHandlerSharedPtr tp_ptr;

        // simulation functions
        void                                                              simulateTree(void);
        void                                                              simulateSSETree(void);
        void                                                              buildSerialSampledRandomBinaryTree(Tree *psi, std::vector<TopologyNode*> &nodes, const std::vector<double> &ages);
        std::vector<double>                                               simulateCoalescentAges() const;
        void                                                              initializeEmptyCharData();

        // virtual methods that may be overwritten, but then the derived class should call this methods
        virtual void                                                      getAffected(RbOrderedSet<DagNode *>& affected, const DagNode* affecter);                                  //!< get affected nodes
        virtual void                                                      keepSpecialization(const DagNode* affecter);
        virtual void                                                      restoreSpecialization(const DagNode *restorer);
        virtual void                                                      touchSpecialization(const DagNode *toucher, bool touchAll);

        // Parameter management functions. You need to override both if you have additional parameters
        virtual void                                                      swapParameterInternal(const DagNode *oldP, const DagNode *newP);                                    //!< Swap a parameter
        void                                                              executeMethod(const std::string &n, const std::vector<const DagNode*> &args, RbVector<double> &rv) const;
        RevLanguage::RevPtr<RevLanguage::RevVariable>                     executeProcedure(const std::string &name, const std::vector<DagNode *> args, bool &found);
        
        // dirty nodes
        void                                                              resizeVectors(size_t num_nodes);
        void                                                              recursivelyFlagNodeDirty(const TopologyNode& n);

        // updating function
        void                                                              prepareParameters(bool force = false);
        void                                                              updateTree(bool force = false);
        void                                                              updateData(bool force = false);
        void                                                              updateTreeData(bool force = false);
        void                                                              updateRootFrequency(bool force = false);
        void                                                              updateLambda(bool force = false);
        void                                                              updateMu(bool force = false);
        void                                                              updatePhi(bool force = false);
        void                                                              updateDelta(bool force = false);
        void                                                              updateUpsilon(bool force = false);
        void                                                              updateGamma(bool force = false);
        void                                                              updateRho(bool force = false);
        void                                                              updateXi(bool force = false);
        void                                                              updateEta(bool force = false);
        void                                                              updateOmega(bool force = false);
        void                                                              updateZeta(bool force = false);

        std::vector<double>                                               RbToStd(const Simplex &obj);
        std::vector<double>                                               RbToStd(const RbVector<double> &obj);
        std::vector< std::vector<double> >                                RbToStd(const RbVector< RbVector<double>> &obj);
        std::vector< std::vector< std::vector<double> > >                 RbToStd(const RbVector< RateGenerator > &obj);
        std::vector< std::vector< std::vector<double> > >                 RbToStd(const RbVector< MatrixReal > &obj);
        std::vector< std::map< std::vector<unsigned>, double > >          RbToStd(const RbVector< CladogeneticProbabilityMatrix > &obj);
        void                                                              checkTimesAreAscending(const RbVector<double> &obj);

        // taxa
        std::vector<Taxon> taxa;

        // revbayes parameters
		const TypedDagNode< double>*                                     age;
		const std::string&                                               condition_type;
		const size_t                                                     num_states;
		bool                                                             use_origin;
		bool                                                             zero_indexed;

		const TypedDagNode< Simplex >*                                   root_frequency;
		bool                                                             root_frequency_dirty = true;

		const TypedDagNode< RbVector< double> >*                         lambda_const = nullptr;
		const TypedDagNode< RbVector< RbVector<double> > >*              lambda_var = nullptr;
		const TypedDagNode< RbVector< double > >*                        lambda_times = nullptr;
		bool                                                             lambda_dirty = true;

		const TypedDagNode< RbVector< double> >*                         mu_const = nullptr;
		const TypedDagNode< RbVector< RbVector<double> > >*              mu_var = nullptr;
		const TypedDagNode< RbVector< double > >*                        mu_times = nullptr;
		bool                                                             mu_dirty = true;

		const TypedDagNode< RbVector< double> >*                         phi_const = nullptr;
		const TypedDagNode< RbVector< RbVector<double> > >*              phi_var = nullptr;
		const TypedDagNode< RbVector< double > >*                        phi_times = nullptr;
		bool                                                             phi_dirty = true;

		const TypedDagNode< RbVector< double> >*                         delta_const = nullptr;
		const TypedDagNode< RbVector< RbVector<double> > >*              delta_var = nullptr;
		const TypedDagNode< RbVector< double > >*                        delta_times = nullptr;
		bool                                                             delta_dirty = true;

		const TypedDagNode< RbVector< RbVector<double> > >*              upsilon = nullptr;
		const TypedDagNode< RbVector< double > >*                        upsilon_times = nullptr;
		bool                                                             upsilon_dirty = true;

		const TypedDagNode< RbVector< RbVector<double> > >*              gamma = nullptr;
		const TypedDagNode< RbVector< double > >*                        gamma_times = nullptr;
		bool                                                             gamma_dirty = true;

		const TypedDagNode< double >*                                    rho_simple = nullptr;
		const TypedDagNode< RbVector< RbVector<double> > >*              rho = nullptr;
		const TypedDagNode< RbVector< double > >*                        rho_times = nullptr;
		bool                                                             rho_dirty = true;

		const TypedDagNode< RbVector< RbVector<double> > >*              xi = nullptr;
		const TypedDagNode< RbVector< double > >*                        xi_times = nullptr;
		bool                                                             xi_dirty = true;

		const TypedDagNode< double >*                                    eta_simple = nullptr;
		const TypedDagNode< RateGenerator >*                             eta_const = nullptr;
		const TypedDagNode< RbVector< RateGenerator > >*                 eta_var = nullptr;
		const TypedDagNode< RbVector< double > >*                        eta_times = nullptr;
		bool                                                             eta_dirty = true;

		const TypedDagNode< CladogeneticProbabilityMatrix >*             omega_const = nullptr;
		const TypedDagNode< RbVector< CladogeneticProbabilityMatrix > >* omega_var = nullptr;
		const TypedDagNode< RbVector< double > >*                        omega_times = nullptr;
		bool                                                             omega_dirty = true;

		const TypedDagNode< RbVector< MatrixReal > >*                    zeta = nullptr;
		bool                                                             zeta_dirty = true;

		// misc. book-keeping
		bool tree_dirty = true;
        std::vector<bool> dirty_nodes;


    };
    
}

#endif
