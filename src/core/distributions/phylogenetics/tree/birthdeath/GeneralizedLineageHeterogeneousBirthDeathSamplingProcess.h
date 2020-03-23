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
																 size_t                        n_proc_);

        // pure virtual member functions
        virtual GeneralizedLineageHeterogeneousBirthDeathSamplingProcess* clone(void) const;
        virtual                                                           ~GeneralizedLineageHeterogeneousBirthDeathSamplingProcess(void);                                                              //!< Virtual destructor
        
        double                                                            computeLnProbability(void);
        void                                                              fireTreeChangeEvent(const TopologyNode &n, const unsigned& m=0);                                                 //!< The tree has changed and we want to know which part.
        const AbstractHomologousDiscreteCharacterData&                    getCharacterData() const;
        virtual void                                                      redrawValue(void);
        virtual void                                                      setValue(Tree *v, bool f=false);                                                                    //!< Set the current value, e.g. attach an observation (clamp)
        
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
        double current_ln_prob;
        double old_ln_prob;
        bool   probability_dirty;

        // tensorphylo interface
        size_t                                    n_proc;
        TensorPhylo::DistributionHandlerSharedPtr tp_ptr;

        // simulation functions
        void                                                              simulateTree(void);
        void                                                              buildSerialSampledRandomBinaryTree(Tree *psi, std::vector<TopologyNode*> &nodes, const std::vector<double> &ages);
        std::vector<double>                                               simulateCoalescentAges() const;
        void                                                              initializeEmptyCharData();

        // virtual methods that may be overwritten, but then the derived class should call this methods
        virtual void                                                      getAffected(RbOrderedSet<DagNode *>& affected, DagNode* affecter);                                  //!< get affected nodes
        virtual void                                                      keepSpecialization(DagNode* affecter);
        virtual void                                                      restoreSpecialization(DagNode *restorer);
        virtual void                                                      touchSpecialization(DagNode *toucher, bool touchAll);

        // Parameter management functions. You need to override both if you have additional parameters
        virtual void                                                      swapParameterInternal(const DagNode *oldP, const DagNode *newP);                                    //!< Swap a parameter
        void                                                              executeMethod(const std::string &n, const std::vector<const DagNode*> &args, RbVector<double> &rv) const;
        RevLanguage::RevPtr<RevLanguage::RevVariable>                     executeProcedure(const std::string &name, const std::vector<DagNode *> args, bool &found);
        
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

        // taxa
        std::vector<Taxon> taxa;

        // revbayes parameters
		const TypedDagNode< double>*                                     age;
		const std::string&                                               condition_type;
		const size_t                                                     num_states;
		const TypedDagNode< Simplex >*                                   root_frequency;

		const TypedDagNode< RbVector< double> >*                         lambda_const;
		const TypedDagNode< RbVector< RbVector<double> > >*              lambda_var;
		const TypedDagNode< RbVector< double > >*                        lambda_times;

		const TypedDagNode< RbVector< double> >*                         mu_const;
		const TypedDagNode< RbVector< RbVector<double> > >*              mu_var;
		const TypedDagNode< RbVector< double > >*                        mu_times;

		const TypedDagNode< RbVector< double> >*                         phi_const;
		const TypedDagNode< RbVector< RbVector<double> > >*              phi_var;
		const TypedDagNode< RbVector< double > >*                        phi_times;

		const TypedDagNode< RbVector< double> >*                         delta_const;
		const TypedDagNode< RbVector< RbVector<double> > >*              delta_var;
		const TypedDagNode< RbVector< double > >*                        delta_times;

		const TypedDagNode< RbVector< RbVector<double> > >*              upsilon;
		const TypedDagNode< RbVector< double > >*                        upsilon_times;

		const TypedDagNode< RbVector< RbVector<double> > >*              gamma;
		const TypedDagNode< RbVector< double > >*                        gamma_times;

		const TypedDagNode< double >*                                    rho_simple;
		const TypedDagNode< RbVector< RbVector<double> > >*              rho;
		const TypedDagNode< RbVector< double > >*                        rho_times;

		const TypedDagNode< RbVector< RbVector<double> > >*              xi;
		const TypedDagNode< RbVector< double > >*                        xi_times;

		const TypedDagNode< double >*                                    eta_simple;
		const TypedDagNode< RateGenerator >*                             eta_const;
		const TypedDagNode< RbVector< RateGenerator > >*                 eta_var;
		const TypedDagNode< RbVector< double > >*                        eta_times;

		const TypedDagNode< CladogeneticProbabilityMatrix >*             omega_const;
		const TypedDagNode< RbVector< CladogeneticProbabilityMatrix > >* omega_var;
		const TypedDagNode< RbVector< double > >*                        omega_times;

		const TypedDagNode< RbVector< MatrixReal > >*                    zeta;
		bool                                                             use_origin;

		// book-keeping
		bool tree_dirty;
		bool root_freq_dirty;
		bool lambda_dirty;
		bool mu_dirty;
		bool phi_dirty;
		bool delta_dirty;
		bool upsilon_dirty;
		bool gamma_dirty;
		bool rho_dirty;
		bool xi_dirty;
		bool eta_dirty;
		bool omega_dirty;
		bool zeta_dirty;

    };
    
}

#endif
