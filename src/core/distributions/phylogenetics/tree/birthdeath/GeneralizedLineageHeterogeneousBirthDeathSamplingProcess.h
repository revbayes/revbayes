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
        GeneralizedLineageHeterogeneousBirthDeathSamplingProcess(
			const TypedDagNode<double>*                                     age_,
			const std::string&                                              condition_type_,
			const TypedDagNode<Simplex >*                                   root_frequency_,
			const TypedDagNode<RbVector< RbVector< double > > >*            lambda_,
			const TypedDagNode<RbVector< double > >*                        lambda_times_,
			const TypedDagNode<RbVector< RbVector< double > > >*            mu_,
			const TypedDagNode<RbVector< double > >*                        mu_times_,
			const TypedDagNode<RbVector< RbVector< double > > >*            phi_,
			const TypedDagNode<RbVector< double > >*                        phi_times_,
			const TypedDagNode<RbVector< RbVector< double > > >*            delta_,
			const TypedDagNode<RbVector< double > >*                        delta_times_,
			const TypedDagNode<RbVector< RbVector< double > > >*            upsilon_,
			const TypedDagNode<RbVector< double > >*                        upsilon_times_,
			const TypedDagNode<RbVector< RbVector< double > > >*            gamma_,
			const TypedDagNode<RbVector< double > >*                        gamma_times_,
			const TypedDagNode<RbVector< RbVector< double > > >*            rho_,
			const TypedDagNode<RbVector< double > >*                        rho_times_,
			const TypedDagNode<RbVector< RbVector< double > > >*            xi_,
			const TypedDagNode<RbVector< double > >*                        xi_times_,
			const TypedDagNode<RbVector< RateGenerator > >*                 eta_,
			const TypedDagNode<RbVector< double > >*                        eta_times_,
			const TypedDagNode<RbVector< CladogeneticProbabilityMatrix > >* omega_,
			const TypedDagNode<RbVector< double > >*                        omega_times_,
			const TypedDagNode<RbVector< MatrixReal > >*                    zeta_,
			bool                                                            use_origin_);
        
        // pure virtual member functions
        virtual GeneralizedLineageHeterogeneousBirthDeathSamplingProcess* clone(void) const;
        virtual                                                           ~GeneralizedLineageHeterogeneousBirthDeathSamplingProcess(void);                                                              //!< Virtual destructor
        
        double                                                            computeLnProbability(void);
        void                                                              fireTreeChangeEvent(const TopologyNode &n, const unsigned& m=0);                                                 //!< The tree has changed and we want to know which part.
        const AbstractHomologousDiscreteCharacterData&                    getCharacterData() const;
        virtual void                                                      redrawValue(void);
        virtual void                                                      setValue(Tree *v, bool f=false);                                                                    //!< Set the current value, e.g. attach an observation (clamp)
        
    protected:
        
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

        // revbayes parameters
		const TypedDagNode< double>*                                     age;
		const std::string&                                               condition_type;
		const TypedDagNode< Simplex >*                                   root_frequency;
		const TypedDagNode< RbVector< RbVector<double> > >*              lambda;
		const TypedDagNode< RbVector< double > >*                        lambda_times;
		const TypedDagNode< RbVector< RbVector<double> > >*              mu;
		const TypedDagNode< RbVector< double > >*                        mu_times;
		const TypedDagNode< RbVector< RbVector<double> > >*              phi;
		const TypedDagNode< RbVector< double > >*                        phi_times;
		const TypedDagNode< RbVector< RbVector<double> > >*              delta;
		const TypedDagNode< RbVector< double > >*                        delta_times;
		const TypedDagNode< RbVector< RbVector<double> > >*              upsilon;
		const TypedDagNode< RbVector< double > >*                        upsilon_times;
		const TypedDagNode< RbVector< RbVector<double> > >*              gamma;
		const TypedDagNode< RbVector< double > >*                        gamma_times;
		const TypedDagNode< RbVector< RbVector<double> > >*              rho;
		const TypedDagNode< RbVector< double > >*                        rho_times;
		const TypedDagNode< RbVector< RbVector<double> > >*              xi;
		const TypedDagNode< RbVector< double > >*                        xi_times;
		const TypedDagNode< RbVector< RateGenerator > >*                 eta;
		const TypedDagNode< RbVector< double > >*                        eta_times;
		const TypedDagNode< RbVector< CladogeneticProbabilityMatrix > >* omega;
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
