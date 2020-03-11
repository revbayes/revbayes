#ifndef GeneralizedLineageHeterogeneousBirthDeathSamplingProcess_H
#define GeneralizedLineageHeterogeneousBirthDeathSamplingProcess_H

#include "AbstractHomologousDiscreteCharacterData.h"
#include "TreeDiscreteCharacterData.h"
#include "CladogeneticSpeciationRateMatrix.h"
#include "RateGenerator.h"
#include "Simplex.h"
#include "Taxon.h"
#include "Tree.h"
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
     * This is the generalized lineage-heterogenous (state-dependent?) birth-death-sampling process.
     * An interface for tensorphylo
     *
     * Michael R. May and Xavier Meyer 2020/3/08
     *
     */
    class GeneralizedLineageHeterogeneousBirthDeathSamplingProcess : public TypedDistribution<Tree>, public TreeChangeEventListener, public MemberObject< RbVector<double> > {
        
    public:
        GeneralizedLineageHeterogeneousBirthDeathSamplingProcess();
        
        // pure virtual member functions
        virtual GeneralizedLineageHeterogeneousBirthDeathSamplingProcess* clone(void) const;
        virtual                                                           ~GeneralizedLineageHeterogeneousBirthDeathSamplingProcess(void);                                                              //!< Virtual destructor
        
        double                                                            computeLnProbability(void);
        void                                                              fireTreeChangeEvent(const TopologyNode &n, const unsigned& m=0);                                                 //!< The tree has changed and we want to know which part.
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
        
    };
    
}

#endif
