#ifndef Distribution_H
#define Distribution_H

#include "Cloneable.h"
#include "Parallelizable.h"
#include "RbOrderedSet.h"
#include "RevPtr.h"
#include "RevVariable.h"

#include <iostream>
#include <set>
#include <vector>

namespace RevBayesCore {
    
    class DagNode;
    
    /**
     * @brief Distribution: interface for all core distributions
     *
     * Distributions are typically associated with stochastic nodes but can also be used as input parameters
     * to other distributions, like the DPP, or to functions.
     *
     * All stochastic nodes hold a distribution pointer. The value of the stochastic node is returned via
     * a call to getValue in the distribution.
     *
     * Some distributions require a distribution as a parameter, e.g. a generating distribution. Thus,
     * we need to implement distributions as objects storable in DAG nodes.
     *
     * Every distribution owns its value and hence this class is templated. Owning the value
     * has the advantage that calls to update can modify the value instead of creating a new object.
     * This is beneficial in distributions on large objects, making it possible to work with partial
     * updates and restores, and partial probability calculations.
     *
     * Each derived distribution is responsible for managing its parameters, swapping them and returning
     * a set of pointers to them.
     */
    class Distribution : public Cloneable, public Parallelizable {
        
    public:
        // destructor
        virtual                                                 ~Distribution(void);
        
        // public methods
        virtual void                                            bootstrap(void);                                                                    //!< Draw a new random value from the distribution
        virtual RevLanguage::RevPtr<RevLanguage::RevVariable>   executeProcedure(const std::string &n, const std::vector<DagNode*> args, bool &f);  //!< execute the procedure
        virtual void                                            getAffected(RbOrderedSet<DagNode *>& affected, const DagNode* affecter);            //!< get affected nodes
        virtual std::vector<double>                             getMixtureProbabilities(void) const;
        virtual size_t                                          getNumberOfMixtureElements(void) const;                                             //!< Get the number of elements for this value
        const std::vector<const DagNode*>&                      getParameters(void) const;                                                          //!< get the parameters of the function
        void                                                    keep(const DagNode* affecter);
        virtual void                                            reInitialized( void );                                                              //!< The model was re-initialized
        void                                                    restore(const DagNode *restorer);
        virtual void                                            setMcmcMode(bool tf);                                                               //!< Change the likelihood computation to or from MCMC mode.
        void                                                    swapParameter(const DagNode *oldP, const DagNode *newP);                            //!< Exchange the parameter
        virtual void                                            touch(const DagNode *toucher, bool touchAll);
        
        // pure virtual public methods
        virtual Distribution*                                   clone(void) const = 0;                                                              //!< Clone the distribution
        virtual double                                          computeLnProbability(void) = 0;                                                     //!< Compute the ln probability
        virtual void                                            redrawValue(SimulationCondition c = SimulationCondition::MCMC) = 0;                 //!< Draw a new random value from the distribution
        
    protected:
        Distribution(void);                                                                                                                         //!< Default constructor
        Distribution(const Distribution &f);                                                                                                        //!< Copy constructor
        Distribution&                                           operator=(const Distribution &f);                                                   //!< Assignment operator
        
        // keep specialization for derived classes
        virtual void                                            keepSpecialization(const DagNode* affecter);
        virtual void                                            restoreSpecialization(const DagNode *restorer);
        virtual void                                            touchSpecialization(const DagNode *toucher, bool touchAll);
        
        // swap parameter methods for internal use of derived classes
        void                                                    addParameter(const DagNode* p);                                                     //!< add a parameter to the function
        void                                                    removeParameter(const DagNode* p);                                                  //!< remove a parameter from the function
        virtual void                                            swapParameterInternal(const DagNode *oldP, const DagNode *newP) = 0;                //!< Exchange the parameter

        
    private:
        
        std::vector<const DagNode*>                             parameters;
        
    };
    
    // Global functions using the class
    std::ostream&                                               operator<<(std::ostream& o, const Distribution& x);                                 //!< Overloaded output operator
    
}

#endif
