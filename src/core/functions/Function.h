#ifndef Function_H
#define Function_H

#include "Cloneable.h"
#include "RbOrderedSet.h"

#include <iostream>
#include <set>
#include <vector>

namespace RevBayesCore {
    
    class DagNode;
    
    class Function : public Cloneable {
        
    public:
        // constructors and destructor
        virtual                                    ~Function(void);
               
        // public methods
        bool                                        forceUpdates(void) const;                                               //!< Does this method forces the DAG node to always call update even if not touched?
        virtual void                                getAffected(RbOrderedSet<DagNode *>& affected, DagNode* affecter);          //!< get affected nodes
        const std::vector<const DagNode*>&          getParameters(void) const;                                              //!< get the parameters of the function
        virtual void                                keep(DagNode* affecter);
        virtual void                                reInitialized( void );                                                  //!< The model was re-initialized
        virtual void                                restore(DagNode *restorer);
        void                                        setForceUpdates(bool tf);                                               //!< Does this method forces the DAG node to always call update even if not touched?
        void                                        swapParameter(const DagNode *oldP, const DagNode *newP);                //!< Exchange the parameter
        virtual void                                touch(DagNode *toucher );

        // pure virtual public methors
        virtual Function*                           clone(void) const = 0;                                                  //!< Clone the function
        virtual void                                update(void) = 0;                                                       //!< Clone the function
        
    protected:
        Function(void);                                                                                             //!< Default constructor
        Function(const Function &f);                                                                                //!< Copy constructor
        Function&                                   operator=(const Function &f);                                           //!< Assignment operator
        
        void                                        addParameter(const DagNode* p);                                         //!< add a parameter to the function
        void                                        removeParameter(const DagNode* p);                                      //!< remove a parameter from the function
        virtual void                                swapParameterInternal(const DagNode *oldP, const DagNode *newP) = 0;    //!< Exchange the parameter
                
    private:
        
        std::vector<const DagNode*>                 parameters;
        bool                                        force_update;
    };
    
    // Global functions using the class
    std::ostream&                                   operator<<(std::ostream& o, const Function& x);                                //!< Overloaded output operator
    
}

#endif
