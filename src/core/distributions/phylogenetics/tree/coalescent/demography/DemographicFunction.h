#ifndef DemographicFunction_H
#define DemographicFunction_H

#include <iostream>
#include <vector>

#include "Cloneable.h"

namespace RevBayesCore {
    
    class DagNode;
    
    /**
     * @brief DemographicFunction: interface for all core demographic functions.
     *
     *  This class contains a vector of variables represented as DAG nodes.
     */
    class DemographicFunction : public Cloneable {
        
    public:
        // destructor
        virtual                                                 ~DemographicFunction(void);                                 //!< Destructor
        
        // public methods
        virtual void                                            addVariable(const DagNode *n);                              //!< Add a variable to the nodes vector
        virtual void                                            swapNode(const DagNode *oldN, const DagNode *newN);         //!< Replacing a DAG node containing a variabe with a different one
        
        // getters and setters
        const std::vector<const DagNode *>&                     getDagNodes(void) const;                                    //!< Get the nodes vector

        // pure virtual public methods
        virtual DemographicFunction*                            clone(void) const = 0;                                      //!< Clone the DemographicFunction
        virtual double                                          getDemographic(double t) const = 0;                         //!< Returns the demographic function N(t) at time t.
        virtual double                                          getIntegral(double start, double finish) const = 0;         //!< Calculates the integral 1/N(x) dx between start and finish.
        virtual double                                          getWaitingTime(double time, double lambda, double ploidy=1.0) const = 0;       //!< Calculates the waiting time until the next coalescent event.

    protected:
        DemographicFunction(void);                                                                                          //!< Default constructor
        DemographicFunction(const DemographicFunction &f);                                                                  //!< Copy constructor
        DemographicFunction&                                    operator=(const DemographicFunction &f);                    //!< Assignment operator
        
        virtual void                                            swapNodeInternal(const DagNode *oldN, const DagNode *newN) {} //!< Internally replacing a DAG node containing a variable with a different one

    private:
        
        std::vector<const DagNode *>                            variables;

    };
    
    // Global functions using the class
    std::ostream&                                               operator<<(std::ostream& o, const DemographicFunction& x);      //!< Overloaded output stream
    
}

std::ostream&                                                   operator<<(std::ostream& o, const RevBayesCore::DemographicFunction& x);  //!< Overloaded output stream

#endif /* DemographicFunction_H */
