#ifndef DemographicFunction_H
#define DemographicFunction_H

#include "Cloneable.h"

#include <iostream>
#include <set>
#include <vector>

namespace RevBayesCore {
    
    class DagNode;
    
    /**
     * @brief DemographicFunction: interface for all core DemographicFunctions
     *
     *
     */
    class DemographicFunction : public Cloneable {
        
    public:
        // destructor
        virtual                                                 ~DemographicFunction(void);
        
        // public methods
        virtual void                                            addVariable(const DagNode *n);
        virtual void                                            swapNode(const DagNode *oldN, const DagNode *newN);
//        virtual void                                            removeVariable(DagNode *n);
        
        // getters and setters
        const std::vector<const DagNode *>&                     getDagNodes(void) const;                                                            //!< Get the nodes vector
//        void                                                    setDagNodes(const std::vector<DagNode *>& args);

        // pure virtual public methods
        virtual DemographicFunction*                            clone(void) const = 0;                                      //!< Clone the DemographicFunction
        virtual double                                          getDemographic(double t) const = 0;                         //!< demographic function N(t) at time t
//        virtual double                                          getLogDemographic(double t) const = 0;                      //!< log-demographic function N(t) at time t
//        virtual double                                          getIntensity(double t) const = 0;                           //!< demographic intensity function at time t (= integral 1/N(x) dx from 0 to t).
//        virtual double                                          getInverseIntensity(double x) const = 0;                    //!< inverse demographic intensity function
        virtual double                                          getIntegral(double start, double finish) const = 0;         //!< Calculates the integral 1/N(x) dx between start and finish.
//        virtual int                                             getNumberOfArguments(void) const = 0;                       //!< number of arguments for this function
//        virtual std::string                                     getArgumentName(int n) const = 0;                           //!<
//        virtual double                                          getArgument(int n) const = 0;                               //!< the n'th argument of this function.
//        virtual void                                            setArgument(int n, double value) const = 0;                 //!< Sets the value of the nth argument of this function
//        virtual double                                          getLowerBound(int n) const = 0;                             //!< the lower bound of the nth argument of this function.
//        virtual double                                          getUpperBound(int n) const = 0;                             //!< the upper bound of the nth argument of this function.
//        virtual double                                          getThreshold(void) const = 0;                               //!< A threshold for underflow on calculation of likelihood of internode intervals.

    protected:
        DemographicFunction(void);                                                                                          //!< Default constructor
        DemographicFunction(const DemographicFunction &f);                                                                  //!< Copy constructor
        DemographicFunction&                                    operator=(const DemographicFunction &f);                    //!< Assignment operator
        
        virtual void                                            swapNodeInternal(const DagNode *oldN, const DagNode *newN) {}

    private:
        
        std::vector<const DagNode *>                            variables;

    };
    
    // Global functions using the class
    std::ostream&                                               operator<<(std::ostream& o, const DemographicFunction& x);      //!< Overloaded output stream
    
}

std::ostream&                                                   operator<<(std::ostream& o, const RevBayesCore::DemographicFunction& x);  //!< Overloaded output stream

#endif /* DemographicFunction_H */
