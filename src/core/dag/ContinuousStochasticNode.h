#ifndef ContinuousStochasticNode_H
#define ContinuousStochasticNode_H

#include <iosfwd>
#include <string>
#include <vector>

#include "StochasticNode.h"
#include "DynamicNode.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore {
class ContinuousDistribution;
    
    class ContinuousStochasticNode : public StochasticNode<double> {
        
    public:
        ContinuousStochasticNode(const std::string &name, ContinuousDistribution *d);
        virtual                                            ~ContinuousStochasticNode(void) {}                                 //!< Virtual destructor
        
        // public member functions
        virtual ContinuousStochasticNode*                   clone(void) const;                                          //!< Create an independent clone
        virtual double                                      getMax(void) const;                                         //!< Get the maximum value the variable can be
        virtual double                                      getMin(void) const;                                         //!< Get the minimum value the variable can be
                
    };
    
}

#endif
