#ifndef RlContinuousStochasticNode_H
#define RlContinuousStochasticNode_H

#include "ContinuousStochasticNode.h"
#include "RevMemberObject.h"
#include "RlContinuousDistribution.h"
#include "RlPositiveContinuousDistribution.h"
#include "Probability.h"

namespace RevLanguage {
    
    class ContinuousStochasticNode : public RevBayesCore::ContinuousStochasticNode, public RevMemberObject {
        
    public:
        ContinuousStochasticNode(const std::string& n, RevBayesCore::ContinuousDistribution* dist, ContinuousDistribution* rlDist);
        ContinuousStochasticNode(const std::string& n, RevBayesCore::ContinuousDistribution* dist, PositiveContinuousDistribution* rlDist);
        ContinuousStochasticNode(const std::string& n, RevBayesCore::ContinuousDistribution* dist, TypedDistribution<Probability>* rlDist);
        ContinuousStochasticNode(const ContinuousStochasticNode &n);
        virtual                            ~ContinuousStochasticNode(void);

        ContinuousStochasticNode&                   operator=(const ContinuousStochasticNode &n);

        
        // public methods
        ContinuousStochasticNode*                   clone(void) const;                                                                              //!< Clone the node
        virtual RevPtr<RevVariable>                 executeMethod(const std::string& name, const std::vector<Argument>& args, bool &found);         //!< Execute member method (if applicable)
        const MethodTable&                          getMethods( void ) const;                                                                       //!< Get the member methods
        Distribution&                               getRlDistribution(void);                                                                        //!< Get the Rev distribution
        const Distribution&                         getRlDistribution(void) const;                                                                  //!< Get the Rev distribution (const)
        void                                        printStructureInfo(std::ostream &o, bool verbose=false) const;                                  //!< Print information on structure
        
    private:
        
        Distribution*                               rlDistribution;                                                                                 //!< Rev distribution
        MethodTable                                 methods;
        
        
    };
    
}

#endif

