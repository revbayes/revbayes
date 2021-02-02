#ifndef RlContinuousStochasticNode_H
#define RlContinuousStochasticNode_H

#include <iosfwd>
#include <vector>

#include "ContinuousStochasticNode.h"
#include "RevMemberObject.h"
#include "MethodTable.h"
#include "RevPtr.h"

namespace RevBayesCore { class ContinuousDistribution; }

namespace RevLanguage {
class Argument;
class ContinuousDistribution;
class Distribution;
class PositiveContinuousDistribution;
class Probability;
class RevVariable;
template <typename rlType> class TypedDistribution;
    
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

