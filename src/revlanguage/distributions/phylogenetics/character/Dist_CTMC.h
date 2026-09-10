#ifndef Dist_CTMC_H
#define Dist_CTMC_H

#include <cmath>
#include <iosfwd>
#include <string>
#include <vector>

#include "AbstractDiscreteTaxonData.h"
#include "RlAbstractDiscreteTaxonData.h"
#include "RlTypedDistribution.h"
#include "DagMemberFunction.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "MethodTable.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlDagMemberFunction.h"
#include "RlDeterministicNode.h"
#include "RlStochasticNode.h"
#include "RlTypedFunction.h"
#include "StochasticNode.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"
#include "TypedFunction.h"

namespace RevLanguage {
class TypeSpec;

    class Dist_CTMC : public TypedDistribution< AbstractDiscreteTaxonData > {

    public:
        Dist_CTMC( void );
        virtual ~Dist_CTMC();

        // Basic utility functions
        Dist_CTMC*                                      clone(void) const;                                                                      //!< Clone the object
        static const std::string&                       getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::string                                     getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        virtual MethodTable                             getDistributionMethods( void ) const;                                                                       //!< Get the member methods
        const TypeSpec&                                 getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                              getParameterRules(void) const;                                                          //!< Get member rules (const)
        void                                            printValue(std::ostream& o) const;                                                      //!< Print the general information on the function ('usage')


        // Distribution functions you have to override
        RevBayesCore::TypedDistribution< RevBayesCore::AbstractDiscreteTaxonData >*      createDistribution(void) const;

    protected:

        void                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable


    private:

        RevPtr<const RevVariable>                       q;
        RevPtr<const RevVariable>                       site_rates;
        RevPtr<const RevVariable>                       site_rates_probs;
        RevPtr<const RevVariable>                       site_matrices;
        RevPtr<const RevVariable>                       site_matrices_probs;
        RevPtr<const RevVariable>                       root_frequencies;
        RevPtr<const RevVariable>                       nSites;
        RevPtr<const RevVariable>                       type;
        RevPtr<const RevVariable>                       coding;

    };

}

#endif
