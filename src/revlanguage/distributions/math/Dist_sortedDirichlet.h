#ifndef Dist_SortedDirichlet_H
#define Dist_SortedDirichlet_H

#include "SortedDirichletDistribution.h"
#include "RlTypedDistribution.h"
#include "RlSimplex.h"

namespace RevLanguage {

    /**
     * The RevLanguage wrapper of the sorted Dirichlet distribution.
     *
     * The RevLanguage wrapper of the Dirichlet distribution simply
     * manages the interactions through the Rev with our core.
     * That is, the internal distribution object can be constructed and hooked up
     * in a model graph.
     * See the DirichletDistributionSorted for more details.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2012-08-08, version 1.0
     *
     */
    class Dist_sortedDirichlet : public TypedDistribution<Simplex> {

    public:
        Dist_sortedDirichlet( void );                                       //!< constructor
        virtual ~Dist_sortedDirichlet();                                    //!< destructor

        // Basic utility functions
        Dist_sortedDirichlet*     clone(void) const;                        //!< Clone the object
        static const std::string& getClassType(void);                       //!< Get Rev type
        static const TypeSpec&    getClassTypeSpec(void);                   //!< Get class type spec
        const TypeSpec&           getTypeSpec(void) const;                  //!< Get the type spec of the instance
        std::string               getDistributionFunctionName(void) const;  //!< Get the Rev-name for this distribution.
        const MemberRules&        getParameterRules(void) const;            //!< Get member rules (const)
        void                      printValue(std::ostream& o) const;        //!< Print the general information on the function ('usage')

        // Distribution functions you have to override
        RevBayesCore::SortedDirichletDistribution* createDistribution(void) const;

    protected:
        void setConstParameter( const std::string& name
                              , const RevPtr<const RevVariable> &var
                              );                                            //!< Set member variable

    private:
        RevPtr<const RevVariable> alpha;

    };
}

#endif
