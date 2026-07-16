#ifndef Dist_FossilRecord_H
#define Dist_FossilRecord_H

#include "FossilRecordProcess.h"
#include "ModelVector.h"
#include "RlTaxon.h"
#include "RlTypedDistribution.h"

namespace RevLanguage {

    /**
     * RevLanguage wrapper of the fossil-record (reporting) model, dnFossilRecord.
     *
     * The observation half of the factored FBD-range model: conditions on a birth-death-range
     * skeleton (a MatrixReal-valued dnFBDRP node) and contributes the fossil-record log-density.
     * Value is the taxon vector (clamped to the observed occurrences).
     *
     *     rec ~ dnFossilRecord( skeleton=bd, reporting="uniform", taxa=taxa )
     *     rec.clamp( taxa )
     */
    class Dist_FossilRecord : public TypedDistribution< ModelVector<Taxon> > {

    public:
        Dist_FossilRecord( void );

        Dist_FossilRecord*                              clone(void) const;
        static const std::string&                       getClassType(void);
        static const TypeSpec&                          getClassTypeSpec(void);
        std::string                                     getDistributionFunctionName(void) const;
        const TypeSpec&                                 getTypeSpec(void) const;
        const MemberRules&                              getParameterRules(void) const;

        RevBayesCore::FossilRecordProcess*              createDistribution(void) const;

    protected:
        void                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);

        RevPtr<const RevVariable>                       skeleton;
        RevPtr<const RevVariable>                       reporting;
        RevPtr<const RevVariable>                       taxa;
    };

}

#endif
