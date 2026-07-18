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
     * The observation half of the factored FBD-range model: conditions on a birth-death range
     * process (a dnFBDRP or dnFBDSP node) and contributes the fossil-record log-density. The
     * occurrences come from that range process; value is the taxon vector.
     *
     *     rec ~ dnFossilRecord( ranges=bd, complete=false )
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

        RevPtr<const RevVariable>                       ranges;
        RevPtr<const RevVariable>                       complete;
    };

}

#endif
