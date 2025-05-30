#ifndef PseudoObservation_H
#define PseudoObservation_H

#include <limits>
#include <cassert>
#include <iostream>
#include "Cloneable.h"


/*
 * The PseudoObservation object is used to handle situations where:
 * - the data is unspecified
 * - the distribution is unspecified
 * - the likelihood function Pr(data|parameter) is specified as a function f(parameter).
 *
 * In RevBayes this would look like:
 *
 *   x ~ dnPseudo(parameter)
 *   x.clamp(f)
 *
 * The likelihood Pr(x|parameter) is then specified to be f(parameter).
 *
 * This only makes sense if x is clamped, so if the user just writes:
 *
 *   x ~ dnPseudo(parameter)
 *
 * then we set x=E where E is the empty pseudodata.  Then Pr(x) = E(x) = 1.
 */

namespace RevBayesCore
{

    class PseudoObservation: public Cloneable
    {
    public:
        virtual ~PseudoObservation() override {};

        PseudoObservation* clone() const override { return new PseudoObservation(*this);}

        // ModelVector< > somehow requires this.
        bool operator==(const PseudoObservation&) const { return false; }
        bool operator!=(const PseudoObservation&) const { return true; }
        bool operator<(const PseudoObservation&) const { return false; }
        bool operator<=(const PseudoObservation&) const { return false; }
    };

    // We need this for TypedDagNode<SiteMixtureModel> for some reason...
    inline std::ostream& operator<<(std::ostream& o, const PseudoObservation&)
    {
        o<<"PseudoObservation";
        return o;
    }
}
#endif
