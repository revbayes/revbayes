#ifndef PseudoData_H
#define PseudoData_H

#include <limits>
#include <cassert>
#include "Cloneable.h"


/*
 * The PseudoData object is used to handle situations where:
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

    template <typename T>
    class PseudoData: public Cloneable
    {
    public:
        typedef std::function<double (const T&)> func_t;

    private:
        func_t log_likelihood;

    public:
        PseudoData(const PseudoData<T>&) = default;
        PseudoData(func_t f):log_likelihood(f) {}
        PseudoData(): PseudoData( [](const T&) { return 0.0; }) {}

        virtual ~PseudoData() override {};

        PseudoData<T>* clone() const override { return new PseudoData<T>(*this);}

        // Call this on the parameter to get the log-likelihood of the unspecified data.
        double operator()(const T& x) const { return log_likelihood(x); }

        // ModelVector< > somehow requires this.
        bool operator==(const PseudoData&) const { return false; }
        bool operator!=(const PseudoData&) const { return true; }
        bool operator<(const PseudoData&) const { return false; }
        bool operator<=(const PseudoData&) const { return false; }
    };

    // We need this for TypedDagNode<SiteMixtureModel> for some reason...
    template <typename T>
    std::ostream&                                       operator<<(std::ostream& o, const PseudoData<T>& x)
    {
        o<<"PseudoData<"<<typeid(T).name()<<">";
        return o;
    }}

#endif
