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
    struct PseudoData: public Cloneable
    {
        virtual double operator()(const T&) = 0;
        virtual ~PseudoData() {};

        PseudoData* clone() const = 0;

        // ModelVector< > somehow requires this.
        bool operator==(const PseudoData&) const { return false; }
        bool operator!=(const PseudoData&) const { return true; }
        bool operator<(const PseudoData&) const { return false; }
        bool operator<=(const PseudoData&) const { return false; }
    };


    // We observed nothing.
    template <typename T>
    struct PseudoDataNothing: public PseudoData<T>
    {
        PseudoDataNothing* clone() const {return new PseudoDataNothing(*this);}
        double operator()(const T&) {return 0;}
    };

    // Whatever we observed, the data has decreasing probability outside [a,b]
    struct PseudoDataInterval: public PseudoData<double>
    {
        double a;
        double b;
        double lambda;
        double operator()(const double& x)
        {
            if (x > a and x < b)
                return 0;
            else
            {
                double d = std::min(std::abs(x-a), std::abs(x-b));
                return -lambda * d;
            }
        }

        PseudoDataInterval* clone() const {return new PseudoDataInterval(*this);}

        PseudoDataInterval(double A, double B, double L):a(A), b(B), lambda(L) {assert(lambda >= 0);}
    };

    // We need this for TypedDagNode<SiteMixtureModel> for some reason...
    template <typename T>
    std::ostream&                                       operator<<(std::ostream& o, const PseudoData<T>& x)
    {
        o<<"PseudoData<"<<typeid(T).name()<<">";
        return o;
    }}

#endif
