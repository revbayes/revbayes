#ifndef LIKELIHOOD_H
#define LIKELIHOOD_H

#include <ostream>
#include <limits>
#include "LogSum.h"

namespace RevBayesCore
{

class PseudoDataLikelihood
{
    double log_value = -std::numeric_limits<double>::infinity();
public:
    constexpr double  log() const {return log_value;}
    constexpr double& log()       {return log_value;}

    PseudoDataLikelihood& operator +=(const PseudoDataLikelihood& y) {loginc(log_value, y.log()); return *this;}
    PseudoDataLikelihood& operator -=(const PseudoDataLikelihood& y) {logdiff(log_value, y.log()); return *this;}

    constexpr PseudoDataLikelihood& operator *=(const PseudoDataLikelihood& y) {log_value += y.log(); return *this;}
    constexpr PseudoDataLikelihood& operator /=(const PseudoDataLikelihood& y) {log_value -= y.log(); return *this;}

    explicit operator double() const {return exp(log_value);}

    PseudoDataLikelihood() = default;

    explicit PseudoDataLikelihood(double l)
        :log_value(std::log(l))
    { }
};

using std::log;
    
constexpr double log(const PseudoDataLikelihood x)
{
    return x.log();
}

inline PseudoDataLikelihood expToPseudoDataLikelihood(double l)
{
    PseudoDataLikelihood L;
    L.log() = l;
    return L;
}

inline PseudoDataLikelihood operator+(PseudoDataLikelihood x, PseudoDataLikelihood y) {
    PseudoDataLikelihood z = x;
    z += y;
    return z;
}

inline PseudoDataLikelihood operator-(PseudoDataLikelihood x, PseudoDataLikelihood y) {
    PseudoDataLikelihood z = x;
    z -= y;
    return z;
}

inline PseudoDataLikelihood operator*(PseudoDataLikelihood x, PseudoDataLikelihood y) {
    PseudoDataLikelihood z = x;
    z *= y;
    return z;
}

inline PseudoDataLikelihood operator/(PseudoDataLikelihood x, PseudoDataLikelihood y) {
    PseudoDataLikelihood z=x;
    z /= y;
    return z;
}

constexpr bool operator==(PseudoDataLikelihood x, PseudoDataLikelihood y) {
    return log(x)==log(y);
}

constexpr bool operator< (PseudoDataLikelihood x, PseudoDataLikelihood y) {
    return log(x)<log(y);
}

constexpr bool operator> (PseudoDataLikelihood x, PseudoDataLikelihood y) {
    return log(x)>log(y);
}

constexpr bool operator<=(PseudoDataLikelihood x, PseudoDataLikelihood y) {
    return log(x)<=log(y);
}

constexpr bool operator>=(PseudoDataLikelihood x, PseudoDataLikelihood y) {
    return log(x)>=log(y);
}

inline std::ostream& operator<<(std::ostream& o, PseudoDataLikelihood x)
{
    o<<"PDL"<<x.log();
    return o;
}
    
}
    
#endif
