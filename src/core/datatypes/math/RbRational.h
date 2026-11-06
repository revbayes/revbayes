#ifndef RbRational_H
#define RbRational_H

/* Exact rational arithmetic for the reversible-jump substitution-model code.
 *
 * Two options, chosen at compile time:
 *
 *   default         mpq_class is the class below, a wrapper around
 *                   boost::multiprecision::cpp_rational. Boost is already a
 *                   RevBayes dependency, so nothing new has to be installed.
 *
 *   -DRB_USE_GMP    mpq_class is GMP's own (gmpxx.h). The executable must then
 *                   be linked against gmpxx and gmp.
 *
 * The wrapper keeps the name mpq_class and provides exactly the part of the
 * gmpxx interface the reversible jump MCMC uses, so the code that uses it compiles
 * under either choice.
 * 
 * Note that the boost default is slower than the GMP code, which is quite fast.
 *
 * @date Last modified: $Date$
 * @author John Huelsenbeck
*/

#if defined(RB_USE_GMP)

#include <gmpxx.h>

#else

#include <boost/multiprecision/cpp_int.hpp>
#include <cfloat>
#include <cmath>
#include <ostream>
#include <stdexcept>
#include <type_traits>
#include <utility>



    class mpq_class {

        public:
            typedef boost::multiprecision::cpp_int      integer_type;
            typedef boost::multiprecision::cpp_rational rational_type;

                                    mpq_class(void) : value(0) {}
                                    mpq_class(double x) { assignDouble(x); }
                                    mpq_class(float x) { assignDouble(x); }
                                    template <typename I, typename std::enable_if<std::is_integral<I>::value && !std::is_same<I,bool>::value, int>::type = 0>
                                    mpq_class(I i) : value(i) {}
                                    template <typename I, typename J, typename std::enable_if<std::is_integral<I>::value && std::is_integral<J>::value, int>::type = 0>
                                    mpq_class(I n, J d) : value(n) { value /= rational_type(d); }
            explicit                mpq_class(const rational_type& r) : value(r) {}
            explicit                mpq_class(rational_type&& r) : value(std::move(r)) {}

            double                  get_d(void) const;
            const rational_type&    get_value(void) const { return value; }

            mpq_class&              operator+=(const mpq_class& rhs) { value += rhs.value; return *this; }
            mpq_class&              operator-=(const mpq_class& rhs) { value -= rhs.value; return *this; }
            mpq_class&              operator*=(const mpq_class& rhs) { value *= rhs.value; return *this; }
            mpq_class&              operator/=(const mpq_class& rhs) { value /= rhs.value; return *this; }
            mpq_class               operator-(void) const { return mpq_class( rational_type(-value) ); }
            mpq_class               operator+(void) const { return *this; }

            friend mpq_class        operator+(const mpq_class& a, const mpq_class& b) { return mpq_class( rational_type(a.value + b.value) ); }
            friend mpq_class        operator-(const mpq_class& a, const mpq_class& b) { return mpq_class( rational_type(a.value - b.value) ); }
            friend mpq_class        operator*(const mpq_class& a, const mpq_class& b) { return mpq_class( rational_type(a.value * b.value) ); }
            friend mpq_class        operator/(const mpq_class& a, const mpq_class& b) { return mpq_class( rational_type(a.value / b.value) ); }

            friend bool             operator==(const mpq_class& a, const mpq_class& b) { return a.value == b.value; }
            friend bool             operator!=(const mpq_class& a, const mpq_class& b) { return a.value != b.value; }
            friend bool             operator< (const mpq_class& a, const mpq_class& b) { return a.value <  b.value; }
            friend bool             operator<=(const mpq_class& a, const mpq_class& b) { return a.value <= b.value; }
            friend bool             operator> (const mpq_class& a, const mpq_class& b) { return a.value >  b.value; }
            friend bool             operator>=(const mpq_class& a, const mpq_class& b) { return a.value >= b.value; }

            friend std::ostream&    operator<<(std::ostream& os, const mpq_class& q);

        private:
            void                    assignDouble(double x);

            rational_type           value;
    };

    inline void mpq_class::assignDouble(double x) {

        if ( std::isfinite(x) == false )
            throw std::domain_error("mpq_class: a non-finite double has no rational value");
        if ( x == 0.0 )
            {
            value = 0;
            return;
            }

        int e = 0;
        double m = std::frexp(x, &e);                   // x = m * 2^e, with 1/2 <= |m| < 1
        m = std::ldexp(m, DBL_MANT_DIG);                // now an integer of at most DBL_MANT_DIG bits
        e -= DBL_MANT_DIG;
        integer_type n = static_cast<long long>(m);
        if ( e >= 0 )
            {
            value = rational_type( integer_type(n << e) );
            }
        else
            {
            value  = rational_type( n );
            value /= rational_type( integer_type(integer_type(1) << -e) );
            }
    }

    inline double mpq_class::get_d(void) const {

        integer_type num = boost::multiprecision::numerator(value);
        if ( num == 0 )
            return 0.0;
        const integer_type den = boost::multiprecision::denominator(value);   // always positive

        bool negative = false;
        if ( num < 0 )
            {
            negative = true;
            num = -num;
            }

        // with this shift the integer quotient has DBL_MANT_DIG + 1 or + 2 bits
        long shift = static_cast<long>(DBL_MANT_DIG) + 1 - ( static_cast<long>(boost::multiprecision::msb(num)) - static_cast<long>(boost::multiprecision::msb(den)) );
        integer_type q;
        if ( shift >= 0 )
            q = integer_type(num << static_cast<unsigned long>(shift)) / den;
        else
            q = num / integer_type(den << static_cast<unsigned long>(-shift));

        // drop the bits a double cannot hold; integer shifts truncate, as mpq_get_d does
        long excess = static_cast<long>(boost::multiprecision::msb(q)) + 1 - static_cast<long>(DBL_MANT_DIG);
        if ( excess > 0 )
            {
            q >>= static_cast<unsigned long>(excess);
            shift -= excess;
            }

        double d = std::ldexp( q.convert_to<double>(), static_cast<int>(-shift) );
        return ( negative == true ? -d : d );
    }

    /* The same format as GMP: "n/d", or just "n" when the denominator is one. */
    inline std::ostream& operator<<(std::ostream& os, const mpq_class& q) {

        os << boost::multiprecision::numerator(q.value);
        if ( boost::multiprecision::denominator(q.value) != 1 )
            os << "/" << boost::multiprecision::denominator(q.value);
        return os;
    }

#endif

#include <stdexcept>
#include <string>

    /* Read a rational written by operator<< ("n/d", or just "n"), exactly. This is
       how checkpoints keep the rate matrix exact across a restart. */
    inline mpq_class mpqFromString(const std::string& s) {

#   if defined(RB_USE_GMP)
        mpq_class q(s, 10);
        if ( q.get_den() == 0 )
            throw std::invalid_argument("mpqFromString: zero denominator in \"" + s + "\"");
        q.canonicalize();
        return q;
#   else
        size_t slash = s.find('/');
        mpq_class::integer_type n( s.substr(0, slash) );
        mpq_class::integer_type d( slash == std::string::npos ? std::string("1") : s.substr(slash + 1) );
        if ( d == 0 )
            throw std::invalid_argument("mpqFromString: zero denominator in \"" + s + "\"");
        return mpq_class( mpq_class::rational_type(n, d) );
#   endif
    }

#endif
