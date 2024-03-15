#ifndef PLANE_H
#define PLANE_H

#include <gmpxx.h>
#include "Vector.h"


namespace RevBayesCore {

    class Plane {
        
        /**
         * A light-weight class to represent a plane in 3D. Planes are represented in the general
         * form, with parameters ax + by + cz + d = 0. These plane parameters are GMP rational
         * numbers.
         *
         * @copyright Copyright 2009-
         * @author The RevBayes Development Core Team (John Huelsenbeck)
         * @since 2014-11-18, version 1.0
         */
        
    public:
                            Plane(void);
                            Plane(Vector pt1, Vector pt2, Vector pt3);
                            Plane(Vector& pt1, Vector& pt2, Vector& pt3);
        bool                operator==(const Plane& rhs) const;
        bool                operator<(const Plane& rhs) const;
        const mpq_class&    getA(void) const { return a; }
        const mpq_class&    getB(void) const { return b; }
        const mpq_class&    getC(void) const { return c; }
        const mpq_class&    getD(void) const { return d; }
        std::string         getStr(void);
        void                set(Vector pt1, Vector pt2, Vector pt3);
        
    private:
        mpq_class           a;
        mpq_class           b;
        mpq_class           c;
        mpq_class           d;
        
        friend std::ostream& operator<<(std::ostream& os, const Plane& pl);
    };
}

inline std::ostream& operator<<(std::ostream& os, const RevBayesCore::Plane& pl) {

    os << "{" << pl.getA() << ", " << pl.getB() << ", " << pl.getC() << ", " << pl.getD() << "}";
    return os;
}

#endif
