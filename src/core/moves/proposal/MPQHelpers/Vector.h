#ifndef Vector_hpp
#define Vector_hpp

#include <gmpxx.h>
#include <iostream>



namespace RevBayesCore {

class Vector {
    
    /**
     * A light-weight class to represent a vector in 3D, with coordinates x, y, z. The
     * coordinates are GMP rational numbers.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (John Huelsenbeck)
     * @since 2014-11-18, version 1.0
     */
    
public:
    Vector(void);
    Vector(const Vector& v);
    Vector(int xI, int yI, int zI);
    Vector(double xD, double yD, double zD);
    Vector(mpq_class& xQ, mpq_class& yQ, mpq_class& zQ);
    bool            operator==(Vector& rhs);
    bool            operator!=(Vector& rhs);
    Vector&         operator+=(Vector& rhs);
    Vector          cross(Vector& rhs);
    mpq_class       distanceSquared(const Vector& vec) const;
    mpq_class&      getX(void) { return x; }
    mpq_class&      getY(void) { return y; }
    mpq_class&      getZ(void) { return z; }
    std::string     getStr(void);
    double          length(void);
    void            normalize(void);
    void            set(mpq_class& xQ, mpq_class& yQ, mpq_class& zQ);
    void            setX(mpq_class& xQ) { x = xQ; }
    void            setY(mpq_class& yQ) { y = yQ; }
    void            setZ(mpq_class& zQ) { z = zQ; }
    mpq_class       x;
    mpq_class       y;
    mpq_class       z;
    
    friend std::ostream& operator<<(std::ostream& os, const Vector& pt);
};

}

inline std::ostream& operator<<(std::ostream& os, const RevBayesCore::Vector& pt) {

    os << "{" << pt.x << ", " << pt.y << ", " << pt.z << "}";
    return os;
}

#endif
