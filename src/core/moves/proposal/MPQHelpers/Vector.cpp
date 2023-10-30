#include "Vector.h"
#include <cmath>


using namespace RevBayesCore;

Vector::Vector(void) {

    this->x = 0;
    this->y = 0;
    this->z = 0;
}

Vector::Vector(const Vector& v) {

    this->x = v.x;
    this->y = v.y;
    this->z = v.z;
}

Vector::Vector(int xI, int yI, int zI) {

    this->x = xI;
    this->y = yI;
    this->z = zI;
}

Vector::Vector(double xD, double yD, double zD) {

    this->x = xD;
    this->y = yD;
    this->z = zD;
}

Vector::Vector(mpq_class& xQ, mpq_class& yQ, mpq_class& zQ) {

    this->x = xQ;
    this->y = yQ;
    this->z = zQ;
}

bool Vector::operator==(Vector& rhs) {

    if (this->x != rhs.x)
        return false;
    if (this->y != rhs.y)
        return false;
    if (this->z != rhs.z)
        return false;
    return true;
}

bool Vector::operator!=(Vector& rhs) {

    return !(*this == rhs);
}

Vector& Vector::operator+=(Vector& rhs) {

    x += rhs.x; y += rhs.y; z += rhs.z;
    return *this;
}

Vector Vector::cross(Vector& rhs) {

    // return Vector(y*rhs.z - z*rhs.y, z*rhs.x - x*rhs.z, x*rhs.y - y*rhs.x)
    mpq_class yz = this->y * rhs.z;
    mpq_class zy = this->z * rhs.y;
    mpq_class zx = this->z * rhs.x;
    mpq_class xz = this->x * rhs.z;
    mpq_class xy = this->x * rhs.y;
    mpq_class yx = this->y * rhs.x;
    mpq_class diff1 = yz-zy;
    mpq_class diff2 = zx-xz;
    mpq_class diff3 = xy-yx;
    Vector v(diff1, diff2, diff3);
    return v;
}

inline mpq_class Vector::distanceSquared(const Vector& vec) const {

    return (vec.x-x)*(vec.x-x) + (vec.y-y)*(vec.y-y) + (vec.z-z)*(vec.z-z);
}

std::string Vector::getStr(void) {

    std::string str = "{";
    str += std::to_string(x.get_d()) + ", " + std::to_string(y.get_d()) + ", " + std::to_string(z.get_d());
    str += "}";
    return str;
}

double Vector::length(void) {

    mpq_class sqQ = (x * x) + (y * y) + (z * z);
    double sqD = sqQ.get_d();
    return std::sqrtf(sqD);
}

void Vector::normalize(void) {

    //@@const double EPSILON = 0.000001f;
    mpq_class dist = x * x + y * y + z * z;
    if (dist == 0)
        return; // or return *this

    mpf_class distInv = 1 / dist;
    mpf_class rDistInv = distInv;
    mpf_class invLength = sqrt(rDistInv);
    mpq_class invLengthQ(invLength);
    
    x *= invLengthQ;
    y *= invLengthQ;
    z *= invLengthQ;
    //return *this;
}

void Vector::set(mpq_class& xQ, mpq_class& yQ, mpq_class& zQ) {

    x = xQ;
    y = yQ;
    z = zQ;
}
