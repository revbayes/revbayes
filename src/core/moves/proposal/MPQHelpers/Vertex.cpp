#include "Vertex.h"

using namespace RevBayesCore;

Vertex::Vertex(void) : Vector() {

    from = nullptr;
    to = nullptr;
}

Vertex::Vertex(const Vector& v) : Vector(v) {

    from = nullptr;
    to = nullptr;
}

Vertex::Vertex(Vertex& v) : Vector(v) {

    from = nullptr;
    to = nullptr;
}

Vertex::Vertex(mpq_class& xq, mpq_class& yq, mpq_class& zq) : Vector(xq, yq, zq) {

    from = nullptr;
    to = nullptr;
}

bool Vertex::operator==(const Vertex& rhs) const {

    if (this->x != rhs.x)
        return false;
    if (this->y != rhs.y)
        return false;
    if (this->z != rhs.z)
        return false;
    return true;
}

void Vertex::clean(void) {

    from = nullptr;
    to = nullptr;
    this->x = 0;
    this->y = 0;
    this->z = 0;
}
