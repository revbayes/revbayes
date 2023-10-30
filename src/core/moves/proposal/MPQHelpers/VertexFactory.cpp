#include "Vertex.h"
#include "VertexFactory.h"

using namespace RevBayesCore;

VertexFactory::VertexFactory(void) {

}

VertexFactory::~VertexFactory(void) {

    for (std::set<Vertex*>::iterator v=allocatedVertices.begin(); v != allocatedVertices.end(); v++)
        delete (*v);
}

void VertexFactory::drainPool(void) {

    for (std::vector<Vertex*>::iterator v=vertexPool.begin(); v != vertexPool.end(); v++)
        {
        allocatedVertices.erase( *v );
        delete (*v);
        }
}

Vertex* VertexFactory::getVertex(Vector& v) {

    Vertex* newV = getVertex();
    newV->setX(v.x);
    newV->setY(v.y);
    newV->setZ(v.z);
    newV->setFrom(nullptr);
    newV->setTo(nullptr);
    return newV;
}

Vertex* VertexFactory::getVertex(Vertex& v) {

    Vertex* newV = getVertex();
    newV->setX(v.x);
    newV->setY(v.y);
    newV->setZ(v.z);
    newV->setFrom(v.getFrom());
    newV->setTo(v.getTo());
    return newV;
}

Vertex* VertexFactory::getVertex(void) {

    if ( vertexPool.empty() == true )
        {
        /* If the pool is empty, we allocate a new vertex and return it. We
           do not need to add it to the pool. */
        Vertex* v = new Vertex;
        allocatedVertices.insert( v );
        onLoan.insert(v);
        return v;
        }
    
    // Return a node from the node pool, remembering to remove it from the pool.
    Vertex* v = vertexPool.back();
    vertexPool.pop_back();
    onLoan.insert(v);
    return v;
}

void VertexFactory::recallAllVertices(void) {

    for (Vertex* v : onLoan)
        {
        v->clean();
        vertexPool.push_back(v);
        }
    onLoan.clear();
}

void VertexFactory::returnToPool(Vertex* v) {

    v->clean();
    vertexPool.push_back( v );
    onLoan.erase(v);
}
