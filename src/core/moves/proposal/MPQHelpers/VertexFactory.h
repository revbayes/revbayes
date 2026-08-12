#ifndef VertexFactory_H
#define VertexFactory_H

#include <set>
#include <vector>
#include "Vector.h"


namespace RevBayesCore {

    // this header uses Vertex only through a pointer; it used to rely on every
    // caller having included Vertex.h first, which happened to be true
    class Vertex;

    class VertexFactory {
        
        /**
         * A class to manage the vertices used in constructing polyhedra. This class
         * hands out and retrieves instances of the Vertex class.
         *
         * This used to be a singleton. It no longer is. The pool and the two sets
         * below are mutable state, and a singleton shares that state across every
         * thread in the process, so two Metropolis-coupled chains proposing a
         * reversible-jump move at the same time would hand each other the same
         * Vertex. A Polyhedron is the only thing that needs vertices, so each
         * Polyhedron now owns a factory of its own and the state is no longer shared.
         *
         * @copyright Copyright 2009-
         * @author The RevBayes Development Core Team (John Huelsenbeck)
         * @since 2014-11-18, version 1.0
         */
        
    public:
                                VertexFactory(void);
                               ~VertexFactory(void);
                                VertexFactory(const VertexFactory&) = delete;
        VertexFactory&          operator=(const VertexFactory&) = delete;
        void                    drainPool(void);
        Vertex*                 getVertex(void);
        Vertex*                 getVertex(Vertex& v);
        Vertex*                 getVertex(Vector& v);
        int                     getNumAllocated(void) { return (int)allocatedVertices.size(); }
        int                     getNumOnLoan(void) { return (int)onLoan.size(); }
        void                    returnToPool(Vertex* nde);
        void                    recallAllVertices(void);
        
    private:
        std::vector<Vertex*>    vertexPool;
        std::set<Vertex*>       allocatedVertices;
        std::set<Vertex*>       onLoan;
    };
}

#endif
