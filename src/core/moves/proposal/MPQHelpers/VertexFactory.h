#ifndef VertexFactory_H
#define VertexFactory_H

#include <set>
#include <vector>
#include "Vector.h"


namespace RevBayesCore {

    class VertexFactory {
        
        /**
         * A singleton class to manage vertices which are used in constructing
         * polyhedra. This class hands out and retrieves instances of the
         * Vertex class.
         *
         * @copyright Copyright 2009-
         * @author The RevBayes Development Core Team (John Huelsenbeck)
         * @since 2014-11-18, version 1.0
         */
        
    public:
        static VertexFactory&   vertexFactoryInstance(void)
                                    {
                                    static VertexFactory singleNodeFactory;
                                    return singleNodeFactory;
                                    }
        void                    drainPool(void);
        Vertex*                 getVertex(void);
        Vertex*                 getVertex(Vertex& v);
        Vertex*                 getVertex(Vector& v);
        int                     getNumAllocated(void) { return (int)allocatedVertices.size(); }
        void                    returnToPool(Vertex* nde);
        void                    recallAllVertices(void);
        
    private:
                                VertexFactory(void);
                                VertexFactory(const VertexFactory&);
                                VertexFactory& operator=(const VertexFactory&);
                               ~VertexFactory(void);
        std::vector<Vertex*>    vertexPool;
        std::set<Vertex*>       allocatedVertices;
        std::set<Vertex*>       onLoan;
    };
}

#endif
