#ifndef Vertex_H
#define Vertex_H

#include "Vector.h"



namespace RevBayesCore {

    class Vertex : public Vector {
        
        /**
         * A light-weight class to represent a vertex in a facet in 3D. This
         * class inherits from Vector, which has the x,y,z coordinates of the
         * Vertex. This class simply adds pointers to the neighboring Vertices.
         *
         * @copyright Copyright 2009-
         * @author The RevBayes Development Core Team (John Huelsenbeck)
         * @since 2014-11-18, version 1.0
         */
        
    public:
                        Vertex(void);
                        Vertex(const Vector& v);
                        Vertex(Vertex& v);
                        Vertex(mpq_class& xq, mpq_class& yq, mpq_class& zq);
        bool            operator==(const Vertex& rhs) const;
        void            clean(void);
        Vertex*         getFrom(void) { return from; }
        Vertex*         getTo(void) { return to; }
        void            setFrom(Vertex* v) { from = v; }
        void            setTo(Vertex* v) { to = v; }
        
    private:
        Vertex*         from;
        Vertex*         to;
    };
}

#endif
