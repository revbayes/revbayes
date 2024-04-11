#ifndef Polyhedron_hpp
#define Polyhedron_hpp

#include <gmpxx.h>   // a C++ wrapper around the GMP C implementation
#include <map>
#include <string>
#include <vector>
#include "Plane.h"


namespace RevBayesCore {

    class RandomNumberGenerator;
    class RateMatrix;
    class Vector;
    class Vertex;

    // a 4 X 4 matrix of GMP rationals used in determining if a point is in a tetrahedron
    class MpqMatrix {
        
        public:
            mpq_class&          operator()(size_t r, size_t c) { return this->m[r * 4 + c]; }
            const mpq_class&    operator()(size_t r, size_t c) const { return this->m[r * 4 + c]; }
            void                print(void);
            void                set(int idx, Vector* v);
            
        private:
            mpq_class           m[16];
    };

    // used to pick a random point in the polyhedron and the value in the vector_volume_map
    struct VectorInfo {
        
        mpq_class           volume;
        double              alphaC;
    };

    // Do you really want to type out the entire map declaration or would you rather
    // use a short-hand version? I thought so.
    typedef std::map<Plane*,std::vector<Vertex*> >                     plane_vertex_map;
    typedef std::map< std::pair<Plane*,Plane*>, std::vector<Vertex*> > line_vertex_map;
    typedef std::map<Vector*,VectorInfo>                               vector_volume_map;



    class Polyhedron {
        
        /**
         * This class is used to calculate the probability of proposing a non-reversible model from
         * a time reversible model for reversible jump MCMC. Proposing a non-reversible model from
         * a time reversible one involves the generation of three uniform(0,1) random variables,
         * u1, u2, and u3. The probability density would appear to be 1 X 1 X 1 = 1, but this is not
         * the case because there are constraints on the values that result in a valid non-reversible
         * model (e.g., one with all q_{ij} >= 0 (i != j). The constraints are (in Latex form),
         *
         * 0 \leq w_{AC} + w_{CG} (2 u_1 - 1) + w_{CT}(2 u_2 - 1) \leq 2 w_{AC}
         * 0 \leq w_{AG} - w_{CG} (2 u_1 - 1) + w_{GT}(2 u_3 - 1) \leq 2 w_{AG}
         * 0 \leq w_{AT} - w_{CT} (2 u_2 - 1) - w_{GT}(2 u_3 - 1) \leq 2 w_{AT}
         *
         * where w_{ij} = \pi_i q_{ij} is the average rate of changing from nucleotide i to j. The
         * class implements the constraints as planes (in 3D space) with the coordinates being (u1,u2,u3).
         * Six of the planes represent the bounds of the uniform(0,1) random variables (i.e., they
         * represent the unit cube). The other six planes represent the six constraints, above. The
         * planes for the unit cube never change and are initialized on instantiation of a Polyhedron
         * object. The other six planes (from the constraints, above) are initialized each time the
         * lnProbability* function is called. In fact, the lnProbability* functions take as a parameter
         * a vector of weights: w_{AC}, w_{AG}, w_{AT}, w_{CG}, w_{CT}, w_{GT}.
         *
         * The constraints form a polyhedron, with the details depending on the weights. The goal
         * is to calculate the volume of this polyhedron because the inverse of the volume is
         * the probability of proposing the non-reversible model from the time reversible model,
         * with the time reversible model being the point in the very center of the cube with
         * coordinates (1/2, 1/2, 1/2). All combinations of three of the 12 planes are checked
         * to see if they intersect (at a point, of course). If they do intersect, the intersecting
         * point is checked to see if it falls on the constraints described above. If the point
         * does satisfy the constraint, then it is a vertex in the polyhedron.
         *
         * After all of the vertices are found, they are organized into facets. The vertices for
         * a facet are ordered in clockwise fashion. This allows easy triangulation of each facet
         * to form tetrahedra (a polyhedron with four facets and four vertices). The tetrahedra
         * that compose the polyhedron are formed such that they all share the center point
         * (1/2, 1/2, 1/2). Calculating the volume of a tetrahedron is pretty easy, just requiring
         * the calculation of the determinant of a matrix containing the vertices.
         *
         * All calculations are done with the GMP rational class. The GMP library is a high-
         * performance C/C++ library for arbitrary-precision arithmetic. I think the rational
         * class from the GMP library uses two of the integer objects to represent the numerator
         * and denominator. Why do I use the GMP class at all? I do it to allow for exact
         * bounds checks and to avoid problems with threshold comparisons (which a normal implementation
         * using doubles would require). I took care to make certain that all of the key calculations
         * are always with rational numbers, never requiring the use of doubles or even
         * mpf_class (real number representation in the GMP library).
         *
         * Besides calculating the volume of the polyhedron, this class will also generate
         * random points from the polyhedron. The algorithm for doing this is pretty cool. You
         * can represent a point in a tetrahedron with vertices V1, V2, V3, and V4 as a mixture
         * of those four vertices:
         *
         * P = V1 * b1 + V2 * b2 + V3 * b3 + V4 * (1-s-t-u)b4
         *
         * Here, b1, b2, b3, and b4 are called the barycentric coordinates, with b1+b2+b3+b4=1.
         * You can generate a point uniformly and at random from a tetrahedron by generating
         * the barycentric coordinates from a Dirichlet(1, 1, 1, 1) probability distribution.
         * This code also allows you to bias the random point to be nearer to the center of
         * the cube (1/2, 1/2, 1/2). Let's denote V1 as the center vertex with coordinates
         * (1/2, 1/2, 1/2). Simply generate the barycentric coordinates from a
         * Dirichlet(alphaT, 1, 1, 1) distribution. If alphaT > 1, then the random points
         * will on average be closer to the center. The degree of clustering of the points
         * around the center will depend on the magnitude of alphaT. Interestingly, you could
         * bias the points to be nearer to the facets of the polyhedron than to the center by
         * setting alphaT less than 1. I don't know why you would want to do that, but there it is.
         *
         * Now that we know how to generate a random point from a tetrahedron, it's easy
         * to draw a point from the polyhedron. Simply pick a tetrahedron from the polyhedron
         * in proportion to its volume (compared to the volume of the polyhedron). Then, pick
         * a point from the tetrahedron, as described in the previous paragraph.
         *
         * This class can also calculate the probability of the random point. This was a rather
         * involved equation to derive, but it's pretty simple in the end. This works for the
         * uniform and biased ways of generating the random point.
         *
         * Again, the main way to interact with this class is through the functions
         * lnProbabilityForward and lnProbabilityReverse. Both functions take as input an
         * array of weights and a Vector (point). lnProbabilityForward generates a random
         * point from the polyhedron, initializing the point passed in as a reference to the
         * function as that random point, and returns the log probability density of that point.
         * The point for lnProbabilityReverse is assumed to be initialized. The code figures
         * out which tetrahedron the point resides in and returns the log probability
         * density of that point.
         *
         * The idea is that lnProbabilityForward is called elsewhere in the program when
         * one proposes a move from the time-reversible model to the non-reversible model.
         * The lnProbabilityReverse function is called when proposing a move from a non-reversible
         * model to a time reversible model, in which case we have to calculate the
         * probability of the imagined reverse move (from a time reversible to a non-reversible
         * model). This is why the lnProbabilityReverse function assumes the point that
         * is passed in as a function parameter is already initialized. The idea is the
         * point is initialized to represent the non-reversible model that one is moving from
         * when proposing the time reversible model.
         *
         * One more thing (really!). This class uses not only the GMP library, but implementations
         * of a Vector, Vertex, and Plane. Also, Vertex objects are managed by a singleton
         * class (called VertexFactory). The Vector, Vertex, and Plane classes all use GMP
         * rational numbers in their innards.
         *
         * @copyright Copyright 2009-
         * @author The RevBayes Development Core Team (John Huelsenbeck)
         * @since 2014-11-18, version 1.0
         */
        
    public:
                            Polyhedron(void);
                            Polyhedron(const Polyhedron& p) = delete;
        void                certify(void);
        void                setAlphaT(double x) { alphaT = x; }
        double              lnProbabilityForward(std::vector<mpq_class>& W, Vector& pt);
        double              lnProbabilityReverse(std::vector<mpq_class>& W, Vector& pt);
        
    private:
        void                calculateFacetVolume(Plane* pln, std::vector<Vertex*>& vertices, mpq_class& vol);
        void                calculateTetrahedronVolume(Vector* v1, Vector* v2, Vector* v3, mpq_class& vol);
        void                clearTetrahedraMap(void);
        void                computeLandU(MpqMatrix& aMat, MpqMatrix& lMat, MpqMatrix& uMat);
        mpq_class           det(MpqMatrix& m);
        Vertex*             findOtherVertex(Vertex* from, Vertex* v, Plane* pln);
        void                initializeFacets(void);
        void                initializePlanes(void);
        void                insertVertex(Plane* p1, Plane* p2, Vertex* v);
        bool                intersect(Plane& plane1, Plane& plane2, Plane& plane3, Vector& intersection);
        bool                isInTetrahedron(Vector* pt, Vector* center, Vector* v1, Vector* v2, Vector* v3, mpq_class& b1, mpq_class& b2, mpq_class& b3, mpq_class& b4);
        bool                isValid(Vector& pt);
        void                sampleTetrahedron(Plane* pln, Vector* center, Vector* v1, Vector* v2, Vector* v3, Vector& pt, VectorInfo& info);
        void                setWeights(std::vector<mpq_class>& W);
        
        mpq_class           wAC;                       // the parameters of the polyhedron are
        mpq_class           wAG;                       // six weights, w_{ij} = pi_i q_{ij}
        mpq_class           wAT;
        mpq_class           wCG;
        mpq_class           wCT;
        mpq_class           wGT;
        
        mpq_class           xzMaxA;                    // used when setting up the constraints and
        mpq_class           xzMaxB;                    // initializing the planes for those constraints
        mpq_class           xzMinA;
        mpq_class           xzMinB;
        mpq_class           xyMaxA;
        mpq_class           xyMaxB;
        mpq_class           xyMinA;
        mpq_class           xyMinB;
        mpq_class           yzMaxA;
        mpq_class           yzMaxB;
        mpq_class           yzMinA;
        mpq_class           yzMinB;
        
        mpq_class           zeroQ;                    // commonly used constants
        mpq_class           oneQ;
        mpq_class           oneHalfQ;
        mpq_class           twoQ;
        
        mpq_class           aQ;                       // used in calculateTetrahedronVolume and
        mpq_class           bQ;                       // are the elements of a 3 X 3 matrix for
        mpq_class           cQ;                       // which the determinant will be calculated
        mpq_class           dQ;
        mpq_class           eQ;
        mpq_class           fQ;
        mpq_class           gQ;
        mpq_class           hQ;
        mpq_class           iQ;
        
        Vector              center;                   // vectors representing the constraints directly used
        Vector              xzMinA_Zero_Zero;         // for initializing the planes for those constraints
        Vector              xzMaxA_Zero_Zero;
        Vector              xzMinA_One_Zero;
        Vector              xzMaxA_One_Zero;
        Vector              xzMinB_Zero_One;
        Vector              xzMaxB_Zero_One;
        Vector              zero_xyMaxA_Zero;
        Vector              zero_xyMinA_Zero;
        Vector              one_xyMaxB_Zero;
        Vector              one_xyMinB_Zero;
        Vector              zero_xyMaxA_One;
        Vector              zero_xyMinA_One;
        Vector              zero_Zero_yzMaxA;
        Vector              zero_Zero_yzMinA;
        Vector              zero_One_yzMaxB;
        Vector              zero_One_yzMinB;
        Vector              one_Zero_yzMaxA;
        Vector              one_Zero_yzMinA;
        
        Plane               front;                    // the planes representing the constraints
        Plane               back;                     // they are persistent until this polyhedron object dies
        Plane               top;
        Plane               bottom;
        Plane               left;
        Plane               right;
        Plane               xz1;
        Plane               xz2;
        Plane               xy1;
        Plane               xy2;
        Plane               yz1;
        Plane               yz2;
        
        std::vector<Plane*> planes;                   // a persistent array holding pointers to the planes
        plane_vertex_map    verticesMap;              // a map holding the vertices on a plane (i.e., the facets of the polyhedron)
        line_vertex_map     linesMap;                 // a map where the key is two planes and the value is an array
                                                      // of vertices for sorting out facets
        bool                randomlySample;           // should a random point be sampled?
        bool                pointFoundInPolyhedron;   // for a sanity check that the point was found only once in the polyhedron
        vector_volume_map   tetrahedra;               // a map where the key is a random vertex and the value is the volume of the tetrahedron
        Vector              randomPoint;              // the point representing the non-reversible model
        double              alphaT;                   // the parameter of the Dirichlet(alphaT,1,1,1) distribution for biasing the random point
        double              alphaC;                   // the barycentric coordinate for the center vertex (for biasing the random point)
            
        mpq_class           sumJacobians;             // the sum of the determinants of the tetrahedra, for volume/probability calculations
    };
}

#endif
