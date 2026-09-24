#include <iomanip>
#include <iostream>
#include <map>

#include "DistributionDirichlet.h"
#include "RateMatrix_MPQ.h"
#include "Polyhedron.h"
#include "RandomNumberGenerator.h"
#include "RandomNumberFactory.h"
#include <cmath>
#include "RbConstants.h"
#include "RbException.h"
#include "RbMathFunctions.h"
#include "Vertex.h"
#include "VertexFactory.h"


using namespace RevBayesCore;



void MpqMatrix::print(void) {

    for (int i=0; i<4; i++)
        {
        for (int j=0; j<4; j++)
            std::cout << m[i*4+j].get_d() << " ";
        std::cout << std::endl;
        }
}

void MpqMatrix::set(int idx, Vector* v) {

    m[idx*4 + 0] = v->getX();
    m[idx*4 + 1] = v->getY();
    m[idx*4 + 2] = v->getZ();
}

Polyhedron::Polyhedron(void) {

    // this polyhedron owns its vertex pool; it is not shared with any other
    // polyhedron, and in particular not with a polyhedron on another chain
    vertexFactory = new VertexFactory;

    // initialize commonly used constants
    zeroQ     = 0;
    oneQ      = 1;
    oneHalfQ  = 1;
    oneHalfQ /= 2;
    twoQ      = 2;
    
    randomlySample = false;
    alphaT = 1.0; // the default value for the Dirichlet(alphaT,1,1,1) used to randomly draw points from tetrahedra

    /* The center vertex is constant and must remain so; see the note in the header.
       (1/2,1/2,1/2) is the time reversible matrix, which is a valid point of the
       polyhedron whatever the weights are, so measuring every tetrahedron from here
       cannot produce a degenerate volume. */
    reversibleCenter.setX(oneHalfQ);
    reversibleCenter.setY(oneHalfQ);
    reversibleCenter.setZ(oneHalfQ);
    center           = reversibleCenter;
    desiredCenter    = reversibleCenter;
    useDesiredCenter = false;
    numCenterClipped = 0;

    // stop this fraction of the way to the boundary, so the centre is never on a facet
    centerShrink  = 99;
    centerShrink /= 100;

    num_degenerate_weights = 0;
    num_point_not_valid    = 0;
    num_point_not_located  = 0;
    num_bad_alpha_c        = 0;

    // set up fixed planes of cube
    front.set( Vector(zeroQ, zeroQ, zeroQ), Vector(zeroQ,  oneQ, zeroQ), Vector( oneQ,  oneQ, zeroQ) );
    back.set( Vector(zeroQ, zeroQ,  oneQ), Vector(zeroQ,  oneQ,  oneQ), Vector( oneQ,  oneQ,  oneQ) );
    top.set( Vector(zeroQ,  oneQ, zeroQ), Vector(zeroQ,  oneQ,  oneQ), Vector( oneQ,  oneQ,  oneQ) );
    bottom.set( Vector(zeroQ, zeroQ, zeroQ), Vector(zeroQ, zeroQ,  oneQ), Vector( oneQ, zeroQ,  oneQ) );
    left.set( Vector(zeroQ, zeroQ, zeroQ), Vector(zeroQ,  oneQ, zeroQ), Vector(zeroQ,  oneQ,  oneQ) );
    right.set( Vector( oneQ, zeroQ, zeroQ), Vector( oneQ,  oneQ, zeroQ), Vector( oneQ,  oneQ,  oneQ) );

    // planes persist and live in the planes vector
    planes.push_back(&front);
    planes.push_back(&back);
    planes.push_back(&top);
    planes.push_back(&bottom);
    planes.push_back(&left);
    planes.push_back(&right);
    planes.push_back(&xz1);
    planes.push_back(&xz2);
    planes.push_back(&xy1);
    planes.push_back(&xy2);
    planes.push_back(&yz1);
    planes.push_back(&yz2);
    
    // initialize the volume
    sumJacobians = 0;
}

/**
 * Find the volume of a polyhedron formed by the facet vertices and the center point, (1/2,1/2,1/2).
 *
 * \param[in]    pln       The plane for the facet
 * \param[in]    vertices  The vertices of the facet
 * \param[in]    vol       A reference to the volume, which will be added to
 */
Polyhedron::~Polyhedron(void) {

    // frees every vertex the factory ever allocated
    delete vertexFactory;
}

void Polyhedron::calculateFacetVolume(Plane* pln, std::vector<Vertex*>& vertices, mpq_class& vol) {

    // loop over triangulations of the facet
    Vector pt;
    Vertex* v1 = vertices[0];
    Vertex* p = v1->getTo();
    do
        {
        // form triangulation, with v1 common to all triangulations of this facet
        Vertex* v2 = p;
        Vertex* v3 = p->getTo();
        
        // calculate volume
        mpq_class tetrahedronVolume;
        calculateTetrahedronVolume(v1, v2, v3, tetrahedronVolume);
        vol += tetrahedronVolume;
        
        if (randomlySample == true)
            {
            // add a random point from tetrahedron to tetrahedra map for later use
            Vector* newV = new Vector(pt);
            VectorInfo info;
            info.volume = tetrahedronVolume;
            sampleTetrahedron(pln, &center, v1, v2, v3, *newV, info);
            tetrahedra.insert( std::make_pair(newV,info) );
            //tetrahedra.insert( std::make_pair(newV,tetrahedronVolume) );
            }
        else
            {
            // check if the point is in this tetrahedron
            mpq_class b1;
            mpq_class b2;
            mpq_class b3;
            mpq_class b4;
            if (isInTetrahedron(&randomPoint, &center, v1, v2, v3, b1, b2, b3, b4) == true)
                {
                //std::cout << b1.get_d() << " " << b2.get_d() << " " << b3.get_d() << " " << b4.get_d() << " " << std::endl;
                alphaC = b1.get_d();
                if (pointFoundInPolyhedron == true)
                    throw(RbException("Polyhedron: Point already found in polyhedron"));
                pointFoundInPolyhedron = true;
                mpq_class x = center.getX() * b1 + v1->getX() * b2 + v2->getX() * b3 + v3->getX() * b4;
                mpq_class y = center.getY() * b1 + v1->getY() * b2 + v2->getY() * b3 + v3->getY() * b4;
                mpq_class z = center.getZ() * b1 + v1->getZ() * b2 + v2->getZ() * b3 + v3->getZ() * b4;
                Vector t;
                t.set(x, y, z); // Vector pt = v0*a + v1*s + v2*t + v3*u;
                if (t != randomPoint)
                    throw(RbException("Polyhedron: Test point does not match random point"));
                }
            }
        
        // on to the next triangulation
        p = v3;
        } while (p->getTo() != v1);
}

/**
 * Calculate the volume (up to a factor of 1/6) of a tetrahedron. I know that a tetrahedron is
 * formed from four vertices, but here we pass in three vectors, which are calculated by subtracting
 * one of the vertices (in this case the center point) from the other three.
 *
 * \param[in]    v1        The first vector
 * \param[in]    v2        The second vector
 * \param[in]    v3        The third vector
 * \param[in]    vol       A reference to the volume of the tetrahedron, passed in as a reference
 *                         to avoid copying of a mpq_class object
 */
void Polyhedron::calculateTetrahedronVolume(Vector* v1, Vector* v2, Vector* v3, mpq_class& vol) {

    /* The fourth vertex of every tetrahedron is the center of the polyhedron, so the
       volume is one sixth of the determinant of the three edge vectors leading away
       from it. This used to subtract oneHalfQ, which is the same thing only as long
       as the center is the time reversible point. It is not any more, and getting
       this wrong would not throw: the tetrahedra would simply be measured about the
       wrong apex, sumJacobians would be wrong, and the proposal density used in the
       Hastings ratio would be quietly wrong with it. */
    aQ = v1->getX() - center.getX();
    bQ = v2->getX() - center.getX();
    cQ = v3->getX() - center.getX();
    dQ = v1->getY() - center.getY();
    eQ = v2->getY() - center.getY();
    fQ = v3->getY() - center.getY();
    gQ = v1->getZ() - center.getZ();
    hQ = v2->getZ() - center.getZ();
    iQ = v3->getZ() - center.getZ();
    
    // volume is 1/6 of the determinant
    vol = (aQ * eQ * iQ) - (aQ * fQ * hQ) - (bQ * dQ * iQ) + (bQ * fQ * gQ) + (cQ * dQ * hQ) - (cQ * eQ * gQ);
    //vol /= 6;                  // this will be taken care of in the probability density with a Gamma factor of 3! = 6
    if (vol < 0)
        vol = -vol;
}

void Polyhedron::certify(void) {

    RandomNumberGenerator* rng = GLOBAL_RNG;
    
    std::vector<double> alphaPi(4, 0.25);
    std::vector<double> alphaR(6, 1.0);
    
    std::vector<double> r = RbStatistics::Dirichlet::rv(alphaR, *rng);
    std::vector<double> pi = RbStatistics::Dirichlet::rv(alphaPi, *rng);
    std::vector<mpq_class> wts(6);
    mpq_class sum = 0;
    for (int i=0, k=0; i<4; i++)
        {
        for (int j=i+1; j<4; j++)
            {
            wts[k] = (pi[i] * r[k] + pi[j] * r[k]) / 2.0;
            sum += wts[k];
            k++;
            }
        }
    sum *= 2;
    for (int i=0; i<6; i++)
        wts[i] /= sum;
    setWeights(wts);
    std::cout << sumJacobians.get_d() << std::endl;
}

void Polyhedron::clearTetrahedraMap(void) {

    for (vector_volume_map::iterator it = tetrahedra.begin(); it != tetrahedra.end(); it++)
        delete it->first;
    tetrahedra.clear();
}

/**
 * Hard-wired determinant of a 4 X 4 matrix of GMP rationals.
 */
mpq_class Polyhedron::det(MpqMatrix& m) {

    mpq_class d = m(0,3) * m(1,2) * m(2,1) * m(3,0) - m(0,2) * m(1,3) * m(2,1) * m(3,0) - m(0,3) * m(1,1) * m(2,2) * m(3,0) +
                  m(0,1) * m(1,3) * m(2,2) * m(3,0) + m(0,2) * m(1,1) * m(2,3) * m(3,0) - m(0,1) * m(1,2) * m(2,3) * m(3,0) -
                  m(0,3) * m(1,2) * m(2,0) * m(3,1) + m(0,2) * m(1,3) * m(2,0) * m(3,1) + m(0,3) * m(1,0) * m(2,2) * m(3,1) -
                  m(0,0) * m(1,3) * m(2,2) * m(3,1) - m(0,2) * m(1,0) * m(2,3) * m(3,1) + m(0,0) * m(1,2) * m(2,3) * m(3,1) +
                  m(0,3) * m(1,1) * m(2,0) * m(3,2) - m(0,1) * m(1,3) * m(2,0) * m(3,2) - m(0,3) * m(1,0) * m(2,1) * m(3,2) +
                  m(0,0) * m(1,3) * m(2,1) * m(3,2) + m(0,1) * m(1,0) * m(2,3) * m(3,2) - m(0,0) * m(1,1) * m(2,3) * m(3,2) -
                  m(0,2) * m(1,1) * m(2,0) * m(3,3) + m(0,1) * m(1,2) * m(2,0) * m(3,3) + m(0,2) * m(1,0) * m(2,1) * m(3,3) -
                  m(0,0) * m(1,2) * m(2,1) * m(3,3) - m(0,1) * m(1,0) * m(2,2) * m(3,3) + m(0,0) * m(1,1) * m(2,2) * m(3,3);
    return d;
}

Vertex* Polyhedron::findOtherVertex(Vertex* from, Vertex* v, Plane* pln) {

    for (auto lne : linesMap)
        {
        if (lne.first.first == pln || lne.first.second == pln)
            {
            bool inList = false;
            for (int i=0, n=(int)lne.second.size(); i<n; i++)
                {
                if (lne.second[i] == v)
                    {
                    inList = true;
                    break;
                    }
                }
                
            if (inList == true)
                {
                if (lne.second[0] != from && lne.second[1] != from)
                    {
                    if (v == lne.second[0])
                        return lne.second[1];
                    else
                        return lne.second[0];
                    }
                }
            
            }
        }
    return nullptr;
}

/**
 * Set up all of the facets by ordering the vertices for each. Here, the
 * volume of each facet (formed by the vertices of the facet and the
 * vertex at the center of the polyhedron at 1/2, 1/2, 1/2) is calculated
 * by triangulating each facet into its component tetrahedra.
 */
void Polyhedron::initializeFacets(void) {
    
    clearTetrahedraMap();

    mpq_class vol;
    sumJacobians = 0.0;
    // loop over all planes/facets of the polyhedron
    for (auto pln : verticesMap)
        {
        // clean vertices
        for (int i=0, n=(int)pln.second.size(); i<n; i++)
            {
            pln.second[i]->setTo(nullptr);
            pln.second[i]->setFrom(nullptr);
            }
            
        // order the vertices
        Vertex* first = pln.second[0];
        Vertex* v = first;
        do {
            Vertex* nextV = findOtherVertex(v->getFrom(), v, pln.first);
            v->setTo(nextV);
            nextV->setFrom(v);
            v = nextV;
            } while (v != first);
            
        // calculate the facet volume and (potentially) randomly sample
        calculateFacetVolume(pln.first, pln.second, vol);
        }
    sumJacobians = vol;
    
    if (randomlySample == true)
        {
        RandomNumberGenerator* rng = GLOBAL_RNG;
        mpq_class u = rng->uniform01();
        u *= sumJacobians;
        mpq_class sumVol;
        for (auto tet : tetrahedra)
            {
            sumVol += tet.second.volume;
            if (u < sumVol)
                {
                randomPoint.set(tet.first->getX(), tet.first->getY(), tet.first->getZ());
                alphaC = tet.second.alphaC;
                break;
                }
            }
        clearTetrahedraMap();
        }
}

/**
 * Find all the vertices of the polyhedron.
 */
void Polyhedron::initializePlanes(void) {

    // make planes that will slice up the cube
    // u1 -> (-wAG + wCG - wGT + 2.0 * u3 * wGT) / (2.0 * wCG)   max in x,z
    // u1 -> ( wAG + wCG - wGT + 2.0 * u3 * wGT) / (2.0 * wCG)   min in x,z
    // u2 -> ( wAC + wCG - 2.0 * u1 * wCG + wCT) / (2.0 * wCT)   max in x,y
    // u2 -> (-wAC + wCG - 2.0 * u1 * wCG + wCT) / (2.0 * wCT)   min in x,y
    // u3 -> (-wAT + wCT - 2.0 * u2 * wCT + wGT) / (2.0 * wGT)   max in y,z
    // u3 -> ( wAT + wCT - 2.0 * u2 * wCT + wGT) / (2.0 * wGT)   min in y,z
        
    xzMaxA = (-wAG + wCG - wGT + twoQ * zeroQ * wGT) / (twoQ * wCG);
    xzMaxB = (-wAG + wCG - wGT + twoQ * oneQ  * wGT) / (twoQ * wCG);
    xzMinA = ( wAG + wCG - wGT + twoQ * zeroQ * wGT) / (twoQ * wCG);
    xzMinB = ( wAG + wCG - wGT + twoQ * oneQ  * wGT) / (twoQ * wCG);

    xyMaxA = (-wAC + wCG - twoQ * zeroQ * wCG + wCT) / (twoQ * wCT);
    xyMaxB = (-wAC + wCG - twoQ * oneQ  * wCG + wCT) / (twoQ * wCT);
    xyMinA = ( wAC + wCG - twoQ * zeroQ * wCG + wCT) / (twoQ * wCT);
    xyMinB = ( wAC + wCG - twoQ * oneQ  * wCG + wCT) / (twoQ * wCT);

    yzMaxA = (-wAT + wCT - twoQ * zeroQ * wCT + wGT) / (twoQ * wGT);
    yzMaxB = (-wAT + wCT - twoQ * oneQ  * wCT + wGT) / (twoQ * wGT);
    yzMinA = ( wAT + wCT - twoQ * zeroQ * wCT + wGT) / (twoQ * wGT);
    yzMinB = ( wAT + wCT - twoQ * oneQ  * wCT + wGT) / (twoQ * wGT);

    xzMaxA_Zero_Zero.set(xzMaxA, zeroQ, zeroQ); // Vector(xzMaxA, zeroQ, zeroQ)
    xzMinA_Zero_Zero.set(xzMinA, zeroQ, zeroQ); // Vector(xzMinA, zeroQ, zeroQ)
    xzMaxA_One_Zero.set(xzMaxA, oneQ, zeroQ);   // Vector(xzMaxA,  oneQ, zeroQ)
    xzMinA_One_Zero.set(xzMinA, oneQ, zeroQ);   // Vector(xzMinA,  oneQ, zeroQ)
    xzMaxB_Zero_One.set(xzMaxB, zeroQ, oneQ);   // Vector(xzMaxB, zeroQ,  oneQ)
    xzMinB_Zero_One.set(xzMinB, zeroQ, oneQ);   // Vector(xzMinB, zeroQ,  oneQ)
    zero_xyMaxA_Zero.set(zeroQ, xyMaxA, zeroQ); // Vector(zeroQ, xyMaxA, zeroQ)
    zero_xyMinA_Zero.set(zeroQ, xyMinA, zeroQ); // Vector(zeroQ, xyMinA, zeroQ)
    one_xyMaxB_Zero.set(oneQ, xyMaxB, zeroQ);   // Vector(oneQ, xyMaxB, zeroQ)
    one_xyMinB_Zero.set(oneQ, xyMinB, zeroQ);   // Vector(oneQ, xyMinB, zeroQ)
    zero_xyMaxA_One.set(zeroQ, xyMaxA, oneQ);   // Vector(zeroQ, xyMaxA, oneQ)
    zero_xyMinA_One.set(zeroQ, xyMinA, oneQ);   // Vector(zeroQ, xyMinA, oneQ)
    zero_Zero_yzMaxA.set(zeroQ, zeroQ, yzMaxA); // Vector(zeroQ, zeroQ, yzMaxA)
    zero_Zero_yzMinA.set(zeroQ, zeroQ, yzMinA); // Vector(zeroQ, zeroQ, yzMinA)
    zero_One_yzMaxB.set(zeroQ, oneQ, yzMaxB);   // Vector(zeroQ, oneQ, yzMaxB)
    zero_One_yzMinB.set(zeroQ, oneQ, yzMinB);   // Vector(zeroQ, oneQ, yzMinB)
    one_Zero_yzMaxA.set(oneQ, zeroQ, yzMaxA);   // Vector(oneQ, zeroQ, yzMaxA)
    one_Zero_yzMinA.set(oneQ, zeroQ, yzMinA);   // Vector(oneQ, zeroQ, yzMinA)

    // set up non-fixed planes
    xz1.set( xzMaxA_Zero_Zero, xzMaxB_Zero_One, xzMaxA_One_Zero );
    xz2.set( xzMinA_Zero_Zero, xzMinB_Zero_One, xzMinA_One_Zero );
    xy1.set( zero_xyMaxA_Zero, one_xyMaxB_Zero, zero_xyMaxA_One );
    xy2.set( zero_xyMinA_Zero, one_xyMinB_Zero, zero_xyMinA_One );
    yz1.set( zero_Zero_yzMaxA, zero_One_yzMaxB, one_Zero_yzMaxA );
    yz2.set( zero_Zero_yzMinA, zero_One_yzMinB, one_Zero_yzMinA );
    
    // empty out the tetrahedra map in preparation for finding random points
    // in triangulations of the polyhedron
    if (randomlySample == true)
        clearTetrahedraMap();
        
    // note that this checks all 12 choose 3 combinations of planes for intersection even though
    // six pairs of the planes are parallel to one another!
    VertexFactory& vf = *vertexFactory;
    verticesMap.clear();
    linesMap.clear();
    for (int i=0, n1 = (int)planes.size(); i<n1; i++)
        {
        for (int j=i+1, n2 = (int)planes.size(); j<n2; j++)
            {
            for (int k=j+1, n3 = (int)planes.size(); k<n3; k++)
                {
                if (i != j && i != k && j != k)
                    {
                    Vertex* intersectionPoint = vf.getVertex();
                    bool planesIntersect = intersect(*planes[i], *planes[j], *planes[k], *intersectionPoint);
                    if (planesIntersect == true && isValid(*intersectionPoint) == true)
                        {
                        // add intersection Vector to planes map
                        plane_vertex_map::iterator it = verticesMap.find(planes[i]);
                        if (it == verticesMap.end())
                            {
                            std::vector<Vertex*> vec;
                            vec.push_back(intersectionPoint);
                            verticesMap.insert( std::make_pair(planes[i],vec) );
                            }
                        else
                            it->second.push_back(intersectionPoint);

                        it = verticesMap.find(planes[j]);
                        if (it == verticesMap.end())
                            {
                            std::vector<Vertex*> vec;
                            vec.push_back(intersectionPoint);
                            verticesMap.insert( std::make_pair(planes[j],vec) );
                            }
                        else
                            it->second.push_back(intersectionPoint);

                        it = verticesMap.find(planes[k]);
                        if (it == verticesMap.end())
                            {
                            std::vector<Vertex*> vec;
                            vec.push_back(intersectionPoint);
                            verticesMap.insert( std::make_pair(planes[k],vec) );
                            }
                        else
                            it->second.push_back(intersectionPoint);
                            
                            
                        insertVertex(planes[i], planes[j], intersectionPoint);
                        insertVertex(planes[i], planes[k], intersectionPoint);
                        insertVertex(planes[j], planes[k], intersectionPoint);
                        }
                        
                    }
                }
            }
        }
        

    // set up facets
    initializeFacets();
    
    // clean up
    vf.recallAllVertices();
}

void Polyhedron::insertVertex(Plane* p1, Plane* p2, Vertex* v) {

    std::pair<Plane*,Plane*> key(p1, p2);
    if (p2 < p1)
        key = std::make_pair(p2, p1);
        
    line_vertex_map::iterator it = linesMap.find(key);
    if (it == linesMap.end())
        {
        std::vector<Vertex*> vec;
        vec.push_back(v);
        linesMap.insert( std::make_pair(key,vec) );
        }
    else
        {
        it->second.push_back(v);
        }
}

/**
 * Test whether three planes intersect and, if so, initialize the intersection point
 *
 * \param[in]    plane1       The first plane
 * \param[in]    plane2       The second plane
 * \param[in]    plane3       The third plane
 * \param[in]    intersection The intersecting point of the three planes
 * \return  A boolean indicating whether or not the planes intersect at a point
 */
bool Polyhedron::intersect(Plane& plane1, Plane& plane2, Plane& plane3, Vector& intersection) {

    const mpq_class& a1 = plane1.getA();
    const mpq_class& b1 = plane1.getB();
    const mpq_class& c1 = plane1.getC();
    const mpq_class& d1 = plane1.getD();
    const mpq_class& a2 = plane2.getA();
    const mpq_class& b2 = plane2.getB();
    const mpq_class& c2 = plane2.getC();
    const mpq_class& d2 = plane2.getD();
    const mpq_class& a3 = plane3.getA();
    const mpq_class& b3 = plane3.getB();
    const mpq_class& c3 = plane3.getC();
    const mpq_class& d3 = plane3.getD();

    mpq_class detA  = a1 * (b2 * c3 - c2 * b3) + b1 * (c2 * a3 - a2 * c3) + c1 * (a2 * b3 - b2 * a3);
    if (detA == 0)
        return false;
    mpq_class detAx = -d1 * (b2 * c3 - c2 * b3) - d2 * (b3 * c1 - c3 * b1) - d3 * (b1 * c2 - c1 * b2);
    mpq_class detAy = -d1 * (c2 * a3 - a2 * c3) - d2 * (c3 * a1 - a3 * c1) - d3 * (c1 * a2 - a1 * c2);
    mpq_class detAz = -d1 * (a2 * b3 - b2 * a3) - d2 * (a3 * b1 - b3 * a1) - d3 * (a1 * b2 - b1 * a2);
    
    mpq_class x = detAx / detA;
    mpq_class y = detAy / detA;
    mpq_class z = detAz / detA;
    
    intersection.setX(x);
    intersection.setY(y);
    intersection.setZ(z);
    
    return true;
}

/**
 * Test whether a point is in a tetrahedron and, if so, initialize its barycentric coordinates.
 *
 * \param[in]    pt       The point to be tested.
 * \param[in]    center   The first vertex (which also is the center vertex of the polyhedron)
 * \param[in]    v1       The first vertex on the facet
 * \param[in]    v2       The second vertex on the facet
 * \param[in]    v3       The third vertex on the facet
 * \param[in]    b1       A reference to the first barycentric coordinate, which might be initialized
 * \param[in]    b2       A reference to the second barycentric coordinate, which might be initialized
 * \param[in]    b3       A reference to the third barycentric coordinate, which might be initialized
 * \param[in]    b4       A reference to the fourth barycentric coordinate, which might be initialized
 * \return  A boolean indicating whether or not the point is in the tetrahedron
 */
bool Polyhedron::isInTetrahedron(Vector* pt, Vector* center, Vector* v1, Vector* v2, Vector* v3, mpq_class& b1, mpq_class& b2, mpq_class& b3, mpq_class& b4) {

    /* Let the tetrahedron have vertices

        V1 = (x1, y1, z1)
        V2 = (x2, y2, z2)
        V3 = (x3, y3, z3)
        V4 = (x4, y4, z4)

      and your test point be

        P = (x, y, z).

      Then the point P is in the tetrahedron if following five determinants all have the same sign.

         |x1 y1 z1 1|
    D0 = |x2 y2 z2 1|
         |x3 y3 z3 1|
         |x4 y4 z4 1|

         |x  y  z  1|
    D1 = |x2 y2 z2 1|
         |x3 y3 z3 1|
         |x4 y4 z4 1|

         |x1 y1 z1 1|
    D2 = |x  y  z  1|
         |x3 y3 z3 1|
         |x4 y4 z4 1|

         |x1 y1 z1 1|
    D3 = |x2 y2 z2 1|
         |x  y  z  1|
         |x4 y4 z4 1|

         |x1 y1 z1 1|
    D4 = |x2 y2 z2 1|
         |x3 y3 z3 1|
         |x  y  z  1|
         
    If the point is in the tetrahedron, the determinants, above, can be used to determine the
    barycentric coordinates of the point as bi = Di / D0, where bi is the i-th barycentric
    coordinate. */
         
    MpqMatrix m;
    
    // d0
    m.set(0, center);
    m.set(1, v1);
    m.set(2, v2);
    m.set(3, v3);
    for (int i=0; i<4; i++)
        m(i,3) = 1;
    mpq_class d0 = det(m);
    bool positiveD0 = false;
    if (d0 > 0)
        positiveD0 = true;
    
    // d1
    m.set(0, pt);
    mpq_class d1 = det(m);
    bool positiveDi = false;
    if (d1 > 0)
        positiveDi = true;
    if (positiveD0 != positiveDi)
        return false;

    // d2
    m.set(0, center);
    m.set(1, pt);
    mpq_class d2 = det(m);
    positiveDi = false;
    if (d2 > 0)
        positiveDi = true;
    if (positiveD0 != positiveDi)
        return false;

    // d3
    m.set(1, v1);
    m.set(2, pt);
    mpq_class d3 = det(m);
    positiveDi = false;
    if (d3 > 0)
        positiveDi = true;
    if (positiveD0 != positiveDi)
        return false;

    // d4
    m.set(2, v2);
    m.set(3, pt);
    mpq_class d4 = det(m);
    positiveDi = false;
    if (d4 > 0)
        positiveDi = true;
    if (positiveD0 != positiveDi)
        return false;

    // We made it! The point, pt, is in the tetrahedron! We can now calculate
    // its barycentric coordinates
    b1 = d1 / d0;
    b2 = d2 / d0;
    b3 = d3 / d0;
    b4 = d4 / d0;
        
    return true;
}

bool Polyhedron::isValid(Vector& pt) {

    mpq_class& x = pt.getX();
    mpq_class& y = pt.getY();
    mpq_class& z = pt.getZ();

    // first, check that the points are part of or inside of the unit cube
    if (x < 0 || x > 1)
        return false;
    if (y < 0 || y > 1)
        return false;
    if (z < 0 || z > 1)
        return false;

    // second, check that the points satisfy the constraints:
    // 0 \leq w_{AC} + w_{CG} (2 u_1 - 1) + w_{CT}(2 u_2 - 1) \leq 2 w_{AC}
    // 0 \leq w_{AG} - w_{CG} (2 u_1 - 1) + w_{GT}(2 u_3 - 1) \leq 2 w_{AG}
    // 0 \leq w_{AT} - w_{CT} (2 u_2 - 1) - w_{GT}(2 u_3 - 1) \leq 2 w_{AT}
    mpq_class v1 = wAC + wCG * (2 * x - 1) + wCT * (2 * y - 1);
    mpq_class v2 = wAG - wCG * (2 * x - 1) + wGT * (2 * z - 1);
    mpq_class v3 = wAT - wCT * (2 * y - 1) - wGT * (2 * z - 1);
    mpq_class v1Max = 2 * wAC;
    mpq_class v2Max = 2 * wAG;
    mpq_class v3Max = 2 * wAT;
    if ( (v1 >= 0) && (v1 <= v1Max) )
        ;
    else
        return false;
    if ( (v2 >= 0) && (v2 <= v2Max) )
        ;
    else
        return false;
    if ( (v3 >= 0) && (v3 <= v3Max) )
        ;
    else
        return false;
        
    return true;
}

/* The log density of drawing the point pt from the polyhedron defined by W.

   Returns negative infinity, rather than throwing, on every path that cannot
   produce a valid density. The caller treats that as a refusal to propose, which
   the MCMC turns into a rejection. Throwing was the wrong behaviour here for two
   reasons: RevBayes rethrows any exception that is not a MATH_ERROR, so a point
   landing on a facet could abort an analysis hours into a run; and an exception
   caught and turned into a rejection is invisible, which is precisely how a
   systematically failing jump can look like overwhelming evidence. */
double Polyhedron::lnProbabilityForward(std::vector<mpq_class>& W, Vector& pt) {

    if (weightsAreUsable(W) == false)
        return RbConstants::Double::neginf;

    // set up the polyhedron, which will also randomly sample and initialize pt
    randomlySample = true;
    setWeights(W);
    
    // initialize pt to the randomly selected point and do a sanity check
    pt.set(randomPoint.getX(), randomPoint.getY(), randomPoint.getZ());
    randomlySample = false;

    if (isValid(pt) == false)
        {
        num_point_not_valid++;
        reportFailure("the randomly drawn point is not inside the polyhedron", num_point_not_valid);
        return RbConstants::Double::neginf;
        }
    /* Note there is deliberately no check on pointFoundInPolyhedron here. That flag
       belongs to the reverse direction, where a supplied point has to be located in
       the triangulation; in the forward direction the point is drawn from a
       tetrahedron that was chosen first, so there is nothing to locate. */
    if (sumJacobians <= 0)
        {
        num_degenerate_weights++;
        reportFailure("the polyhedron has no volume", num_degenerate_weights);
        return RbConstants::Double::neginf;
        }
    
    // calculate the probability of the randomly proposed point, pt
    double lnProb = -log(sumJacobians.get_d());
    lnProb += RbMath::lnGamma(alphaT + 3.0) - RbMath::lnGamma(alphaT);

    /* The barycentric weight on the center vertex is zero for any point lying on a
       facet, and log(0) is negative infinity. With the default alphaT of one the
       coefficient is zero, and in IEEE arithmetic zero times negative infinity is
       NaN, not zero: the term would silently poison the density rather than drop
       out of it. So the term is only formed when it is actually needed. */
    if (alphaT != 1.0)
        {
        if (alphaC <= 0.0)
            {
            num_bad_alpha_c++;
            reportFailure("the point lies on a facet, so the concentration term of the proposal density is undefined", num_bad_alpha_c);
            return RbConstants::Double::neginf;
            }
        lnProb += (alphaT - 1.0) * log(alphaC);
        }

    if (std::isfinite(lnProb) == false)
        {
        num_bad_alpha_c++;
        reportFailure("the forward proposal density is not finite", num_bad_alpha_c);
        return RbConstants::Double::neginf;
        }
    return lnProb;
}

/* The log density that the forward move would have drawn the point pt. See the note
   on lnProbabilityForward for why this returns negative infinity rather than
   throwing.

   This is the direction that matters most. A failure here refuses the jump from the
   non-reversible model back to the time-reversible one; if it fails systematically,
   the chain can enter the non-reversible model and never leave, the tuning pushes
   the prior odds further and further toward reversibility trying to balance a chain
   that cannot be balanced, and the result looks exactly like evidence for
   non-reversibility of unbounded strength. */
double Polyhedron::lnProbabilityReverse(std::vector<mpq_class>& W, Vector& pt) {

    if (weightsAreUsable(W) == false)
        return RbConstants::Double::neginf;

    // set up the polyhedron, which will also locate the point
    randomlySample = false;
    pointFoundInPolyhedron = false;
    randomPoint = pt;                 // randomPoint now represents the point passed in to this function and
    setWeights(W);                    // must be set before setWeights() is called
    
    if (isValid(pt) == false)
        {
        num_point_not_valid++;
        reportFailure("the point representing the reverse move is not inside the polyhedron", num_point_not_valid);
        return RbConstants::Double::neginf;
        }
    if (pointFoundInPolyhedron == false)
        {
        num_point_not_located++;
        reportFailure("the point representing the reverse move could not be located in any tetrahedron, which happens when it lies on a facet", num_point_not_located);
        return RbConstants::Double::neginf;
        }
    if (sumJacobians <= 0)
        {
        num_degenerate_weights++;
        reportFailure("the polyhedron has no volume", num_degenerate_weights);
        return RbConstants::Double::neginf;
        }
    
    // calculate the probability of proposing the point, pt, passed in as a parameter
    double lnProb = -log(sumJacobians.get_d());
    lnProb += RbMath::lnGamma(alphaT + 3.0) - RbMath::lnGamma(alphaT);

    // see the note on the same term in lnProbabilityForward
    if (alphaT != 1.0)
        {
        if (alphaC <= 0.0)
            {
            num_bad_alpha_c++;
            reportFailure("the reverse point lies on a facet, so the concentration term of the proposal density is undefined", num_bad_alpha_c);
            return RbConstants::Double::neginf;
            }
        lnProb += (alphaT - 1.0) * log(alphaC);
        }

    if (std::isfinite(lnProb) == false)
        {
        num_point_not_located++;
        reportFailure("the reverse proposal density is not finite", num_point_not_located);
        return RbConstants::Double::neginf;
        }
    return lnProb;
}

void Polyhedron::sampleTetrahedron(Plane* pln, Vector* center, Vector* v1, Vector* v2, Vector* v3, Vector& pt, VectorInfo& info) {
    
    RandomNumberGenerator* rng = GLOBAL_RNG;
    
    std::vector<double> alpha(4, 1.0);
    alpha[0] = alphaT;
    std::vector<double> dirRv = RbStatistics::Dirichlet::rv(alpha, *rng);
    mpq_class a = dirRv[0];
    mpq_class s = dirRv[1];
    mpq_class t = dirRv[2];
    mpq_class u = 1 - a - s - t;
    info.alphaC = a.get_d();
    mpq_class x = center->getX() * a + v1->getX() * s + v2->getX() * t + v3->getX() * u;
    mpq_class y = center->getY() * a + v1->getY() * s + v2->getY() * t + v3->getY() * u;
    mpq_class z = center->getZ() * a + v1->getZ() * s + v2->getZ() * t + v3->getZ() * u;
    pt.set(x, y, z); // pt = vC*a + v1*s + v2*t + v3*u;
}

/* Print a failure the first time it happens and then on every power of ten.

   An MCMC can call these paths millions of times, so printing each one would bury
   the run in output; printing none would hide a systematic failure completely,
   which is the situation this whole mechanism exists to expose. */
void Polyhedron::reportFailure(const char* what, long count) {

    long p = 1;
    while (p < count)
        p *= 10;
    if (p != count)
        return;

    std::cerr << "Polyhedron: " << what;
    if (count == 1)
        std::cerr << " (first occurrence)";
    else
        std::cerr << " (occurrence " << count << ")";
    std::cerr << std::endl;
}

/* Are the current weights capable of defining a non-degenerate polyhedron?

   Every backbone weight has to be strictly positive. If one is zero, say w_CG, then
   both w_CG and w_GC are zero for every value of u1, the corresponding coordinate
   drops out of the map entirely, and the polyhedron is degenerate: it has no volume
   in that direction and the point that generated it cannot be located. The
   reversible jump is not defined in that case and the proposal has to be refused. */
bool Polyhedron::weightsAreUsable(std::vector<mpq_class>& W) {

    if (W.size() != 6)
        {
        num_degenerate_weights++;
        reportFailure("the weight vector does not have six elements", num_degenerate_weights);
        return false;
        }

    /* Test the weights that were passed in, not the members: the members are only
       assigned inside setWeights, which has not run yet when this is called. */
    if (W[0] <= 0 || W[1] <= 0 || W[2] <= 0 || W[3] <= 0 || W[4] <= 0 || W[5] <= 0)
        {
        num_degenerate_weights++;
        reportFailure("a backbone weight is zero, so the polyhedron is degenerate and the reversible-jump move cannot be made", num_degenerate_weights);
        return false;
        }
    return true;
}

void Polyhedron::setCenter(const Vector& v) {

    desiredCenter    = v;
    useDesiredCenter = true;
}

void Polyhedron::setCenter(double x, double y, double z) {

    Vector v(x, y, z);
    setCenter(v);
}

void Polyhedron::useReversibleCenter(void) {

    desiredCenter    = reversibleCenter;
    useDesiredCenter = false;
}

/* The twelve constraints that define the polyhedron, as values that must all be at
   or above zero. Six are the unit cube and six are the two-sided bounds on the three
   weight sums; see isValid, which tests the same quantities. Each is affine in
   (u1,u2,u3), which is what makes the centre clipping below exact. */
void Polyhedron::constraintValues(const Vector& pt, std::vector<mpq_class>& g) const {

    if (g.size() != 12)
        g.resize(12);

    const mpq_class& x = pt.getX();
    const mpq_class& y = pt.getY();
    const mpq_class& z = pt.getZ();

    mpq_class v1 = wAC + wCG * (2 * x - 1) + wCT * (2 * y - 1);
    mpq_class v2 = wAG - wCG * (2 * x - 1) + wGT * (2 * z - 1);
    mpq_class v3 = wAT - wCT * (2 * y - 1) - wGT * (2 * z - 1);

    g[0]  = x;
    g[1]  = 1 - x;
    g[2]  = y;
    g[3]  = 1 - y;
    g[4]  = z;
    g[5]  = 1 - z;
    g[6]  = v1;
    g[7]  = 2 * wAC - v1;
    g[8]  = v2;
    g[9]  = 2 * wAG - v2;
    g[10] = v3;
    g[11] = 2 * wAT - v3;
}

/* Decide where the centre of the triangulation goes for the polyhedron the current
   weights describe.

   The requested centre may lie outside it, because the shape depends on the weights
   and they move with every proposal. Rather than refuse, walk the point back along
   the straight line joining it to (1/2,1/2,1/2), which is always inside, and stop
   just short of where that line leaves the polyhedron.

   Writing C for the reversible centre and P for the requested one, the line is
   X(t) = C + t (P - C), and every constraint g is affine, so along the line

       g(t) = g(C) + t ( g(P) - g(C) )

   Each constraint that decreases along the line gives an upper bound on t of
   g(C) / (g(C) - g(P)), and the smallest of those bounds is where the line leaves.
   All of it is rational arithmetic, so the crossing is found exactly rather than by
   bisection, and multiplying by centerShrink then keeps the centre off the facet,
   where the barycentric weight would be zero and the proposal density undefined.

   This has to be a deterministic function of the weights, and it is: both directions
   of the reversible jump call setWeights with the same backbone, so both build the
   same triangulation about the same centre and their densities remain comparable. */
void Polyhedron::chooseCenter(void) {

    if (useDesiredCenter == false)
        {
        center = reversibleCenter;
        return;
        }

    std::vector<mpq_class> gC(12), gP(12);
    constraintValues(reversibleCenter, gC);
    constraintValues(desiredCenter, gP);

    mpq_class t = 1;
    bool clipped = false;
    for (int i=0; i<12; i++)
        {
        if (gC[i] <= 0)
            {
            /* The reversible centre is itself on the boundary, which means a backbone
               weight is zero and the polyhedron is degenerate. weightsAreUsable
               refuses the move before this is reached, so this is only a safety net. */
            center = reversibleCenter;
            return;
            }
        if (gP[i] < gC[i])
            {
            mpq_class bound = gC[i] / (gC[i] - gP[i]);
            if (bound < t)
                {
                t = bound;
                clipped = true;
                }
            }
        }

    if (clipped == true)
        {
        t *= centerShrink;
        numCenterClipped++;
        }

    mpq_class cx = reversibleCenter.getX() + t * (desiredCenter.getX() - reversibleCenter.getX());
    mpq_class cy = reversibleCenter.getY() + t * (desiredCenter.getY() - reversibleCenter.getY());
    mpq_class cz = reversibleCenter.getZ() + t * (desiredCenter.getZ() - reversibleCenter.getZ());
    center.set(cx, cy, cz);
}

void Polyhedron::setWeights(std::vector<mpq_class>& W) {

    // assign instance variables representing the
    // parameter of the polyhedron from W
    this->wAC = W[0];
    this->wAG = W[1];
    this->wAT = W[2];
    this->wCG = W[3];
    this->wCT = W[4];
    this->wGT = W[5];
    
    // the centre depends on the weights just assigned, so it is chosen here
    chooseCenter();
    
    // construct the polyhedron
    initializePlanes();
}
