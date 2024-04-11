#include <iomanip>
#include <iostream>
#include <map>

#include "DistributionDirichlet.h"
#include "RateMatrix_MPQ.h"
#include "Polyhedron.h"
#include "RandomNumberGenerator.h"
#include "RandomNumberFactory.h"
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
    
    // initialize commonly used constants
    zeroQ     = 0;
    oneQ      = 1;
    oneHalfQ  = 1;
    oneHalfQ /= 2;
    twoQ      = 2;
    
    randomlySample = false;
    alphaT = 1.0; // the default value for the Dirichlet(alphaT,1,1,1) used to randomly draw points from tetrahedra

    // the center vertex is constant, so initialize it here
    center.setX(oneHalfQ);
    center.setY(oneHalfQ);
    center.setZ(oneHalfQ);

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

    aQ = v1->getX() - oneHalfQ;
    bQ = v2->getX() - oneHalfQ;
    cQ = v3->getX() - oneHalfQ;
    dQ = v1->getY() - oneHalfQ;
    eQ = v2->getY() - oneHalfQ;
    fQ = v3->getY() - oneHalfQ;
    gQ = v1->getZ() - oneHalfQ;
    hQ = v2->getZ() - oneHalfQ;
    iQ = v3->getZ() - oneHalfQ;
    
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
    VertexFactory& vf = VertexFactory::vertexFactoryInstance();
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

double Polyhedron::lnProbabilityForward(std::vector<mpq_class>& W, Vector& pt) {

    // set up the polyhedron, which will also randomly sample and initialize pt
    randomlySample = true;
    setWeights(W);
    
    // initialize pt to the randomly selected point and do a sanity check
    pt.set(randomPoint.getX(), randomPoint.getY(), randomPoint.getZ());
    if (isValid(pt) == false)
        throw(RbException("Polyhedron: Random point is not in polyhedron"));
    randomlySample = false;
    
    // calculate the probability of the randomly proposed point, pt
    double lnProb = -log(sumJacobians.get_d());
    lnProb += RbMath::lnGamma(alphaT + 3.0) - RbMath::lnGamma(alphaT);
    lnProb += (alphaT - 1.0) * log(alphaC);
    return lnProb;
}

double Polyhedron::lnProbabilityReverse(std::vector<mpq_class>& W, Vector& pt) {

    // set up the polyhedron, which will also locate the point
    randomlySample = false;
    pointFoundInPolyhedron = false;
    randomPoint = pt;                 // randomPoint now represents the point passed in to this function and
    setWeights(W);                    // must be set before setWeights() is called
    
    // some sanity checks
    if (isValid(pt) == false)
        throw(RbException("Polyhedron: Point representing reverse move not in polyhedron"));
    if (pointFoundInPolyhedron == false)
        throw(RbException("Polyhedron: Did not find point in polyhedron"));
    
    // calculate the probability of proposing the point, pt, passed in as a parameter
    double lnProb = -log(sumJacobians.get_d());
    lnProb += RbMath::lnGamma(alphaT + 3.0) - RbMath::lnGamma(alphaT);
    lnProb += (alphaT - 1.0) * log(alphaC);
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

void Polyhedron::setWeights(std::vector<mpq_class>& W) {

    // assign instance variables representing the
    // parameter of the polyhedron from W
    this->wAC = W[0];
    this->wAG = W[1];
    this->wAT = W[2];
    this->wCG = W[3];
    this->wCT = W[4];
    this->wGT = W[5];
    
    // construct the polyhedron
    initializePlanes();
}
