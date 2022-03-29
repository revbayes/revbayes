/**
 * @file
 * This file contains the declaration of Matrix,
 * a container type used to hold value matrices for the inference machinery.
 *
 * @brief Declaration of Matrix
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date: 2012-03-10 12:55:11 +0100 (Sat, 10 Mar 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2012-05-08, version 1.0
 *
 * $Id: Matrix.h 1327 2012-03-10 11:55:11Z hoehna $
 */

#ifndef Simplex_H
#define Simplex_H

#include <cstddef>
#include <iostream>
#include <vector>

#include "MemberObject.h"
#include "RbVector.h"

namespace RevBayesCore {
class DagNode;
    
    class Simplex : public RbVector<double> {
        
    public:
        Simplex(void);                       //!< Default constructor required by revlanguage use of this class
        Simplex(size_t n);
        Simplex(size_t n, double v);
        Simplex(const std::vector<double> &v);
        Simplex(const Simplex& m);
        Simplex(Simplex&& m) = default;
        virtual                                ~Simplex(void);

        Simplex& operator=(const Simplex&) = default;
        Simplex& operator=(Simplex&&) = default;
        
        // global operators
        RbVector<double>                        operator+(double b) const;                                          //!< operator + for matrix + scalar
        RbVector<double>                        operator-(double b) const;                                          //!< operator - for scalar
        RbVector<double>                        operator*(double b) const;                                          //!< operator * for scalar
//        std::vector<double>                     operator*(const std::vector<double> &b) const;
        
        Simplex*                                clone(void) const;
        
        // utility funcions
        virtual void                            initFromString(const std::string &s);
        
    protected:
        
        void                                    normalize(void);
        
        
    };
    
    // Global functions using the class
    std::ostream&                       operator<<(std::ostream& o, const Simplex& x);                                           //!< Overloaded output operator
    
    RbVector<double>                        operator*(const double &a, const Simplex& B);                            //!< operator * for scalar * simplex

    
}

#endif

