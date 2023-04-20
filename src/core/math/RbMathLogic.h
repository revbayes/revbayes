/**
 * @file RbMathLogic
 * This file contains some logic function on numbers.
 *
 * @brief Implementation of simple RlBoolean algebra on numbers.
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date$
 * @author The RevBayes core development team
 * @license GPL version 3
 * @version 1.0
 * @since 2011-03-17, version 1.0
 *
 * $Id$
 */


#ifndef RbMathLogic_H
#define RbMathLogic_H

#include <vector>

namespace RevBayesCore {
    
    namespace RbMath {
    
        bool                        compApproximatelyEqual(double a, double b);                                         //!< Compares two doubles for equality
        bool                        compApproximatelyEqual(double a, double b, double epsilon);                         //!< Compares two doubles for equality
        bool                        compEssentiallyEqual(double a, double b);                                           //!< Compares two doubles for equality
        bool                        compEssentiallyEqual(double a, double b, double epsilon);                           //!< Compares two doubles for equality
        bool                        compDefinitelyGreaterThan(double a, double b);                                      //!< Tests whether one number is greater than another
        bool                        compDefinitelyGreaterThan(double a, double b, double epsilon);                      //!< Tests whether one number is greater than another
        bool                        compDefinitelyLessThan(double a, double b);                                         //!< Tests whether one number is less than another
        bool                        compDefinitelyLessThan(double a, double b, double epsilon);                         //!< Tests whether one number is less than another
        bool                        isAComputableNumber(double x);
        bool                        isFinite(double x);                                                                 //!< Tests whether a double is finite
        //bool                        isFinite(size_t x); // size_t is always finite (WP)
        bool                        isInt(double x);                                                                    //!< Tests whether a double is actually an interger
        bool                        isNan(double x);                                                                    //!< Tests whether a double is NAN (not a number)
        double                      max(double a, double b);                                                            //!< Returns the max of two real numbers
        double                      min(double a, double b);                                                            //!< Returns the min of two real numbers
        int                         min(int a, int b);                                                                  //!< Returns the min of two real numbers
        int                         sgn(double x);                                                                      //!< Returns the sign of a number
    }
    
}

#endif
