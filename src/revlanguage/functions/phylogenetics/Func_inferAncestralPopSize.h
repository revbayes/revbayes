#ifndef Func_inferAncestralPopSize_H
#define Func_inferAncestralPopSize_H

// any other includes? tt
#include "RlMatrixReal.h"
#include "MatrixReal.h"
#include "RlTypedFunction.h"

#include <string>

namespace RevLanguage {

  /*
  * Add description
  */

  class Func_inferAncestralPopSize : public TypedFunction<MatrixReal> {

  public:
    Func_inferAncestralPopSize ( void );

    // Basic utility functions
    Func_inferAncestralPopSize*                                       clone(void) const; // Clone the object
    static const std::string&                                         getClassType(void); // Get the Rev type
    static const TypeSpec&                                            getClassTypeSpec(void); // Get the class type
    std::string                                                       getFunctionName(void) const; // Get the primary name
    const TypeSpec&                                                   getTypeSpec(void) const; // Get the type spec

    // Function functions you have to override
    RevBayesCore::TypedFunction< RevBayesCore::MatrixReal >*          createFunction(void) const; // Create internal function object
    const ArgumentRules&                                              getArgumentRules(void) const; // Get argument rules
    //const MemberRules&                                                 getParameterRules(void) const; //Get member rules

  protected:

      void                                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);

  private:

      RevPtr<const RevVariable>                               start_age;                                                                              //!< The time of origin
      RevPtr<const RevVariable>                               lambda;                                                                                 //!< The speciation rate
      RevPtr<const RevVariable>                               mu;                                                                                     //!< The extinction rate
      RevPtr<const RevVariable>                               psi;                                                                                    //!< The serial sampling rate
      RevPtr<const RevVariable>                               omega;                                                                                  //!< The occurrence sampling rate
      RevPtr<const RevVariable>                               rho;                                                                                    //!< The taxon sampling fraction
      RevPtr<const RevVariable>                               removalPr;                                                                              //!< The removal probability after sampling
      RevPtr<const RevVariable>                               maxHiddenLin;                                                                           //!< The number of hidden lineages (algorithm accuracy)
      std::string                                             start_condition;                                                                        //!< The start condition of the process (rootAge/originAge)
      RevPtr<const RevVariable>                               condition;                                                                              //!< The conditioning of the process ("time" or "survaval")
      RevPtr<const RevVariable>                               occurrence_ages;                                                                        //!< The occurrence ages of incomplete fossils
      RevPtr<const RevVariable>                               extant;                                                                                 //!< The number of extant taxa
      RevPtr<const RevVariable>                               time_points;                                                                            //!< The times at which density is computed
      RevPtr<const RevVariable>                               timeTree;                                                                               //!< The tree for ancestral pop. size inference
      RevPtr<const RevVariable>                               timeline;                                                                               //!< The rate shifts timeline


  };

}

#endif
