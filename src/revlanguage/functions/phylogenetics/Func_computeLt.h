#ifndef Func_computeLt_H
#define Func_computeLt_H

// any other includes? tt
#include "RealPos.h"
#include "RlTypedFunction.h"

#include <string>

namespace RevLanguage {

  /*
  * Add description
  */

  class Func_computeLt : public TypedFunction<RealPos> {

  public:
    Func_computeLt ( void );

    // Basic utility functions
    Func_computeLt*                 clone(void) const; // Clone the object
    static const std::string&       getClassType(void); // Get the Rev type
    static const TypeSpec&          getClassTypeSpec(void); // Get the class type
    std::string                     getFunctionName(void) const; // Get the primary name
    const TypeSpec&                 getTypeSpec(void) const; // Get the type spec

    // Function functions you have to override
    RevBayesCore::TypedFunction<double>*      createFunction(void) const; // Create internal function object
    const ArgumentRules&                      getArgumentRules(void) const; // Get argument rules
    //const MemberRules&                          getParameterRules(void) const; //Get member rules

  protected:

      void                                     setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);

  private:

      RevPtr<const RevVariable>                               origin;
      RevPtr<const RevVariable>                               lambda;
      RevPtr<const RevVariable>                               mu;
      RevPtr<const RevVariable>                               psi;
      RevPtr<const RevVariable>                               omega;
      RevPtr<const RevVariable>                               removal;
      RevPtr<const RevVariable>                               rho;
      RevPtr<const RevVariable>                               listA;
      RevPtr<const RevVariable>                               listB;
      RevPtr<const RevVariable>                               listC;
      RevPtr<const RevVariable>                               listD;
      RevPtr<const RevVariable>                               listE;
      RevPtr<const RevVariable>                               listF;
      RevPtr<const RevVariable>                               listG;


  };

}

#endif
