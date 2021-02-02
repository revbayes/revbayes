#ifndef RlAverageDistanceMatrix_H
#define RlAverageDistanceMatrix_H

#include "ModelObject.h"
#include "AverageDistanceMatrix.h"
#include "TypedDagNode.h"

#include <ostream>
#include <string>

namespace RevLanguage {
    
   /**
    * The RevLanguage wrapper of the average distance matrix.
    *
    * The RevLanguage wrapper of the average distance matrix simply
    * manages the interactions through the Rev with our core.
    * That is, the internal move object can be constructed and hooked up
    * in a DAG-nove (variable) that it works on.
    * See AverageDistanceMatrix.h for more details.
    *
    *
    * @author David Cerny
    * @since 2019-10-24, version 1.0
    *
    */
    class AverageDistanceMatrix : public ModelObject<RevBayesCore::AverageDistanceMatrix>  {
            
        public:
            
                                                AverageDistanceMatrix(void);                                                            //!< Default constructor
                                                AverageDistanceMatrix(RevBayesCore::AverageDistanceMatrix *a);                                 //!< Default constructor
                                                AverageDistanceMatrix(const RevBayesCore::AverageDistanceMatrix &v);                           //!< Default constructor
                                                AverageDistanceMatrix(RevBayesCore::TypedDagNode<RevBayesCore::AverageDistanceMatrix> * mat);     //!< Default constructor

            // the value type definition
            AverageDistanceMatrix*              clone(void) const;                                                                      //!< Clone object
            static const std::string&           getClassType(void);                                                                     //!< Get Rev type
            static const TypeSpec&              getClassTypeSpec(void);                                                                 //!< Get class type spec
            const TypeSpec&                     getTypeSpec(void) const;                                                                //!< Get language type of the object

            // A STL-like vector function
            virtual size_t                      size(void) const;                                                                       //!< Get number of elements in container

            // Member method inits
            RevPtr<RevVariable>                 executeMethod(const std::string& name, const std::vector<Argument>& args, bool &f);     //!< Override to map member methods to internal functions

            std::string                         getGuiName(void) { return "Average Distance Matrix"; }
            std::string                         getGuiUnicodeSymbol(void) { return ""; }
            std::string                         getGuiInfo(void) { return ""; }
            
        private:
            void                                initializeMethods(void);
    };

}

#endif /* defined(__revbayes__AverageDistanceMatrix__) */
