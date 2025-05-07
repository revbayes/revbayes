/**
 * @file
 * This file contains the Workspace function that adds types and functions
 * to the global workspace, registering them with the interpreter/compiler
 * during the process.
 *
 * @brief Function registering language objects
 *
 * Instructions
 *
 * This is the central registry of Rev objects. It is a large file and needs
 * to be properly organized to facilitate maintenance. Follow these simple
 * guidelines to ensure that your additions follow the existing structure.
 *
 * 1. All headers are added in groups corresponding to directories in the
 *    revlanguage code base.
 * 2. All objects (types, distributions, and functions) are registered in
 *    groups corresponding to directories in the revlanguage code base.
 * 3. All entries in each group are listed in alphabetical order.
 *
 * Some explanation of the directory structure is provided in the comments
 * in this file. Consult these comments if you are uncertain about where
 * to add your objects in the code.
 */

#include <sstream>
#include <vector>
#include <cstdlib>
#include <cstdio>
#include <string>

/* Files including helper classes */
#include "RbException.h"
#include "RlUserInterface.h"
#include "Workspace.h"

/// Miscellaneous types ///

#include "BinaryAddition.h"                               // for BinaryAddition
#include "BinaryDivision.h"                               // for BinaryDivision
#include "BinaryMultiplication.h"                         // for BinaryMulti...
#include "BinarySubtraction.h"                            // for BinarySubtr...
#include "ConstantNode.h"                                 // for ConstantNode
#include "DagNode.h"                                      // for DagNode
#include "DeterministicNode.h"                            // for Determinist...
#include "DynamicNode.h"                                  // for DynamicNode
#include "EquationFunction.h"                             // for EquationFun...
#include "Func__conversion.h"                             // for Func__conve...
#include "GreaterEqualFunction.h"                         // for GreaterEqua...
#include "GreaterThanFunction.h"                          // for GreaterThan...
#include "IfElseFunction.h"                               // for IfElseFunction
#include "IndirectReferenceFunction.h"                    // for IndirectRef...
#include "LessEqualFunction.h"                            // for LessEqualFu...
#include "LessThanFunction.h"                             // for LessThanFun...
#include "MatrixReal.h"                                   // for MatrixReal
#include "ModelObject.h"                                  // for ModelObject
#include "NotEqualFunction.h"                             // for NotEqualFun...
#include "RbBoolean.h"                                    // for Boolean
#include "RbUtil.h"                                       // for operator*
#include "RbVector.h"                                     // for RbVector
#include "RbVectorImpl.h"                                 // for RbVectorImpl
#include "ReplicateFunction.h"                            // for ReplicateFu...
#include "ReplicateEventsFunction.h"                      // for ReplicateFu...
#include "RevPtr.h"                                       // for RevPtr
#include "RlBranchLengthTree.h"                           // for BranchLengt...
#include "RlConstantNode.h"                               // for ConstantNode
#include "RlDeterministicNode.h"                          // for Determinist...
#include "RlRateGenerator.h"                              // for RateGenerator
#include "RlSiteMixtureModel.h"                               // for MixtureModel
#include "RlTimeTree.h"                                   // for TimeTree
#include "RlTypedFunction.h"                              // for TypedFunction
#include "ScalarMatrixMultiplication.h"                   // for ScalarMatri...
#include "ScalarVectorAddition.h"                         // for ScalarVecto...
#include "ScalarVectorDivision.h"                         // for ScalarVecto...
#include "ScalarVectorMultiplication.h"                   // for ScalarVecto...
#include "ScalarVectorSubtraction.h"                      // for ScalarVecto...
#include "TypeConversionFunction.h"                       // for TypeConvers...
#include "TypedDagNode.h"                                 // for TypedDagNode
#include "TypedFunction.h"                                // for TypedFunction
#include "UnaryMinus.h"                                   // for UnaryMinus
#include "UserFunctionNode.h"                             // for UserFunctio...
#include "VectorAppendElement.h"                          // for VectorAppen...
#include "VectorAppendVector.h"                           // for VectorAppen...
#include "VectorReverse.h"                                // for VectorRever...
#include "VectorScalarDivision.h"                         // for VectorScala...
#include "VectorScalarSubtraction.h"                      // for VectorScala...

/* Base types (in folder "datatypes") */

/* Primitive types (in folder "datatypes/basic") */
#include "Integer.h"
#include "IntegerPos.h"
#include "Natural.h"
#include "Probability.h"
#include "RlBoolean.h"
#include "RlString.h"
#include "Real.h"
#include "RealPos.h"

/* Container types (in folder "datatypes/container") */
#include "ModelVector.h"

/* Evolution types (in folder "datatypes/evolution") */

/* Character state types (in folder "datatypes/evolution/character") */
#include "RlDiscreteCharacterState.h"

/* Character data types (in folder "datatypes/evolution/datamatrix") */

/* Tree types (in folder "datatypes/evolution/trees") */

/* Taxon types (in folder "datatypes/evolution") */

/* Math types (in folder "datatypes/math") */
#include "RlMatrixReal.h"
#include "RlMatrixRealSymmetric.h"
#include "RlSimplex.h"


/// Distributions ///

/* Distribution types (in folder "distributions") */

/* Character evolution models (in folder "distributions/evolution/character") */

/* Argument rules (in folder "functions/argumentrules") */

/* Basic functions (in folder "functions/basic"). */

/* These are core functions for the Rev environment, providing user help
   and other essential services. */

#include "Func_append.h"
#include "Func_appendVector.h"
#include "Func_clear.h"
#include "Func_exists.h"
#include "Func_getOption.h"
#include "Func_getwd.h"
#include "Func_help.h"
#include "Func_ifelse.h"
#include "Func_license.h"
#include "Func_listOptions.h"
#include "Func_ls.h"
#include "Func_printSeed.h"
#include "Func_profile.h"
#include "Func_quit.h"
#include "Func_range.h"
#include "Func_replicate.h"
#include "Func_reverse.h"
#include "Func_replicateEvents.h"
#include "Func_seed.h"
#include "Func_seq.h"
#include "Func_setOption.h"
#include "Func_setwd.h"
#include "Func_structure.h"
#include "Func_system.h"
#include "Func_time.h"
#include "Func_type.h"

#if defined( RB_BEAGLE )
#include "Func_beagleInfo.h"
#endif /* RB_BEAGLE */

/* Internal functions (in folder ("functions/internal") */

/* These are functions that are typically not called explicitly but implicitly
   through parsing of a Rev statement. Examples include a statement like '1 + 2',
   which results in the builtin '_add' function being called.

   Exceptions include Func_range and Func_vector, which are both used for implicit
   and explicit calls. They are therefore considered basic functions instead of
   internal functions.

   All internal functions have function calls that start with an underscore character,
   and therefore their class names have two underscore characters. They are typically
   templated. */

#include "Func__add.h"
#include "Func__and.h"
#include "Func__div.h"
#include "Func__eq.h"
#include "Func__ge.h"
#include "Func__gt.h"
#include "Func__le.h"
#include "Func__lt.h"
#include "Func__mult.h"
#include "Func__mod.h"
#include "Func__ne.h"
#include "Func__or.h"
#include "Func__unot.h"
#include "Func__scalarMatrixMult.h"
#include "Func__scalarVectorAdd.h"
#include "Func__scalarVectorDiv.h"
#include "Func__scalarVectorMult.h"
#include "Func__scalarVectorSub.h"
#include "Func__sub.h"
#include "Func__uminus.h"
#include "Func__vectorIndexOperator.h"
#include "Func__simplexIndexOperator.h"
#include "Func__vectorScalarDiv.h"
#include "Func__vectorScalarSub.h"


/* Input/output functions (in folder "functions/io") */
#include "Func_ancestralStateTree.h"
#include "Func_annotateTree.h"
#include "Func_characterMapTree.h"
#include "Func_consensusTree.h"
#include "Func_convertToPhylowood.h"
#include "Func_fileExists.h"
#include "Func_listFiles.h"
#include "Func_loadPlugin.h"
#include "Func_readRegionalFeatures.h"
#include "Func_maxdiff.h"
#include "Func_mapTree.h"
#include "Func_mccTree.h"
#include "Func_module.h"
#include "Func_readAncestralStateTreeTrace.h"
#include "Func_readAncestralStateTrace.h"
#include "Func_readAtlas.h"
#include "Func_readBranchLengthTrees.h"
#include "Func_readDelimitedCharacterData.h"
#include "Func_readCharacterData.h"
#include "Func_readContinuousCharacterData.h"
#include "Func_readDelimitedDataFile.h"
#include "Func_readDiscreteCharacterData.h"
#include "Func_readDistanceMatrix.h"
#include "Func_readMatrix.h"
#include "Func_readRelativeNodeAgeConstraints.h"
#include "Func_readRelativeNodeAgeWeightedConstraints.h"
#include "Func_readStochasticVariableTrace.h"
#include "Func_readTrace.h"
#include "Func_readTrees.h"
#include "Func_readTreeTrace.h"
#include "Func_readVCF.h"
#include "Func_source.h"
#include "Func_summarizeCharacterMaps.h"
#include "Func_TaxonReader.h"
#include "Func_treeTrace.h"
#include "Func_write.h"
#include "Func_writeDelimitedCharacterData.h"
#include "Func_writeFasta.h"
#include "Func_writeNexus.h"


/* Math functions (in folder "functions/math") */
//#include "Func_abs.h"
//#include "Func_ceil.h"
//#include "Func_diagonalMatrix.h"
//#include "Func_exp.h"
//#include "Func_floor.h"
//#include "Func_lnProbability.h"
//#include "Func_hyperbolicTangent.h"
//#include "Func_ln.h"
//#include "Func_log.h"
//#include "Func_max.h"
//#include "Func_mean.h"
//#include "Func_min.h"
//#include "Func_normalize.h"
#include "Func_power.h"
#include "Func_powerVector.h"
//#include "Func_round.h"
//#include "Func_simplex.h"
//#include "Func_simplexFromVector.h"
//#include "Func_sum.h"
//#include "Func_sumPositive.h"
//#include "Func_standardDeviation.h"
//#include "Func_sqrt.h"
//#include "Func_trunc.h"
//#include "Func_variance.h"


/* Statistics functions (in folder "functions/statistics") */
/* These are functions related to statistical distributions */
//#include "Func_discretizeBeta.h"
//#include "Func_discretizeBetaQuadrature.h"
//#include "Func_discretizeGamma.h"
//#include "Func_discretizeGammaQuadrature.h"
//#include "Func_discretizeLognormalQuadrature.h"
//#include "Func_discretizeDistribution.h"
//#include "Func_discretizePositiveDistribution.h"
//#include "Func_dppConcFromMean.h"
//#include "Func_dppMeanFromConc.h"
//#include "Func_fnNormalizedQuantile.h"
//#include "Func_numUniqueInVector.h"
//#include "Func_stirling.h"
//#include "Func_varianceCovarianceMatrix.h"
//#include "Func_decomposedVarianceCovarianceMatrix.h"


/** Initialize global workspace */
void RevLanguage::Workspace::initializeBasicGlobalWorkspace(void)
{

    try
    {
        ///////////////////////////////////////////
        /* Add functions (in "functions" folder) */
        ///////////////////////////////////////////

        /* Basic functions (in folder "functions/basic") */

        // regular functions
        addFunction( new Func_append<Integer>()             );
        addFunction( new Func_append<Real>()                );
        addFunction( new Func_append<Natural>()             );
        addFunction( new Func_append<RealPos>()             );
        addFunction( new Func_append<RlString>()            );
        addFunction( new Func_append<RlBoolean>()           );
        addFunction( new Func_appendVector<Integer>()       );
        addFunction( new Func_appendVector<Real>()          );
        addFunction( new Func_appendVector<Natural>()       );
        addFunction( new Func_appendVector<RealPos>()       );
        addFunction( new Func_appendVector<RlString>()      );
        addFunction( new Func_appendVector<RlBoolean>()     );
        addFunction( new Func_reverse<Integer>()            );
        addFunction( new Func_reverse<Real>()               );
        addFunction( new Func_reverse<Natural>()            );
        addFunction( new Func_reverse<RealPos>()            );
        addFunction( new Func_reverse<RlString>()           );
        addFunction( new Func_reverse<RlBoolean>()          );
        addFunction( new Func_clear()                       );
        addFunction( new Func_exists()                      );
        addFunction( new Func_getwd()                       );
        addFunction( new Func_getOption()                   );
        addFunction( new Func_help()                        );
        addFunction( new Func_ifelse<Natural>()             );
        addFunction( new Func_ifelse<Integer>()             );
        addFunction( new Func_ifelse<Real>()                );
        addFunction( new Func_ifelse<RealPos>()             );
        addFunction( new Func_ifelse<RlString>()            );
        addFunction( new Func_license()                     );
        addFunction( new Func_listOptions()                 );
        addFunction( new Func_ls()                          );
        addFunction( new Func_printSeed()                   );
        addFunction( new Func_profile()                     );
        addFunction( new Func_quit()                        );
        addFunction( new Func_replicate<Integer>()          );
        addFunction( new Func_replicate<Real>()             );
        addFunction( new Func_replicate<Natural>()          );
        addFunction( new Func_replicate<RealPos>()          );
        addFunction( new Func_replicate<RlString>()         );
        addFunction( new Func_replicate<RlBoolean>()        );
        addFunction( new Func_replicate<SiteMixtureModel>() );
        addFunction( new Func_seed()                        );
        addFunction( new Func_replicateEvents<Real>()       );
        addFunction( new Func_replicateEvents<RealPos>()    );
        addFunction( new Func_replicateEvents<Probability>());
        addFunction( new Func_seq<Integer>()                );
        addFunction( new Func_seq<Real>()                   );
        addFunction( new Func_setOption()                   );
        addFunction( new Func_setwd()                       );
        addFunction( new Func_structure()                   );
        addFunction( new Func_system()                      );
        addFunction( new Func_time()                        );
        addFunction( new Func_type()                        );

#if defined( RB_BEAGLE )
        addFunction( new Func_beagleInfo()                  );
#endif /* RB_BEAGLE */


        /* Internal functions (in folder "functions/internal") */
        /* Note: These are functions that are called implicitly, and the name of which, if
         called explicitly, starts with an underscore character. */

        // not templated logical functions
        addFunction( new Func__and()   );
        addFunction( new Func__or()    );
        addFunction( new Func__unot()  );

        // range function (x:y)
        addFunction( new Func_range()  );

        // logical templated functions
        addFunction( new Func__eq<             Integer,          Integer >()             );
        addFunction( new Func__eq<                Real,             Real >()             );
        addFunction( new Func__eq<             Integer,             Real >()             );
        addFunction( new Func__eq<                Real,          Integer >()             );
        addFunction( new Func__eq<           RlBoolean,        RlBoolean >()             );
        addFunction( new Func__eq<            RlString,         RlString >()             );
        addFunction( new Func__eq<             Simplex,          Simplex >()             );
        addFunction( new Func__eq<       RateGenerator,    RateGenerator >()             );
        addFunction( new Func__eq<            TimeTree,         TimeTree >()             );
        addFunction( new Func__eq<    BranchLengthTree, BranchLengthTree >()             );
        addFunction( new Func__eq<DiscreteCharacterState, DiscreteCharacterState >()     );
        addFunction( new Func__ge<             Integer,          Integer >()             );
        addFunction( new Func__ge<                Real,             Real >()             );
        addFunction( new Func__ge<             Integer,             Real >()             );
        addFunction( new Func__ge<                Real,          Integer >()             );
        addFunction( new Func__ge<           RlBoolean,        RlBoolean >()             );
        addFunction( new Func__gt<             Integer,          Integer >()             );
        addFunction( new Func__gt<                Real,             Real >()             );
        addFunction( new Func__gt<           RlBoolean,        RlBoolean >()             );
        addFunction( new Func__le<             Integer,          Integer >()             );
        addFunction( new Func__le<                Real,             Real >()             );
        addFunction( new Func__le<             Integer,             Real >()             );
        addFunction( new Func__le<                Real,          Integer >()             );
        addFunction( new Func__le<           RlBoolean,        RlBoolean >()             );
        addFunction( new Func__lt<             Integer,          Integer >()             );
        addFunction( new Func__lt<                Real,             Real >()             );
        addFunction( new Func__lt<           RlBoolean,        RlBoolean >()             );
        addFunction( new Func__ne<             Integer,          Integer >()             );
        addFunction( new Func__ne<                Real,             Real >()             );
        addFunction( new Func__ne<             Integer,             Real >()             );
        addFunction( new Func__ne<                Real,          Integer >()             );
        addFunction( new Func__ne<           RlBoolean,        RlBoolean >()             );
        addFunction( new Func__ne<            RlString,         RlString >()             );
        addFunction( new Func__ne<       RateGenerator,    RateGenerator >()             );
        addFunction( new Func__ne<             Simplex,          Simplex >()             );
        addFunction( new Func__ne<            TimeTree,         TimeTree >()             );
        addFunction( new Func__ne<    BranchLengthTree, BranchLengthTree >()             );

        // unary minus (e.g. -a)
        addFunction( new Func__uminus<Integer,             Integer>()  );
        addFunction( new Func__uminus<Natural,             Integer>()  );
        addFunction( new Func__uminus<Real,                Real>()        );
        addFunction( new Func__uminus<RealPos,             Real>()     );
        addFunction( new Func__uminus<MatrixReal,          MatrixReal>()     );
        addFunction( new Func__uminus<MatrixRealSymmetric, MatrixRealSymmetric>()     );

        // addition (e.g. a+b )
        addFunction( new Func__add< Natural                , Natural                          , Natural                            >(  )   );
        addFunction( new Func__add< Integer                , Integer                          , Integer                            >(  )   );
        addFunction( new Func__add< Real                   , Real                             , Real                               >(  )   );
        addFunction( new Func__add< RealPos                , RealPos                          , RealPos                            >(  )   );
        addFunction( new Func__add< RlString               , RlString                         , RlString                           >(  )   );
        addFunction( new Func__add< RlString               , Real                             , RlString                           >(  )   );
        addFunction( new Func__add< RlString               , Integer                          , RlString                           >(  )   );
        addFunction( new Func__add< RlString               , ModelVector<Natural>             , RlString                           >(  )   );
        addFunction( new Func__add< RlString               , ModelVector<Integer>             , RlString                           >(  )   );
        addFunction( new Func__add< RlString               , ModelVector<RealPos>             , RlString                           >(  )   );
        addFunction( new Func__add< RlString               , ModelVector<Real>                , RlString                           >(  )   );
        addFunction( new Func__add< MatrixReal             , MatrixReal                       , MatrixReal                         >(  )   );
        addFunction( new Func__add< MatrixRealSymmetric    , MatrixRealSymmetric              , MatrixRealSymmetric                >(  )   );
        addFunction( new Func__add< ModelVector<Natural>   , ModelVector<Natural>             , ModelVector<Natural>               >(  )   );
        addFunction( new Func__add< ModelVector<Integer>   , ModelVector<Integer>             , ModelVector<Integer>               >(  )   );
        addFunction( new Func__add< ModelVector<RealPos>   , ModelVector<RealPos>             , ModelVector<RealPos>               >(  )   );
        addFunction( new Func__add< ModelVector<Real>      , ModelVector<Real>                , ModelVector<Real>                  >(  )   );
        addFunction( new Func__scalarVectorAdd<Natural     , ModelVector<Natural>             , ModelVector<Natural>               >(  )   );
        addFunction( new Func__scalarVectorAdd<Integer     , ModelVector<Integer>             , ModelVector<Integer>               >(  )   );
        addFunction( new Func__scalarVectorAdd<Real        , ModelVector<Real>                , ModelVector<Real>                  >(  )   );
        addFunction( new Func__scalarVectorAdd<RealPos     , ModelVector<RealPos>             , ModelVector<RealPos>               >(  )   );

        // division
        addFunction( new Func__div< RealPos                            , Natural               , RealPos                   >(  )  );
        addFunction( new Func__div< Real                               , Integer               , Real                      >(  )  );
        addFunction( new Func__div< Real                               , Real                  , Real                      >(  )  );
        addFunction( new Func__div< RealPos                            , RealPos               , RealPos                   >(  )  );
        addFunction( new Func__div< ModelVector<RealPos>               , ModelVector<Natural>  , ModelVector<RealPos>      >(  )  );
        addFunction( new Func__div< ModelVector<Real>                  , ModelVector<Integer>  , ModelVector<Real>         >(  )  );
        addFunction( new Func__div< ModelVector<RealPos>               , ModelVector<RealPos>  , ModelVector<RealPos>      >(  )  );
        addFunction( new Func__div< ModelVector<Real>                  , ModelVector<Real>     , ModelVector<RealPos>      >(  )  );
        addFunction( new Func__vectorScalarDiv<Natural                 , RealPos               , RealPos                   >(  )   );
        addFunction( new Func__vectorScalarDiv<RealPos                 , Natural               , RealPos                   >(  )   );
        addFunction( new Func__vectorScalarDiv<Integer                 , Real                  , Real                      >(  )   );
        addFunction( new Func__vectorScalarDiv<Real                    , Integer               , Real                      >(  )   );
        addFunction( new Func__vectorScalarDiv<Real                    , Real                  , Real                      >(  )   );
        addFunction( new Func__vectorScalarDiv<RealPos                 , RealPos               , RealPos                   >(  )   );
        addFunction( new Func__scalarVectorDiv<RealPos                 , Natural               , RealPos                   >(  )   );
        addFunction( new Func__scalarVectorDiv<Natural                 , RealPos               , RealPos                   >(  )   );
        addFunction( new Func__scalarVectorDiv<Real                    , Integer               , Real                      >(  )   );
        addFunction( new Func__scalarVectorDiv<Integer                 , Real                  , Real                      >(  )   );
        addFunction( new Func__scalarVectorDiv<Real                    , Real                  , Real                      >(  )   );
        addFunction( new Func__scalarVectorDiv<RealPos                 , RealPos               , RealPos                   >(  )   );

        // multiplication
        addFunction( new Func__mult< Natural               , Natural               , Natural               >(  )  );
        addFunction( new Func__mult< Integer               , Integer               , Integer               >(  )  );
        addFunction( new Func__mult< Real                  , Real                  , Real                  >(  )  );
        addFunction( new Func__mult< RealPos               , RealPos               , RealPos               >(  )  );
        addFunction( new Func__mult< Simplex               , RealPos               , ModelVector<RealPos>  >(  )  );
        addFunction( new Func__mult< Simplex               , Real                  , ModelVector<Real>     >(  )  );
        addFunction( new Func__mult< ModelVector<Natural>  , ModelVector<Natural>  , ModelVector<Natural>  >(  )  );
        addFunction( new Func__mult< ModelVector<Integer>  , ModelVector<Integer>  , ModelVector<Integer>  >(  )  );
        addFunction( new Func__mult< ModelVector<RealPos>  , ModelVector<RealPos>  , ModelVector<RealPos>  >(  )  );
        addFunction( new Func__mult< ModelVector<Real>     , ModelVector<Real>     , ModelVector<Real>     >(  )  );
        addFunction( new Func__scalarVectorMult<Natural    , ModelVector<Natural>  , ModelVector<Natural>  >(  )   );
        addFunction( new Func__scalarVectorMult<Integer    , ModelVector<Integer>  , ModelVector<Integer>  >(  )   );
        addFunction( new Func__scalarVectorMult<Real       , ModelVector<Real>     , ModelVector<Real>     >(  )   );
        addFunction( new Func__scalarVectorMult<RealPos    , ModelVector<RealPos>  , ModelVector<RealPos>  >(  )   );
        addFunction( new Func__scalarMatrixMult<Natural    , MatrixReal  , MatrixReal  >(  )   );
        addFunction( new Func__scalarMatrixMult<Integer    , MatrixReal  , MatrixReal  >(  )   );
        addFunction( new Func__scalarMatrixMult<Real       , MatrixReal  , MatrixReal  >(  )   );
        addFunction( new Func__scalarMatrixMult<RealPos    , MatrixReal  , MatrixReal  >(  )   );
        addFunction( new Func__scalarMatrixMult<Natural    , MatrixRealSymmetric  , MatrixRealSymmetric  >(  )   );
        addFunction( new Func__scalarMatrixMult<Integer    , MatrixRealSymmetric  , MatrixRealSymmetric  >(  )   );
        addFunction( new Func__scalarMatrixMult<Real       , MatrixRealSymmetric  , MatrixRealSymmetric  >(  )   );
        addFunction( new Func__scalarMatrixMult<RealPos    , MatrixRealSymmetric  , MatrixRealSymmetric  >(  )   );

        // subtraction
        addFunction( new Func__sub< Integer                            , Integer               , Integer               >(  )  );
        addFunction( new Func__sub< Real                               , Real                  , Real                  >(  )  );
        addFunction( new Func__sub< MatrixReal                         , MatrixReal            , MatrixReal            >(  )  );
        addFunction( new Func__sub< MatrixRealSymmetric                , MatrixRealSymmetric   , MatrixRealSymmetric     >(  )  );
        addFunction( new Func__sub< ModelVector<Integer>               , ModelVector<Integer>  , ModelVector<Integer>  >(  )  );
        addFunction( new Func__sub< ModelVector<Real>                  , ModelVector<Real>     , ModelVector<Real>     >(  )  );
        addFunction( new Func__vectorScalarSub<Integer                 , Integer               , Integer                   >(  )   );
        addFunction( new Func__vectorScalarSub<Real                    , Real                  , Real                      >(  )   );
        addFunction( new Func__scalarVectorSub<Integer                 , Integer               , Integer                   >(  )   );
        addFunction( new Func__scalarVectorSub<Real                    , Real                  , Real                      >(  )   );

        // modulo
        addFunction( new Func__mod() );

        // exponentiation
        addFunction( new Func_power() );
        addFunction( new Func_powerVector() );


        // type conversion
//        addFunction( "_Natural2Integer",            new Func__conversion<Natural, Integer>()        );
        addFunction( new Func__conversion<Natural, Real>()           );
        addFunction( new Func__conversion<Natural, RealPos>()        );
        addFunction( new Func__conversion<Integer, Real>()           );
//        addFunction( "_Probability2RealPos",        new Func__conversion<Probability, RealPos>()    );
//        addFunction( "_Probability2Real",           new Func__conversion<Probability, Real>()       );
//        addFunction( "_RealPos2Real",               new Func__conversion<RealPos, Real>()           );


        addFunction( new Func__conversion<ModelVector<Natural>, ModelVector<Integer> >()         );
        addFunction( new Func__conversion<ModelVector<Natural>, ModelVector<Real> >()            );
        addFunction( new Func__conversion<ModelVector<Natural>, ModelVector<RealPos> >()         );
        addFunction( new Func__conversion<ModelVector<Integer>, ModelVector<Real> >()            );
        addFunction( new Func__conversion<ModelVector<RealPos>, ModelVector<Real> >()            );
        addFunction( new Func__conversion<ModelVector<Probability>, ModelVector<RealPos> >()     );
        addFunction( new Func__conversion<ModelVector<Probability>, ModelVector<Real> >()        );
        addFunction( new Func__conversion<Simplex, ModelVector<Real> >()                         );
        addFunction( new Func__conversion<Simplex, ModelVector<RealPos> >()                      );
        addFunction( new Func__conversion<ModelVector<ModelVector<Natural> >, ModelVector<ModelVector<Integer> > >()         );
//        addFunction( new Func__conversion<ModelVector<ModelVector<Natural> >, ModelVector<ModelVector<Real> > >()            );
//        addFunction( new Func__conversion<ModelVector<ModelVector<Natural> >, ModelVector<ModelVector<RealPos> > >()         );
//        addFunction( new Func__conversion<ModelVector<ModelVector<Integer> >, ModelVector<ModelVector<Real> > >()            );
        addFunction( new Func__conversion<ModelVector<ModelVector<RealPos> >, ModelVector<ModelVector<Real> > >()            );
        addFunction( new Func__conversion<ModelVector<ModelVector<Probability> >, ModelVector<ModelVector<RealPos> > >()     );
        addFunction( new Func__conversion<ModelVector<ModelVector<Probability> >, ModelVector<ModelVector<Real> > >()        );


        // index operation
        addFunction( new Func__simplexIndexOperator() );



        /* Input/output functions (in folder "functions/io") */
        addFunction( new Func_ancestralStateTree()                      );
        addFunction( new Func_annotateTree()                            );
		addFunction( new Func_characterMapTree()                        );
        addFunction( new Func_consensusTree()                           );
        addFunction( new Func_convertToPhylowood()                      );
        addFunction( new Func_fileExists()                              );
        addFunction( new Func_listFiles()                               );
        addFunction( new Func_readRegionalFeatures()                    );
        addFunction( new Func_maxdiff()                                 );
        addFunction( new Func_mapTree()                                 );
        addFunction( new Func_mccTree()                                 );
        addFunction( new Func_loadPlugin()                                 );
        addFunction( new Func_module()                                  );
        addFunction( new Func_readAncestralStateTreeTrace()             );
		addFunction( new Func_readAncestralStateTrace()                 );
        addFunction( new Func_readAtlas()                               );
		addFunction( new Func_readBranchLengthTrees()                   );
        addFunction( new Func_readContinuousCharacterData()             );
        addFunction( new Func_readDelimitedCharacterData()              );
        addFunction( new Func_readDelimitedDataFile()                   );
        addFunction( new Func_readDiscreteCharacterData()               );
		addFunction( new Func_readDistanceMatrix()                      );
        addFunction( new Func_readCharacterData()              );
        addFunction( new Func_readMatrix()                              );
        addFunction( new Func_readRelativeNodeAgeConstraints()          );
        addFunction( new Func_readRelativeNodeAgeWeightedConstraints()  );
        addFunction( new Func_TaxonReader()                             );
        addFunction( new Func_readStochasticVariableTrace()             );
        addFunction( new Func_readTrace()                               );
        addFunction( new Func_readTrees()                               );
        addFunction( new Func_readTreeTrace()                           );
		addFunction( new Func_readVCF()                                 );
        addFunction( new Func_source()                                  );
        addFunction( new Func_summarizeCharacterMaps()                  );
        addFunction( new Func_treeTrace()                               );
        addFunction( new Func_write()                                   );
        addFunction( new Func_writeDelimitedCharacterData()             );
        addFunction( new Func_writeFasta()                              );
        addFunction( new Func_writeNexus()                              );

    }
    catch(RbException& rbException)
    {

        RBOUT("Caught an exception while initializing the workspace\n");
        std::ostringstream msg;
        rbException.print(msg);
        msg << std::endl;
        RBOUT(msg.str());

        RBOUT("Please report this bug to the RevBayes Development Core Team");

        RBOUT("Press any character to exit the program.");
        getchar();
        exit(1);
    }

}
