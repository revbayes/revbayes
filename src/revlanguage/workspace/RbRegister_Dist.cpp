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
#include <cmath>
#include <cstdio>
#include <string>

/* Files including helper classes */
#include "AddContinuousDistribution.h"
#include "AddDistribution.h"
#include "RbException.h"
#include "RlUserInterface.h"
#include "Workspace.h"

/// Miscellaneous types ///

#include "AbstractHomologousDiscreteCharacterData.h"
#include "AverageDistanceMatrix.h"
#include "ConstantNode.h"
#include "ContinuousCharacterData.h"
#include "DagMemberFunction.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DirichletProcessPriorDistribution.h"
#include "DistanceMatrix.h"
#include "DistributionMemberFunction.h"
#include "DynamicNode.h"
#include "EmpiricalSampleDistribution.h"
#include "EventDistribution.h"
#include "IndirectReferenceFunction.h"
#include "MatrixBoolean.h"
#include "MatrixReal.h"
#include "MixtureDistribution.h"
#include "ModelObject.h"
#include "MultiValueEvent.h"
#include "ProbabilityDensityFunction.h"
#include "RateGenerator.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "RevPtr.h"
#include "ReversibleJumpMixtureConstantDistribution.h"
#include "RlAbstractHomologousDiscreteCharacterData.h"
#include "RlAverageDistanceMatrix.h"
#include "RlConstantNode.h"
#include "RlContinuousCharacterData.h"
#include "RlDagMemberFunction.h"
#include "RlDeterministicNode.h"
#include "RlDistanceMatrix.h"
#include "RlDistributionMemberFunction.h"
#include "RlMultiValueEvent.h"
#include "RlOrderedEventTimes.h"
#include "RlOrderedEvents.h"
#include "RlStochasticNode.h"
#include "RlTimeTree.h"
#include "RlTree.h"
#include "RlTypedDistribution.h"
#include "RlTypedFunction.h"
#include "Simplex.h"
#include "StochasticNode.h"
#include "Tree.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"
#include "TypedFunction.h"
#include "UniformPartitioningDistribution.h"
#include "UserFunctionNode.h"
#include "WeightedSampleDistribution.h"
#include "WorkspaceToCoreWrapperObject.h"

/* Base types (in folder "datatypes") */

/* Primitive types (in folder "datatypes/basic") */
#include "Integer.h"
#include "IntegerPos.h"
#include "Natural.h"
#include "Probability.h"
#include "Real.h"
#include "RealPos.h"

/* Container types (in folder "datatypes/container") */
#include "ModelVector.h"
#include "WorkspaceVector.h"

/* Evolution types (in folder "datatypes/phylogenetics") */

/* Character state types (in folder "datatypes/phylogenetics/character") */

/* Character data types (in folder "datatypes/phylogenetics/characterdata") */

/* Tree types (in folder "datatypes/phylogenetics/trees") */


/// Types ///

/* These types are needed as template types for the moves */
#include "RlBranchLengthTree.h"
#include "RlRateGenerator.h"

/* Math types (in folder "datatypes/math") */
#include "RlMatrixBoolean.h"
#include "RlMatrixReal.h"
#include "RlMatrixRealSymmetric.h"
#include "RlSimplex.h"

/// Distributions ///

/* Distribution types (in folder "distributions") */

#include "Dist_EmpiricalSample.h"
#include "Dist_WeightedSample.h"

/* Character evolution models (in folder "distributions/phylogenetics/character") */
#include "Dist_CTMC.h"
#include "Dist_phyloCTMC.h"
#include "Dist_phyloCTMCDASequence.h"
#include "Dist_phyloCTMCDASiteIID.h"
#include "Dist_phyloCTMCClado.h"
#include "Dist_phyloCTMCDollo.h"

/* Branch rate priors (in folder "distributions/phylogenetics/tree") */

/* Trait evolution models (in folder "distributions/phylogenetics/branchrates") */
#include "Dist_PhyloBranchRateBM.h"
#include "Dist_PhyloBrownian.h"
#include "Dist_PhyloBrownianMVN.h"
#include "Dist_PhyloBrownianREML.h"
#include "Dist_PhyloBrownianMultiSampleREML.h"
#include "Dist_PhyloBrownianProcessStateDependent.h"
#include "Dist_PhyloMvtBrownian.h"
#include "Dist_PhyloMultiSampleOrnsteinUhlenbeck.h"
#include "Dist_PhyloMultiSampleOrnsteinUhlenbeckREML.h"
#include "Dist_PhyloMultivariateBrownianREML.h"
#include "Dist_PhyloMultivariateBrownianMultiSampleREML.h"
#include "Dist_PhyloOrnsteinUhlenbeck.h"
#include "Dist_PhyloOrnsteinUhlenbeckMVN.h"
#include "Dist_PhyloOrnsteinUhlenbeckPruning.h"
#include "Dist_PhyloOrnsteinUhlenbeckThreePoint.h"
#include "Dist_PhyloOrnsteinUhlenbeckStateDependent.h"
#include "Dist_PhyloBrownianProcessStateDependentTrend.h"
#include "Dist_PhyloWhiteNoise.h"

/* Tree priors (in folder "distributions/phylogenetics/tree") */
#include "Dist_bdp.h"
#include "Dist_bdp_complete.h"
#include "Dist_BDSTP.h"
#include "Dist_FBDP.h"
#include "Dist_BirthDeathBurstProcess.h"
#include "Dist_BranchRateTree.h"
#include "Dist_CharacterDependentBirthDeathProcess.h"
#include "Dist_Coalescent.h"
#include "Dist_CoalescentDemography.h"
#include "Dist_CoalescentSkyline.h"
#include "Dist_conditionedBirthDeathShiftProcessContinuous.h"
#include "Dist_ConstrainedTopology.h"
#include "Dist_ConstrainedUnrootedTopology.h"
#include "Dist_ConstrainedNodeAge.h"
#include "Dist_ConstrainedNodeOrder.h"
#include "Dist_WeightedConstrainedNodeOrder.h"
#include "Dist_DuplicationLoss.h"
#include "Dist_FBDRP.h"
#include "Dist_FBDSP.h"
#include "Dist_GLHBDSP.h"
#include "Dist_constPopMultispCoal.h"
#include "Dist_divDepYuleProcess.h"
#include "Dist_empiricalTree.h"
#include "Dist_episodicBirthDeath.h"
#include "Dist_heterogeneousRateBirthDeath.h"
#include "Dist_multispeciesCoalescentInverseGammaPrior.h"
#include "Dist_multispeciesCoalescentUniformPrior.h"
#include "Dist_MultispeciesCoalescentMigration.h"
#include "Dist_outgroupBirthDeath.h"
#include "Dist_PhylodynamicBDP.h"
#include "Dist_phyloDistanceGamma.h"
#include "Dist_sampledSpeciationBirthDeathProcess.h"
#include "Dist_occurrenceBirthDeathProcess.h"
#include "Dist_TimeVaryingStateDependentSpeciationExtinctionProcess.h"
#include "Dist_UltrametricTree.h"
#include "Dist_uniformTimeTree.h"
#include "Dist_uniformSerialSampledTimeTree.h"
#include "Dist_uniformTopology.h"
#include "Dist_uniformTopologyBranchLength.h"

/* Distributions on simple variables (in folder "distributions/math") */
#include "Dist_bernoulli.h"
#include "Dist_beta.h"
#include "Dist_bimodalLnorm.h"
#include "Dist_bimodalNorm.h"
#include "Dist_binomial.h"
#include "Dist_bivariatePoisson.h"
#include "Dist_categorical.h"
#include "Dist_Cauchy.h"
#include "Dist_chisq.h"
#include "Dist_cppNormal.h"
#include "Dist_decomposedInverseWishart.h"
#include "Dist_dirichlet.h"
#include "Dist_exponential.h"
#include "Dist_exponentialError.h"
#include "Dist_exponentialNegativeOffset.h"
#include "Dist_gamma.h"
#include "Dist_geom.h"
#include "Dist_GilbertGraph.h"
#include "Dist_halfCauchy.h"
#include "Dist_halfCauchyPositive.h"
#include "Dist_halfNormal.h"
#include "Dist_halfNormalPositive.h"
#include "Dist_inverseGamma.h"
#include "Dist_inverseWishart.h"
#include "Dist_Laplace.h"
#include "Dist_LKJ.h"
#include "Dist_LKJPartial.h"
#include "Dist_lnorm.h"
#include "Dist_lnormNegativeOffset.h"
#include "Dist_logExponential.h"
#include "Dist_logUniform.h"
#include "Dist_multinomial.h"
#include "Dist_multivariateNorm.h"
#include "Dist_nbinomial.h"
#include "Dist_norm.h"
#include "Dist_normTruncated.h"
#include "Dist_normTruncatedPositive.h"
#include "Dist_pointMass.h"
#include "Dist_pointMassPositive.h"
#include "Dist_poisson.h"
#include "Dist_scaledDirichlet.h"
#include "Dist_softBoundUniformNormal.h"
#include "Dist_studentT.h"
#include "Dist_unif.h"
#include "Dist_unifPositive.h"
#include "Dist_unifProbability.h"
#include "Dist_UniformInteger.h"
#include "Dist_UniformNatural.h"
#include "Dist_varianceGamma.h"
#include "Dist_whiteNoise.h"
#include "Dist_wishart.h"
#include "Process_OrnsteinUhlenbeck.h"

/* Mixture distributions (in folder "distributions/mixture") */
#include "Dist_AutocorrelatedEvent.h"
#include "Dist_dpp.h"
#include "Dist_event.h"
#include "Dist_IID.h"
#include "Dist_Log.h"
#include "Dist_MultivariateLog.h"
#include "Dist_markovTimes.h"
#include "Dist_markovEvents.h"
#include "Dist_mixture.h"
#include "Dist_mixtureAnalytical.h"
#include "Dist_mixtureVector.h"
#include "Dist_MultiValueEvent.h"
#include "Dist_reversibleJumpMixtureConstant.h"
#include "Dist_upp.h"

#include "Transform_Exp.h"
#include "Transform_Log.h"
#include "Transform_Logit.h"
#include "Transform_InvLogit.h"

#include "Transform_Add.h"
#include "Transform_Sub1.h"
#include "Transform_Sub2.h"
#include "Transform_Mul.h"

#include "Transform_Vector_Exp.h"
#include "Transform_Vector_Log.h"
#include "Transform_Vector_Logit.h"
#include "Transform_Vector_Invlogit.h"

/// Functions ///

/* Helper functions for creating functions (in folder "functions") */
#include "DistributionFunctionPdf.h"
#include "DistributionFunctionRv.h"

/* Argument rules (in folder "functions/argumentrules") */


/** Initialize global workspace */
void RevLanguage::Workspace::initializeDistGlobalWorkspace(void)
{

    try
    {
        ///////////////////////////////////////////////////
        /* Add distributions (in folder "distributions") */
        ///////////////////////////////////////////////////
        addType( new WorkspaceVector<Distribution>( ) );


        /* Evolutionary processes (in folder "distributions/phylogenetics") */

        /* Branch rate processes (in folder "distributions/phylogenetics/branchrate") */

        // white noise process
        AddDistribution< ModelVector<RealPos>       >( new Dist_PhyloWhiteNoise()          );

        /* trait evolution (in folder "distributions/phylogenetics/branchrate") */

        AddDistribution< ModelVector<RealPos>       >( new Dist_PhyloBranchRateBM()                             );

        // brownian motion
        AddDistribution< ModelVector<Real>          >( new Dist_PhyloBrownian()                                 );
        AddDistribution< ContinuousCharacterData    >( new Dist_PhyloBrownianREML()                             );
        AddDistribution< ContinuousCharacterData    >( new Dist_PhyloBrownianMVN()                              );
        AddDistribution< ContinuousCharacterData    >( new Dist_PhyloBrownianMultiSampleREML()                  );
        AddDistribution< ContinuousCharacterData    >( new Dist_PhyloBrownianProcessStateDependent()            );
        AddDistribution< ContinuousCharacterData    >( new Dist_PhyloMultiSampleOrnsteinUhlenbeck()             );
        AddDistribution< ContinuousCharacterData    >( new Dist_PhyloMultiSampleOrnsteinUhlenbeckREML()         );
        AddDistribution< ContinuousCharacterData    >( new Dist_PhyloMultivariateBrownianREML()                 );
        AddDistribution< ContinuousCharacterData    >( new Dist_PhyloMultivariateBrownianMultiSampleREML()      );
        AddDistribution< ModelVector<Real>          >( new Dist_PhyloOrnsteinUhlenbeck()                        );
        AddDistribution< ContinuousCharacterData    >( new Dist_PhyloOrnsteinUhlenbeckMVN()                     );
        AddDistribution< ContinuousCharacterData    >( new Dist_PhyloOrnsteinUhlenbeckPruning()                 );
        AddDistribution< ContinuousCharacterData    >( new Dist_PhyloOrnsteinUhlenbeckThreePoint()              );
        AddDistribution< ContinuousCharacterData    >( new Dist_PhyloOrnsteinUhlenbeckStateDependent()          );
        AddDistribution< ContinuousCharacterData    >( new Dist_PhyloBrownianProcessStateDependentTrend()       );

        // multivariate brownian motion
        AddDistribution< ModelVector< ModelVector<Real> > >( new Dist_PhyloMvtBrownian() );

        /* Character state evolution processes (in folder "distributions/phylogenetics/character") */

        // simple phylogenetic CTMC on fixed number of discrete states
//        AddDistribution< AbstractHomologousDiscreteCharacterData >( new Dist_phyloCTMC() );
//        AddDistribution< AbstractHomologousDiscreteCharacterData >( new Dist_phyloDACTMC() );
//        AddDistribution< AbstractHomologousDiscreteCharacterData >( new Dist_phyloCTMCClado() );
        addDistribution( new Dist_CTMC() );
        addDistribution( new Dist_phyloCTMC() );
        addDistribution( new Dist_phyloCTMCDASequence() );
        addDistribution( new Dist_phyloCTMCDASiteIID() );
        addDistribution( new Dist_phyloCTMCClado() );
        addDistribution( new Dist_phyloCTMCDollo() );

        /* Tree distributions (in folder "distributions/phylogenetics/tree") */

        // constant rate birth-death process
        AddDistribution< TimeTree                   >( new Dist_bdp());
        AddDistribution< TimeTree                   >( new Dist_bdp_complete());

        AddDistribution< TimeTree                   >( new Dist_BirthDeathBurstProcess());

        AddDistribution< TimeTree                   >( new Dist_CharacterDependentBirthDeathProcess() );
        AddDistribution< TimeTree                   >( new Dist_heterogeneousRateBirthDeath() );
        AddDistribution< TimeTree                   >( new Dist_conditionedBirthDeathShiftProcessContinuous() );
        AddDistribution< TimeTree                   >( new Dist_outgroupBirthDeath() );
        AddDistribution< TimeTree                   >( new Dist_sampledSpeciationBirthDeathProcess() );
        AddDistribution< TimeTree                   >( new Dist_TimeVaryingStateDependentSpeciationExtinctionProcess() );
        AddDistribution< TimeTree                   >( new Dist_GLHBDSP() );

        // fossilized-birth-death range processes
        AddDistribution< MatrixReal                 >( new Dist_FBDRP());
        AddDistribution< TimeTree                   >( new Dist_FBDSP());

        // birth-death-sampling-treatment processes and submodels
        AddDistribution< TimeTree                   >( new Dist_BDSTP());
        AddDistribution< TimeTree                   >( new Dist_FBDP());
        AddDistribution< TimeTree                   >( new Dist_PhylodynamicBDP());

        // diversity-dependent pure-birth process
        AddDistribution< TimeTree                   >( new Dist_divDepYuleProcess() );

        // episodic birth-death process
        AddDistribution< TimeTree                   >( new Dist_episodicBirthDeath() );

        // coalescent (constant population sizes)
        AddDistribution< TimeTree                   >( new Dist_Coalescent() );
        
        // coalescent (population sizes via demography functions)
        AddDistribution< TimeTree                   >( new Dist_CoalescentDemography() );

        // coalescent (skyline population sizes)
        AddDistribution< TimeTree                   >( new Dist_CoalescentSkyline() );

        // duplication loss process
        AddDistribution< TimeTree                   >( new Dist_DuplicationLoss() );

        // multispecies coalescent (per branch constant population sizes)
        AddDistribution< TimeTree                   >( new Dist_constPopMultispCoal() );
        AddDistribution< ModelVector<TimeTree>      >( new Dist_multispeciesCoalescentInverseGammaPrior() );
        AddDistribution< TimeTree                   >( new Dist_multispeciesCoalescentUniformPrior() );
        AddDistribution< TimeTree                   >( new Dist_MultispeciesCoalescentMigration() );

        // constrained node age distribution
        AddDistribution< TimeTree                   >( new Dist_ConstrainedNodeAge() );

        // constrained node order distribution
        AddDistribution< TimeTree                   >( new Dist_ConstrainedNodeOrder() );

        // constrained node order distribution
        AddDistribution< TimeTree                   >( new Dist_WeightedConstrainedNodeOrder() );

        // constrained topology distribution
        AddDistribution< TimeTree                   >( new Dist_ConstrainedTopology() );

        // constrained topology distribution
        AddDistribution< BranchLengthTree           >( new Dist_ConstrainedUnrootedTopology() );

        // uniform time tree distribution
        AddDistribution< TimeTree                   >( new Dist_uniformTimeTree() );

        // uniform serial-sampled time tree distribution
        AddDistribution< TimeTree                   >( new Dist_uniformSerialSampledTimeTree() );

        // occurrence birth death process tree distribution
        AddDistribution< TimeTree                   >( new Dist_occurrenceBirthDeathProcess() );

        // uniform topology distribution
        AddDistribution< BranchLengthTree           >( new Dist_uniformTopology() );

        // uniform topology with branch lengths distribution
        AddDistribution< BranchLengthTree           >( new Dist_uniformTopologyBranchLength() );

        // empirical tree distributions
        AddDistribution< Tree                       >( new Dist_empiricalTree() );

        // ultrametric tree distributions
        AddDistribution< TimeTree                   >( new Dist_UltrametricTree() );

        // branch rate tree distributions
        AddDistribution< BranchLengthTree           >( new Dist_BranchRateTree() );

        // Distance Matrix Gamma distribution
        AddDistribution< DistanceMatrix             >( new Dist_phyloDistanceGamma() );


        /* Statistical distributions on simple variables (in folder "distributions/math") */

        // bernoulli distribution
        AddDistribution< Natural                    >( new Dist_bernoulli() );

        // binomial distribution
        AddDistribution< Natural                    >( new Dist_binomial() );

        // bivariate poisson distribution
        AddDistribution< ModelVector<Natural>       >( new Dist_bivariatePoisson() );

        // negative binomial distribution
        AddDistribution< Natural                    >( new Dist_nbinomial() );

        // beta distribution
//        AddContinuousDistribution< Probability >( new Dist_beta() );
        AddDistribution< Probability                >( new Dist_beta() );

        // bimodal normal distribution
        AddContinuousDistribution< Real             >( new Dist_bimodalNorm() );

        // bimodal lognormal distribution
        AddContinuousDistribution< RealPos          >( new Dist_bimodalLnorm() );

        // categorical distribution
        AddDistribution< Natural                    >( new Dist_categorical() );

        // Cauchy distribution
        AddContinuousDistribution< Real             >( new Dist_Cauchy() );

        // chi-square distribution
        AddContinuousDistribution< RealPos          >( new Dist_chisq() );

        // Student's t distribution
        AddContinuousDistribution< Real             >(new Dist_studentT() );

        // compound Poisson w/ normal kernel
        AddDistribution< Real                       >( new Dist_cppNormal() );

        // dirichlet distribution
        AddDistribution< Simplex                    >( new Dist_dirichlet() );

        // scaled Dirichlet distribution
        AddDistribution< Simplex                    >( new Dist_scaledDirichlet() );

        // gamma distribution
        AddContinuousDistribution< RealPos          >( new Dist_gamma() );

        // geometric distribution
        AddDistribution< Natural                    >( new Dist_geom() );

        // half-Cauchy distribution
        AddContinuousDistribution< Real             >( new Dist_halfCauchy() );
        AddContinuousDistribution< RealPos          >( new Dist_halfCauchyPositive() );

        // half-Normal distribution
        AddContinuousDistribution< Real             >( new Dist_halfNormal() );
        AddContinuousDistribution< RealPos          >( new Dist_halfNormalPositive() );

        // inverse-gamma distribution
        AddContinuousDistribution< RealPos          >( new Dist_inverseGamma() );

        // point mass distribution
        AddDistribution< Real                       >( new Dist_pointMass() );
        AddDistribution< RealPos                    >( new Dist_pointMassPositive() );

        // poisson distribution
        AddDistribution< Natural                    >( new Dist_poisson() );

        // exponential distribution
        AddContinuousDistribution< RealPos          >( new Dist_exponential() );
        AddContinuousDistribution< Real             >( new Dist_exponentialNegativeOffset() );

        // Laplace distribution
        AddContinuousDistribution< Real             >( new Dist_Laplace() );

        // LKJ distribution
        AddDistribution< MatrixRealSymmetric        >( new Dist_LKJ() );
        AddDistribution< MatrixRealSymmetric        >( new Dist_LKJPartial() );

        // random graph distributions
        AddDistribution< MatrixRealSymmetric        >( new Dist_GilbertGraph() );

        // lognormal distribution
        AddContinuousDistribution< RealPos          >( new Dist_lnorm() );
        AddContinuousDistribution< Real             >( new Dist_lnormNegativeOffset() );

        // LogExponential distribution
        AddContinuousDistribution< Real             >( new Dist_logExponential() );

        // LogUniform distribution
        AddContinuousDistribution< RealPos          >( new Dist_logUniform() );

        // multinomial distribution
        AddDistribution< ModelVector<Natural>       >( new Dist_multinomial() );

        // multivariate normal distribution
        AddDistribution< ModelVector<Real>          >( new Dist_multivariateNorm());

        // normal distribution
        AddContinuousDistribution< Real             >( new Dist_norm() );
        AddContinuousDistribution< Real             >( new Dist_normTruncated() );
        AddContinuousDistribution< RealPos          >( new Dist_normTruncatedPositive() );

        // Uniform distribution with normal distributed bounds
        AddContinuousDistribution< Real             >( new Dist_SoftBoundUniformNormal() );

        // uniform distribution
        AddContinuousDistribution< Real             >( new Dist_unif() );
        AddContinuousDistribution< RealPos          >( new Dist_unifPositive() );
//        AddContinuousDistribution< Probability      >( new Dist_unifProbability() );
        AddDistribution< Probability                >( new Dist_unifProbability() );
        AddDistribution< Integer                    >( new Dist_UniformInteger() );
        AddDistribution< Natural                    >( new Dist_UniformNatural() );
        AddContinuousDistribution< Real             >( new Dist_varianceGamma() );

        // White-Noise process
        AddContinuousDistribution< RealPos          >( new Dist_whiteNoise() );

        // Wishart distribution
        AddDistribution< MatrixRealSymmetric        >( new Dist_wishart() );

        // inverse Wishart distribution
        AddDistribution< MatrixRealSymmetric        >( new Dist_inverseWishart() );

        // and the so-called "decomposed" Inverse Wishart
        AddDistribution< MatrixReal                 >( new Dist_decomposedInverseWishart() );

        // Exponential error distribution for matrix distance from average distance matrix
        AddDistribution< AverageDistanceMatrix      >( new Dist_exponentialError());

        /* Empirical sample distributions (in folder "distributions/mixture") */
        AddDistribution< ModelVector<Natural>       >( new Dist_EmpiricalSample<Natural>());
        AddDistribution< ModelVector<Real>          >( new Dist_EmpiricalSample<Real>());
        AddDistribution< ModelVector<RealPos>       >( new Dist_EmpiricalSample<RealPos>());
        AddDistribution< ModelVector<TimeTree>      >( new Dist_EmpiricalSample<TimeTree>());
        AddDistribution< ModelVector< ModelVector<TimeTree> >       >( new Dist_EmpiricalSample< ModelVector<TimeTree> >());
        AddDistribution< ModelVector<BranchLengthTree>              >( new Dist_EmpiricalSample<BranchLengthTree>());
        AddDistribution< ModelVector<TimeTree>      >( new Dist_WeightedSample<TimeTree>());
        AddDistribution< ModelVector<AbstractHomologousDiscreteCharacterData>      >( new Dist_WeightedSample<AbstractHomologousDiscreteCharacterData>());



        // dirichlet process prior distribution
        AddDistribution< ModelVector<Real>          >( new Dist_dpp<Real>()         );
        AddDistribution< ModelVector<RealPos>       >( new Dist_dpp<RealPos>()      );
        AddDistribution< ModelVector<Natural>       >( new Dist_dpp<Natural>()      );
        AddDistribution< ModelVector<Integer>       >( new Dist_dpp<Integer>()      );
        AddDistribution< ModelVector<Probability>   >( new Dist_dpp<Probability>()  );
        AddDistribution< ModelVector<Simplex>       >( new Dist_dpp<Simplex>()      );

        // event distribution
        AddDistribution< MultiValueEvent            >( new Dist_AutocorrelatedEvent() );
        AddDistribution< ModelVector<Real>          >( new Dist_event<Real>()         );
        AddDistribution< ModelVector<RealPos>       >( new Dist_event<RealPos>()      );
        AddDistribution< ModelVector<Natural>       >( new Dist_event<Natural>()      );
        AddDistribution< ModelVector<Integer>       >( new Dist_event<Integer>()      );
        AddDistribution< ModelVector<Probability>   >( new Dist_event<Probability>()  );
        AddDistribution< MultiValueEvent            >( new Dist_MultiValueEvent()     );

        // IID distribution
        AddDistribution< ModelVector<Real>          >( new Dist_IID<Real>()         );
        AddDistribution< ModelVector<RealPos>       >( new Dist_IID<RealPos>()      );
        AddDistribution< ModelVector<Natural>       >( new Dist_IID<Natural>()      );
        AddDistribution< ModelVector<Integer>       >( new Dist_IID<Integer>()      );
        AddDistribution< ModelVector<Probability>   >( new Dist_IID<Probability>()  );

        AddDistribution< RealPos                    >( new Dist_Log()               );
        AddDistribution< ModelVector<RealPos>       >( new Dist_MultivariateLog()   );
        AddDistribution< RealPos                    >( new Transform_Exp()          );
        AddDistribution< Real                       >( new Transform_Log()          );
        AddDistribution< Real                       >( new Transform_Logit()        );
        AddDistribution< Probability                >( new Transform_InvLogit()     );
        AddDistribution< RealPos                    >( new Transform_Add<RealPos    , false>() );
        AddDistribution< Real                       >( new Transform_Add<Real       , false>() );
        AddDistribution< Probability                >( new Transform_Mul<Probability, false>() );
        AddDistribution< RealPos                    >( new Transform_Mul<RealPos    , false>() );
        AddDistribution< Real                       >( new Transform_Mul<Real       , false>() );

        AddDistribution< Real                       >( new Transform_Sub1()        );
        AddDistribution< Real                       >( new Transform_Sub2()        );
        AddDistribution< RealPos                    >( new Transform_Add<RealPos>()     );
        AddDistribution< Real                       >( new Transform_Add<Real>()        );
        AddDistribution< Probability                >( new Transform_Mul<Probability>() );
        AddDistribution< RealPos                    >( new Transform_Mul<RealPos>()     );
        AddDistribution< Real                       >( new Transform_Mul<Real>()        );

        AddDistribution< ModelVector<RealPos>       >( new Transform_Vector_Exp()   );
        AddDistribution< ModelVector<Real>          >( new Transform_Vector_Log()   );
        AddDistribution< ModelVector<Real>          >( new Transform_Vector_Logit() );
        AddDistribution< ModelVector<Probability>   >( new Transform_Vector_InvLogit() );

        // uniform partitions prior
        AddDistribution< ModelVector<RealPos>       >( new Dist_upp<RealPos>() );

        // mixture distribution
        AddDistribution< Real                       >( new Dist_mixture<Real>() );
        AddDistribution< RealPos                    >( new Dist_mixture<RealPos>() );
        AddDistribution< Natural                    >( new Dist_mixture<Natural>() );
        AddDistribution< Integer                    >( new Dist_mixture<Integer>() );
        AddDistribution< Probability                >( new Dist_mixture<Probability>() );
        AddDistribution< Simplex                    >( new Dist_mixture<Simplex>() );
        AddDistribution< ModelVector<Real>          >( new Dist_mixture< ModelVector<Real> >() );
        AddDistribution< ModelVector<RealPos>       >( new Dist_mixture< ModelVector<RealPos> >() );
//        AddDistribution< RateGenerator              >( new Dist_mixture<RateGenerator>() );
        addDistribution( new Dist_mixture<RateGenerator>() );
        AddDistribution< TimeTree                   >( new Dist_mixture<TimeTree>() );


        // analytical mixture distribution
        AddDistribution< Real                       >( new Dist_mixtureAnalytical<Real>() );
        AddDistribution< RealPos                    >( new Dist_mixtureAnalytical<RealPos>() );
        AddDistribution< Natural                    >( new Dist_mixtureAnalytical<Natural>() );
        AddDistribution< Integer                    >( new Dist_mixtureAnalytical<Integer>() );
        AddDistribution< Probability                >( new Dist_mixtureAnalytical<Probability>() );
        AddDistribution< Simplex                    >( new Dist_mixtureAnalytical<Simplex>() );
        AddDistribution< ModelVector<Real>          >( new Dist_mixtureAnalytical< ModelVector<Real> >() );
        AddDistribution< ModelVector<RealPos>       >( new Dist_mixtureAnalytical< ModelVector<RealPos> >() );
        AddDistribution< TimeTree                   >( new Dist_mixtureAnalytical<TimeTree>() );

        AddDistribution< ModelVector<Real>          >( new Dist_mixtureVector<Real>() );
        AddDistribution< ModelVector<RealPos>       >( new Dist_mixtureVector<RealPos>() );

        // Ornstein-Uhlenbeck process
        AddDistribution< Real                       >( new OrnsteinUhlenbeckProcess() );

        // mixture distribution
        AddDistribution< Real                       >( new Dist_reversibleJumpMixtureConstant<Real>() );
        AddDistribution< RealPos                    >( new Dist_reversibleJumpMixtureConstant<RealPos>() );
        AddDistribution< Natural                    >( new Dist_reversibleJumpMixtureConstant<Natural>() );
        AddDistribution< Integer                    >( new Dist_reversibleJumpMixtureConstant<Integer>() );
        AddDistribution< Probability                >( new Dist_reversibleJumpMixtureConstant<Probability>() );
        AddDistribution< Simplex                    >( new Dist_reversibleJumpMixtureConstant<Simplex>() );
        AddDistribution< ModelVector<Natural>       >( new Dist_reversibleJumpMixtureConstant<ModelVector<Natural> >() );
        AddDistribution< TimeTree                   >( new Dist_reversibleJumpMixtureConstant<TimeTree>() );
        AddDistribution< BranchLengthTree           >( new Dist_reversibleJumpMixtureConstant<BranchLengthTree>() );

        // markov events
        AddDistribution< RlOrderedEventTimes          >( new Dist_markovTimes() );
        AddDistribution< RlOrderedEvents<Real>        >( new Dist_markovEvents<Real>() );
        AddDistribution< RlOrderedEvents<RealPos>     >( new Dist_markovEvents<RealPos>() );
        AddDistribution< RlOrderedEvents<Probability> >( new Dist_markovEvents<Probability>() );

        AddDistribution< RlOrderedEvents<ModelVector<Real> >        >( new Dist_markovEvents<ModelVector<Real>        >() );
        AddDistribution< RlOrderedEvents<ModelVector<RealPos> >     >( new Dist_markovEvents<ModelVector<RealPos>     >() );
        AddDistribution< RlOrderedEvents<ModelVector<Probability> > >( new Dist_markovEvents<ModelVector<Probability> >() );

        /* Now we have added all primitive and complex data types and can start type checking */
        Workspace::globalWorkspace().typesInitialized = true;
        Workspace::userWorkspace().typesInitialized   = true;

    }
    catch(RbException& rbException)
    {

        RBOUT("Caught an exception while initializing distributions in the workspace\n");
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
