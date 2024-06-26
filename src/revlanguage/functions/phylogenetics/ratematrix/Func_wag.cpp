#include <cstddef>
#include <iostream>
#include <string>
#include <vector>

#include "Func_wag.h"

#include "RlDeterministicNode.h"
#include "RlRateMatrix.h"
#include "RlSimplex.h"
#include "GenericFunction.h"
#include "Simplex.h"
#include "TypedDagNode.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "RateGenerator.h"
#include "RbException.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "RlTypedFunction.h"
#include "StringUtilities.h"
#include "TypeSpec.h"
#include "TypedFunction.h"

#include "ConcreteTimeReversibleRateMatrix.h"

namespace RevLanguage { class RevObject; }

using namespace RevLanguage;

using std::vector;

RevBayesCore::MatrixReal wag_exchange_rates()
{
    RevBayesCore::MatrixReal m(20);
    m[ 0][ 0] = 0.0000000; m[ 1][ 0] = 0.5515710; m[ 2][ 0] = 0.5098480; m[ 3][ 0] = 0.7389980; m[ 4][ 0] = 1.0270400;
    m[ 5][ 0] = 0.9085980; m[ 6][ 0] = 1.5828500; m[ 7][ 0] = 1.4167200; m[ 8][ 0] = 0.3169540; m[ 9][ 0] = 0.1933350;
    m[10][ 0] = 0.3979150; m[11][ 0] = 0.9062650; m[12][ 0] = 0.8934960; m[13][ 0] = 0.2104940; m[14][ 0] = 1.4385500;
    m[15][ 0] = 3.3707900; m[16][ 0] = 2.1211100; m[17][ 0] = 0.1131330; m[18][ 0] = 0.2407350; m[19][ 0] = 2.0060100;
    m[ 0][ 1] = 0.5515710; m[ 1][ 1] = 0.0000000; m[ 2][ 1] = 0.6353460; m[ 3][ 1] = 0.1473040; m[ 4][ 1] = 0.5281910;
    m[ 5][ 1] = 3.0355000; m[ 6][ 1] = 0.4391570; m[ 7][ 1] = 0.5846650; m[ 8][ 1] = 2.1371500; m[ 9][ 1] = 0.1869790;
    m[10][ 1] = 0.4976710; m[11][ 1] = 5.3514200; m[12][ 1] = 0.6831620; m[13][ 1] = 0.1027110; m[14][ 1] = 0.6794890;
    m[15][ 1] = 1.2241900; m[16][ 1] = 0.5544130; m[17][ 1] = 1.1639200; m[18][ 1] = 0.3815330; m[19][ 1] = 0.2518490;
    m[ 0][ 2] = 0.5098480; m[ 1][ 2] = 0.6353460; m[ 2][ 2] = 0.0000000; m[ 3][ 2] = 5.4294200; m[ 4][ 2] = 0.2652560;
    m[ 5][ 2] = 1.5436400; m[ 6][ 2] = 0.9471980; m[ 7][ 2] = 1.1255600; m[ 8][ 2] = 3.9562900; m[ 9][ 2] = 0.5542360;
    m[10][ 2] = 0.1315280; m[11][ 2] = 3.0120100; m[12][ 2] = 0.1982210; m[13][ 2] = 0.0961621; m[14][ 2] = 0.1950810;
    m[15][ 2] = 3.9742300; m[16][ 2] = 2.0300600; m[17][ 2] = 0.0719167; m[18][ 2] = 1.0860000; m[19][ 2] = 0.1962460;
    m[ 0][ 3] = 0.7389980; m[ 1][ 3] = 0.1473040; m[ 2][ 3] = 5.4294200; m[ 3][ 3] = 0.0000000; m[ 4][ 3] = 0.0302949;
    m[ 5][ 3] = 0.6167830; m[ 6][ 3] = 6.1741600; m[ 7][ 3] = 0.8655840; m[ 8][ 3] = 0.9306760; m[ 9][ 3] = 0.0394370;
    m[10][ 3] = 0.0848047; m[11][ 3] = 0.4798550; m[12][ 3] = 0.1037540; m[13][ 3] = 0.0467304; m[14][ 3] = 0.4239840;
    m[15][ 3] = 1.0717600; m[16][ 3] = 0.3748660; m[17][ 3] = 0.1297670; m[18][ 3] = 0.3257110; m[19][ 3] = 0.1523350;
    m[ 0][ 4] = 1.0270400; m[ 1][ 4] = 0.5281910; m[ 2][ 4] = 0.2652560; m[ 3][ 4] = 0.0302949; m[ 4][ 4] = 0.0000000;
    m[ 5][ 4] = 0.0988179; m[ 6][ 4] = 0.0213520; m[ 7][ 4] = 0.3066740; m[ 8][ 4] = 0.2489720; m[ 9][ 4] = 0.1701350;
    m[10][ 4] = 0.3842870; m[11][ 4] = 0.0740339; m[12][ 4] = 0.3904820; m[13][ 4] = 0.3980200; m[14][ 4] = 0.1094040;
    m[15][ 4] = 1.4076600; m[16][ 4] = 0.5129840; m[17][ 4] = 0.7170700; m[18][ 4] = 0.5438330; m[19][ 4] = 1.0021400;
    m[ 0][ 5] = 0.9085980; m[ 1][ 5] = 3.0355000; m[ 2][ 5] = 1.5436400; m[ 3][ 5] = 0.6167830; m[ 4][ 5] = 0.0988179;
    m[ 5][ 5] = 0.0000000; m[ 6][ 5] = 5.4694700; m[ 7][ 5] = 0.3300520; m[ 8][ 5] = 4.2941100; m[ 9][ 5] = 0.1139170;
    m[10][ 5] = 0.8694890; m[11][ 5] = 3.8949000; m[12][ 5] = 1.5452600; m[13][ 5] = 0.0999208; m[14][ 5] = 0.9333720;
    m[15][ 5] = 1.0288700; m[16][ 5] = 0.8579280; m[17][ 5] = 0.2157370; m[18][ 5] = 0.2277100; m[19][ 5] = 0.3012810;
    m[ 0][ 6] = 1.5828500; m[ 1][ 6] = 0.4391570; m[ 2][ 6] = 0.9471980; m[ 3][ 6] = 6.1741600; m[ 4][ 6] = 0.0213520;
    m[ 5][ 6] = 5.4694700; m[ 6][ 6] = 0.0000000; m[ 7][ 6] = 0.5677170; m[ 8][ 6] = 0.5700250; m[ 9][ 6] = 0.1273950;
    m[10][ 6] = 0.1542630; m[11][ 6] = 2.5844300; m[12][ 6] = 0.3151240; m[13][ 6] = 0.0811339; m[14][ 6] = 0.6823550;
    m[15][ 6] = 0.7049390; m[16][ 6] = 0.8227650; m[17][ 6] = 0.1565570; m[18][ 6] = 0.1963030; m[19][ 6] = 0.5887310;
    m[ 0][ 7] = 1.4167200; m[ 1][ 7] = 0.5846650; m[ 2][ 7] = 1.1255600; m[ 3][ 7] = 0.8655840; m[ 4][ 7] = 0.3066740;
    m[ 5][ 7] = 0.3300520; m[ 6][ 7] = 0.5677170; m[ 7][ 7] = 0.0000000; m[ 8][ 7] = 0.2494100; m[ 9][ 7] = 0.0304501;
    m[10][ 7] = 0.0613037; m[11][ 7] = 0.3735580; m[12][ 7] = 0.1741000; m[13][ 7] = 0.0499310; m[14][ 7] = 0.2435700;
    m[15][ 7] = 1.3418200; m[16][ 7] = 0.2258330; m[17][ 7] = 0.3369830; m[18][ 7] = 0.1036040; m[19][ 7] = 0.1872470;
    m[ 0][ 8] = 0.3169540; m[ 1][ 8] = 2.1371500; m[ 2][ 8] = 3.9562900; m[ 3][ 8] = 0.9306760; m[ 4][ 8] = 0.2489720;
    m[ 5][ 8] = 4.2941100; m[ 6][ 8] = 0.5700250; m[ 7][ 8] = 0.2494100; m[ 8][ 8] = 0.0000000; m[ 9][ 8] = 0.1381900;
    m[10][ 8] = 0.4994620; m[11][ 8] = 0.8904320; m[12][ 8] = 0.4041410; m[13][ 8] = 0.6793710; m[14][ 8] = 0.6961980;
    m[15][ 8] = 0.7401690; m[16][ 8] = 0.4733070; m[17][ 8] = 0.2625690; m[18][ 8] = 3.8734400; m[19][ 8] = 0.1183580;
    m[ 0][ 9] = 0.1933350; m[ 1][ 9] = 0.1869790; m[ 2][ 9] = 0.5542360; m[ 3][ 9] = 0.0394370; m[ 4][ 9] = 0.1701350;
    m[ 5][ 9] = 0.1139170; m[ 6][ 9] = 0.1273950; m[ 7][ 9] = 0.0304501; m[ 8][ 9] = 0.1381900; m[ 9][ 9] = 0.0000000;
    m[10][ 9] = 3.1709700; m[11][ 9] = 0.3238320; m[12][ 9] = 4.2574600; m[13][ 9] = 1.0594700; m[14][ 9] = 0.0999288;
    m[15][ 9] = 0.3194400; m[16][ 9] = 1.4581600; m[17][ 9] = 0.2124830; m[18][ 9] = 0.4201700; m[19][ 9] = 7.8213000;
    m[ 0][10] = 0.3979150; m[ 1][10] = 0.4976710; m[ 2][10] = 0.1315280; m[ 3][10] = 0.0848047; m[ 4][10] = 0.3842870;
    m[ 5][10] = 0.8694890; m[ 6][10] = 0.1542630; m[ 7][10] = 0.0613037; m[ 8][10] = 0.4994620; m[ 9][10] = 3.1709700;
    m[10][10] = 0.0000000; m[11][10] = 0.2575550; m[12][10] = 4.8540200; m[13][10] = 2.1151700; m[14][10] = 0.4158440;
    m[15][10] = 0.3447390; m[16][10] = 0.3266220; m[17][10] = 0.6653090; m[18][10] = 0.3986180; m[19][10] = 1.8003400;
    m[ 0][11] = 0.9062650; m[ 1][11] = 5.3514200; m[ 2][11] = 3.0120100; m[ 3][11] = 0.4798550; m[ 4][11] = 0.0740339;
    m[ 5][11] = 3.8949000; m[ 6][11] = 2.5844300; m[ 7][11] = 0.3735580; m[ 8][11] = 0.8904320; m[ 9][11] = 0.3238320;
    m[10][11] = 0.2575550; m[11][11] = 0.0000000; m[12][11] = 0.9342760; m[13][11] = 0.0888360; m[14][11] = 0.5568960;
    m[15][11] = 0.9671300; m[16][11] = 1.3869800; m[17][11] = 0.1375050; m[18][11] = 0.1332640; m[19][11] = 0.3054340;
    m[ 0][12] = 0.8934960; m[ 1][12] = 0.6831620; m[ 2][12] = 0.1982210; m[ 3][12] = 0.1037540; m[ 4][12] = 0.3904820;
    m[ 5][12] = 1.5452600; m[ 6][12] = 0.3151240; m[ 7][12] = 0.1741000; m[ 8][12] = 0.4041410; m[ 9][12] = 4.2574600;
    m[10][12] = 4.8540200; m[11][12] = 0.9342760; m[12][12] = 0.0000000; m[13][12] = 1.1906300; m[14][12] = 0.1713290;
    m[15][12] = 0.4939050; m[16][12] = 1.5161200; m[17][12] = 0.5157060; m[18][12] = 0.4284370; m[19][12] = 2.0584500;
    m[ 0][13] = 0.2104940; m[ 1][13] = 0.1027110; m[ 2][13] = 0.0961621; m[ 3][13] = 0.0467304; m[ 4][13] = 0.3980200;
    m[ 5][13] = 0.0999208; m[ 6][13] = 0.0811339; m[ 7][13] = 0.0499310; m[ 8][13] = 0.6793710; m[ 9][13] = 1.0594700;
    m[10][13] = 2.1151700; m[11][13] = 0.0888360; m[12][13] = 1.1906300; m[13][13] = 0.0000000; m[14][13] = 0.1614440;
    m[15][13] = 0.5459310; m[16][13] = 0.1719030; m[17][13] = 1.5296400; m[18][13] = 6.4542800; m[19][13] = 0.6498920;
    m[ 0][14] = 1.4385500; m[ 1][14] = 0.6794890; m[ 2][14] = 0.1950810; m[ 3][14] = 0.4239840; m[ 4][14] = 0.1094040;
    m[ 5][14] = 0.9333720; m[ 6][14] = 0.6823550; m[ 7][14] = 0.2435700; m[ 8][14] = 0.6961980; m[ 9][14] = 0.0999288;
    m[10][14] = 0.4158440; m[11][14] = 0.5568960; m[12][14] = 0.1713290; m[13][14] = 0.1614440; m[14][14] = 0.0000000;
    m[15][14] = 1.6132800; m[16][14] = 0.7953840; m[17][14] = 0.1394050; m[18][14] = 0.2160460; m[19][14] = 0.3148870;
    m[ 0][15] = 3.3707900; m[ 1][15] = 1.2241900; m[ 2][15] = 3.9742300; m[ 3][15] = 1.0717600; m[ 4][15] = 1.4076600;
    m[ 5][15] = 1.0288700; m[ 6][15] = 0.7049390; m[ 7][15] = 1.3418200; m[ 8][15] = 0.7401690; m[ 9][15] = 0.3194400;
    m[10][15] = 0.3447390; m[11][15] = 0.9671300; m[12][15] = 0.4939050; m[13][15] = 0.5459310; m[14][15] = 1.6132800;
    m[15][15] = 0.0000000; m[16][15] = 4.3780200; m[17][15] = 0.5237420; m[18][15] = 0.7869930; m[19][15] = 0.2327390;
    m[ 0][16] = 2.1211100; m[ 1][16] = 0.5544130; m[ 2][16] = 2.0300600; m[ 3][16] = 0.3748660; m[ 4][16] = 0.5129840;
    m[ 5][16] = 0.8579280; m[ 6][16] = 0.8227650; m[ 7][16] = 0.2258330; m[ 8][16] = 0.4733070; m[ 9][16] = 1.4581600;
    m[10][16] = 0.3266220; m[11][16] = 1.3869800; m[12][16] = 1.5161200; m[13][16] = 0.1719030; m[14][16] = 0.7953840;
    m[15][16] = 4.3780200; m[16][16] = 0.0000000; m[17][16] = 0.1108640; m[18][16] = 0.2911480; m[19][16] = 1.3882300;
    m[ 0][17] = 0.1131330; m[ 1][17] = 1.1639200; m[ 2][17] = 0.0719167; m[ 3][17] = 0.1297670; m[ 4][17] = 0.7170700;
    m[ 5][17] = 0.2157370; m[ 6][17] = 0.1565570; m[ 7][17] = 0.3369830; m[ 8][17] = 0.2625690; m[ 9][17] = 0.2124830;
    m[10][17] = 0.6653090; m[11][17] = 0.1375050; m[12][17] = 0.5157060; m[13][17] = 1.5296400; m[14][17] = 0.1394050;
    m[15][17] = 0.5237420; m[16][17] = 0.1108640; m[17][17] = 0.0000000; m[18][17] = 2.4853900; m[19][17] = 0.3653690;
    m[ 0][18] = 0.2407350; m[ 1][18] = 0.3815330; m[ 2][18] = 1.0860000; m[ 3][18] = 0.3257110; m[ 4][18] = 0.5438330;
    m[ 5][18] = 0.2277100; m[ 6][18] = 0.1963030; m[ 7][18] = 0.1036040; m[ 8][18] = 3.8734400; m[ 9][18] = 0.4201700;
    m[10][18] = 0.3986180; m[11][18] = 0.1332640; m[12][18] = 0.4284370; m[13][18] = 6.4542800; m[14][18] = 0.2160460;
    m[15][18] = 0.7869930; m[16][18] = 0.2911480; m[17][18] = 2.4853900; m[18][18] = 0.0000000; m[19][18] = 0.3147300;
    m[ 0][19] = 2.0060100; m[ 1][19] = 0.2518490; m[ 2][19] = 0.1962460; m[ 3][19] = 0.1523350; m[ 4][19] = 1.0021400;
    m[ 5][19] = 0.3012810; m[ 6][19] = 0.5887310; m[ 7][19] = 0.1872470; m[ 8][19] = 0.1183580; m[ 9][19] = 7.8213000;
    m[10][19] = 1.8003400; m[11][19] = 0.3054340; m[12][19] = 2.0584500; m[13][19] = 0.6498920; m[14][19] = 0.3148870;
    m[15][19] = 0.2327390; m[16][19] = 1.3882300; m[17][19] = 0.3653690; m[18][19] = 0.3147300; m[19][19] = 0.0000000;

    return m;
}

vector<double> wag_empirical_freqs()
{
    std::vector<double> empirical_freqs(20);

    empirical_freqs[ 0] = 0.0866279;
    empirical_freqs[ 1] = 0.0439720;
    empirical_freqs[ 2] = 0.0390894;
    empirical_freqs[ 3] = 0.0570451;
    empirical_freqs[ 4] = 0.0193078;
    empirical_freqs[ 5] = 0.0367281;
    empirical_freqs[ 6] = 0.0580589;
    empirical_freqs[ 7] = 0.0832518;
    empirical_freqs[ 8] = 0.0244313;
    empirical_freqs[ 9] = 0.048466;
    empirical_freqs[10] = 0.086209;
    empirical_freqs[11] = 0.0620286;
    empirical_freqs[12] = 0.0195027;
    empirical_freqs[13] = 0.0384319;
    empirical_freqs[14] = 0.0457631;
    empirical_freqs[15] = 0.0695179;
    empirical_freqs[16] = 0.0610127;
    empirical_freqs[17] = 0.0143859;
    empirical_freqs[18] = 0.0352742;
    empirical_freqs[19] = 0.0708956;

    return empirical_freqs;
}

typedef RevBayesCore::ConcreteTimeReversibleRateMatrix CGTR;

CGTR* WAG(const RevBayesCore::Simplex& pi)
{
    if (pi.size() != 20)
        throw RbException()<<"fnWAG: need exactly 20 amino acid frequencies, but got "<<pi.size();

    return new CGTR( wag_exchange_rates(), pi);
}



/** default constructor */
Func_wag::Func_wag( void ) : TypedFunction<RateMatrix>( )
{
    bf_empirical = new Simplex( RevBayesCore::Simplex( wag_empirical_freqs() ));
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_wag* Func_wag::clone( void ) const
{
    
    return new Func_wag( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::RateGenerator >* Func_wag::createFunction( void ) const
{
    const RevObject& obj = args[0].getVariable()->getRevObject();

    RevBayesCore::TypedDagNode<RevBayesCore::Simplex>* bf = static_cast<const Simplex &>( obj ).getDagNode();

    return RevBayesCore::generic_function_ptr< RevBayesCore::RateGenerator >( WAG, bf );
}


/* Get argument rules */
const ArgumentRules& Func_wag::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "aaFrequencies", Simplex::getClassTypeSpec(), "The stationary frequencies of the states.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, bf_empirical ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_wag::getClassType(void)
{
    
    static std::string rev_type = "Func_wag";
    
    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_wag::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_wag::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnWAG";
    
    return f_name;
}


const TypeSpec& Func_wag::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
