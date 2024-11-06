%{
/**
 * @file
 * Grammar specification in bison format for RevBayes, a computing environment 
 * for evolutionary analysis, particularly Bayesian phylogenetic inference. The
 * language used by RevBayes is referred to as Rev.
 *
 * The grammar borrows heavily from the R grammar specification in the gram.y 
 * file of the R source code but deviates significantly in many respects, being
 * more similar to object-oriented languages like C++ or Java. The model description
 * syntax is inspired by the language used originally by BUGS to describe complex 
 * stochastic models. Unlike BUGS and similar programs, such as JAGS, RevBayes
 * allows models to be built in an interpreted (interactive) environment.
 *
 * @brief Grammar specification in bison format
 *
 * @author Fredrik Ronquist and Sebastian Hoehna
 */

/* The following statements go into the resulting C code */

#include "Environment.h"
#include "Integer.h"
#include "Natural.h"
#include "Probability.h"
#include "Parser.h"
#include "RlBoolean.h"
#include "RbException.h"
#include "RlString.h"
#include "Real.h"
#include "RealPos.h"
#include "SyntaxElement.h"
#include "SyntaxAdditionAssignment.h"
#include "SyntaxBinaryExpr.h"
#include "SyntaxClassDef.h"
#include "SyntaxConstant.h"
#include "SyntaxConstantAssignment.h"
#include "SyntaxDecrement.h"
#include "SyntaxDeterministicAssignment.h"
#include "SyntaxDivisionAssignment.h"
#include "SyntaxForLoop.h"
#include "SyntaxFormal.h"
#include "SyntaxFunctionCall.h"
#include "SyntaxFunctionDef.h"
#include "SyntaxIncrement.h"
#include "SyntaxIndexOperation.h"
#include "SyntaxLabeledExpr.h"
#include "SyntaxMultiplicationAssignment.h"
#include "SyntaxReferenceAssignment.h"
#include "SyntaxStatement.h"
#include "SyntaxStochasticAssignment.h"
#include "SyntaxSubtractionAssignment.h"
#include "SyntaxUnaryExpr.h"
#include "SyntaxVariableDecl.h"
#include "SyntaxVariable.h"
#include "SyntaxPipePlaceholder.h"
#include "SyntaxWorkspaceVariableAssignment.h"
#include "Workspace.h"

#include <iostream>
#include <cstdint>
#include <list>
#include <sstream>
#include <string>

using namespace RevLanguage;

extern int yylex(void);
extern char *yytext;
extern Environment *executionEnvironment;

/* The function yyerror handles errors. It is defined below. */
int yyerror(const char *);
RevLanguage::SyntaxElement* xxpipe(RevLanguage::SyntaxElement* arg, RevLanguage::SyntaxElement* fxnCallE);

/* We use the global parser to execute the syntax tree */
Parser& parser = Parser::getParser();


#define YY_NEVER_INTERACTIVE
%}

/* Bison needs a union of the types handled by the parser */
%union {
    char*                                           c_string;
    std::string*                                    string;
    double                                          realValue;
    std::int64_t                                             longIntValue;
    bool                                            boolValue;
    RevLanguage::SyntaxElement*                     syntaxElement;
    RevLanguage::SyntaxVariable*                    syntaxVariable;
    RevLanguage::SyntaxFunctionCall*                syntaxFunctionCall;
    RevLanguage::SyntaxLabeledExpr*                 syntaxLabeledExpr;
    RevLanguage::SyntaxFormal*                      syntaxFormal;
    std::list<RevLanguage::SyntaxElement*>*         syntaxElementList;
    std::list<RevLanguage::SyntaxLabeledExpr*>*     argumentList;
    std::list<RevLanguage::SyntaxFormal*>*          formalList;
};

/* Return types of the elements handled by the parser */
%type <c_string> NAME STRING
%type <realValue> REAL RBINF
%type <longIntValue> INT RBNULL
%type <boolValue> FALSE TRUE
%type <string> RBTAB
%type <string> identifier typeSpec optDims dimList 
%type <syntaxElement> variable
%type <syntaxFunctionCall> functionCall fxnCall
%type <syntaxLabeledExpr> argument
%type <syntaxFormal> formal
%type <syntaxElement> constant
%type <syntaxElement> statement expression stmt_or_expr
%type <syntaxElement> arrowAssign tildeAssign equationAssign workspaceAssign
%type <syntaxElement> referenceAssign
%type <syntaxElement> additionAssign subtractionAssign multiplicationAssign divisionAssign
%type <syntaxElement> declaration classDef memberDef
%type <syntaxElement> functionDef procedureDef
%type <syntaxElement> forStatement ifStatement whileStatement
%type <syntaxElement> forCond cond returnStatement
%type <syntaxElement> nextStatement breakStatement
%type <syntaxElementList> elementList optElements
%type <syntaxElementList> stmts stmtList
%type <syntaxElementList> memberDefs
%type <argumentList> argumentList optArguments vectorList vector
%type <formalList> formalList optFormals

/* Tokens returned by the lexer and handled by the parser */
%token REAL INT NAME STRING RBNULL RBTAB FALSE TRUE RBINF
%token FUNCTION PROCEDURE CLASS FOR IN IF ELSE WHILE NEXT BREAK RETURN
%token MOD_CONST MOD_DYNAMIC MOD_STOCHASTIC MOD_DETERMINISTIC PROTECTED
%token ARROW_ASSIGN TILDE_ASSIGN EQUATION_ASSIGN WORKSPACE_ASSIGN
%token REFERENCE_ASSIGN
%token ADDITION_ASSIGN SUBTRACTION_ASSIGN MULTIPLICATION_ASSIGN DIVISION_ASSIGN 
%token DECREMENT INCREMENT
%token EQUAL 
%token AND OR AND2 OR2 GT GE LT LE EQ NE
%token PIPE
%token PIPE_PLACEHOLDER
%token END_OF_INPUT

/* Destructors */
%destructor { for (std::list<SyntaxElement*>::iterator it=$$->begin(); it != $$->end(); it++) { SyntaxElement* theElement = *it; delete theElement; }; delete ($$); } elementList optElements stmts stmtList memberDefs 
%destructor { for (std::list<SyntaxLabeledExpr*>::iterator it=$$->begin(); it != $$->end(); it++) { SyntaxLabeledExpr* theElement = *it; delete theElement; }; delete ($$); } argumentList optArguments vectorList vector
%destructor { for (std::list<SyntaxFormal*>::iterator it=$$->begin(); it != $$->end(); it++) { SyntaxFormal* theElement = *it; delete theElement; }; delete ($$); } formalList optFormals
%destructor { delete ($$); } identifier typeSpec optDims dimList 
%destructor { delete ($$); } variable functionCall fxnCall argument formal constant
%destructor { delete ($$); } statement expression stmt_or_expr 
%destructor { delete ($$); } arrowAssign tildeAssign equationAssign workspaceAssign
%destructor { delete ($$); } referenceAssign
%destructor { delete ($$); } additionAssign subtractionAssign multiplicationAssign divisionAssign
%destructor { delete ($$); } declaration classDef memberDef 
%destructor { delete ($$); } functionDef procedureDef 
%destructor { delete ($$); } forStatement ifStatement whileStatement 
%destructor { delete ($$); } forCond cond returnStatement
%destructor { delete ($$); } nextStatement breakStatement

/* Use location tracking */
%locations

/*
 * Precedence table, low to high, with order of evaluation.
 */
%left       '?'
%left       WHILE FOR
%right      IF
%left       ELSE
%right      ARROW_ASSIGN
%right      TILDE_ASSIGN
%right      EQUATION_ASSIGN
%right      REFERENCE_ASSIGN
%right      ADDITION_ASSIGN
%right      SUBTRACTION_ASSIGN
%right      MULTIPLICATION_ASSIGN
%right      DIVISION_ASSIGN
%right      EQUAL
%left       OR OR2
%left       AND AND2
%left       UNOT
%nonassoc   GT GE LT LE EQ NE
%left       '+' '-'
%left       '*' '/'
%left       PIPE
%left       ':' '%'
%right      DECREMENT INCREMENT
%left       UMINUS UPLUS
%right      '^'
%left       '.' UAND
%nonassoc   '(' '['

%%

/*
 * The grammar is based both on R and on more modern object-
 * oriented programming languages like C++ and Java. Like R,
 * it is a language intended for interactive use. It does
 * not support implicit looping or vector indexing like R
 * does. This means that explicit loops will have to be used
 * in more situations than in R. On the other hand, this
 * makes the language less ambiguous and more intuitive for
 * users familiar with Java or C++. All language objects
 * inherit from "RevObject". RevBayes allows types of variables
 * and function arguments, as well as return types of func-
 * tions, to be specified. However, types default to "RevObject",
 * allowing apparently untyped expressions to be evaluated.
 * Global variables can be declared but are more typically
 * created by assigning to them using various assignment
 * operators.
 *
 * RevBayes allows a natural description of graphical represen-
 * tations of complex probability models. The most essential
 * components are tilde assignment ('~') for creating stochastic
 * nodes and equation assignment (':=') for creating deterministic
 * nodes in model DAGs. Constant nodes are created using
 * arrow assignment ('<-'), as in R.
 *
 * The Rev language distinguishes between three ways of passing
 * arguments: constant, dynamic and mutable. The "constant" method
 * is effectively equivalent to passing by value. The "dynamic"
 * method is different in dynamic expressions, occurring on the
 * right-hand side of deterministic or stochastic assignments.
 * In those contexts, functions receive dynamically update versions
 * of the variables instead of constant values. This includes
 * the constructor functions of distributions. Finally, the
 * "mutable" method is similar to "dynamic", but allows the
 * function to modify the object passed in. Functions with
 * mutable arguments can only be used in the workspace, they
 * are not allowed in model DAGs because they potentially
 * violate the logic of DAGs, creating implicit loops.
 * 
 * Unlike R, right assignment ('->') is not supported. Also,
 * a number of constructs that are interpreted as expressions
 * in R are treated as statements or declarations in RevBayes
 * and have no return value (or the return value is discarded).
 * These constructs include control statements and function def-
 * initions, among others.
 *
 * For more discussion of the language, see the RevBayes
 * documentation at http://revbayes.net.
 */

prog    :       END_OF_INPUT
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Bison encountered end_of_input; ignored\n");
#endif
                    return 0;
                }
        |       '\n'
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Bison encountered newline; ignored\n");
#endif
                    return 0;
                }
        |       stmt_or_expr '\n'
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Bison trying to execute statement or expression\n");
#endif
                    int rv = parser.execute($1, *executionEnvironment);
                    delete $1;
                    return rv;
                }
        |       stmt_or_expr ';'
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Bison trying to execute statement or expression\n");
#endif
                    int rv =  parser.execute($1, *executionEnvironment);
                    delete $1;
                    return rv;
                }
        |       declaration '\n'
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Bison trying to execute declaration\n");
#endif
                    int rv =  parser.execute($1, *executionEnvironment);
                    delete $1;
                    return rv;
                }
        |       declaration ';'
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Bison trying to execute declaration\n");
#endif
                    int rv =  parser.execute($1, *executionEnvironment);
                    delete $1;
                    return rv;
                }
        |       '?' identifier '\n'
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Bison trying to get help for symbol\n");
#endif
                    int rv =  parser.help(*$2);
                    delete $2;
                    return rv;
                }
        |       '?' identifier ';'
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Bison trying to get help for symbol\n");
#endif
                    int rv =  parser.help(*$2);
                    delete $2;
                    return rv;
                }
        |       '?' identifier '.' identifier '\n'
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Bison trying to get help for symbol\n");
#endif
                    int rv =  parser.help(*$2, *$4);
                    delete $2;
                    delete $4;
                    return rv;
                }
        |       '?' identifier '.' identifier ';'
                {
#ifdef DEBUG_BISON_FLEX
                printf("Bison trying to get help for symbol\n");
#endif
                int rv =  parser.help(*$2, *$4);
                delete $2;
                delete $4;
                return rv;
                }
        |       '?' functionCall '\n'
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Bison trying to get help for function call\n");
#endif
                    int rv =  parser.help($2);
                    delete $2;
                    return rv;
                }
        |       '?' functionCall ';'
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Bison trying to get help for function call\n");
#endif
                    int rv =  parser.help($2);
                    delete $2;
                    return rv;
                }
        |       error '\n'
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Bison error when reading line %d position %d to line %d position %d\n", @$.first_line, @$.first_column, @$.last_line, @$.last_column);
#endif
                    YYABORT;
                }
        |       error ';'
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Bison error when reading line %d position %d to line %d position %d\n", @$.first_line, @$.first_column, @$.last_line, @$.last_column);
#endif
                    YYABORT;
                }
        ;

expression  :   constant                    { $$ = $1; }
            |   PIPE_PLACEHOLDER            { $$ = new SyntaxPipePlaceholder(); }

            |   vector                      { $$ = new SyntaxFunctionCall("v", $1); }

            |   '(' expression ')'          { $$ = $2; }

            |   '-' expression %prec UMINUS { $$ = new SyntaxUnaryExpr(SyntaxUnaryExpr::UMinus, $2); }
            |   '+' expression %prec UPLUS  { $$ = new SyntaxUnaryExpr(SyntaxUnaryExpr::UPlus, $2); }
            |   '!' expression %prec UNOT   { $$ = new SyntaxUnaryExpr(SyntaxUnaryExpr::UNot, $2); }
            |   AND expression %prec UAND   { $$ = new SyntaxUnaryExpr(SyntaxUnaryExpr::UAnd, $2); }

            |   DECREMENT variable          { $$ = new SyntaxDecrement( $2, false ); }
            |   variable DECREMENT          { $$ = new SyntaxDecrement( $1, true ); }
            |   INCREMENT variable          { $$ = new SyntaxIncrement( $2, false ); }
            |   variable INCREMENT          { $$ = new SyntaxIncrement( $1, true ); }

            |   expression ':' expression   { $$ = new SyntaxBinaryExpr(SyntaxBinaryExpr::Range, $1, $3); }

            |   expression '+' expression   { $$ = new SyntaxBinaryExpr(SyntaxBinaryExpr::Add, $1, $3); }
            |   expression '-' expression   { $$ = new SyntaxBinaryExpr(SyntaxBinaryExpr::Sub, $1, $3); }
            |   expression '*' expression   { $$ = new SyntaxBinaryExpr(SyntaxBinaryExpr::Mul, $1, $3); }
            |   expression '/' expression   { $$ = new SyntaxBinaryExpr(SyntaxBinaryExpr::Div, $1, $3); }
            |   expression '^' expression   { $$ = new SyntaxBinaryExpr(SyntaxBinaryExpr::Exp, $1, $3); }
            |   expression '%' expression   { $$ = new SyntaxBinaryExpr(SyntaxBinaryExpr::Mod, $1, $3); }

            |   expression LT expression    { $$ = new SyntaxBinaryExpr(SyntaxBinaryExpr::Lt, $1, $3); }
            |   expression LE expression    { $$ = new SyntaxBinaryExpr(SyntaxBinaryExpr::Le, $1, $3); }
            |   expression EQ expression    { $$ = new SyntaxBinaryExpr(SyntaxBinaryExpr::Eq, $1, $3); }
            |   expression NE expression    { $$ = new SyntaxBinaryExpr(SyntaxBinaryExpr::Ne, $1, $3); }
            |   expression GE expression    { $$ = new SyntaxBinaryExpr(SyntaxBinaryExpr::Ge, $1, $3); }
            |   expression GT expression    { $$ = new SyntaxBinaryExpr(SyntaxBinaryExpr::Gt, $1, $3); }

            |   expression AND expression   { $$ = new SyntaxBinaryExpr(SyntaxBinaryExpr::And, $1, $3); }
            |   expression OR expression    { $$ = new SyntaxBinaryExpr(SyntaxBinaryExpr::Or, $1, $3); }
            |   expression AND2 expression  { $$ = new SyntaxBinaryExpr(SyntaxBinaryExpr::And2, $1, $3); }
            |   expression OR2 expression   { $$ = new SyntaxBinaryExpr(SyntaxBinaryExpr::Or2, $1, $3); }

            |   expression PIPE expression  { $$ = xxpipe($1, $3); }

            |   arrowAssign                 { $$ = $1; }
            |   equationAssign              { $$ = $1; }
            |   tildeAssign                 { $$ = $1; }
            |   workspaceAssign             { $$ = $1; }
            |   referenceAssign             { $$ = $1; }

            |   additionAssign              { $$ = $1; }
            |   subtractionAssign           { $$ = $1; }
            |   multiplicationAssign        { $$ = $1; }
            |   divisionAssign              { $$ = $1; }

            |   functionCall                { $$ = $1; }

            |   variable                    { $$ = $1; }

            ;

arrowAssign     :   expression ARROW_ASSIGN expression
                    { 
#ifdef DEBUG_BISON_FLEX
                        printf("Parser inserting arrow assignment (ARROW_ASSIGN) in syntax tree\n");
#endif
                        $$ = new SyntaxConstantAssignment($1, $3);
                    }
                ;

tildeAssign     :   expression TILDE_ASSIGN expression
                    {
#ifdef DEBUG_BISON_FLEX
                        printf("Parser inserting tilde assignment (TILDE_ASSIGN) in syntax tree\n");
#endif
                        $$ = new SyntaxStochasticAssignment($1, $3);
                    }
                ;

equationAssign  :   expression EQUATION_ASSIGN expression
                    {
#ifdef DEBUG_BISON_FLEX
                        printf("Parser inserting equation assignment (EQUATION_ASSIGN) in syntax tree\n");
#endif
                        $$ = new SyntaxDeterministicAssignment($1, $3); 
                    }
                ;

workspaceAssign :   expression EQUAL expression
                    {
#ifdef DEBUG_BISON_FLEX
                        printf("Parser inserting workspace assignment (WORKSPACE_ASSIGN) in syntax tree\n");
#endif
                        $$ = new SyntaxWorkspaceVariableAssignment($1, $3);
                    }
                ;

referenceAssign :   expression REFERENCE_ASSIGN expression
                    {
#ifdef DEBUG_BISON_FLEX
                        printf("Parser inserting reference assignment (REFERENCE_ASSIGN) in syntax tree\n");
#endif
                        $$ = new SyntaxReferenceAssignment($1, $3);
                    }
                ;

additionAssign  :   expression ADDITION_ASSIGN expression
                    {
#ifdef DEBUG_BISON_FLEX
                        printf("Parser inserting addition assignment (ADDITION_ASSIGN) in syntax tree\n");
#endif
                        $$ = new SyntaxAdditionAssignment($1, $3); 
                    }
                ;

subtractionAssign   :   expression SUBTRACTION_ASSIGN expression
                        {
#ifdef DEBUG_BISON_FLEX
                            printf("Parser inserting subtraction assignment (SUBTRACTION_ASSIGN) in syntax tree\n");
#endif
                            $$ = new SyntaxSubtractionAssignment($1, $3); 
                        }
                    ;

multiplicationAssign    :   expression MULTIPLICATION_ASSIGN expression
                            {
#ifdef DEBUG_BISON_FLEX
                                printf("Parser inserting multiplication assignment (MULTIPLICATION_ASSIGN) in syntax tree\n");
#endif
                                $$ = new SyntaxMultiplicationAssignment($1, $3); 
                            }
                        ;

divisionAssign  :   expression DIVISION_ASSIGN expression
                    {
#ifdef DEBUG_BISON_FLEX
                        printf("Parser inserting division assignment (DIVISION_ASSIGN) in syntax tree\n");
#endif
                        $$ = new SyntaxDivisionAssignment($1, $3); 
                    }
                ;

variable    :   identifier optElements
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Parser inserting variable (NAMED_VAR) in syntax tree\n");
#endif
                    $$ = new SyntaxVariable(*$1);
                    for (std::list<SyntaxElement*>::iterator it=$2->begin(); it!=$2->end(); ++it)
                    {
                        $$ = new SyntaxIndexOperation($$,*it);
                    }
                    delete $1;
                    delete $2;
                }
            |   fxnCall elementList
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Parser inserting variable (FUNCTION_VAR) in syntax tree\n");
#endif
                    $$ = $1;
                    for (auto& element: *$2)
                    {
                        $$ = new SyntaxIndexOperation($$, element);
                    }
                    delete $2;
                }
            |   '(' expression ')' elementList
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Parser inserting variable (EXPRESSION_VAR) in syntax tree\n");
#endif
                    $$ = $2;
                    for (auto& element: *$4)
                    {
                        $$ = new SyntaxIndexOperation($$, element);
                    }
                    delete $4;
                }
            |   variable '.' fxnCall elementList
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Parser inserting member variable (FUNCTION_VAR) in syntax tree\n");
#endif
                    $3->setBaseVariable($1);
                    $$ = $3;
                    for (auto& element: *$4)
                    {
                        $$ = new SyntaxIndexOperation($$, element);
                    }
                    delete $4;
                }
            |   PIPE_PLACEHOLDER
                {
                    $$ = new SyntaxPipePlaceholder;
                }
            ;

optElements :   /* empty */                     { $$ = new std::list<SyntaxElement*>(); }
            |   elementList                     { $$ = $1; }
            ;

elementList :   '[' expression ']'              { $$ = new std::list<SyntaxElement*>(1, $2); }
            |   '[' ']'                         { $$ = new std::list<SyntaxElement*>(); }
            |   elementList '[' expression ']'  { $1->push_back($3); $$ = $1; }
            |   elementList '[' ']'             { $1->push_back( NULL ); $$ = $1; }
            ;

fxnCall     :   identifier '(' optArguments ')'
                {
                    $$ = new SyntaxFunctionCall(*$1, $3);
                    delete $1;
                }
            ;

functionCall    :   fxnCall
                    {
#ifdef DEBUG_BISON_FLEX
                        printf("Parser inserting function call in syntax tree\n");
#endif
                        $$ = $1;
                    }
                |   variable '.' fxnCall
                    {
#ifdef DEBUG_BISON_FLEX
                        printf("Parser inserting member call in syntax tree\n");
#endif
                        $3->setBaseVariable($1);
                        $$ = $3;
                    }
                |   functionCall '.' fxnCall
                    {
#ifdef DEBUG_BISON_FLEX
                        printf("Parser inserting member call in syntax tree\n");
#endif
                        $3->setBaseVariable($1);
                        $$ = $3;
                    }
                ;


optArguments    :   /* empty */             { $$ = new std::list<SyntaxLabeledExpr*>(); }
                |   argumentList            { $$ = $1; }
                ;

argumentList    :   argument                    { $$ = new std::list<SyntaxLabeledExpr*>(1,$1); }
                |   argumentList ',' argument   { $1->push_back($3); $$ = $1; }
                ;

argument   :    expression
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Parser inserting unlabeled argument in syntax tree\n");
#endif
                    $$ = new SyntaxLabeledExpr( "", $1);
                }
            |   identifier EQUAL expression
                { 
#ifdef DEBUG_BISON_FLEX
                    printf("Parser inserting labeled argument in syntax tree\n");
#endif
                    $$ = new SyntaxLabeledExpr(*$1, $3);
                    delete $1;
                }
            ;

functionDef :   FUNCTION identifier '(' optFormals ')' stmts
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Parser inserting function definition in syntax tree\n");
#endif
                    $$ = new SyntaxFunctionDef("", *$2, $4, $6);
                    delete $2;
                }
 
            |   FUNCTION identifier optDims identifier '(' optFormals ')' stmts
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Parser inserting typed function definition in syntax tree\n");
#endif
                    $2->append(*$3);
                    $$ = new SyntaxFunctionDef(*$2, *$4, $6, $8);
                    delete $2;
                    delete $3;
                    delete $4;
                }
            ;

procedureDef :   PROCEDURE identifier '(' optFormals ')' stmts
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Parser inserting procedure definition in syntax tree\n");
#endif
                    $$ = new SyntaxFunctionDef("", *$2, $4, $6, true);
                    delete $2;
                }
 
            |   PROCEDURE identifier optDims identifier '(' optFormals ')' stmts
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Parser inserting typed procedure definition in syntax tree\n");
#endif
                    $2->append(*$3);
                    $$ = new SyntaxFunctionDef(*$2, *$4, $6, $8, true);
                    delete $2;
                    delete $3;
                    delete $4;
                }
            ;

optFormals  :   /* empty */             { $$ = new std::list<SyntaxFormal*>(); }
            |   formalList              { $$ = $1; }
            ;

formalList  :   formal                  { $$ = new std::list<SyntaxFormal*>(1, $1); }
            |   formalList ',' formal   { $1->push_back($3); $$ = $1; }
            ;

formal      :   identifier
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Inserting labeled formal argument without default in syntax tree\n");
#endif
                    $$ = new SyntaxFormal(*$1, NULL );
                    delete $1;
                }
            |   identifier EQUAL expression
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Inserting labeled formal argument with default in syntax tree\n");
#endif
                    $$ = new SyntaxFormal(*$1, $3);
                    delete $1;
                }
            |   typeSpec identifier
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Inserting typed labeled formal argument without default in syntax tree\n");
#endif
                    $$ = new SyntaxFormal(*$1, *$2, NULL);
                    delete $1;
                    delete $2;
                }
            |   typeSpec identifier EQUAL expression
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Inserting typed labeled formal argument with default in syntax tree\n");
#endif
                    $$ = new SyntaxFormal(*$1, *$2, $4);
                    delete $1;
                    delete $2;
                }
            ;

typeSpec    :   identifier optDims                      { $1->append(*($2)); delete $2; $$ = $1; }
            |   MOD_CONST identifier optDims            { $2->append(*($3)); $2->insert(0, "const ");           delete $3; $$ = $2; }
            |   MOD_DYNAMIC identifier optDims          { $2->append(*($3)); $2->insert(0, "dynamic ");         delete $3; $$ = $2; }
            |   MOD_STOCHASTIC identifier optDims       { $2->append(*($3)); $2->insert(0, "stochastic ");      delete $3; $$ = $2; }
            |   MOD_DETERMINISTIC identifier optDims    { $2->append(*($3)); $2->insert(0, "deterministic ");   delete $3; $$ = $2; }
            ;

optDims     :   /* empty */                 { $$ = new std::string(""); }
            |   dimList                     { $$ = $1; }
            ;

dimList     :   '[' ']'                     { $$ = new std::string("[]"); }
            |   dimList '[' ']'             { $1->append("[]"); $$ = $1; }
            ;

stmts       :   '{' stmtList '}'                { $$ = $2; }
            |   stmt_or_expr
                {
                    std::list<SyntaxElement*>* stmts = new std::list<SyntaxElement*>();
                    stmts->push_back($1);
                    $$ = stmts;
                }
            ;

stmtList    :   /* empty */                 { $$ = new std::list<SyntaxElement*>(); }
            |   stmt_or_expr                { $$ = new std::list<SyntaxElement*>(1, $1); }
            |   stmtList ';' stmt_or_expr   { $1->push_back($3); $$ = $1; }
            |   stmtList ';'                { $$ = $1; }
            |   stmtList '\n' stmt_or_expr  { $1->push_back($3); $$ = $1; }
            |   stmtList '\n'               { $$ = $1; }
            ;

statement   :   ifStatement             { $$ = $1; }
            |   forStatement            { $$ = $1; }
            |   whileStatement          { $$ = $1; }
            |   nextStatement           { $$ = $1; }
            |   breakStatement          { $$ = $1; }
            |   returnStatement         { $$ = $1; }
            ;

stmt_or_expr    :   statement           { $$ = $1; }
                |   expression          { $$ = $1; }
                ;

declaration     :   classDef            { $$ = $1; }
                |   functionDef         { $$ = $1; }
                |   procedureDef        { $$ = $1; }
                |   identifier optElements identifier
                    {
#ifdef DEBUG_BISON_FLEX
                        printf("Parser inserting variable declaration in syntax tree\n");
#endif
                        $$ = new SyntaxVariableDecl(*$1, $2, *$3);
                        delete $1;
                        delete $3;
                    }
                ;

memberDefs      :   /* empty */                 { $$ = new std::list<SyntaxElement*>(); }
                |   memberDef                   { $$ = new std::list<SyntaxElement*>(1, $1); }
                |   memberDefs ';' memberDef    { $1->push_back($3); $$ = $1; }
                |   memberDefs ';'              { $$ = $1; }
                |   memberDefs '\n' memberDef   { $1->push_back($3); $$ = $1; }
                |   memberDefs '\n'             { $$ = $1; }
                ;

memberDef   :   formal              { $$ = $1; }
            |   PROTECTED formal    { $2->setIsProtected(); $$ = $2; }
            |   functionDef         { $$ = $1; }
            |   procedureDef        { $$ = $1; }
            ;

classDef    :   CLASS identifier ':' identifier '{' memberDefs '}'
                {
#ifdef DEBUG_BISON_FLEX
                printf("Parser inserting class definition (CLASS_DEF) in syntax tree\n");
#endif
                    $$ = new SyntaxClassDef(*$2, *$4, $6);
                    delete $2;
                    delete $4;
                }
            ;

ifStatement :   IF cond stmts               { $$ = new SyntaxStatement(SyntaxStatement::If, $2, $3); }
            |   IF cond stmts ELSE stmts    { $$ = new SyntaxStatement(SyntaxStatement::IfElse, $2, $3, $5); }

cond    :   '(' expression ')'    { $$ = $2; }
        ;

forStatement    :   FOR forCond stmts   { $$ = new SyntaxStatement(SyntaxStatement::For, $2, $3); }
                ;

forCond     :   '(' identifier IN expression ')'    { $$ = new SyntaxForLoop(*$2, $4); delete $2; }
            ;

whileStatement  :   WHILE cond stmts    { $$ = new SyntaxStatement(SyntaxStatement::While, $2, $3); }
                ;

nextStatement   :   NEXT    { $$ = new SyntaxStatement(SyntaxStatement::Next); }
                ;

breakStatement  :   BREAK   { $$ = new SyntaxStatement(SyntaxStatement::Break); }
                ;

returnStatement :   RETURN              { $$ = new SyntaxStatement(SyntaxStatement::Return); }
                |   RETURN expression   { $$ = new SyntaxStatement(SyntaxStatement::Return, $2); }
                ;

identifier  :   NAME    { $$ = new std::string($1); }
            ;


vector      :   '[' vectorList ']'      { $$ = $2; }
            ;

vectorList  :   vectorList ',' expression   { $1->push_back(new SyntaxLabeledExpr( "", $3) ); $$ = $1; }
            |   expression
                {
                $$ = new std::list<SyntaxLabeledExpr*>(1, new SyntaxLabeledExpr("", $1) );
                }
            ;

constant    :   FALSE
                {
                    $$ = new SyntaxConstant(new RlBoolean(false) );
                }
            |   TRUE
                {
                    $$ = new SyntaxConstant(new RlBoolean(true) );
                }
            |   RBNULL
                {
                    $$ = new SyntaxConstant( NULL );
                }
            |   RBTAB
                {
                    $$ = new SyntaxConstant( new RlString("\t") );
                }
            |   RBINF
                {
                    $$ = new SyntaxConstant( new RealPos( RbConstants::Double::inf ) );
                }
            |   INT
                {
                    if ( $1 < 0 ) {
                        $$ = new SyntaxConstant(new Integer($1) );
                    }
                    else { 
                        $$ = new SyntaxConstant(new Natural($1) );
                    }
                }
            |   STRING
                {
                    $$ = new SyntaxConstant(new RlString($1) );
                }
            |   REAL
                {
                    /* This code records and preserves input format of the real */
                    /*
                    int prec;
                    bool sci;
                    int i=0;
                    while (yytext[i]!='\0') {
                    if (yytext[i]=='E' || yytext[i]=='e')
                        break;
                        i++;
                    }
                    if (yytext[i]!='\0') {
                        sci = true;
                        prec = i;
                        for (i=0; yytext[i]!='\0'; i++) {
                            if (yytext[i]=='.') {
                                prec--;
                            break;
                            }
                        }
                    }
                    else {
                        sci = false;
                        for (i=0; yytext[i]!='\0'; i++) {
                            if (yytext[i]=='.') {
                                break;
                            }
                        }
                        prec = (int)(strlen(yytext)) - 1 - i;
                    }
                    Real* real;
                    if ($1 >= 0.0)
                        real = new RealPos($1);
                    else
                        real = new Real($1);
                    real->setPrecision(prec);
                    real->setScientific(sci);
                    */
                    
                    if ($1 >= 0.0 && $1 <= 1.0) {
                        #ifdef DEBUG_BISON_FLEX
                        printf("Parser inserting RealPos constant in syntax tree \n");
                        #endif
                        $$ = new SyntaxConstant(new Probability($1) );
                    }
                    else if ($1 >= 0.0) {
#ifdef DEBUG_BISON_FLEX
                        printf("Parser inserting RealPos constant in syntax tree \n");
#endif
                        $$ = new SyntaxConstant(new RealPos($1) );
                    }
                    else {
#ifdef DEBUG_BISON_FLEX
                        printf("Parser inserting Real constant in syntax tree \n");
#endif
                        $$ = new SyntaxConstant(new Real($1) );
                    }
                }
            ;

%%


/* Definition of yyerror. */
int yyerror(const char *msg)
{
#ifdef DEBUG_BISON_FLEX
    printf("Bison code said: %s\n", msg);
#endif
    if ( foundNewline == true )
        foundErrorBeforeEnd = false;
    else
        foundErrorBeforeEnd = true;

    yylloc.first_column = yycolumn - yyleng;
    yylloc.last_column  = yycolumn - 1;

    return 1;
}

RevLanguage::SyntaxElement* xxpipe(RevLanguage::SyntaxElement* arg,
                                   RevLanguage::SyntaxElement* fxnCallE)
{
    if (auto fxnCall = dynamic_cast<RevLanguage::SyntaxFunctionCall*>(fxnCallE))
    {
        fxnCall->pipeAddArg(arg);
    }
    else
    {
        delete arg;
        delete fxnCallE;
        throw RbException("The pipe operator requires a function call as RHS.");
    }

    return fxnCallE;
}
