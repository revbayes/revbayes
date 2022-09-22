/* A Bison parser, made by GNU Bison 3.8.2.  */

/* Bison implementation for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2015, 2018-2021 Free Software Foundation,
   Inc.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* DO NOT RELY ON FEATURES THAT ARE NOT DOCUMENTED in the manual,
   especially those whose name start with YY_ or yy_.  They are
   private implementation details that can be changed or removed.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output, and Bison version.  */
#define YYBISON 30802

/* Bison version string.  */
#define YYBISON_VERSION "3.8.2"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1




/* First part of user prologue.  */
#line 1 "./grammar.y"

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
#include <list>
#include <sstream>
#include <string>

using namespace RevLanguage;

extern int yylex(void);
extern char *yytext;
extern Environment *executionEnvironment;

/* The function yyerror handles errors. It is defined below. */
int yyerror(const char *);

/* We use the global parser to execute the syntax tree */
Parser& parser = Parser::getParser();


#define YY_NEVER_INTERACTIVE

#line 150 "./grammar.tab.c"

# ifndef YY_CAST
#  ifdef __cplusplus
#   define YY_CAST(Type, Val) static_cast<Type> (Val)
#   define YY_REINTERPRET_CAST(Type, Val) reinterpret_cast<Type> (Val)
#  else
#   define YY_CAST(Type, Val) ((Type) (Val))
#   define YY_REINTERPRET_CAST(Type, Val) ((Type) (Val))
#  endif
# endif
# ifndef YY_NULLPTR
#  if defined __cplusplus
#   if 201103L <= __cplusplus
#    define YY_NULLPTR nullptr
#   else
#    define YY_NULLPTR 0
#   endif
#  else
#   define YY_NULLPTR ((void*)0)
#  endif
# endif

#include "grammar.tab.h"
/* Symbol kind.  */
enum yysymbol_kind_t
{
  YYSYMBOL_YYEMPTY = -2,
  YYSYMBOL_YYEOF = 0,                      /* "end of file"  */
  YYSYMBOL_YYerror = 1,                    /* error  */
  YYSYMBOL_YYUNDEF = 2,                    /* "invalid token"  */
  YYSYMBOL_REAL = 3,                       /* REAL  */
  YYSYMBOL_INT = 4,                        /* INT  */
  YYSYMBOL_NAME = 5,                       /* NAME  */
  YYSYMBOL_STRING = 6,                     /* STRING  */
  YYSYMBOL_RBNULL = 7,                     /* RBNULL  */
  YYSYMBOL_RBTAB = 8,                      /* RBTAB  */
  YYSYMBOL_FALSE = 9,                      /* FALSE  */
  YYSYMBOL_TRUE = 10,                      /* TRUE  */
  YYSYMBOL_RBINF = 11,                     /* RBINF  */
  YYSYMBOL_FUNCTION = 12,                  /* FUNCTION  */
  YYSYMBOL_PROCEDURE = 13,                 /* PROCEDURE  */
  YYSYMBOL_CLASS = 14,                     /* CLASS  */
  YYSYMBOL_FOR = 15,                       /* FOR  */
  YYSYMBOL_IN = 16,                        /* IN  */
  YYSYMBOL_IF = 17,                        /* IF  */
  YYSYMBOL_ELSE = 18,                      /* ELSE  */
  YYSYMBOL_WHILE = 19,                     /* WHILE  */
  YYSYMBOL_NEXT = 20,                      /* NEXT  */
  YYSYMBOL_BREAK = 21,                     /* BREAK  */
  YYSYMBOL_RETURN = 22,                    /* RETURN  */
  YYSYMBOL_MOD_CONST = 23,                 /* MOD_CONST  */
  YYSYMBOL_MOD_DYNAMIC = 24,               /* MOD_DYNAMIC  */
  YYSYMBOL_MOD_STOCHASTIC = 25,            /* MOD_STOCHASTIC  */
  YYSYMBOL_MOD_DETERMINISTIC = 26,         /* MOD_DETERMINISTIC  */
  YYSYMBOL_PROTECTED = 27,                 /* PROTECTED  */
  YYSYMBOL_ARROW_ASSIGN = 28,              /* ARROW_ASSIGN  */
  YYSYMBOL_TILDE_ASSIGN = 29,              /* TILDE_ASSIGN  */
  YYSYMBOL_EQUATION_ASSIGN = 30,           /* EQUATION_ASSIGN  */
  YYSYMBOL_WORKSPACE_ASSIGN = 31,          /* WORKSPACE_ASSIGN  */
  YYSYMBOL_REFERENCE_ASSIGN = 32,          /* REFERENCE_ASSIGN  */
  YYSYMBOL_ADDITION_ASSIGN = 33,           /* ADDITION_ASSIGN  */
  YYSYMBOL_SUBTRACTION_ASSIGN = 34,        /* SUBTRACTION_ASSIGN  */
  YYSYMBOL_MULTIPLICATION_ASSIGN = 35,     /* MULTIPLICATION_ASSIGN  */
  YYSYMBOL_DIVISION_ASSIGN = 36,           /* DIVISION_ASSIGN  */
  YYSYMBOL_DECREMENT = 37,                 /* DECREMENT  */
  YYSYMBOL_INCREMENT = 38,                 /* INCREMENT  */
  YYSYMBOL_EQUAL = 39,                     /* EQUAL  */
  YYSYMBOL_AND = 40,                       /* AND  */
  YYSYMBOL_OR = 41,                        /* OR  */
  YYSYMBOL_AND2 = 42,                      /* AND2  */
  YYSYMBOL_OR2 = 43,                       /* OR2  */
  YYSYMBOL_GT = 44,                        /* GT  */
  YYSYMBOL_GE = 45,                        /* GE  */
  YYSYMBOL_LT = 46,                        /* LT  */
  YYSYMBOL_LE = 47,                        /* LE  */
  YYSYMBOL_EQ = 48,                        /* EQ  */
  YYSYMBOL_NE = 49,                        /* NE  */
  YYSYMBOL_PIPE = 50,                      /* PIPE  */
  YYSYMBOL_PIPE_PLACEHOLDER = 51,          /* PIPE_PLACEHOLDER  */
  YYSYMBOL_END_OF_INPUT = 52,              /* END_OF_INPUT  */
  YYSYMBOL_53_ = 53,                       /* '?'  */
  YYSYMBOL_UNOT = 54,                      /* UNOT  */
  YYSYMBOL_55_ = 55,                       /* '+'  */
  YYSYMBOL_56_ = 56,                       /* '-'  */
  YYSYMBOL_57_ = 57,                       /* '*'  */
  YYSYMBOL_58_ = 58,                       /* '/'  */
  YYSYMBOL_59_ = 59,                       /* ':'  */
  YYSYMBOL_60_ = 60,                       /* '%'  */
  YYSYMBOL_UMINUS = 61,                    /* UMINUS  */
  YYSYMBOL_UPLUS = 62,                     /* UPLUS  */
  YYSYMBOL_63_ = 63,                       /* '^'  */
  YYSYMBOL_64_ = 64,                       /* '.'  */
  YYSYMBOL_UAND = 65,                      /* UAND  */
  YYSYMBOL_66_ = 66,                       /* '('  */
  YYSYMBOL_67_ = 67,                       /* '['  */
  YYSYMBOL_68_n_ = 68,                     /* '\n'  */
  YYSYMBOL_69_ = 69,                       /* ';'  */
  YYSYMBOL_70_ = 70,                       /* ')'  */
  YYSYMBOL_71_ = 71,                       /* '!'  */
  YYSYMBOL_72_ = 72,                       /* ']'  */
  YYSYMBOL_73_ = 73,                       /* ','  */
  YYSYMBOL_74_ = 74,                       /* '{'  */
  YYSYMBOL_75_ = 75,                       /* '}'  */
  YYSYMBOL_YYACCEPT = 76,                  /* $accept  */
  YYSYMBOL_prog = 77,                      /* prog  */
  YYSYMBOL_expression = 78,                /* expression  */
  YYSYMBOL_arrowAssign = 79,               /* arrowAssign  */
  YYSYMBOL_tildeAssign = 80,               /* tildeAssign  */
  YYSYMBOL_equationAssign = 81,            /* equationAssign  */
  YYSYMBOL_workspaceAssign = 82,           /* workspaceAssign  */
  YYSYMBOL_referenceAssign = 83,           /* referenceAssign  */
  YYSYMBOL_additionAssign = 84,            /* additionAssign  */
  YYSYMBOL_subtractionAssign = 85,         /* subtractionAssign  */
  YYSYMBOL_multiplicationAssign = 86,      /* multiplicationAssign  */
  YYSYMBOL_divisionAssign = 87,            /* divisionAssign  */
  YYSYMBOL_variable = 88,                  /* variable  */
  YYSYMBOL_optElements = 89,               /* optElements  */
  YYSYMBOL_elementList = 90,               /* elementList  */
  YYSYMBOL_fxnCall = 91,                   /* fxnCall  */
  YYSYMBOL_functionCall = 92,              /* functionCall  */
  YYSYMBOL_optArguments = 93,              /* optArguments  */
  YYSYMBOL_argumentList = 94,              /* argumentList  */
  YYSYMBOL_argument = 95,                  /* argument  */
  YYSYMBOL_functionDef = 96,               /* functionDef  */
  YYSYMBOL_procedureDef = 97,              /* procedureDef  */
  YYSYMBOL_optFormals = 98,                /* optFormals  */
  YYSYMBOL_formalList = 99,                /* formalList  */
  YYSYMBOL_formal = 100,                   /* formal  */
  YYSYMBOL_typeSpec = 101,                 /* typeSpec  */
  YYSYMBOL_optDims = 102,                  /* optDims  */
  YYSYMBOL_dimList = 103,                  /* dimList  */
  YYSYMBOL_stmts = 104,                    /* stmts  */
  YYSYMBOL_stmtList = 105,                 /* stmtList  */
  YYSYMBOL_statement = 106,                /* statement  */
  YYSYMBOL_stmt_or_expr = 107,             /* stmt_or_expr  */
  YYSYMBOL_declaration = 108,              /* declaration  */
  YYSYMBOL_memberDefs = 109,               /* memberDefs  */
  YYSYMBOL_memberDef = 110,                /* memberDef  */
  YYSYMBOL_classDef = 111,                 /* classDef  */
  YYSYMBOL_ifStatement = 112,              /* ifStatement  */
  YYSYMBOL_cond = 113,                     /* cond  */
  YYSYMBOL_forStatement = 114,             /* forStatement  */
  YYSYMBOL_forCond = 115,                  /* forCond  */
  YYSYMBOL_whileStatement = 116,           /* whileStatement  */
  YYSYMBOL_nextStatement = 117,            /* nextStatement  */
  YYSYMBOL_breakStatement = 118,           /* breakStatement  */
  YYSYMBOL_returnStatement = 119,          /* returnStatement  */
  YYSYMBOL_identifier = 120,               /* identifier  */
  YYSYMBOL_vector = 121,                   /* vector  */
  YYSYMBOL_vectorList = 122,               /* vectorList  */
  YYSYMBOL_constant = 123                  /* constant  */
};
typedef enum yysymbol_kind_t yysymbol_kind_t;




#ifdef short
# undef short
#endif

/* On compilers that do not define __PTRDIFF_MAX__ etc., make sure
   <limits.h> and (if available) <stdint.h> are included
   so that the code can choose integer types of a good width.  */

#ifndef __PTRDIFF_MAX__
# include <limits.h> /* INFRINGES ON USER NAME SPACE */
# if defined __STDC_VERSION__ && 199901 <= __STDC_VERSION__
#  include <stdint.h> /* INFRINGES ON USER NAME SPACE */
#  define YY_STDINT_H
# endif
#endif

/* Narrow types that promote to a signed type and that can represent a
   signed or unsigned integer of at least N bits.  In tables they can
   save space and decrease cache pressure.  Promoting to a signed type
   helps avoid bugs in integer arithmetic.  */

#ifdef __INT_LEAST8_MAX__
typedef __INT_LEAST8_TYPE__ yytype_int8;
#elif defined YY_STDINT_H
typedef int_least8_t yytype_int8;
#else
typedef signed char yytype_int8;
#endif

#ifdef __INT_LEAST16_MAX__
typedef __INT_LEAST16_TYPE__ yytype_int16;
#elif defined YY_STDINT_H
typedef int_least16_t yytype_int16;
#else
typedef short yytype_int16;
#endif

/* Work around bug in HP-UX 11.23, which defines these macros
   incorrectly for preprocessor constants.  This workaround can likely
   be removed in 2023, as HPE has promised support for HP-UX 11.23
   (aka HP-UX 11i v2) only through the end of 2022; see Table 2 of
   <https://h20195.www2.hpe.com/V2/getpdf.aspx/4AA4-7673ENW.pdf>.  */
#ifdef __hpux
# undef UINT_LEAST8_MAX
# undef UINT_LEAST16_MAX
# define UINT_LEAST8_MAX 255
# define UINT_LEAST16_MAX 65535
#endif

#if defined __UINT_LEAST8_MAX__ && __UINT_LEAST8_MAX__ <= __INT_MAX__
typedef __UINT_LEAST8_TYPE__ yytype_uint8;
#elif (!defined __UINT_LEAST8_MAX__ && defined YY_STDINT_H \
       && UINT_LEAST8_MAX <= INT_MAX)
typedef uint_least8_t yytype_uint8;
#elif !defined __UINT_LEAST8_MAX__ && UCHAR_MAX <= INT_MAX
typedef unsigned char yytype_uint8;
#else
typedef short yytype_uint8;
#endif

#if defined __UINT_LEAST16_MAX__ && __UINT_LEAST16_MAX__ <= __INT_MAX__
typedef __UINT_LEAST16_TYPE__ yytype_uint16;
#elif (!defined __UINT_LEAST16_MAX__ && defined YY_STDINT_H \
       && UINT_LEAST16_MAX <= INT_MAX)
typedef uint_least16_t yytype_uint16;
#elif !defined __UINT_LEAST16_MAX__ && USHRT_MAX <= INT_MAX
typedef unsigned short yytype_uint16;
#else
typedef int yytype_uint16;
#endif

#ifndef YYPTRDIFF_T
# if defined __PTRDIFF_TYPE__ && defined __PTRDIFF_MAX__
#  define YYPTRDIFF_T __PTRDIFF_TYPE__
#  define YYPTRDIFF_MAXIMUM __PTRDIFF_MAX__
# elif defined PTRDIFF_MAX
#  ifndef ptrdiff_t
#   include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  endif
#  define YYPTRDIFF_T ptrdiff_t
#  define YYPTRDIFF_MAXIMUM PTRDIFF_MAX
# else
#  define YYPTRDIFF_T long
#  define YYPTRDIFF_MAXIMUM LONG_MAX
# endif
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif defined __STDC_VERSION__ && 199901 <= __STDC_VERSION__
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned
# endif
#endif

#define YYSIZE_MAXIMUM                                  \
  YY_CAST (YYPTRDIFF_T,                                 \
           (YYPTRDIFF_MAXIMUM < YY_CAST (YYSIZE_T, -1)  \
            ? YYPTRDIFF_MAXIMUM                         \
            : YY_CAST (YYSIZE_T, -1)))

#define YYSIZEOF(X) YY_CAST (YYPTRDIFF_T, sizeof (X))


/* Stored state numbers (used for stacks). */
typedef yytype_int16 yy_state_t;

/* State numbers in computations.  */
typedef int yy_state_fast_t;

#ifndef YY_
# if defined YYENABLE_NLS && YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(Msgid) dgettext ("bison-runtime", Msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(Msgid) Msgid
# endif
#endif


#ifndef YY_ATTRIBUTE_PURE
# if defined __GNUC__ && 2 < __GNUC__ + (96 <= __GNUC_MINOR__)
#  define YY_ATTRIBUTE_PURE __attribute__ ((__pure__))
# else
#  define YY_ATTRIBUTE_PURE
# endif
#endif

#ifndef YY_ATTRIBUTE_UNUSED
# if defined __GNUC__ && 2 < __GNUC__ + (7 <= __GNUC_MINOR__)
#  define YY_ATTRIBUTE_UNUSED __attribute__ ((__unused__))
# else
#  define YY_ATTRIBUTE_UNUSED
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YY_USE(E) ((void) (E))
#else
# define YY_USE(E) /* empty */
#endif

/* Suppress an incorrect diagnostic about yylval being uninitialized.  */
#if defined __GNUC__ && ! defined __ICC && 406 <= __GNUC__ * 100 + __GNUC_MINOR__
# if __GNUC__ * 100 + __GNUC_MINOR__ < 407
#  define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN                           \
    _Pragma ("GCC diagnostic push")                                     \
    _Pragma ("GCC diagnostic ignored \"-Wuninitialized\"")
# else
#  define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN                           \
    _Pragma ("GCC diagnostic push")                                     \
    _Pragma ("GCC diagnostic ignored \"-Wuninitialized\"")              \
    _Pragma ("GCC diagnostic ignored \"-Wmaybe-uninitialized\"")
# endif
# define YY_IGNORE_MAYBE_UNINITIALIZED_END      \
    _Pragma ("GCC diagnostic pop")
#else
# define YY_INITIAL_VALUE(Value) Value
#endif
#ifndef YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_END
#endif
#ifndef YY_INITIAL_VALUE
# define YY_INITIAL_VALUE(Value) /* Nothing. */
#endif

#if defined __cplusplus && defined __GNUC__ && ! defined __ICC && 6 <= __GNUC__
# define YY_IGNORE_USELESS_CAST_BEGIN                          \
    _Pragma ("GCC diagnostic push")                            \
    _Pragma ("GCC diagnostic ignored \"-Wuseless-cast\"")
# define YY_IGNORE_USELESS_CAST_END            \
    _Pragma ("GCC diagnostic pop")
#endif
#ifndef YY_IGNORE_USELESS_CAST_BEGIN
# define YY_IGNORE_USELESS_CAST_BEGIN
# define YY_IGNORE_USELESS_CAST_END
#endif


#define YY_ASSERT(E) ((void) (0 && (E)))

#if !defined yyoverflow

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined EXIT_SUCCESS
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
      /* Use EXIT_SUCCESS as a witness for stdlib.h.  */
#     ifndef EXIT_SUCCESS
#      define EXIT_SUCCESS 0
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's 'empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (0)
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined EXIT_SUCCESS \
       && ! ((defined YYMALLOC || defined malloc) \
             && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef EXIT_SUCCESS
#    define EXIT_SUCCESS 0
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined EXIT_SUCCESS
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined EXIT_SUCCESS
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* !defined yyoverflow */

#if (! defined yyoverflow \
     && (! defined __cplusplus \
         || (defined YYLTYPE_IS_TRIVIAL && YYLTYPE_IS_TRIVIAL \
             && defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yy_state_t yyss_alloc;
  YYSTYPE yyvs_alloc;
  YYLTYPE yyls_alloc;
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (YYSIZEOF (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (YYSIZEOF (yy_state_t) + YYSIZEOF (YYSTYPE) \
             + YYSIZEOF (YYLTYPE)) \
      + 2 * YYSTACK_GAP_MAXIMUM)

# define YYCOPY_NEEDED 1

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack_alloc, Stack)                           \
    do                                                                  \
      {                                                                 \
        YYPTRDIFF_T yynewbytes;                                         \
        YYCOPY (&yyptr->Stack_alloc, Stack, yysize);                    \
        Stack = &yyptr->Stack_alloc;                                    \
        yynewbytes = yystacksize * YYSIZEOF (*Stack) + YYSTACK_GAP_MAXIMUM; \
        yyptr += yynewbytes / YYSIZEOF (*yyptr);                        \
      }                                                                 \
    while (0)

#endif

#if defined YYCOPY_NEEDED && YYCOPY_NEEDED
/* Copy COUNT objects from SRC to DST.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(Dst, Src, Count) \
      __builtin_memcpy (Dst, Src, YY_CAST (YYSIZE_T, (Count)) * sizeof (*(Src)))
#  else
#   define YYCOPY(Dst, Src, Count)              \
      do                                        \
        {                                       \
          YYPTRDIFF_T yyi;                      \
          for (yyi = 0; yyi < (Count); yyi++)   \
            (Dst)[yyi] = (Src)[yyi];            \
        }                                       \
      while (0)
#  endif
# endif
#endif /* !YYCOPY_NEEDED */

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  87
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   1070

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  76
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  48
/* YYNRULES -- Number of rules.  */
#define YYNRULES  159
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  279

/* YYMAXUTOK -- Last valid token kind.  */
#define YYMAXUTOK   311


/* YYTRANSLATE(TOKEN-NUM) -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex, with out-of-bounds checking.  */
#define YYTRANSLATE(YYX)                                \
  (0 <= (YYX) && (YYX) <= YYMAXUTOK                     \
   ? YY_CAST (yysymbol_kind_t, yytranslate[YYX])        \
   : YYSYMBOL_YYUNDEF)

/* YYTRANSLATE[TOKEN-NUM] -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex.  */
static const yytype_int8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
      68,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    71,     2,     2,     2,    60,     2,     2,
      66,    70,    57,    55,    73,    56,    64,    58,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,    59,    69,
       2,     2,     2,    53,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,    67,     2,    72,    63,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    74,     2,    75,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    33,    34,
      35,    36,    37,    38,    39,    40,    41,    42,    43,    44,
      45,    46,    47,    48,    49,    50,    51,    52,    54,    61,
      62,    65
};

#if YYDEBUG
/* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
static const yytype_int16 yyrline[] =
{
       0,   237,   237,   244,   251,   260,   269,   278,   287,   296,
     305,   315,   325,   334,   343,   350,   359,   361,   363,   365,
     366,   367,   368,   370,   371,   372,   373,   375,   377,   378,
     379,   380,   381,   382,   384,   385,   386,   387,   388,   389,
     391,   392,   393,   394,   396,   398,   399,   400,   401,   402,
     404,   405,   406,   407,   409,   411,   415,   424,   433,   442,
     451,   460,   469,   478,   487,   496,   509,   521,   533,   548,
     549,   552,   553,   554,   555,   558,   565,   572,   580,   588,
     596,   597,   600,   601,   604,   611,   621,   630,   643,   652,
     665,   666,   669,   670,   673,   681,   689,   698,   709,   710,
     711,   712,   713,   716,   717,   720,   721,   724,   725,   733,
     734,   735,   736,   737,   738,   741,   742,   743,   744,   745,
     746,   749,   750,   753,   754,   755,   756,   767,   768,   769,
     770,   771,   772,   775,   776,   777,   778,   781,   792,   793,
     795,   798,   801,   804,   807,   810,   813,   814,   817,   821,
     824,   825,   831,   835,   839,   843,   847,   851,   860,   864
};
#endif

/** Accessing symbol of state STATE.  */
#define YY_ACCESSING_SYMBOL(State) YY_CAST (yysymbol_kind_t, yystos[State])

#if YYDEBUG || 0
/* The user-facing name of the symbol whose (internal) number is
   YYSYMBOL.  No bounds checking.  */
static const char *yysymbol_name (yysymbol_kind_t yysymbol) YY_ATTRIBUTE_UNUSED;

/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "\"end of file\"", "error", "\"invalid token\"", "REAL", "INT", "NAME",
  "STRING", "RBNULL", "RBTAB", "FALSE", "TRUE", "RBINF", "FUNCTION",
  "PROCEDURE", "CLASS", "FOR", "IN", "IF", "ELSE", "WHILE", "NEXT",
  "BREAK", "RETURN", "MOD_CONST", "MOD_DYNAMIC", "MOD_STOCHASTIC",
  "MOD_DETERMINISTIC", "PROTECTED", "ARROW_ASSIGN", "TILDE_ASSIGN",
  "EQUATION_ASSIGN", "WORKSPACE_ASSIGN", "REFERENCE_ASSIGN",
  "ADDITION_ASSIGN", "SUBTRACTION_ASSIGN", "MULTIPLICATION_ASSIGN",
  "DIVISION_ASSIGN", "DECREMENT", "INCREMENT", "EQUAL", "AND", "OR",
  "AND2", "OR2", "GT", "GE", "LT", "LE", "EQ", "NE", "PIPE",
  "PIPE_PLACEHOLDER", "END_OF_INPUT", "'?'", "UNOT", "'+'", "'-'", "'*'",
  "'/'", "':'", "'%'", "UMINUS", "UPLUS", "'^'", "'.'", "UAND", "'('",
  "'['", "'\\n'", "';'", "')'", "'!'", "']'", "','", "'{'", "'}'",
  "$accept", "prog", "expression", "arrowAssign", "tildeAssign",
  "equationAssign", "workspaceAssign", "referenceAssign", "additionAssign",
  "subtractionAssign", "multiplicationAssign", "divisionAssign",
  "variable", "optElements", "elementList", "fxnCall", "functionCall",
  "optArguments", "argumentList", "argument", "functionDef",
  "procedureDef", "optFormals", "formalList", "formal", "typeSpec",
  "optDims", "dimList", "stmts", "stmtList", "statement", "stmt_or_expr",
  "declaration", "memberDefs", "memberDef", "classDef", "ifStatement",
  "cond", "forStatement", "forCond", "whileStatement", "nextStatement",
  "breakStatement", "returnStatement", "identifier", "vector",
  "vectorList", "constant", YY_NULLPTR
};

static const char *
yysymbol_name (yysymbol_kind_t yysymbol)
{
  return yytname[yysymbol];
}
#endif

#define YYPACT_NINF (-220)

#define yypact_value_is_default(Yyn) \
  ((Yyn) == YYPACT_NINF)

#define YYTABLE_NINF (-104)

#define yytable_value_is_error(Yyn) \
  ((Yyn) == YYTABLE_NINF)

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
static const yytype_int16 yypact[] =
{
     364,   -13,  -220,  -220,  -220,  -220,  -220,  -220,  -220,  -220,
    -220,    11,    11,    11,   -61,   -46,   -46,  -220,  -220,   517,
       2,     2,   517,  -220,  -220,   517,   517,   517,   517,   517,
    -220,   517,    23,   754,  -220,  -220,  -220,  -220,  -220,  -220,
    -220,  -220,  -220,   -29,   -38,   -30,  -220,  -220,  -220,     3,
      18,  -220,  -220,  -220,  -220,  -220,  -220,  -220,    33,  -220,
    -220,  -220,  -220,    38,    42,   -23,    11,   292,   517,   292,
     292,   754,    33,   517,   -14,   -38,   -14,     7,   754,    15,
      29,   -44,   -44,   598,   754,    -3,   278,  -220,   517,   517,
     517,   517,   517,   517,   517,   517,   517,   517,   517,   517,
     517,   517,   517,   517,   517,   517,   517,    11,   517,   517,
     517,   517,   517,   517,   517,  -220,  -220,    11,   385,   -20,
      11,  -220,  -220,  -220,  -220,   517,    11,   -20,    19,     6,
      11,    14,    19,    11,    11,    57,   463,  -220,  -220,   637,
      84,  -220,  -220,   676,    11,  -220,  -220,    11,  -220,  -220,
     -38,  -220,   517,   754,   786,   818,   850,   881,   911,   937,
     962,   987,   278,   206,   278,   206,  1007,  1007,  1007,  1007,
    1007,  1007,  -220,    55,    17,    17,   -33,   -33,   -44,   -44,
     -44,   -38,  -220,   503,   439,  -220,   754,    40,    49,  -220,
     -35,  -220,    11,    11,    11,    11,    53,    51,  -220,    11,
      -2,  -220,    61,    68,    59,    69,    67,   517,   -47,  -220,
    -220,   292,   -38,   -38,    44,   -20,   754,   -20,  -220,  -220,
     557,  -220,   517,   517,    75,    75,    75,    75,   292,    19,
     104,   517,  -220,    19,  -220,   292,    19,    36,   715,   463,
     463,  -220,  -220,  -220,  -220,  -220,  -220,   754,  -220,  -220,
    -220,  -220,  -220,  -220,   517,   754,    76,  -220,    78,    19,
    -220,  -220,  -220,   -17,  -220,  -220,  -220,  -220,   754,   292,
     292,  -220,    36,    36,  -220,  -220,  -220,  -220,  -220
};

/* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
   Performed when YYTABLE does not specify something else to do.  Zero
   means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       0,     0,   159,   157,   148,   158,   154,   155,   152,   153,
     156,     0,     0,     0,     0,     0,     0,   144,   145,   146,
       0,     0,     0,    44,     2,     0,     0,     0,     0,     0,
       3,     0,     0,   122,    45,    47,    46,    48,    49,    50,
      51,    52,    53,    55,    76,    54,   124,   125,   121,     0,
       0,   123,   115,   116,   117,   118,   119,   120,    69,    17,
      16,    14,    15,   103,   103,     0,     0,     0,     0,     0,
       0,   147,    69,     0,    23,     0,    25,    22,     0,    54,
      69,    20,    19,     0,   151,     0,    21,     1,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    24,    26,     0,     0,    66,
       0,     4,     5,     6,     7,    80,    65,    70,    90,     0,
       0,   104,    90,     0,     0,     0,   109,   141,   108,     0,
     138,   143,    65,     0,     0,    12,    13,     0,     8,     9,
      18,   149,     0,    56,    57,    58,    60,    61,    62,    63,
      64,    59,    40,    41,    42,    43,    39,    38,    34,    35,
      36,    37,    79,     0,    28,    29,    30,    31,    27,    33,
      32,    77,    72,     0,     0,    78,    84,     0,    81,    82,
      69,   126,     0,     0,     0,     0,     0,    91,    92,     0,
      94,   105,     0,     0,     0,     0,     0,     0,     0,   110,
     140,     0,     0,     0,     0,    67,   150,    68,    71,    74,
       0,    75,     0,     0,   103,   103,   103,   103,     0,     0,
      96,     0,    98,    90,   106,     0,    90,   127,     0,   114,
     112,   107,   139,    10,    11,    73,    83,    85,    99,   100,
     101,   102,    86,    93,     0,    95,     0,    88,     0,     0,
     135,   136,   133,     0,   128,   142,   113,   111,    97,     0,
       0,   134,   132,   130,   137,    87,    89,   131,   129
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
    -220,  -220,    63,  -220,  -220,  -220,  -220,  -220,  -220,  -220,
    -220,  -220,    95,    87,   -11,    -6,   124,  -220,  -220,   -72,
     178,   179,  -130,  -220,  -219,  -220,   -18,  -220,   -31,  -220,
    -220,     1,  -220,  -220,  -154,  -220,  -220,   167,  -220,  -220,
    -220,  -220,  -220,  -220,     0,  -220,  -220,  -220
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
       0,    32,    33,    34,    35,    36,    37,    38,    39,    40,
      41,    42,    43,   142,   127,    44,    45,   187,   188,   189,
     260,   261,   196,   197,   198,   199,   130,   131,   137,   208,
      48,   138,    50,   263,   264,    51,    52,    69,    53,    67,
      54,    55,    56,    57,    72,    59,    85,    60
};

/* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule whose
   number is the opposite.  If YYTABLE_NINF, syntax error.  */
static const yytype_int16 yytable[] =
{
      58,    49,   204,  -103,   223,    66,   107,     4,   115,   116,
     253,    63,    64,    65,    75,    75,     4,   107,   262,   114,
      68,   239,   240,    87,     4,    80,   112,   113,   241,   118,
     114,   125,   118,   119,   120,   117,   134,   231,   140,   141,
     271,     4,   192,   193,   194,   195,   133,   184,    11,    12,
     144,   272,   273,   262,   262,    61,    62,   107,   274,   192,
     193,   194,   195,   259,   119,   129,   135,   107,    73,   151,
     152,   121,   122,   207,   110,   111,   112,   113,   201,   120,
     114,   203,    71,   145,   146,    77,   123,   124,    78,    81,
      82,    83,    84,   147,    86,   125,   118,   148,   149,   125,
     118,   172,   211,   256,   128,   129,   258,   173,   132,   129,
     221,   181,   243,   244,   185,    74,    76,   173,   277,   278,
     173,   125,   222,   228,   229,   190,   191,   233,   200,   235,
     202,   139,   200,   205,   206,   236,   143,   209,   213,   215,
     234,   237,   129,   254,   173,   126,   269,   214,   270,    79,
     246,   153,   154,   155,   156,   157,   158,   159,   160,   161,
     162,   163,   164,   165,   166,   167,   168,   169,   170,   171,
     217,   174,   175,   176,   177,   178,   179,   180,    46,    47,
     242,   183,   232,    70,     0,     0,     0,     0,   186,     0,
       0,     0,   224,   225,   226,   227,     0,   252,     0,   230,
       0,   215,   217,     0,   257,     0,   248,   249,   250,   251,
       0,     0,     0,     0,     0,   216,     0,     0,     0,     0,
       0,     0,   190,     0,     0,     0,     0,     0,     0,   200,
       0,     0,     0,   200,     0,     0,   200,   200,   275,   276,
     266,   267,     0,     0,     0,     0,    97,   220,    99,     0,
     101,   102,   103,   104,   105,   106,   107,     0,     0,   200,
       0,   108,   109,   110,   111,   112,   113,     0,     0,   114,
     238,     0,   200,   200,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   186,   247,     0,     0,     0,
       0,     0,     0,     0,   255,     2,     3,     4,     5,     6,
       7,     8,     9,    10,     0,     0,     0,    14,     0,    15,
       0,    16,    17,    18,    19,     0,     0,   268,     0,     0,
       0,     0,   101,   102,   103,   104,   105,   106,   107,    20,
      21,     0,    22,   108,   109,   110,   111,   112,   113,     0,
       0,   114,     0,    23,     0,     0,     0,    26,    27,     0,
       0,     0,     0,     0,     0,     0,     0,     0,    28,    29,
       0,     0,     0,    31,     0,     1,   136,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
       0,    15,     0,    16,    17,    18,    19,     0,     2,     3,
       4,     5,     6,     7,     8,     9,    10,     0,     0,     0,
       0,    20,    21,     0,    22,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    23,    24,    25,     0,    26,
      27,     0,    20,    21,     0,    22,     0,     0,     0,     0,
      28,    29,    30,     0,     0,    31,    23,     0,     0,     0,
      26,    27,     2,     3,     4,     5,     6,     7,     8,     9,
      10,    28,    29,     0,     0,     0,    31,   182,     0,     0,
       0,     0,     0,     0,     0,     0,     2,     3,     4,     5,
       6,     7,     8,     9,    10,     0,    20,    21,    14,    22,
      15,     0,    16,    17,    18,    19,     0,     0,     0,     0,
      23,     0,     0,     0,    26,    27,     0,     0,     0,     0,
      20,    21,     0,    22,     0,    28,    29,     0,     0,     0,
      31,   219,     0,     0,    23,     0,     0,     0,    26,    27,
       2,     3,     4,     5,     6,     7,     8,     9,    10,    28,
      29,    88,    89,    90,    31,    91,    92,    93,    94,    95,
       0,     0,    96,    97,    98,    99,   100,   101,   102,   103,
     104,   105,   106,   107,    20,    21,     0,    22,   108,   109,
     110,   111,   112,   113,     0,     0,   114,     0,    23,     0,
       0,     0,    26,    27,     0,   218,     0,     0,     0,     0,
       0,     0,     0,    28,    29,    88,    89,    90,    31,    91,
      92,    93,    94,    95,     0,     0,    96,    97,    98,    99,
     100,   101,   102,   103,   104,   105,   106,   107,     0,     0,
       0,     0,   108,   109,   110,   111,   112,   113,     0,     0,
     114,     0,     0,     0,     0,     0,    88,    89,    90,   245,
      91,    92,    93,    94,    95,     0,     0,    96,    97,    98,
      99,   100,   101,   102,   103,   104,   105,   106,   107,     0,
       0,     0,     0,   108,   109,   110,   111,   112,   113,     0,
       0,   114,     0,     0,     0,    88,    89,    90,   150,    91,
      92,    93,    94,    95,     0,     0,    96,    97,    98,    99,
     100,   101,   102,   103,   104,   105,   106,   107,     0,     0,
       0,     0,   108,   109,   110,   111,   112,   113,     0,     0,
     114,     0,     0,     0,    88,    89,    90,   210,    91,    92,
      93,    94,    95,     0,     0,    96,    97,    98,    99,   100,
     101,   102,   103,   104,   105,   106,   107,     0,     0,     0,
       0,   108,   109,   110,   111,   112,   113,     0,     0,   114,
       0,     0,     0,    88,    89,    90,   212,    91,    92,    93,
      94,    95,     0,     0,    96,    97,    98,    99,   100,   101,
     102,   103,   104,   105,   106,   107,     0,     0,     0,     0,
     108,   109,   110,   111,   112,   113,     0,     0,   114,     0,
       0,     0,    88,    89,    90,   265,    91,    92,    93,    94,
      95,     0,     0,    96,    97,    98,    99,   100,   101,   102,
     103,   104,   105,   106,   107,     0,     0,     0,     0,   108,
     109,   110,   111,   112,   113,    89,    90,   114,    91,    92,
      93,    94,    95,     0,     0,    96,    97,    98,    99,   100,
     101,   102,   103,   104,   105,   106,   107,     0,     0,     0,
       0,   108,   109,   110,   111,   112,   113,     0,    90,   114,
      91,    92,    93,    94,    95,     0,     0,    96,    97,    98,
      99,   100,   101,   102,   103,   104,   105,   106,   107,     0,
       0,     0,     0,   108,   109,   110,   111,   112,   113,     0,
       0,   114,    91,    92,    93,    94,    95,     0,     0,    96,
      97,    98,    99,   100,   101,   102,   103,   104,   105,   106,
     107,     0,     0,     0,     0,   108,   109,   110,   111,   112,
     113,     0,     0,   114,    92,    93,    94,    95,     0,     0,
      96,    97,    98,    99,   100,   101,   102,   103,   104,   105,
     106,   107,     0,     0,     0,     0,   108,   109,   110,   111,
     112,   113,     0,     0,   114,    93,    94,    95,     0,     0,
      96,    97,    98,    99,   100,   101,   102,   103,   104,   105,
     106,   107,     0,     0,     0,     0,   108,   109,   110,   111,
     112,   113,    94,    95,   114,     0,    96,    97,    98,    99,
     100,   101,   102,   103,   104,   105,   106,   107,     0,     0,
       0,     0,   108,   109,   110,   111,   112,   113,    95,     0,
     114,    96,    97,    98,    99,   100,   101,   102,   103,   104,
     105,   106,   107,     0,     0,     0,     0,   108,   109,   110,
     111,   112,   113,     0,     0,   114,    96,    97,    98,    99,
     100,   101,   102,   103,   104,   105,   106,   107,     0,     0,
       0,     0,   108,   109,   110,   111,   112,   113,     0,     0,
     114,  -104,  -104,  -104,  -104,  -104,  -104,   107,     0,     0,
       0,     0,   108,   109,   110,   111,   112,   113,     0,     0,
     114
};

static const yytype_int16 yycheck[] =
{
       0,     0,   132,     5,    39,    66,    50,     5,    37,    38,
     229,    11,    12,    13,    20,    21,     5,    50,   237,    63,
      66,    68,    69,     0,     5,    25,    59,    60,    75,    67,
      63,    66,    67,    44,    64,    64,    59,    39,    69,    70,
     259,     5,    23,    24,    25,    26,    64,    67,    12,    13,
      64,    68,    69,   272,   273,    68,    69,    50,    75,    23,
      24,    25,    26,    27,    75,    67,    66,    50,    66,    72,
      73,    68,    69,    16,    57,    58,    59,    60,    72,    64,
      63,    67,    19,    68,    69,    22,    68,    69,    25,    26,
      27,    28,    29,    64,    31,    66,    67,    68,    69,    66,
      67,   107,    18,   233,    66,    67,   236,   107,    66,    67,
      70,   117,    68,    69,   120,    20,    21,   117,   272,   273,
     120,    66,    73,    70,    73,   125,   126,    66,   128,    70,
     130,    68,   132,   133,   134,    66,    73,   136,   144,   150,
      72,    74,    67,    39,   144,    58,    70,   147,    70,    25,
     222,    88,    89,    90,    91,    92,    93,    94,    95,    96,
      97,    98,    99,   100,   101,   102,   103,   104,   105,   106,
     181,   108,   109,   110,   111,   112,   113,   114,     0,     0,
     211,   118,   200,    16,    -1,    -1,    -1,    -1,   125,    -1,
      -1,    -1,   192,   193,   194,   195,    -1,   228,    -1,   199,
      -1,   212,   213,    -1,   235,    -1,   224,   225,   226,   227,
      -1,    -1,    -1,    -1,    -1,   152,    -1,    -1,    -1,    -1,
      -1,    -1,   222,    -1,    -1,    -1,    -1,    -1,    -1,   229,
      -1,    -1,    -1,   233,    -1,    -1,   236,   237,   269,   270,
     239,   240,    -1,    -1,    -1,    -1,    40,   184,    42,    -1,
      44,    45,    46,    47,    48,    49,    50,    -1,    -1,   259,
      -1,    55,    56,    57,    58,    59,    60,    -1,    -1,    63,
     207,    -1,   272,   273,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   222,   223,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   231,     3,     4,     5,     6,     7,
       8,     9,    10,    11,    -1,    -1,    -1,    15,    -1,    17,
      -1,    19,    20,    21,    22,    -1,    -1,   254,    -1,    -1,
      -1,    -1,    44,    45,    46,    47,    48,    49,    50,    37,
      38,    -1,    40,    55,    56,    57,    58,    59,    60,    -1,
      -1,    63,    -1,    51,    -1,    -1,    -1,    55,    56,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    66,    67,
      -1,    -1,    -1,    71,    -1,     1,    74,     3,     4,     5,
       6,     7,     8,     9,    10,    11,    12,    13,    14,    15,
      -1,    17,    -1,    19,    20,    21,    22,    -1,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    -1,    -1,    -1,
      -1,    37,    38,    -1,    40,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    51,    52,    53,    -1,    55,
      56,    -1,    37,    38,    -1,    40,    -1,    -1,    -1,    -1,
      66,    67,    68,    -1,    -1,    71,    51,    -1,    -1,    -1,
      55,    56,     3,     4,     5,     6,     7,     8,     9,    10,
      11,    66,    67,    -1,    -1,    -1,    71,    72,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,     3,     4,     5,     6,
       7,     8,     9,    10,    11,    -1,    37,    38,    15,    40,
      17,    -1,    19,    20,    21,    22,    -1,    -1,    -1,    -1,
      51,    -1,    -1,    -1,    55,    56,    -1,    -1,    -1,    -1,
      37,    38,    -1,    40,    -1,    66,    67,    -1,    -1,    -1,
      71,    72,    -1,    -1,    51,    -1,    -1,    -1,    55,    56,
       3,     4,     5,     6,     7,     8,     9,    10,    11,    66,
      67,    28,    29,    30,    71,    32,    33,    34,    35,    36,
      -1,    -1,    39,    40,    41,    42,    43,    44,    45,    46,
      47,    48,    49,    50,    37,    38,    -1,    40,    55,    56,
      57,    58,    59,    60,    -1,    -1,    63,    -1,    51,    -1,
      -1,    -1,    55,    56,    -1,    72,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    66,    67,    28,    29,    30,    71,    32,
      33,    34,    35,    36,    -1,    -1,    39,    40,    41,    42,
      43,    44,    45,    46,    47,    48,    49,    50,    -1,    -1,
      -1,    -1,    55,    56,    57,    58,    59,    60,    -1,    -1,
      63,    -1,    -1,    -1,    -1,    -1,    28,    29,    30,    72,
      32,    33,    34,    35,    36,    -1,    -1,    39,    40,    41,
      42,    43,    44,    45,    46,    47,    48,    49,    50,    -1,
      -1,    -1,    -1,    55,    56,    57,    58,    59,    60,    -1,
      -1,    63,    -1,    -1,    -1,    28,    29,    30,    70,    32,
      33,    34,    35,    36,    -1,    -1,    39,    40,    41,    42,
      43,    44,    45,    46,    47,    48,    49,    50,    -1,    -1,
      -1,    -1,    55,    56,    57,    58,    59,    60,    -1,    -1,
      63,    -1,    -1,    -1,    28,    29,    30,    70,    32,    33,
      34,    35,    36,    -1,    -1,    39,    40,    41,    42,    43,
      44,    45,    46,    47,    48,    49,    50,    -1,    -1,    -1,
      -1,    55,    56,    57,    58,    59,    60,    -1,    -1,    63,
      -1,    -1,    -1,    28,    29,    30,    70,    32,    33,    34,
      35,    36,    -1,    -1,    39,    40,    41,    42,    43,    44,
      45,    46,    47,    48,    49,    50,    -1,    -1,    -1,    -1,
      55,    56,    57,    58,    59,    60,    -1,    -1,    63,    -1,
      -1,    -1,    28,    29,    30,    70,    32,    33,    34,    35,
      36,    -1,    -1,    39,    40,    41,    42,    43,    44,    45,
      46,    47,    48,    49,    50,    -1,    -1,    -1,    -1,    55,
      56,    57,    58,    59,    60,    29,    30,    63,    32,    33,
      34,    35,    36,    -1,    -1,    39,    40,    41,    42,    43,
      44,    45,    46,    47,    48,    49,    50,    -1,    -1,    -1,
      -1,    55,    56,    57,    58,    59,    60,    -1,    30,    63,
      32,    33,    34,    35,    36,    -1,    -1,    39,    40,    41,
      42,    43,    44,    45,    46,    47,    48,    49,    50,    -1,
      -1,    -1,    -1,    55,    56,    57,    58,    59,    60,    -1,
      -1,    63,    32,    33,    34,    35,    36,    -1,    -1,    39,
      40,    41,    42,    43,    44,    45,    46,    47,    48,    49,
      50,    -1,    -1,    -1,    -1,    55,    56,    57,    58,    59,
      60,    -1,    -1,    63,    33,    34,    35,    36,    -1,    -1,
      39,    40,    41,    42,    43,    44,    45,    46,    47,    48,
      49,    50,    -1,    -1,    -1,    -1,    55,    56,    57,    58,
      59,    60,    -1,    -1,    63,    34,    35,    36,    -1,    -1,
      39,    40,    41,    42,    43,    44,    45,    46,    47,    48,
      49,    50,    -1,    -1,    -1,    -1,    55,    56,    57,    58,
      59,    60,    35,    36,    63,    -1,    39,    40,    41,    42,
      43,    44,    45,    46,    47,    48,    49,    50,    -1,    -1,
      -1,    -1,    55,    56,    57,    58,    59,    60,    36,    -1,
      63,    39,    40,    41,    42,    43,    44,    45,    46,    47,
      48,    49,    50,    -1,    -1,    -1,    -1,    55,    56,    57,
      58,    59,    60,    -1,    -1,    63,    39,    40,    41,    42,
      43,    44,    45,    46,    47,    48,    49,    50,    -1,    -1,
      -1,    -1,    55,    56,    57,    58,    59,    60,    -1,    -1,
      63,    44,    45,    46,    47,    48,    49,    50,    -1,    -1,
      -1,    -1,    55,    56,    57,    58,    59,    60,    -1,    -1,
      63
};

/* YYSTOS[STATE-NUM] -- The symbol kind of the accessing symbol of
   state STATE-NUM.  */
static const yytype_int8 yystos[] =
{
       0,     1,     3,     4,     5,     6,     7,     8,     9,    10,
      11,    12,    13,    14,    15,    17,    19,    20,    21,    22,
      37,    38,    40,    51,    52,    53,    55,    56,    66,    67,
      68,    71,    77,    78,    79,    80,    81,    82,    83,    84,
      85,    86,    87,    88,    91,    92,    96,    97,   106,   107,
     108,   111,   112,   114,   116,   117,   118,   119,   120,   121,
     123,    68,    69,   120,   120,   120,    66,   115,    66,   113,
     113,    78,   120,    66,    88,    91,    88,    78,    78,    92,
     120,    78,    78,    78,    78,   122,    78,     0,    28,    29,
      30,    32,    33,    34,    35,    36,    39,    40,    41,    42,
      43,    44,    45,    46,    47,    48,    49,    50,    55,    56,
      57,    58,    59,    60,    63,    37,    38,    64,    67,    90,
      64,    68,    69,    68,    69,    66,    89,    90,    66,    67,
     102,   103,    66,   102,    59,   120,    74,   104,   107,    78,
     104,   104,    89,    78,    64,    68,    69,    64,    68,    69,
      70,    72,    73,    78,    78,    78,    78,    78,    78,    78,
      78,    78,    78,    78,    78,    78,    78,    78,    78,    78,
      78,    78,    91,   120,    78,    78,    78,    78,    78,    78,
      78,    91,    72,    78,    67,    91,    78,    93,    94,    95,
     120,   120,    23,    24,    25,    26,    98,    99,   100,   101,
     120,    72,   120,    67,    98,   120,   120,    16,   105,   107,
      70,    18,    70,    91,   120,    90,    78,    90,    72,    72,
      78,    70,    73,    39,   120,   120,   120,   120,    70,    73,
     120,    39,   102,    66,    72,    70,    66,    74,    78,    68,
      69,    75,   104,    68,    69,    72,    95,    78,   102,   102,
     102,   102,   104,   100,    39,    78,    98,   104,    98,    27,
      96,    97,   100,   109,   110,    70,   107,   107,    78,    70,
      70,   100,    68,    69,    75,   104,   104,   110,   110
};

/* YYR1[RULE-NUM] -- Symbol kind of the left-hand side of rule RULE-NUM.  */
static const yytype_int8 yyr1[] =
{
       0,    76,    77,    77,    77,    77,    77,    77,    77,    77,
      77,    77,    77,    77,    77,    77,    78,    78,    78,    78,
      78,    78,    78,    78,    78,    78,    78,    78,    78,    78,
      78,    78,    78,    78,    78,    78,    78,    78,    78,    78,
      78,    78,    78,    78,    78,    78,    78,    78,    78,    78,
      78,    78,    78,    78,    78,    78,    79,    80,    81,    82,
      83,    84,    85,    86,    87,    88,    88,    88,    88,    89,
      89,    90,    90,    90,    90,    91,    92,    92,    92,    92,
      93,    93,    94,    94,    95,    95,    96,    96,    97,    97,
      98,    98,    99,    99,   100,   100,   100,   100,   101,   101,
     101,   101,   101,   102,   102,   103,   103,   104,   104,   105,
     105,   105,   105,   105,   105,   106,   106,   106,   106,   106,
     106,   107,   107,   108,   108,   108,   108,   109,   109,   109,
     109,   109,   109,   110,   110,   110,   110,   111,   112,   112,
     113,   114,   115,   116,   117,   118,   119,   119,   120,   121,
     122,   122,   123,   123,   123,   123,   123,   123,   123,   123
};

/* YYR2[RULE-NUM] -- Number of symbols on the right-hand side of rule RULE-NUM.  */
static const yytype_int8 yyr2[] =
{
       0,     2,     1,     1,     2,     2,     2,     2,     3,     3,
       5,     5,     3,     3,     2,     2,     1,     1,     3,     2,
       2,     2,     2,     2,     2,     2,     2,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     2,     2,     4,     4,     0,
       1,     3,     2,     4,     3,     4,     1,     3,     3,     3,
       0,     1,     1,     3,     1,     3,     6,     8,     6,     8,
       0,     1,     1,     3,     1,     3,     2,     4,     2,     3,
       3,     3,     3,     0,     1,     2,     3,     3,     1,     0,
       1,     3,     2,     3,     2,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     3,     0,     1,     3,
       2,     3,     2,     1,     2,     1,     1,     7,     3,     5,
       3,     3,     5,     3,     1,     1,     1,     2,     1,     3,
       3,     1,     1,     1,     1,     1,     1,     1,     1,     1
};


enum { YYENOMEM = -2 };

#define yyerrok         (yyerrstatus = 0)
#define yyclearin       (yychar = YYEMPTY)

#define YYACCEPT        goto yyacceptlab
#define YYABORT         goto yyabortlab
#define YYERROR         goto yyerrorlab
#define YYNOMEM         goto yyexhaustedlab


#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)                                    \
  do                                                              \
    if (yychar == YYEMPTY)                                        \
      {                                                           \
        yychar = (Token);                                         \
        yylval = (Value);                                         \
        YYPOPSTACK (yylen);                                       \
        yystate = *yyssp;                                         \
        goto yybackup;                                            \
      }                                                           \
    else                                                          \
      {                                                           \
        yyerror (YY_("syntax error: cannot back up")); \
        YYERROR;                                                  \
      }                                                           \
  while (0)

/* Backward compatibility with an undocumented macro.
   Use YYerror or YYUNDEF. */
#define YYERRCODE YYUNDEF

/* YYLLOC_DEFAULT -- Set CURRENT to span from RHS[1] to RHS[N].
   If N is 0, then set CURRENT to the empty location which ends
   the previous symbol: RHS[0] (always defined).  */

#ifndef YYLLOC_DEFAULT
# define YYLLOC_DEFAULT(Current, Rhs, N)                                \
    do                                                                  \
      if (N)                                                            \
        {                                                               \
          (Current).first_line   = YYRHSLOC (Rhs, 1).first_line;        \
          (Current).first_column = YYRHSLOC (Rhs, 1).first_column;      \
          (Current).last_line    = YYRHSLOC (Rhs, N).last_line;         \
          (Current).last_column  = YYRHSLOC (Rhs, N).last_column;       \
        }                                                               \
      else                                                              \
        {                                                               \
          (Current).first_line   = (Current).last_line   =              \
            YYRHSLOC (Rhs, 0).last_line;                                \
          (Current).first_column = (Current).last_column =              \
            YYRHSLOC (Rhs, 0).last_column;                              \
        }                                                               \
    while (0)
#endif

#define YYRHSLOC(Rhs, K) ((Rhs)[K])


/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)                        \
do {                                            \
  if (yydebug)                                  \
    YYFPRINTF Args;                             \
} while (0)


/* YYLOCATION_PRINT -- Print the location on the stream.
   This macro was not mandated originally: define only if we know
   we won't break user code: when these are the locations we know.  */

# ifndef YYLOCATION_PRINT

#  if defined YY_LOCATION_PRINT

   /* Temporary convenience wrapper in case some people defined the
      undocumented and private YY_LOCATION_PRINT macros.  */
#   define YYLOCATION_PRINT(File, Loc)  YY_LOCATION_PRINT(File, *(Loc))

#  elif defined YYLTYPE_IS_TRIVIAL && YYLTYPE_IS_TRIVIAL

/* Print *YYLOCP on YYO.  Private, do not rely on its existence. */

YY_ATTRIBUTE_UNUSED
static int
yy_location_print_ (FILE *yyo, YYLTYPE const * const yylocp)
{
  int res = 0;
  int end_col = 0 != yylocp->last_column ? yylocp->last_column - 1 : 0;
  if (0 <= yylocp->first_line)
    {
      res += YYFPRINTF (yyo, "%d", yylocp->first_line);
      if (0 <= yylocp->first_column)
        res += YYFPRINTF (yyo, ".%d", yylocp->first_column);
    }
  if (0 <= yylocp->last_line)
    {
      if (yylocp->first_line < yylocp->last_line)
        {
          res += YYFPRINTF (yyo, "-%d", yylocp->last_line);
          if (0 <= end_col)
            res += YYFPRINTF (yyo, ".%d", end_col);
        }
      else if (0 <= end_col && yylocp->first_column < end_col)
        res += YYFPRINTF (yyo, "-%d", end_col);
    }
  return res;
}

#   define YYLOCATION_PRINT  yy_location_print_

    /* Temporary convenience wrapper in case some people defined the
       undocumented and private YY_LOCATION_PRINT macros.  */
#   define YY_LOCATION_PRINT(File, Loc)  YYLOCATION_PRINT(File, &(Loc))

#  else

#   define YYLOCATION_PRINT(File, Loc) ((void) 0)
    /* Temporary convenience wrapper in case some people defined the
       undocumented and private YY_LOCATION_PRINT macros.  */
#   define YY_LOCATION_PRINT  YYLOCATION_PRINT

#  endif
# endif /* !defined YYLOCATION_PRINT */


# define YY_SYMBOL_PRINT(Title, Kind, Value, Location)                    \
do {                                                                      \
  if (yydebug)                                                            \
    {                                                                     \
      YYFPRINTF (stderr, "%s ", Title);                                   \
      yy_symbol_print (stderr,                                            \
                  Kind, Value, Location); \
      YYFPRINTF (stderr, "\n");                                           \
    }                                                                     \
} while (0)


/*-----------------------------------.
| Print this symbol's value on YYO.  |
`-----------------------------------*/

static void
yy_symbol_value_print (FILE *yyo,
                       yysymbol_kind_t yykind, YYSTYPE const * const yyvaluep, YYLTYPE const * const yylocationp)
{
  FILE *yyoutput = yyo;
  YY_USE (yyoutput);
  YY_USE (yylocationp);
  if (!yyvaluep)
    return;
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  YY_USE (yykind);
  YY_IGNORE_MAYBE_UNINITIALIZED_END
}


/*---------------------------.
| Print this symbol on YYO.  |
`---------------------------*/

static void
yy_symbol_print (FILE *yyo,
                 yysymbol_kind_t yykind, YYSTYPE const * const yyvaluep, YYLTYPE const * const yylocationp)
{
  YYFPRINTF (yyo, "%s %s (",
             yykind < YYNTOKENS ? "token" : "nterm", yysymbol_name (yykind));

  YYLOCATION_PRINT (yyo, yylocationp);
  YYFPRINTF (yyo, ": ");
  yy_symbol_value_print (yyo, yykind, yyvaluep, yylocationp);
  YYFPRINTF (yyo, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

static void
yy_stack_print (yy_state_t *yybottom, yy_state_t *yytop)
{
  YYFPRINTF (stderr, "Stack now");
  for (; yybottom <= yytop; yybottom++)
    {
      int yybot = *yybottom;
      YYFPRINTF (stderr, " %d", yybot);
    }
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)                            \
do {                                                            \
  if (yydebug)                                                  \
    yy_stack_print ((Bottom), (Top));                           \
} while (0)


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

static void
yy_reduce_print (yy_state_t *yyssp, YYSTYPE *yyvsp, YYLTYPE *yylsp,
                 int yyrule)
{
  int yylno = yyrline[yyrule];
  int yynrhs = yyr2[yyrule];
  int yyi;
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %d):\n",
             yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      YYFPRINTF (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr,
                       YY_ACCESSING_SYMBOL (+yyssp[yyi + 1 - yynrhs]),
                       &yyvsp[(yyi + 1) - (yynrhs)],
                       &(yylsp[(yyi + 1) - (yynrhs)]));
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)          \
do {                                    \
  if (yydebug)                          \
    yy_reduce_print (yyssp, yyvsp, yylsp, Rule); \
} while (0)

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args) ((void) 0)
# define YY_SYMBOL_PRINT(Title, Kind, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif






/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

static void
yydestruct (const char *yymsg,
            yysymbol_kind_t yykind, YYSTYPE *yyvaluep, YYLTYPE *yylocationp)
{
  YY_USE (yyvaluep);
  YY_USE (yylocationp);
  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yykind, yyvaluep, yylocationp);

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  switch (yykind)
    {
    case YYSYMBOL_expression: /* expression  */
#line 144 "./grammar.y"
            { delete (((*yyvaluep).syntaxElement)); }
#line 1452 "./grammar.tab.c"
        break;

    case YYSYMBOL_arrowAssign: /* arrowAssign  */
#line 145 "./grammar.y"
            { delete (((*yyvaluep).syntaxElement)); }
#line 1458 "./grammar.tab.c"
        break;

    case YYSYMBOL_tildeAssign: /* tildeAssign  */
#line 145 "./grammar.y"
            { delete (((*yyvaluep).syntaxElement)); }
#line 1464 "./grammar.tab.c"
        break;

    case YYSYMBOL_equationAssign: /* equationAssign  */
#line 145 "./grammar.y"
            { delete (((*yyvaluep).syntaxElement)); }
#line 1470 "./grammar.tab.c"
        break;

    case YYSYMBOL_workspaceAssign: /* workspaceAssign  */
#line 145 "./grammar.y"
            { delete (((*yyvaluep).syntaxElement)); }
#line 1476 "./grammar.tab.c"
        break;

    case YYSYMBOL_referenceAssign: /* referenceAssign  */
#line 146 "./grammar.y"
            { delete (((*yyvaluep).syntaxElement)); }
#line 1482 "./grammar.tab.c"
        break;

    case YYSYMBOL_additionAssign: /* additionAssign  */
#line 147 "./grammar.y"
            { delete (((*yyvaluep).syntaxElement)); }
#line 1488 "./grammar.tab.c"
        break;

    case YYSYMBOL_subtractionAssign: /* subtractionAssign  */
#line 147 "./grammar.y"
            { delete (((*yyvaluep).syntaxElement)); }
#line 1494 "./grammar.tab.c"
        break;

    case YYSYMBOL_multiplicationAssign: /* multiplicationAssign  */
#line 147 "./grammar.y"
            { delete (((*yyvaluep).syntaxElement)); }
#line 1500 "./grammar.tab.c"
        break;

    case YYSYMBOL_divisionAssign: /* divisionAssign  */
#line 147 "./grammar.y"
            { delete (((*yyvaluep).syntaxElement)); }
#line 1506 "./grammar.tab.c"
        break;

    case YYSYMBOL_variable: /* variable  */
#line 143 "./grammar.y"
            { delete (((*yyvaluep).syntaxElement)); }
#line 1512 "./grammar.tab.c"
        break;

    case YYSYMBOL_optElements: /* optElements  */
#line 139 "./grammar.y"
            { for (std::list<SyntaxElement*>::iterator it=((*yyvaluep).syntaxElementList)->begin(); it != ((*yyvaluep).syntaxElementList)->end(); it++) { SyntaxElement* theElement = *it; delete theElement; }; delete (((*yyvaluep).syntaxElementList)); }
#line 1518 "./grammar.tab.c"
        break;

    case YYSYMBOL_elementList: /* elementList  */
#line 139 "./grammar.y"
            { for (std::list<SyntaxElement*>::iterator it=((*yyvaluep).syntaxElementList)->begin(); it != ((*yyvaluep).syntaxElementList)->end(); it++) { SyntaxElement* theElement = *it; delete theElement; }; delete (((*yyvaluep).syntaxElementList)); }
#line 1524 "./grammar.tab.c"
        break;

    case YYSYMBOL_fxnCall: /* fxnCall  */
#line 143 "./grammar.y"
            { delete (((*yyvaluep).syntaxFunctionCall)); }
#line 1530 "./grammar.tab.c"
        break;

    case YYSYMBOL_functionCall: /* functionCall  */
#line 143 "./grammar.y"
            { delete (((*yyvaluep).syntaxFunctionCall)); }
#line 1536 "./grammar.tab.c"
        break;

    case YYSYMBOL_optArguments: /* optArguments  */
#line 140 "./grammar.y"
            { for (std::list<SyntaxLabeledExpr*>::iterator it=((*yyvaluep).argumentList)->begin(); it != ((*yyvaluep).argumentList)->end(); it++) { SyntaxLabeledExpr* theElement = *it; delete theElement; }; delete (((*yyvaluep).argumentList)); }
#line 1542 "./grammar.tab.c"
        break;

    case YYSYMBOL_argumentList: /* argumentList  */
#line 140 "./grammar.y"
            { for (std::list<SyntaxLabeledExpr*>::iterator it=((*yyvaluep).argumentList)->begin(); it != ((*yyvaluep).argumentList)->end(); it++) { SyntaxLabeledExpr* theElement = *it; delete theElement; }; delete (((*yyvaluep).argumentList)); }
#line 1548 "./grammar.tab.c"
        break;

    case YYSYMBOL_argument: /* argument  */
#line 143 "./grammar.y"
            { delete (((*yyvaluep).syntaxLabeledExpr)); }
#line 1554 "./grammar.tab.c"
        break;

    case YYSYMBOL_functionDef: /* functionDef  */
#line 149 "./grammar.y"
            { delete (((*yyvaluep).syntaxElement)); }
#line 1560 "./grammar.tab.c"
        break;

    case YYSYMBOL_procedureDef: /* procedureDef  */
#line 149 "./grammar.y"
            { delete (((*yyvaluep).syntaxElement)); }
#line 1566 "./grammar.tab.c"
        break;

    case YYSYMBOL_optFormals: /* optFormals  */
#line 141 "./grammar.y"
            { for (std::list<SyntaxFormal*>::iterator it=((*yyvaluep).formalList)->begin(); it != ((*yyvaluep).formalList)->end(); it++) { SyntaxFormal* theElement = *it; delete theElement; }; delete (((*yyvaluep).formalList)); }
#line 1572 "./grammar.tab.c"
        break;

    case YYSYMBOL_formalList: /* formalList  */
#line 141 "./grammar.y"
            { for (std::list<SyntaxFormal*>::iterator it=((*yyvaluep).formalList)->begin(); it != ((*yyvaluep).formalList)->end(); it++) { SyntaxFormal* theElement = *it; delete theElement; }; delete (((*yyvaluep).formalList)); }
#line 1578 "./grammar.tab.c"
        break;

    case YYSYMBOL_formal: /* formal  */
#line 143 "./grammar.y"
            { delete (((*yyvaluep).syntaxFormal)); }
#line 1584 "./grammar.tab.c"
        break;

    case YYSYMBOL_typeSpec: /* typeSpec  */
#line 142 "./grammar.y"
            { delete (((*yyvaluep).string)); }
#line 1590 "./grammar.tab.c"
        break;

    case YYSYMBOL_optDims: /* optDims  */
#line 142 "./grammar.y"
            { delete (((*yyvaluep).string)); }
#line 1596 "./grammar.tab.c"
        break;

    case YYSYMBOL_dimList: /* dimList  */
#line 142 "./grammar.y"
            { delete (((*yyvaluep).string)); }
#line 1602 "./grammar.tab.c"
        break;

    case YYSYMBOL_stmts: /* stmts  */
#line 139 "./grammar.y"
            { for (std::list<SyntaxElement*>::iterator it=((*yyvaluep).syntaxElementList)->begin(); it != ((*yyvaluep).syntaxElementList)->end(); it++) { SyntaxElement* theElement = *it; delete theElement; }; delete (((*yyvaluep).syntaxElementList)); }
#line 1608 "./grammar.tab.c"
        break;

    case YYSYMBOL_stmtList: /* stmtList  */
#line 139 "./grammar.y"
            { for (std::list<SyntaxElement*>::iterator it=((*yyvaluep).syntaxElementList)->begin(); it != ((*yyvaluep).syntaxElementList)->end(); it++) { SyntaxElement* theElement = *it; delete theElement; }; delete (((*yyvaluep).syntaxElementList)); }
#line 1614 "./grammar.tab.c"
        break;

    case YYSYMBOL_statement: /* statement  */
#line 144 "./grammar.y"
            { delete (((*yyvaluep).syntaxElement)); }
#line 1620 "./grammar.tab.c"
        break;

    case YYSYMBOL_stmt_or_expr: /* stmt_or_expr  */
#line 144 "./grammar.y"
            { delete (((*yyvaluep).syntaxElement)); }
#line 1626 "./grammar.tab.c"
        break;

    case YYSYMBOL_declaration: /* declaration  */
#line 148 "./grammar.y"
            { delete (((*yyvaluep).syntaxElement)); }
#line 1632 "./grammar.tab.c"
        break;

    case YYSYMBOL_memberDefs: /* memberDefs  */
#line 139 "./grammar.y"
            { for (std::list<SyntaxElement*>::iterator it=((*yyvaluep).syntaxElementList)->begin(); it != ((*yyvaluep).syntaxElementList)->end(); it++) { SyntaxElement* theElement = *it; delete theElement; }; delete (((*yyvaluep).syntaxElementList)); }
#line 1638 "./grammar.tab.c"
        break;

    case YYSYMBOL_memberDef: /* memberDef  */
#line 148 "./grammar.y"
            { delete (((*yyvaluep).syntaxElement)); }
#line 1644 "./grammar.tab.c"
        break;

    case YYSYMBOL_classDef: /* classDef  */
#line 148 "./grammar.y"
            { delete (((*yyvaluep).syntaxElement)); }
#line 1650 "./grammar.tab.c"
        break;

    case YYSYMBOL_ifStatement: /* ifStatement  */
#line 150 "./grammar.y"
            { delete (((*yyvaluep).syntaxElement)); }
#line 1656 "./grammar.tab.c"
        break;

    case YYSYMBOL_cond: /* cond  */
#line 151 "./grammar.y"
            { delete (((*yyvaluep).syntaxElement)); }
#line 1662 "./grammar.tab.c"
        break;

    case YYSYMBOL_forStatement: /* forStatement  */
#line 150 "./grammar.y"
            { delete (((*yyvaluep).syntaxElement)); }
#line 1668 "./grammar.tab.c"
        break;

    case YYSYMBOL_forCond: /* forCond  */
#line 151 "./grammar.y"
            { delete (((*yyvaluep).syntaxElement)); }
#line 1674 "./grammar.tab.c"
        break;

    case YYSYMBOL_whileStatement: /* whileStatement  */
#line 150 "./grammar.y"
            { delete (((*yyvaluep).syntaxElement)); }
#line 1680 "./grammar.tab.c"
        break;

    case YYSYMBOL_nextStatement: /* nextStatement  */
#line 152 "./grammar.y"
            { delete (((*yyvaluep).syntaxElement)); }
#line 1686 "./grammar.tab.c"
        break;

    case YYSYMBOL_breakStatement: /* breakStatement  */
#line 152 "./grammar.y"
            { delete (((*yyvaluep).syntaxElement)); }
#line 1692 "./grammar.tab.c"
        break;

    case YYSYMBOL_returnStatement: /* returnStatement  */
#line 151 "./grammar.y"
            { delete (((*yyvaluep).syntaxElement)); }
#line 1698 "./grammar.tab.c"
        break;

    case YYSYMBOL_identifier: /* identifier  */
#line 142 "./grammar.y"
            { delete (((*yyvaluep).string)); }
#line 1704 "./grammar.tab.c"
        break;

    case YYSYMBOL_vector: /* vector  */
#line 140 "./grammar.y"
            { for (std::list<SyntaxLabeledExpr*>::iterator it=((*yyvaluep).argumentList)->begin(); it != ((*yyvaluep).argumentList)->end(); it++) { SyntaxLabeledExpr* theElement = *it; delete theElement; }; delete (((*yyvaluep).argumentList)); }
#line 1710 "./grammar.tab.c"
        break;

    case YYSYMBOL_vectorList: /* vectorList  */
#line 140 "./grammar.y"
            { for (std::list<SyntaxLabeledExpr*>::iterator it=((*yyvaluep).argumentList)->begin(); it != ((*yyvaluep).argumentList)->end(); it++) { SyntaxLabeledExpr* theElement = *it; delete theElement; }; delete (((*yyvaluep).argumentList)); }
#line 1716 "./grammar.tab.c"
        break;

    case YYSYMBOL_constant: /* constant  */
#line 143 "./grammar.y"
            { delete (((*yyvaluep).syntaxElement)); }
#line 1722 "./grammar.tab.c"
        break;

      default:
        break;
    }
  YY_IGNORE_MAYBE_UNINITIALIZED_END
}


/* Lookahead token kind.  */
int yychar;

/* The semantic value of the lookahead symbol.  */
YYSTYPE yylval;
/* Location data for the lookahead symbol.  */
YYLTYPE yylloc
# if defined YYLTYPE_IS_TRIVIAL && YYLTYPE_IS_TRIVIAL
  = { 1, 1, 1, 1 }
# endif
;
/* Number of syntax errors so far.  */
int yynerrs;




/*----------.
| yyparse.  |
`----------*/

int
yyparse (void)
{
    yy_state_fast_t yystate = 0;
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus = 0;

    /* Refer to the stacks through separate pointers, to allow yyoverflow
       to reallocate them elsewhere.  */

    /* Their size.  */
    YYPTRDIFF_T yystacksize = YYINITDEPTH;

    /* The state stack: array, bottom, top.  */
    yy_state_t yyssa[YYINITDEPTH];
    yy_state_t *yyss = yyssa;
    yy_state_t *yyssp = yyss;

    /* The semantic value stack: array, bottom, top.  */
    YYSTYPE yyvsa[YYINITDEPTH];
    YYSTYPE *yyvs = yyvsa;
    YYSTYPE *yyvsp = yyvs;

    /* The location stack: array, bottom, top.  */
    YYLTYPE yylsa[YYINITDEPTH];
    YYLTYPE *yyls = yylsa;
    YYLTYPE *yylsp = yyls;

  int yyn;
  /* The return value of yyparse.  */
  int yyresult;
  /* Lookahead symbol kind.  */
  yysymbol_kind_t yytoken = YYSYMBOL_YYEMPTY;
  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;
  YYLTYPE yyloc;

  /* The locations where the error started and ended.  */
  YYLTYPE yyerror_range[3];



#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N), yylsp -= (N))

  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yychar = YYEMPTY; /* Cause a token to be read.  */

  yylsp[0] = yylloc;
  goto yysetstate;


/*------------------------------------------------------------.
| yynewstate -- push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;


/*--------------------------------------------------------------------.
| yysetstate -- set current state (the top of the stack) to yystate.  |
`--------------------------------------------------------------------*/
yysetstate:
  YYDPRINTF ((stderr, "Entering state %d\n", yystate));
  YY_ASSERT (0 <= yystate && yystate < YYNSTATES);
  YY_IGNORE_USELESS_CAST_BEGIN
  *yyssp = YY_CAST (yy_state_t, yystate);
  YY_IGNORE_USELESS_CAST_END
  YY_STACK_PRINT (yyss, yyssp);

  if (yyss + yystacksize - 1 <= yyssp)
#if !defined yyoverflow && !defined YYSTACK_RELOCATE
    YYNOMEM;
#else
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYPTRDIFF_T yysize = yyssp - yyss + 1;

# if defined yyoverflow
      {
        /* Give user a chance to reallocate the stack.  Use copies of
           these so that the &'s don't force the real ones into
           memory.  */
        yy_state_t *yyss1 = yyss;
        YYSTYPE *yyvs1 = yyvs;
        YYLTYPE *yyls1 = yyls;

        /* Each stack pointer address is followed by the size of the
           data in use in that stack, in bytes.  This used to be a
           conditional around just the two extra args, but that might
           be undefined if yyoverflow is a macro.  */
        yyoverflow (YY_("memory exhausted"),
                    &yyss1, yysize * YYSIZEOF (*yyssp),
                    &yyvs1, yysize * YYSIZEOF (*yyvsp),
                    &yyls1, yysize * YYSIZEOF (*yylsp),
                    &yystacksize);
        yyss = yyss1;
        yyvs = yyvs1;
        yyls = yyls1;
      }
# else /* defined YYSTACK_RELOCATE */
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
        YYNOMEM;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
        yystacksize = YYMAXDEPTH;

      {
        yy_state_t *yyss1 = yyss;
        union yyalloc *yyptr =
          YY_CAST (union yyalloc *,
                   YYSTACK_ALLOC (YY_CAST (YYSIZE_T, YYSTACK_BYTES (yystacksize))));
        if (! yyptr)
          YYNOMEM;
        YYSTACK_RELOCATE (yyss_alloc, yyss);
        YYSTACK_RELOCATE (yyvs_alloc, yyvs);
        YYSTACK_RELOCATE (yyls_alloc, yyls);
#  undef YYSTACK_RELOCATE
        if (yyss1 != yyssa)
          YYSTACK_FREE (yyss1);
      }
# endif

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;
      yylsp = yyls + yysize - 1;

      YY_IGNORE_USELESS_CAST_BEGIN
      YYDPRINTF ((stderr, "Stack size increased to %ld\n",
                  YY_CAST (long, yystacksize)));
      YY_IGNORE_USELESS_CAST_END

      if (yyss + yystacksize - 1 <= yyssp)
        YYABORT;
    }
#endif /* !defined yyoverflow && !defined YYSTACK_RELOCATE */


  if (yystate == YYFINAL)
    YYACCEPT;

  goto yybackup;


/*-----------.
| yybackup.  |
`-----------*/
yybackup:
  /* Do appropriate processing given the current state.  Read a
     lookahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to lookahead token.  */
  yyn = yypact[yystate];
  if (yypact_value_is_default (yyn))
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either empty, or end-of-input, or a valid lookahead.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token\n"));
      yychar = yylex ();
    }

  if (yychar <= YYEOF)
    {
      yychar = YYEOF;
      yytoken = YYSYMBOL_YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else if (yychar == YYerror)
    {
      /* The scanner already issued an error message, process directly
         to error recovery.  But do not keep the error token as
         lookahead, it is too special and may lead us to an endless
         loop in error recovery. */
      yychar = YYUNDEF;
      yytoken = YYSYMBOL_YYerror;
      yyerror_range[1] = yylloc;
      goto yyerrlab1;
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yytable_value_is_error (yyn))
        goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the lookahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);
  yystate = yyn;
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END
  *++yylsp = yylloc;

  /* Discard the shifted token.  */
  yychar = YYEMPTY;
  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     '$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];

  /* Default location. */
  YYLLOC_DEFAULT (yyloc, (yylsp - yylen), yylen);
  yyerror_range[1] = yyloc;
  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
  case 2: /* prog: END_OF_INPUT  */
#line 238 "./grammar.y"
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Bison encountered end_of_input; ignored\n");
#endif
                    return 0;
                }
#line 2022 "./grammar.tab.c"
    break;

  case 3: /* prog: '\n'  */
#line 245 "./grammar.y"
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Bison encountered newline; ignored\n");
#endif
                    return 0;
                }
#line 2033 "./grammar.tab.c"
    break;

  case 4: /* prog: stmt_or_expr '\n'  */
#line 252 "./grammar.y"
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Bison trying to execute statement or expression\n");
#endif
                    int rv = parser.execute((yyvsp[-1].syntaxElement), *executionEnvironment);
                    delete (yyvsp[-1].syntaxElement);
                    return rv;
                }
#line 2046 "./grammar.tab.c"
    break;

  case 5: /* prog: stmt_or_expr ';'  */
#line 261 "./grammar.y"
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Bison trying to execute statement or expression\n");
#endif
                    int rv =  parser.execute((yyvsp[-1].syntaxElement), *executionEnvironment);
                    delete (yyvsp[-1].syntaxElement);
                    return rv;
                }
#line 2059 "./grammar.tab.c"
    break;

  case 6: /* prog: declaration '\n'  */
#line 270 "./grammar.y"
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Bison trying to execute declaration\n");
#endif
                    int rv =  parser.execute((yyvsp[-1].syntaxElement), *executionEnvironment);
                    delete (yyvsp[-1].syntaxElement);
                    return rv;
                }
#line 2072 "./grammar.tab.c"
    break;

  case 7: /* prog: declaration ';'  */
#line 279 "./grammar.y"
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Bison trying to execute declaration\n");
#endif
                    int rv =  parser.execute((yyvsp[-1].syntaxElement), *executionEnvironment);
                    delete (yyvsp[-1].syntaxElement);
                    return rv;
                }
#line 2085 "./grammar.tab.c"
    break;

  case 8: /* prog: '?' identifier '\n'  */
#line 288 "./grammar.y"
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Bison trying to get help for symbol\n");
#endif
                    int rv =  parser.help(*(yyvsp[-1].string));
                    delete (yyvsp[-1].string);
                    return rv;
                }
#line 2098 "./grammar.tab.c"
    break;

  case 9: /* prog: '?' identifier ';'  */
#line 297 "./grammar.y"
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Bison trying to get help for symbol\n");
#endif
                    int rv =  parser.help(*(yyvsp[-1].string));
                    delete (yyvsp[-1].string);
                    return rv;
                }
#line 2111 "./grammar.tab.c"
    break;

  case 10: /* prog: '?' identifier '.' identifier '\n'  */
#line 306 "./grammar.y"
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Bison trying to get help for symbol\n");
#endif
                    int rv =  parser.help(*(yyvsp[-3].string), *(yyvsp[-1].string));
                    delete (yyvsp[-3].string);
                    delete (yyvsp[-1].string);
                    return rv;
                }
#line 2125 "./grammar.tab.c"
    break;

  case 11: /* prog: '?' identifier '.' identifier ';'  */
#line 316 "./grammar.y"
                {
#ifdef DEBUG_BISON_FLEX
                printf("Bison trying to get help for symbol\n");
#endif
                int rv =  parser.help(*(yyvsp[-3].string), *(yyvsp[-1].string));
                delete (yyvsp[-3].string);
                delete (yyvsp[-1].string);
                return rv;
                }
#line 2139 "./grammar.tab.c"
    break;

  case 12: /* prog: '?' functionCall '\n'  */
#line 326 "./grammar.y"
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Bison trying to get help for function call\n");
#endif
                    int rv =  parser.help((yyvsp[-1].syntaxFunctionCall));
                    delete (yyvsp[-1].syntaxFunctionCall);
                    return rv;
                }
#line 2152 "./grammar.tab.c"
    break;

  case 13: /* prog: '?' functionCall ';'  */
#line 335 "./grammar.y"
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Bison trying to get help for function call\n");
#endif
                    int rv =  parser.help((yyvsp[-1].syntaxFunctionCall));
                    delete (yyvsp[-1].syntaxFunctionCall);
                    return rv;
                }
#line 2165 "./grammar.tab.c"
    break;

  case 14: /* prog: error '\n'  */
#line 344 "./grammar.y"
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Bison error when reading line %d position %d to line %d position %d\n", (yyloc).first_line, (yyloc).first_column, (yyloc).last_line, (yyloc).last_column);
#endif
                    YYABORT;
                }
#line 2176 "./grammar.tab.c"
    break;

  case 15: /* prog: error ';'  */
#line 351 "./grammar.y"
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Bison error when reading line %d position %d to line %d position %d\n", (yyloc).first_line, (yyloc).first_column, (yyloc).last_line, (yyloc).last_column);
#endif
                    YYABORT;
                }
#line 2187 "./grammar.tab.c"
    break;

  case 16: /* expression: constant  */
#line 359 "./grammar.y"
                                            { (yyval.syntaxElement) = (yyvsp[0].syntaxElement); }
#line 2193 "./grammar.tab.c"
    break;

  case 17: /* expression: vector  */
#line 361 "./grammar.y"
                                            { (yyval.syntaxElement) = new SyntaxFunctionCall("v", (yyvsp[0].argumentList)); }
#line 2199 "./grammar.tab.c"
    break;

  case 18: /* expression: '(' expression ')'  */
#line 363 "./grammar.y"
                                            { (yyval.syntaxElement) = (yyvsp[-1].syntaxElement); }
#line 2205 "./grammar.tab.c"
    break;

  case 19: /* expression: '-' expression  */
#line 365 "./grammar.y"
                                            { (yyval.syntaxElement) = new SyntaxUnaryExpr(SyntaxUnaryExpr::UMinus, (yyvsp[0].syntaxElement)); }
#line 2211 "./grammar.tab.c"
    break;

  case 20: /* expression: '+' expression  */
#line 366 "./grammar.y"
                                            { (yyval.syntaxElement) = new SyntaxUnaryExpr(SyntaxUnaryExpr::UPlus, (yyvsp[0].syntaxElement)); }
#line 2217 "./grammar.tab.c"
    break;

  case 21: /* expression: '!' expression  */
#line 367 "./grammar.y"
                                            { (yyval.syntaxElement) = new SyntaxUnaryExpr(SyntaxUnaryExpr::UNot, (yyvsp[0].syntaxElement)); }
#line 2223 "./grammar.tab.c"
    break;

  case 22: /* expression: AND expression  */
#line 368 "./grammar.y"
                                            { (yyval.syntaxElement) = new SyntaxUnaryExpr(SyntaxUnaryExpr::UAnd, (yyvsp[0].syntaxElement)); }
#line 2229 "./grammar.tab.c"
    break;

  case 23: /* expression: DECREMENT variable  */
#line 370 "./grammar.y"
                                            { (yyval.syntaxElement) = new SyntaxDecrement( (yyvsp[0].syntaxElement), false ); }
#line 2235 "./grammar.tab.c"
    break;

  case 24: /* expression: variable DECREMENT  */
#line 371 "./grammar.y"
                                            { (yyval.syntaxElement) = new SyntaxDecrement( (yyvsp[-1].syntaxElement), true ); }
#line 2241 "./grammar.tab.c"
    break;

  case 25: /* expression: INCREMENT variable  */
#line 372 "./grammar.y"
                                            { (yyval.syntaxElement) = new SyntaxIncrement( (yyvsp[0].syntaxElement), false ); }
#line 2247 "./grammar.tab.c"
    break;

  case 26: /* expression: variable INCREMENT  */
#line 373 "./grammar.y"
                                            { (yyval.syntaxElement) = new SyntaxIncrement( (yyvsp[-1].syntaxElement), true ); }
#line 2253 "./grammar.tab.c"
    break;

  case 27: /* expression: expression ':' expression  */
#line 375 "./grammar.y"
                                            { (yyval.syntaxElement) = new SyntaxBinaryExpr(SyntaxBinaryExpr::Range, (yyvsp[-2].syntaxElement), (yyvsp[0].syntaxElement)); }
#line 2259 "./grammar.tab.c"
    break;

  case 28: /* expression: expression '+' expression  */
#line 377 "./grammar.y"
                                            { (yyval.syntaxElement) = new SyntaxBinaryExpr(SyntaxBinaryExpr::Add, (yyvsp[-2].syntaxElement), (yyvsp[0].syntaxElement)); }
#line 2265 "./grammar.tab.c"
    break;

  case 29: /* expression: expression '-' expression  */
#line 378 "./grammar.y"
                                            { (yyval.syntaxElement) = new SyntaxBinaryExpr(SyntaxBinaryExpr::Sub, (yyvsp[-2].syntaxElement), (yyvsp[0].syntaxElement)); }
#line 2271 "./grammar.tab.c"
    break;

  case 30: /* expression: expression '*' expression  */
#line 379 "./grammar.y"
                                            { (yyval.syntaxElement) = new SyntaxBinaryExpr(SyntaxBinaryExpr::Mul, (yyvsp[-2].syntaxElement), (yyvsp[0].syntaxElement)); }
#line 2277 "./grammar.tab.c"
    break;

  case 31: /* expression: expression '/' expression  */
#line 380 "./grammar.y"
                                            { (yyval.syntaxElement) = new SyntaxBinaryExpr(SyntaxBinaryExpr::Div, (yyvsp[-2].syntaxElement), (yyvsp[0].syntaxElement)); }
#line 2283 "./grammar.tab.c"
    break;

  case 32: /* expression: expression '^' expression  */
#line 381 "./grammar.y"
                                            { (yyval.syntaxElement) = new SyntaxBinaryExpr(SyntaxBinaryExpr::Exp, (yyvsp[-2].syntaxElement), (yyvsp[0].syntaxElement)); }
#line 2289 "./grammar.tab.c"
    break;

  case 33: /* expression: expression '%' expression  */
#line 382 "./grammar.y"
                                            { (yyval.syntaxElement) = new SyntaxBinaryExpr(SyntaxBinaryExpr::Mod, (yyvsp[-2].syntaxElement), (yyvsp[0].syntaxElement)); }
#line 2295 "./grammar.tab.c"
    break;

  case 34: /* expression: expression LT expression  */
#line 384 "./grammar.y"
                                            { (yyval.syntaxElement) = new SyntaxBinaryExpr(SyntaxBinaryExpr::Lt, (yyvsp[-2].syntaxElement), (yyvsp[0].syntaxElement)); }
#line 2301 "./grammar.tab.c"
    break;

  case 35: /* expression: expression LE expression  */
#line 385 "./grammar.y"
                                            { (yyval.syntaxElement) = new SyntaxBinaryExpr(SyntaxBinaryExpr::Le, (yyvsp[-2].syntaxElement), (yyvsp[0].syntaxElement)); }
#line 2307 "./grammar.tab.c"
    break;

  case 36: /* expression: expression EQ expression  */
#line 386 "./grammar.y"
                                            { (yyval.syntaxElement) = new SyntaxBinaryExpr(SyntaxBinaryExpr::Eq, (yyvsp[-2].syntaxElement), (yyvsp[0].syntaxElement)); }
#line 2313 "./grammar.tab.c"
    break;

  case 37: /* expression: expression NE expression  */
#line 387 "./grammar.y"
                                            { (yyval.syntaxElement) = new SyntaxBinaryExpr(SyntaxBinaryExpr::Ne, (yyvsp[-2].syntaxElement), (yyvsp[0].syntaxElement)); }
#line 2319 "./grammar.tab.c"
    break;

  case 38: /* expression: expression GE expression  */
#line 388 "./grammar.y"
                                            { (yyval.syntaxElement) = new SyntaxBinaryExpr(SyntaxBinaryExpr::Ge, (yyvsp[-2].syntaxElement), (yyvsp[0].syntaxElement)); }
#line 2325 "./grammar.tab.c"
    break;

  case 39: /* expression: expression GT expression  */
#line 389 "./grammar.y"
                                            { (yyval.syntaxElement) = new SyntaxBinaryExpr(SyntaxBinaryExpr::Gt, (yyvsp[-2].syntaxElement), (yyvsp[0].syntaxElement)); }
#line 2331 "./grammar.tab.c"
    break;

  case 40: /* expression: expression AND expression  */
#line 391 "./grammar.y"
                                            { (yyval.syntaxElement) = new SyntaxBinaryExpr(SyntaxBinaryExpr::And, (yyvsp[-2].syntaxElement), (yyvsp[0].syntaxElement)); }
#line 2337 "./grammar.tab.c"
    break;

  case 41: /* expression: expression OR expression  */
#line 392 "./grammar.y"
                                            { (yyval.syntaxElement) = new SyntaxBinaryExpr(SyntaxBinaryExpr::Or, (yyvsp[-2].syntaxElement), (yyvsp[0].syntaxElement)); }
#line 2343 "./grammar.tab.c"
    break;

  case 42: /* expression: expression AND2 expression  */
#line 393 "./grammar.y"
                                            { (yyval.syntaxElement) = new SyntaxBinaryExpr(SyntaxBinaryExpr::And2, (yyvsp[-2].syntaxElement), (yyvsp[0].syntaxElement)); }
#line 2349 "./grammar.tab.c"
    break;

  case 43: /* expression: expression OR2 expression  */
#line 394 "./grammar.y"
                                            { (yyval.syntaxElement) = new SyntaxBinaryExpr(SyntaxBinaryExpr::Or2, (yyvsp[-2].syntaxElement), (yyvsp[0].syntaxElement)); }
#line 2355 "./grammar.tab.c"
    break;

  case 44: /* expression: PIPE_PLACEHOLDER  */
#line 396 "./grammar.y"
                                            { (yyval.syntaxElement) = new SyntaxPipePlaceholder(); }
#line 2361 "./grammar.tab.c"
    break;

  case 45: /* expression: arrowAssign  */
#line 398 "./grammar.y"
                                            { (yyval.syntaxElement) = (yyvsp[0].syntaxElement); }
#line 2367 "./grammar.tab.c"
    break;

  case 46: /* expression: equationAssign  */
#line 399 "./grammar.y"
                                            { (yyval.syntaxElement) = (yyvsp[0].syntaxElement); }
#line 2373 "./grammar.tab.c"
    break;

  case 47: /* expression: tildeAssign  */
#line 400 "./grammar.y"
                                            { (yyval.syntaxElement) = (yyvsp[0].syntaxElement); }
#line 2379 "./grammar.tab.c"
    break;

  case 48: /* expression: workspaceAssign  */
#line 401 "./grammar.y"
                                            { (yyval.syntaxElement) = (yyvsp[0].syntaxElement); }
#line 2385 "./grammar.tab.c"
    break;

  case 49: /* expression: referenceAssign  */
#line 402 "./grammar.y"
                                            { (yyval.syntaxElement) = (yyvsp[0].syntaxElement); }
#line 2391 "./grammar.tab.c"
    break;

  case 50: /* expression: additionAssign  */
#line 404 "./grammar.y"
                                            { (yyval.syntaxElement) = (yyvsp[0].syntaxElement); }
#line 2397 "./grammar.tab.c"
    break;

  case 51: /* expression: subtractionAssign  */
#line 405 "./grammar.y"
                                            { (yyval.syntaxElement) = (yyvsp[0].syntaxElement); }
#line 2403 "./grammar.tab.c"
    break;

  case 52: /* expression: multiplicationAssign  */
#line 406 "./grammar.y"
                                            { (yyval.syntaxElement) = (yyvsp[0].syntaxElement); }
#line 2409 "./grammar.tab.c"
    break;

  case 53: /* expression: divisionAssign  */
#line 407 "./grammar.y"
                                            { (yyval.syntaxElement) = (yyvsp[0].syntaxElement); }
#line 2415 "./grammar.tab.c"
    break;

  case 54: /* expression: functionCall  */
#line 409 "./grammar.y"
                                            { (yyval.syntaxElement) = (yyvsp[0].syntaxFunctionCall); }
#line 2421 "./grammar.tab.c"
    break;

  case 55: /* expression: variable  */
#line 411 "./grammar.y"
                                            { (yyval.syntaxElement) = (yyvsp[0].syntaxElement); }
#line 2427 "./grammar.tab.c"
    break;

  case 56: /* arrowAssign: expression ARROW_ASSIGN expression  */
#line 416 "./grammar.y"
                    { 
#ifdef DEBUG_BISON_FLEX
                        printf("Parser inserting arrow assignment (ARROW_ASSIGN) in syntax tree\n");
#endif
                        (yyval.syntaxElement) = new SyntaxConstantAssignment((yyvsp[-2].syntaxElement), (yyvsp[0].syntaxElement));
                    }
#line 2438 "./grammar.tab.c"
    break;

  case 57: /* tildeAssign: expression TILDE_ASSIGN expression  */
#line 425 "./grammar.y"
                    {
#ifdef DEBUG_BISON_FLEX
                        printf("Parser inserting tilde assignment (TILDE_ASSIGN) in syntax tree\n");
#endif
                        (yyval.syntaxElement) = new SyntaxStochasticAssignment((yyvsp[-2].syntaxElement), (yyvsp[0].syntaxElement));
                    }
#line 2449 "./grammar.tab.c"
    break;

  case 58: /* equationAssign: expression EQUATION_ASSIGN expression  */
#line 434 "./grammar.y"
                    {
#ifdef DEBUG_BISON_FLEX
                        printf("Parser inserting equation assignment (EQUATION_ASSIGN) in syntax tree\n");
#endif
                        (yyval.syntaxElement) = new SyntaxDeterministicAssignment((yyvsp[-2].syntaxElement), (yyvsp[0].syntaxElement)); 
                    }
#line 2460 "./grammar.tab.c"
    break;

  case 59: /* workspaceAssign: expression EQUAL expression  */
#line 443 "./grammar.y"
                    {
#ifdef DEBUG_BISON_FLEX
                        printf("Parser inserting workspace assignment (WORKSPACE_ASSIGN) in syntax tree\n");
#endif
                        (yyval.syntaxElement) = new SyntaxWorkspaceVariableAssignment((yyvsp[-2].syntaxElement), (yyvsp[0].syntaxElement));
                    }
#line 2471 "./grammar.tab.c"
    break;

  case 60: /* referenceAssign: expression REFERENCE_ASSIGN expression  */
#line 452 "./grammar.y"
                    {
#ifdef DEBUG_BISON_FLEX
                        printf("Parser inserting reference assignment (REFERENCE_ASSIGN) in syntax tree\n");
#endif
                        (yyval.syntaxElement) = new SyntaxReferenceAssignment((yyvsp[-2].syntaxElement), (yyvsp[0].syntaxElement));
                    }
#line 2482 "./grammar.tab.c"
    break;

  case 61: /* additionAssign: expression ADDITION_ASSIGN expression  */
#line 461 "./grammar.y"
                    {
#ifdef DEBUG_BISON_FLEX
                        printf("Parser inserting addition assignment (ADDITION_ASSIGN) in syntax tree\n");
#endif
                        (yyval.syntaxElement) = new SyntaxAdditionAssignment((yyvsp[-2].syntaxElement), (yyvsp[0].syntaxElement)); 
                    }
#line 2493 "./grammar.tab.c"
    break;

  case 62: /* subtractionAssign: expression SUBTRACTION_ASSIGN expression  */
#line 470 "./grammar.y"
                        {
#ifdef DEBUG_BISON_FLEX
                            printf("Parser inserting subtraction assignment (SUBTRACTION_ASSIGN) in syntax tree\n");
#endif
                            (yyval.syntaxElement) = new SyntaxSubtractionAssignment((yyvsp[-2].syntaxElement), (yyvsp[0].syntaxElement)); 
                        }
#line 2504 "./grammar.tab.c"
    break;

  case 63: /* multiplicationAssign: expression MULTIPLICATION_ASSIGN expression  */
#line 479 "./grammar.y"
                            {
#ifdef DEBUG_BISON_FLEX
                                printf("Parser inserting multiplication assignment (MULTIPLICATION_ASSIGN) in syntax tree\n");
#endif
                                (yyval.syntaxElement) = new SyntaxMultiplicationAssignment((yyvsp[-2].syntaxElement), (yyvsp[0].syntaxElement)); 
                            }
#line 2515 "./grammar.tab.c"
    break;

  case 64: /* divisionAssign: expression DIVISION_ASSIGN expression  */
#line 488 "./grammar.y"
                    {
#ifdef DEBUG_BISON_FLEX
                        printf("Parser inserting division assignment (DIVISION_ASSIGN) in syntax tree\n");
#endif
                        (yyval.syntaxElement) = new SyntaxDivisionAssignment((yyvsp[-2].syntaxElement), (yyvsp[0].syntaxElement)); 
                    }
#line 2526 "./grammar.tab.c"
    break;

  case 65: /* variable: identifier optElements  */
#line 497 "./grammar.y"
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Parser inserting variable (NAMED_VAR) in syntax tree\n");
#endif
                    (yyval.syntaxElement) = new SyntaxVariable(*(yyvsp[-1].string));
                    for (std::list<SyntaxElement*>::iterator it=(yyvsp[0].syntaxElementList)->begin(); it!=(yyvsp[0].syntaxElementList)->end(); ++it)
                    {
                        (yyval.syntaxElement) = new SyntaxIndexOperation((yyval.syntaxElement),*it);
                    }
                    delete (yyvsp[-1].string);
                    delete (yyvsp[0].syntaxElementList);
                }
#line 2543 "./grammar.tab.c"
    break;

  case 66: /* variable: fxnCall elementList  */
#line 510 "./grammar.y"
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Parser inserting variable (FUNCTION_VAR) in syntax tree\n");
#endif
                    (yyval.syntaxElement) = (yyvsp[-1].syntaxFunctionCall);
                    for (auto& element: *(yyvsp[0].syntaxElementList))
                    {
                        (yyval.syntaxElement) = new SyntaxIndexOperation((yyval.syntaxElement), element);
                    }
                    delete (yyvsp[0].syntaxElementList);
                }
#line 2559 "./grammar.tab.c"
    break;

  case 67: /* variable: '(' expression ')' elementList  */
#line 522 "./grammar.y"
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Parser inserting variable (EXPRESSION_VAR) in syntax tree\n");
#endif
                    (yyval.syntaxElement) = (yyvsp[-2].syntaxElement);
                    for (auto& element: *(yyvsp[0].syntaxElementList))
                    {
                        (yyval.syntaxElement) = new SyntaxIndexOperation((yyval.syntaxElement), element);
                    }
                    delete (yyvsp[0].syntaxElementList);
                }
#line 2575 "./grammar.tab.c"
    break;

  case 68: /* variable: variable '.' fxnCall elementList  */
#line 534 "./grammar.y"
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Parser inserting member variable (FUNCTION_VAR) in syntax tree\n");
#endif
                    (yyvsp[-1].syntaxFunctionCall)->setBaseVariable((yyvsp[-3].syntaxElement));
                    (yyval.syntaxElement) = (yyvsp[-1].syntaxFunctionCall);
                    for (auto& element: *(yyvsp[0].syntaxElementList))
                    {
                        (yyval.syntaxElement) = new SyntaxIndexOperation((yyval.syntaxElement), element);
                    }
                    delete (yyvsp[0].syntaxElementList);
                }
#line 2592 "./grammar.tab.c"
    break;

  case 69: /* optElements: %empty  */
#line 548 "./grammar.y"
                                                { (yyval.syntaxElementList) = new std::list<SyntaxElement*>(); }
#line 2598 "./grammar.tab.c"
    break;

  case 70: /* optElements: elementList  */
#line 549 "./grammar.y"
                                                { (yyval.syntaxElementList) = (yyvsp[0].syntaxElementList); }
#line 2604 "./grammar.tab.c"
    break;

  case 71: /* elementList: '[' expression ']'  */
#line 552 "./grammar.y"
                                                { (yyval.syntaxElementList) = new std::list<SyntaxElement*>(1, (yyvsp[-1].syntaxElement)); }
#line 2610 "./grammar.tab.c"
    break;

  case 72: /* elementList: '[' ']'  */
#line 553 "./grammar.y"
                                                { (yyval.syntaxElementList) = new std::list<SyntaxElement*>(); }
#line 2616 "./grammar.tab.c"
    break;

  case 73: /* elementList: elementList '[' expression ']'  */
#line 554 "./grammar.y"
                                                { (yyvsp[-3].syntaxElementList)->push_back((yyvsp[-1].syntaxElement)); (yyval.syntaxElementList) = (yyvsp[-3].syntaxElementList); }
#line 2622 "./grammar.tab.c"
    break;

  case 74: /* elementList: elementList '[' ']'  */
#line 555 "./grammar.y"
                                                { (yyvsp[-2].syntaxElementList)->push_back( NULL ); (yyval.syntaxElementList) = (yyvsp[-2].syntaxElementList); }
#line 2628 "./grammar.tab.c"
    break;

  case 75: /* fxnCall: identifier '(' optArguments ')'  */
#line 559 "./grammar.y"
                {
                    (yyval.syntaxFunctionCall) = new SyntaxFunctionCall(*(yyvsp[-3].string), (yyvsp[-1].argumentList));
                    delete (yyvsp[-3].string);
                }
#line 2637 "./grammar.tab.c"
    break;

  case 76: /* functionCall: fxnCall  */
#line 566 "./grammar.y"
                    {
#ifdef DEBUG_BISON_FLEX
                        printf("Parser inserting function call in syntax tree\n");
#endif
                        (yyval.syntaxFunctionCall) = (yyvsp[0].syntaxFunctionCall);
                    }
#line 2648 "./grammar.tab.c"
    break;

  case 77: /* functionCall: variable '.' fxnCall  */
#line 573 "./grammar.y"
                    {
#ifdef DEBUG_BISON_FLEX
                        printf("Parser inserting member call in syntax tree\n");
#endif
                        (yyvsp[0].syntaxFunctionCall)->setBaseVariable((yyvsp[-2].syntaxElement));
                        (yyval.syntaxFunctionCall) = (yyvsp[0].syntaxFunctionCall);
                    }
#line 2660 "./grammar.tab.c"
    break;

  case 78: /* functionCall: functionCall '.' fxnCall  */
#line 581 "./grammar.y"
                    {
#ifdef DEBUG_BISON_FLEX
                        printf("Parser inserting member call in syntax tree\n");
#endif
                        (yyvsp[0].syntaxFunctionCall)->setBaseVariable((yyvsp[-2].syntaxFunctionCall));
                        (yyval.syntaxFunctionCall) = (yyvsp[0].syntaxFunctionCall);
                    }
#line 2672 "./grammar.tab.c"
    break;

  case 79: /* functionCall: expression PIPE fxnCall  */
#line 589 "./grammar.y"
                    {
                        (yyval.syntaxFunctionCall) = (yyvsp[0].syntaxFunctionCall);
                        (yyval.syntaxFunctionCall)->pipeAddArg((yyvsp[-2].syntaxElement));
                    }
#line 2681 "./grammar.tab.c"
    break;

  case 80: /* optArguments: %empty  */
#line 596 "./grammar.y"
                                            { (yyval.argumentList) = new std::list<SyntaxLabeledExpr*>(); }
#line 2687 "./grammar.tab.c"
    break;

  case 81: /* optArguments: argumentList  */
#line 597 "./grammar.y"
                                            { (yyval.argumentList) = (yyvsp[0].argumentList); }
#line 2693 "./grammar.tab.c"
    break;

  case 82: /* argumentList: argument  */
#line 600 "./grammar.y"
                                                { (yyval.argumentList) = new std::list<SyntaxLabeledExpr*>(1,(yyvsp[0].syntaxLabeledExpr)); }
#line 2699 "./grammar.tab.c"
    break;

  case 83: /* argumentList: argumentList ',' argument  */
#line 601 "./grammar.y"
                                                { (yyvsp[-2].argumentList)->push_back((yyvsp[0].syntaxLabeledExpr)); (yyval.argumentList) = (yyvsp[-2].argumentList); }
#line 2705 "./grammar.tab.c"
    break;

  case 84: /* argument: expression  */
#line 605 "./grammar.y"
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Parser inserting unlabeled argument in syntax tree\n");
#endif
                    (yyval.syntaxLabeledExpr) = new SyntaxLabeledExpr( "", (yyvsp[0].syntaxElement));
                }
#line 2716 "./grammar.tab.c"
    break;

  case 85: /* argument: identifier EQUAL expression  */
#line 612 "./grammar.y"
                { 
#ifdef DEBUG_BISON_FLEX
                    printf("Parser inserting labeled argument in syntax tree\n");
#endif
                    (yyval.syntaxLabeledExpr) = new SyntaxLabeledExpr(*(yyvsp[-2].string), (yyvsp[0].syntaxElement));
                    delete (yyvsp[-2].string);
                }
#line 2728 "./grammar.tab.c"
    break;

  case 86: /* functionDef: FUNCTION identifier '(' optFormals ')' stmts  */
#line 622 "./grammar.y"
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Parser inserting function definition in syntax tree\n");
#endif
                    (yyval.syntaxElement) = new SyntaxFunctionDef("", *(yyvsp[-4].string), (yyvsp[-2].formalList), (yyvsp[0].syntaxElementList));
                    delete (yyvsp[-4].string);
                }
#line 2740 "./grammar.tab.c"
    break;

  case 87: /* functionDef: FUNCTION identifier optDims identifier '(' optFormals ')' stmts  */
#line 631 "./grammar.y"
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Parser inserting typed function definition in syntax tree\n");
#endif
                    (yyvsp[-6].string)->append(*(yyvsp[-5].string));
                    (yyval.syntaxElement) = new SyntaxFunctionDef(*(yyvsp[-6].string), *(yyvsp[-4].string), (yyvsp[-2].formalList), (yyvsp[0].syntaxElementList));
                    delete (yyvsp[-6].string);
                    delete (yyvsp[-5].string);
                    delete (yyvsp[-4].string);
                }
#line 2755 "./grammar.tab.c"
    break;

  case 88: /* procedureDef: PROCEDURE identifier '(' optFormals ')' stmts  */
#line 644 "./grammar.y"
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Parser inserting procedure definition in syntax tree\n");
#endif
                    (yyval.syntaxElement) = new SyntaxFunctionDef("", *(yyvsp[-4].string), (yyvsp[-2].formalList), (yyvsp[0].syntaxElementList), true);
                    delete (yyvsp[-4].string);
                }
#line 2767 "./grammar.tab.c"
    break;

  case 89: /* procedureDef: PROCEDURE identifier optDims identifier '(' optFormals ')' stmts  */
#line 653 "./grammar.y"
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Parser inserting typed procedure definition in syntax tree\n");
#endif
                    (yyvsp[-6].string)->append(*(yyvsp[-5].string));
                    (yyval.syntaxElement) = new SyntaxFunctionDef(*(yyvsp[-6].string), *(yyvsp[-4].string), (yyvsp[-2].formalList), (yyvsp[0].syntaxElementList), true);
                    delete (yyvsp[-6].string);
                    delete (yyvsp[-5].string);
                    delete (yyvsp[-4].string);
                }
#line 2782 "./grammar.tab.c"
    break;

  case 90: /* optFormals: %empty  */
#line 665 "./grammar.y"
                                        { (yyval.formalList) = new std::list<SyntaxFormal*>(); }
#line 2788 "./grammar.tab.c"
    break;

  case 91: /* optFormals: formalList  */
#line 666 "./grammar.y"
                                        { (yyval.formalList) = (yyvsp[0].formalList); }
#line 2794 "./grammar.tab.c"
    break;

  case 92: /* formalList: formal  */
#line 669 "./grammar.y"
                                        { (yyval.formalList) = new std::list<SyntaxFormal*>(1, (yyvsp[0].syntaxFormal)); }
#line 2800 "./grammar.tab.c"
    break;

  case 93: /* formalList: formalList ',' formal  */
#line 670 "./grammar.y"
                                        { (yyvsp[-2].formalList)->push_back((yyvsp[0].syntaxFormal)); (yyval.formalList) = (yyvsp[-2].formalList); }
#line 2806 "./grammar.tab.c"
    break;

  case 94: /* formal: identifier  */
#line 674 "./grammar.y"
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Inserting labeled formal argument without default in syntax tree\n");
#endif
                    (yyval.syntaxFormal) = new SyntaxFormal(*(yyvsp[0].string), NULL );
                    delete (yyvsp[0].string);
                }
#line 2818 "./grammar.tab.c"
    break;

  case 95: /* formal: identifier EQUAL expression  */
#line 682 "./grammar.y"
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Inserting labeled formal argument with default in syntax tree\n");
#endif
                    (yyval.syntaxFormal) = new SyntaxFormal(*(yyvsp[-2].string), (yyvsp[0].syntaxElement));
                    delete (yyvsp[-2].string);
                }
#line 2830 "./grammar.tab.c"
    break;

  case 96: /* formal: typeSpec identifier  */
#line 690 "./grammar.y"
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Inserting typed labeled formal argument without default in syntax tree\n");
#endif
                    (yyval.syntaxFormal) = new SyntaxFormal(*(yyvsp[-1].string), *(yyvsp[0].string), NULL);
                    delete (yyvsp[-1].string);
                    delete (yyvsp[0].string);
                }
#line 2843 "./grammar.tab.c"
    break;

  case 97: /* formal: typeSpec identifier EQUAL expression  */
#line 699 "./grammar.y"
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Inserting typed labeled formal argument with default in syntax tree\n");
#endif
                    (yyval.syntaxFormal) = new SyntaxFormal(*(yyvsp[-3].string), *(yyvsp[-2].string), (yyvsp[0].syntaxElement));
                    delete (yyvsp[-3].string);
                    delete (yyvsp[-2].string);
                }
#line 2856 "./grammar.tab.c"
    break;

  case 98: /* typeSpec: identifier optDims  */
#line 709 "./grammar.y"
                                                        { (yyvsp[-1].string)->append(*((yyvsp[0].string))); delete (yyvsp[0].string); (yyval.string) = (yyvsp[-1].string); }
#line 2862 "./grammar.tab.c"
    break;

  case 99: /* typeSpec: MOD_CONST identifier optDims  */
#line 710 "./grammar.y"
                                                        { (yyvsp[-1].string)->append(*((yyvsp[0].string))); (yyvsp[-1].string)->insert(0, "const ");           delete (yyvsp[0].string); (yyval.string) = (yyvsp[-1].string); }
#line 2868 "./grammar.tab.c"
    break;

  case 100: /* typeSpec: MOD_DYNAMIC identifier optDims  */
#line 711 "./grammar.y"
                                                        { (yyvsp[-1].string)->append(*((yyvsp[0].string))); (yyvsp[-1].string)->insert(0, "dynamic ");         delete (yyvsp[0].string); (yyval.string) = (yyvsp[-1].string); }
#line 2874 "./grammar.tab.c"
    break;

  case 101: /* typeSpec: MOD_STOCHASTIC identifier optDims  */
#line 712 "./grammar.y"
                                                        { (yyvsp[-1].string)->append(*((yyvsp[0].string))); (yyvsp[-1].string)->insert(0, "stochastic ");      delete (yyvsp[0].string); (yyval.string) = (yyvsp[-1].string); }
#line 2880 "./grammar.tab.c"
    break;

  case 102: /* typeSpec: MOD_DETERMINISTIC identifier optDims  */
#line 713 "./grammar.y"
                                                        { (yyvsp[-1].string)->append(*((yyvsp[0].string))); (yyvsp[-1].string)->insert(0, "deterministic ");   delete (yyvsp[0].string); (yyval.string) = (yyvsp[-1].string); }
#line 2886 "./grammar.tab.c"
    break;

  case 103: /* optDims: %empty  */
#line 716 "./grammar.y"
                                            { (yyval.string) = new std::string(""); }
#line 2892 "./grammar.tab.c"
    break;

  case 104: /* optDims: dimList  */
#line 717 "./grammar.y"
                                            { (yyval.string) = (yyvsp[0].string); }
#line 2898 "./grammar.tab.c"
    break;

  case 105: /* dimList: '[' ']'  */
#line 720 "./grammar.y"
                                            { (yyval.string) = new std::string("[]"); }
#line 2904 "./grammar.tab.c"
    break;

  case 106: /* dimList: dimList '[' ']'  */
#line 721 "./grammar.y"
                                            { (yyvsp[-2].string)->append("[]"); (yyval.string) = (yyvsp[-2].string); }
#line 2910 "./grammar.tab.c"
    break;

  case 107: /* stmts: '{' stmtList '}'  */
#line 724 "./grammar.y"
                                                { (yyval.syntaxElementList) = (yyvsp[-1].syntaxElementList); }
#line 2916 "./grammar.tab.c"
    break;

  case 108: /* stmts: stmt_or_expr  */
#line 726 "./grammar.y"
                {
                    std::list<SyntaxElement*>* stmts = new std::list<SyntaxElement*>();
                    stmts->push_back((yyvsp[0].syntaxElement));
                    (yyval.syntaxElementList) = stmts;
                }
#line 2926 "./grammar.tab.c"
    break;

  case 109: /* stmtList: %empty  */
#line 733 "./grammar.y"
                                            { (yyval.syntaxElementList) = new std::list<SyntaxElement*>(); }
#line 2932 "./grammar.tab.c"
    break;

  case 110: /* stmtList: stmt_or_expr  */
#line 734 "./grammar.y"
                                            { (yyval.syntaxElementList) = new std::list<SyntaxElement*>(1, (yyvsp[0].syntaxElement)); }
#line 2938 "./grammar.tab.c"
    break;

  case 111: /* stmtList: stmtList ';' stmt_or_expr  */
#line 735 "./grammar.y"
                                            { (yyvsp[-2].syntaxElementList)->push_back((yyvsp[0].syntaxElement)); (yyval.syntaxElementList) = (yyvsp[-2].syntaxElementList); }
#line 2944 "./grammar.tab.c"
    break;

  case 112: /* stmtList: stmtList ';'  */
#line 736 "./grammar.y"
                                            { (yyval.syntaxElementList) = (yyvsp[-1].syntaxElementList); }
#line 2950 "./grammar.tab.c"
    break;

  case 113: /* stmtList: stmtList '\n' stmt_or_expr  */
#line 737 "./grammar.y"
                                            { (yyvsp[-2].syntaxElementList)->push_back((yyvsp[0].syntaxElement)); (yyval.syntaxElementList) = (yyvsp[-2].syntaxElementList); }
#line 2956 "./grammar.tab.c"
    break;

  case 114: /* stmtList: stmtList '\n'  */
#line 738 "./grammar.y"
                                            { (yyval.syntaxElementList) = (yyvsp[-1].syntaxElementList); }
#line 2962 "./grammar.tab.c"
    break;

  case 115: /* statement: ifStatement  */
#line 741 "./grammar.y"
                                        { (yyval.syntaxElement) = (yyvsp[0].syntaxElement); }
#line 2968 "./grammar.tab.c"
    break;

  case 116: /* statement: forStatement  */
#line 742 "./grammar.y"
                                        { (yyval.syntaxElement) = (yyvsp[0].syntaxElement); }
#line 2974 "./grammar.tab.c"
    break;

  case 117: /* statement: whileStatement  */
#line 743 "./grammar.y"
                                        { (yyval.syntaxElement) = (yyvsp[0].syntaxElement); }
#line 2980 "./grammar.tab.c"
    break;

  case 118: /* statement: nextStatement  */
#line 744 "./grammar.y"
                                        { (yyval.syntaxElement) = (yyvsp[0].syntaxElement); }
#line 2986 "./grammar.tab.c"
    break;

  case 119: /* statement: breakStatement  */
#line 745 "./grammar.y"
                                        { (yyval.syntaxElement) = (yyvsp[0].syntaxElement); }
#line 2992 "./grammar.tab.c"
    break;

  case 120: /* statement: returnStatement  */
#line 746 "./grammar.y"
                                        { (yyval.syntaxElement) = (yyvsp[0].syntaxElement); }
#line 2998 "./grammar.tab.c"
    break;

  case 121: /* stmt_or_expr: statement  */
#line 749 "./grammar.y"
                                        { (yyval.syntaxElement) = (yyvsp[0].syntaxElement); }
#line 3004 "./grammar.tab.c"
    break;

  case 122: /* stmt_or_expr: expression  */
#line 750 "./grammar.y"
                                        { (yyval.syntaxElement) = (yyvsp[0].syntaxElement); }
#line 3010 "./grammar.tab.c"
    break;

  case 123: /* declaration: classDef  */
#line 753 "./grammar.y"
                                        { (yyval.syntaxElement) = (yyvsp[0].syntaxElement); }
#line 3016 "./grammar.tab.c"
    break;

  case 124: /* declaration: functionDef  */
#line 754 "./grammar.y"
                                        { (yyval.syntaxElement) = (yyvsp[0].syntaxElement); }
#line 3022 "./grammar.tab.c"
    break;

  case 125: /* declaration: procedureDef  */
#line 755 "./grammar.y"
                                        { (yyval.syntaxElement) = (yyvsp[0].syntaxElement); }
#line 3028 "./grammar.tab.c"
    break;

  case 126: /* declaration: identifier optElements identifier  */
#line 757 "./grammar.y"
                    {
#ifdef DEBUG_BISON_FLEX
                        printf("Parser inserting variable declaration in syntax tree\n");
#endif
                        (yyval.syntaxElement) = new SyntaxVariableDecl(*(yyvsp[-2].string), (yyvsp[-1].syntaxElementList), *(yyvsp[0].string));
                        delete (yyvsp[-2].string);
                        delete (yyvsp[0].string);
                    }
#line 3041 "./grammar.tab.c"
    break;

  case 127: /* memberDefs: %empty  */
#line 767 "./grammar.y"
                                                { (yyval.syntaxElementList) = new std::list<SyntaxElement*>(); }
#line 3047 "./grammar.tab.c"
    break;

  case 128: /* memberDefs: memberDef  */
#line 768 "./grammar.y"
                                                { (yyval.syntaxElementList) = new std::list<SyntaxElement*>(1, (yyvsp[0].syntaxElement)); }
#line 3053 "./grammar.tab.c"
    break;

  case 129: /* memberDefs: memberDefs ';' memberDef  */
#line 769 "./grammar.y"
                                                { (yyvsp[-2].syntaxElementList)->push_back((yyvsp[0].syntaxElement)); (yyval.syntaxElementList) = (yyvsp[-2].syntaxElementList); }
#line 3059 "./grammar.tab.c"
    break;

  case 130: /* memberDefs: memberDefs ';'  */
#line 770 "./grammar.y"
                                                { (yyval.syntaxElementList) = (yyvsp[-1].syntaxElementList); }
#line 3065 "./grammar.tab.c"
    break;

  case 131: /* memberDefs: memberDefs '\n' memberDef  */
#line 771 "./grammar.y"
                                                { (yyvsp[-2].syntaxElementList)->push_back((yyvsp[0].syntaxElement)); (yyval.syntaxElementList) = (yyvsp[-2].syntaxElementList); }
#line 3071 "./grammar.tab.c"
    break;

  case 132: /* memberDefs: memberDefs '\n'  */
#line 772 "./grammar.y"
                                                { (yyval.syntaxElementList) = (yyvsp[-1].syntaxElementList); }
#line 3077 "./grammar.tab.c"
    break;

  case 133: /* memberDef: formal  */
#line 775 "./grammar.y"
                                    { (yyval.syntaxElement) = (yyvsp[0].syntaxFormal); }
#line 3083 "./grammar.tab.c"
    break;

  case 134: /* memberDef: PROTECTED formal  */
#line 776 "./grammar.y"
                                    { (yyvsp[0].syntaxFormal)->setIsProtected(); (yyval.syntaxElement) = (yyvsp[0].syntaxFormal); }
#line 3089 "./grammar.tab.c"
    break;

  case 135: /* memberDef: functionDef  */
#line 777 "./grammar.y"
                                    { (yyval.syntaxElement) = (yyvsp[0].syntaxElement); }
#line 3095 "./grammar.tab.c"
    break;

  case 136: /* memberDef: procedureDef  */
#line 778 "./grammar.y"
                                    { (yyval.syntaxElement) = (yyvsp[0].syntaxElement); }
#line 3101 "./grammar.tab.c"
    break;

  case 137: /* classDef: CLASS identifier ':' identifier '{' memberDefs '}'  */
#line 782 "./grammar.y"
                {
#ifdef DEBUG_BISON_FLEX
                printf("Parser inserting class definition (CLASS_DEF) in syntax tree\n");
#endif
                    (yyval.syntaxElement) = new SyntaxClassDef(*(yyvsp[-5].string), *(yyvsp[-3].string), (yyvsp[-1].syntaxElementList));
                    delete (yyvsp[-5].string);
                    delete (yyvsp[-3].string);
                }
#line 3114 "./grammar.tab.c"
    break;

  case 138: /* ifStatement: IF cond stmts  */
#line 792 "./grammar.y"
                                            { (yyval.syntaxElement) = new SyntaxStatement(SyntaxStatement::If, (yyvsp[-1].syntaxElement), (yyvsp[0].syntaxElementList)); }
#line 3120 "./grammar.tab.c"
    break;

  case 139: /* ifStatement: IF cond stmts ELSE stmts  */
#line 793 "./grammar.y"
                                            { (yyval.syntaxElement) = new SyntaxStatement(SyntaxStatement::IfElse, (yyvsp[-3].syntaxElement), (yyvsp[-2].syntaxElementList), (yyvsp[0].syntaxElementList)); }
#line 3126 "./grammar.tab.c"
    break;

  case 140: /* cond: '(' expression ')'  */
#line 795 "./grammar.y"
                                  { (yyval.syntaxElement) = (yyvsp[-1].syntaxElement); }
#line 3132 "./grammar.tab.c"
    break;

  case 141: /* forStatement: FOR forCond stmts  */
#line 798 "./grammar.y"
                                        { (yyval.syntaxElement) = new SyntaxStatement(SyntaxStatement::For, (yyvsp[-1].syntaxElement), (yyvsp[0].syntaxElementList)); }
#line 3138 "./grammar.tab.c"
    break;

  case 142: /* forCond: '(' identifier IN expression ')'  */
#line 801 "./grammar.y"
                                                    { (yyval.syntaxElement) = new SyntaxForLoop(*(yyvsp[-3].string), (yyvsp[-1].syntaxElement)); delete (yyvsp[-3].string); }
#line 3144 "./grammar.tab.c"
    break;

  case 143: /* whileStatement: WHILE cond stmts  */
#line 804 "./grammar.y"
                                        { (yyval.syntaxElement) = new SyntaxStatement(SyntaxStatement::While, (yyvsp[-1].syntaxElement), (yyvsp[0].syntaxElementList)); }
#line 3150 "./grammar.tab.c"
    break;

  case 144: /* nextStatement: NEXT  */
#line 807 "./grammar.y"
                            { (yyval.syntaxElement) = new SyntaxStatement(SyntaxStatement::Next); }
#line 3156 "./grammar.tab.c"
    break;

  case 145: /* breakStatement: BREAK  */
#line 810 "./grammar.y"
                            { (yyval.syntaxElement) = new SyntaxStatement(SyntaxStatement::Break); }
#line 3162 "./grammar.tab.c"
    break;

  case 146: /* returnStatement: RETURN  */
#line 813 "./grammar.y"
                                        { (yyval.syntaxElement) = new SyntaxStatement(SyntaxStatement::Return); }
#line 3168 "./grammar.tab.c"
    break;

  case 147: /* returnStatement: RETURN expression  */
#line 814 "./grammar.y"
                                        { (yyval.syntaxElement) = new SyntaxStatement(SyntaxStatement::Return, (yyvsp[0].syntaxElement)); }
#line 3174 "./grammar.tab.c"
    break;

  case 148: /* identifier: NAME  */
#line 817 "./grammar.y"
                        { (yyval.string) = new std::string((yyvsp[0].c_string)); }
#line 3180 "./grammar.tab.c"
    break;

  case 149: /* vector: '[' vectorList ']'  */
#line 821 "./grammar.y"
                                        { (yyval.argumentList) = (yyvsp[-1].argumentList); }
#line 3186 "./grammar.tab.c"
    break;

  case 150: /* vectorList: vectorList ',' expression  */
#line 824 "./grammar.y"
                                            { (yyvsp[-2].argumentList)->push_back(new SyntaxLabeledExpr( "", (yyvsp[0].syntaxElement)) ); (yyval.argumentList) = (yyvsp[-2].argumentList); }
#line 3192 "./grammar.tab.c"
    break;

  case 151: /* vectorList: expression  */
#line 826 "./grammar.y"
                {
                (yyval.argumentList) = new std::list<SyntaxLabeledExpr*>(1, new SyntaxLabeledExpr("", (yyvsp[0].syntaxElement)) );
                }
#line 3200 "./grammar.tab.c"
    break;

  case 152: /* constant: FALSE  */
#line 832 "./grammar.y"
                {
                    (yyval.syntaxElement) = new SyntaxConstant(new RlBoolean(false) );
                }
#line 3208 "./grammar.tab.c"
    break;

  case 153: /* constant: TRUE  */
#line 836 "./grammar.y"
                {
                    (yyval.syntaxElement) = new SyntaxConstant(new RlBoolean(true) );
                }
#line 3216 "./grammar.tab.c"
    break;

  case 154: /* constant: RBNULL  */
#line 840 "./grammar.y"
                {
                    (yyval.syntaxElement) = new SyntaxConstant( NULL );
                }
#line 3224 "./grammar.tab.c"
    break;

  case 155: /* constant: RBTAB  */
#line 844 "./grammar.y"
                {
                    (yyval.syntaxElement) = new SyntaxConstant( new RlString("\t") );
                }
#line 3232 "./grammar.tab.c"
    break;

  case 156: /* constant: RBINF  */
#line 848 "./grammar.y"
                {
                    (yyval.syntaxElement) = new SyntaxConstant( new RealPos( RbConstants::Double::inf ) );
                }
#line 3240 "./grammar.tab.c"
    break;

  case 157: /* constant: INT  */
#line 852 "./grammar.y"
                {
                    if ( (yyvsp[0].longIntValue) < 0 ) {
                        (yyval.syntaxElement) = new SyntaxConstant(new Integer((yyvsp[0].longIntValue)) );
                    }
                    else { 
                        (yyval.syntaxElement) = new SyntaxConstant(new Natural((yyvsp[0].longIntValue)) );
                    }
                }
#line 3253 "./grammar.tab.c"
    break;

  case 158: /* constant: STRING  */
#line 861 "./grammar.y"
                {
                    (yyval.syntaxElement) = new SyntaxConstant(new RlString((yyvsp[0].c_string)) );
                }
#line 3261 "./grammar.tab.c"
    break;

  case 159: /* constant: REAL  */
#line 865 "./grammar.y"
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
                    
                    if ((yyvsp[0].realValue) >= 0.0 && (yyvsp[0].realValue) <= 1.0) {
                        #ifdef DEBUG_BISON_FLEX
                        printf("Parser inserting RealPos constant in syntax tree \n");
                        #endif
                        (yyval.syntaxElement) = new SyntaxConstant(new Probability((yyvsp[0].realValue)) );
                    }
                    else if ((yyvsp[0].realValue) >= 0.0) {
#ifdef DEBUG_BISON_FLEX
                        printf("Parser inserting RealPos constant in syntax tree \n");
#endif
                        (yyval.syntaxElement) = new SyntaxConstant(new RealPos((yyvsp[0].realValue)) );
                    }
                    else {
#ifdef DEBUG_BISON_FLEX
                        printf("Parser inserting Real constant in syntax tree \n");
#endif
                        (yyval.syntaxElement) = new SyntaxConstant(new Real((yyvsp[0].realValue)) );
                    }
                }
#line 3324 "./grammar.tab.c"
    break;


#line 3328 "./grammar.tab.c"

      default: break;
    }
  /* User semantic actions sometimes alter yychar, and that requires
     that yytoken be updated with the new translation.  We take the
     approach of translating immediately before every use of yytoken.
     One alternative is translating here after every semantic action,
     but that translation would be missed if the semantic action invokes
     YYABORT, YYACCEPT, or YYERROR immediately after altering yychar or
     if it invokes YYBACKUP.  In the case of YYABORT or YYACCEPT, an
     incorrect destructor might then be invoked immediately.  In the
     case of YYERROR or YYBACKUP, subsequent parser actions might lead
     to an incorrect destructor call or verbose syntax error message
     before the lookahead is translated.  */
  YY_SYMBOL_PRINT ("-> $$ =", YY_CAST (yysymbol_kind_t, yyr1[yyn]), &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;

  *++yyvsp = yyval;
  *++yylsp = yyloc;

  /* Now 'shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */
  {
    const int yylhs = yyr1[yyn] - YYNTOKENS;
    const int yyi = yypgoto[yylhs] + *yyssp;
    yystate = (0 <= yyi && yyi <= YYLAST && yycheck[yyi] == *yyssp
               ? yytable[yyi]
               : yydefgoto[yylhs]);
  }

  goto yynewstate;


/*--------------------------------------.
| yyerrlab -- here on detecting error.  |
`--------------------------------------*/
yyerrlab:
  /* Make sure we have latest lookahead translation.  See comments at
     user semantic actions for why this is necessary.  */
  yytoken = yychar == YYEMPTY ? YYSYMBOL_YYEMPTY : YYTRANSLATE (yychar);
  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
      yyerror (YY_("syntax error"));
    }

  yyerror_range[1] = yylloc;
  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
         error, discard it.  */

      if (yychar <= YYEOF)
        {
          /* Return failure if at end of input.  */
          if (yychar == YYEOF)
            YYABORT;
        }
      else
        {
          yydestruct ("Error: discarding",
                      yytoken, &yylval, &yylloc);
          yychar = YYEMPTY;
        }
    }

  /* Else will try to reuse lookahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:
  /* Pacify compilers when the user code never invokes YYERROR and the
     label yyerrorlab therefore never appears in user code.  */
  if (0)
    YYERROR;
  ++yynerrs;

  /* Do not reclaim the symbols of the rule whose action triggered
     this YYERROR.  */
  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;      /* Each real token shifted decrements this.  */

  /* Pop stack until we find a state that shifts the error token.  */
  for (;;)
    {
      yyn = yypact[yystate];
      if (!yypact_value_is_default (yyn))
        {
          yyn += YYSYMBOL_YYerror;
          if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYSYMBOL_YYerror)
            {
              yyn = yytable[yyn];
              if (0 < yyn)
                break;
            }
        }

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
        YYABORT;

      yyerror_range[1] = *yylsp;
      yydestruct ("Error: popping",
                  YY_ACCESSING_SYMBOL (yystate), yyvsp, yylsp);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END

  yyerror_range[2] = yylloc;
  ++yylsp;
  YYLLOC_DEFAULT (*yylsp, yyerror_range, 2);

  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", YY_ACCESSING_SYMBOL (yyn), yyvsp, yylsp);

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturnlab;


/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturnlab;


/*-----------------------------------------------------------.
| yyexhaustedlab -- YYNOMEM (memory exhaustion) comes here.  |
`-----------------------------------------------------------*/
yyexhaustedlab:
  yyerror (YY_("memory exhausted"));
  yyresult = 2;
  goto yyreturnlab;


/*----------------------------------------------------------.
| yyreturnlab -- parsing is finished, clean up and return.  |
`----------------------------------------------------------*/
yyreturnlab:
  if (yychar != YYEMPTY)
    {
      /* Make sure we have latest lookahead translation.  See comments at
         user semantic actions for why this is necessary.  */
      yytoken = YYTRANSLATE (yychar);
      yydestruct ("Cleanup: discarding lookahead",
                  yytoken, &yylval, &yylloc);
    }
  /* Do not reclaim the symbols of the rule whose action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
                  YY_ACCESSING_SYMBOL (+*yyssp), yyvsp, yylsp);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif

  return yyresult;
}

#line 925 "./grammar.y"



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


