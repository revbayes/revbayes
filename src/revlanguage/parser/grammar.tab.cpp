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

#line 149 "./grammar.tab.c"

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
  YYSYMBOL_END_OF_INPUT = 51,              /* END_OF_INPUT  */
  YYSYMBOL_52_ = 52,                       /* '?'  */
  YYSYMBOL_UNOT = 53,                      /* UNOT  */
  YYSYMBOL_54_ = 54,                       /* '+'  */
  YYSYMBOL_55_ = 55,                       /* '-'  */
  YYSYMBOL_56_ = 56,                       /* '*'  */
  YYSYMBOL_57_ = 57,                       /* '/'  */
  YYSYMBOL_58_ = 58,                       /* ':'  */
  YYSYMBOL_59_ = 59,                       /* '%'  */
  YYSYMBOL_UMINUS = 60,                    /* UMINUS  */
  YYSYMBOL_UPLUS = 61,                     /* UPLUS  */
  YYSYMBOL_62_ = 62,                       /* '^'  */
  YYSYMBOL_63_ = 63,                       /* '.'  */
  YYSYMBOL_UAND = 64,                      /* UAND  */
  YYSYMBOL_65_ = 65,                       /* '('  */
  YYSYMBOL_66_ = 66,                       /* '['  */
  YYSYMBOL_67_n_ = 67,                     /* '\n'  */
  YYSYMBOL_68_ = 68,                       /* ';'  */
  YYSYMBOL_69_ = 69,                       /* ')'  */
  YYSYMBOL_70_ = 70,                       /* '!'  */
  YYSYMBOL_71_ = 71,                       /* ']'  */
  YYSYMBOL_72_ = 72,                       /* ','  */
  YYSYMBOL_73_ = 73,                       /* '{'  */
  YYSYMBOL_74_ = 74,                       /* '}'  */
  YYSYMBOL_YYACCEPT = 75,                  /* $accept  */
  YYSYMBOL_prog = 76,                      /* prog  */
  YYSYMBOL_expression = 77,                /* expression  */
  YYSYMBOL_arrowAssign = 78,               /* arrowAssign  */
  YYSYMBOL_tildeAssign = 79,               /* tildeAssign  */
  YYSYMBOL_equationAssign = 80,            /* equationAssign  */
  YYSYMBOL_workspaceAssign = 81,           /* workspaceAssign  */
  YYSYMBOL_referenceAssign = 82,           /* referenceAssign  */
  YYSYMBOL_additionAssign = 83,            /* additionAssign  */
  YYSYMBOL_subtractionAssign = 84,         /* subtractionAssign  */
  YYSYMBOL_multiplicationAssign = 85,      /* multiplicationAssign  */
  YYSYMBOL_divisionAssign = 86,            /* divisionAssign  */
  YYSYMBOL_variable = 87,                  /* variable  */
  YYSYMBOL_optElements = 88,               /* optElements  */
  YYSYMBOL_elementList = 89,               /* elementList  */
  YYSYMBOL_fxnCall = 90,                   /* fxnCall  */
  YYSYMBOL_functionCall = 91,              /* functionCall  */
  YYSYMBOL_optArguments = 92,              /* optArguments  */
  YYSYMBOL_argumentList = 93,              /* argumentList  */
  YYSYMBOL_argument = 94,                  /* argument  */
  YYSYMBOL_functionDef = 95,               /* functionDef  */
  YYSYMBOL_procedureDef = 96,              /* procedureDef  */
  YYSYMBOL_optFormals = 97,                /* optFormals  */
  YYSYMBOL_formalList = 98,                /* formalList  */
  YYSYMBOL_formal = 99,                    /* formal  */
  YYSYMBOL_typeSpec = 100,                 /* typeSpec  */
  YYSYMBOL_optDims = 101,                  /* optDims  */
  YYSYMBOL_dimList = 102,                  /* dimList  */
  YYSYMBOL_stmts = 103,                    /* stmts  */
  YYSYMBOL_stmtList = 104,                 /* stmtList  */
  YYSYMBOL_statement = 105,                /* statement  */
  YYSYMBOL_stmt_or_expr = 106,             /* stmt_or_expr  */
  YYSYMBOL_declaration = 107,              /* declaration  */
  YYSYMBOL_memberDefs = 108,               /* memberDefs  */
  YYSYMBOL_memberDef = 109,                /* memberDef  */
  YYSYMBOL_classDef = 110,                 /* classDef  */
  YYSYMBOL_ifStatement = 111,              /* ifStatement  */
  YYSYMBOL_cond = 112,                     /* cond  */
  YYSYMBOL_forStatement = 113,             /* forStatement  */
  YYSYMBOL_forCond = 114,                  /* forCond  */
  YYSYMBOL_whileStatement = 115,           /* whileStatement  */
  YYSYMBOL_nextStatement = 116,            /* nextStatement  */
  YYSYMBOL_breakStatement = 117,           /* breakStatement  */
  YYSYMBOL_returnStatement = 118,          /* returnStatement  */
  YYSYMBOL_identifier = 119,               /* identifier  */
  YYSYMBOL_vector = 120,                   /* vector  */
  YYSYMBOL_vectorList = 121,               /* vectorList  */
  YYSYMBOL_constant = 122                  /* constant  */
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
#define YYFINAL  86
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   1177

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  75
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  48
/* YYNRULES -- Number of rules.  */
#define YYNRULES  158
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  287

/* YYMAXUTOK -- Last valid token kind.  */
#define YYMAXUTOK   310


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
      67,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    70,     2,     2,     2,    59,     2,     2,
      65,    69,    56,    54,    72,    55,    63,    57,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,    58,    68,
       2,     2,     2,    52,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,    66,     2,    71,    62,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    73,     2,    74,     2,     2,     2,     2,
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
      45,    46,    47,    48,    49,    50,    51,    53,    60,    61,
      64
};

#if YYDEBUG
/* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
static const yytype_int16 yyrline[] =
{
       0,   235,   235,   242,   249,   258,   267,   276,   285,   294,
     303,   313,   323,   332,   341,   348,   357,   359,   361,   363,
     364,   365,   366,   368,   369,   370,   371,   373,   375,   376,
     377,   378,   379,   380,   382,   383,   384,   385,   386,   387,
     389,   390,   391,   392,   394,   395,   396,   397,   398,   400,
     401,   402,   403,   405,   407,   410,   419,   428,   437,   446,
     455,   464,   473,   482,   491,   504,   516,   528,   543,   544,
     547,   548,   549,   550,   553,   560,   567,   575,   583,   591,
     592,   595,   596,   599,   606,   616,   625,   638,   647,   660,
     661,   664,   665,   668,   676,   684,   693,   704,   705,   706,
     707,   708,   711,   712,   715,   716,   719,   720,   728,   729,
     730,   731,   732,   733,   736,   737,   738,   739,   740,   741,
     744,   745,   748,   749,   750,   751,   762,   763,   764,   765,
     766,   767,   770,   771,   772,   773,   776,   787,   788,   790,
     793,   796,   799,   802,   805,   808,   809,   812,   816,   819,
     820,   826,   830,   834,   838,   842,   846,   855,   859
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
  "END_OF_INPUT", "'?'", "UNOT", "'+'", "'-'", "'*'", "'/'", "':'", "'%'",
  "UMINUS", "UPLUS", "'^'", "'.'", "UAND", "'('", "'['", "'\\n'", "';'",
  "')'", "'!'", "']'", "','", "'{'", "'}'", "$accept", "prog",
  "expression", "arrowAssign", "tildeAssign", "equationAssign",
  "workspaceAssign", "referenceAssign", "additionAssign",
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

#define YYPACT_NINF (-223)

#define yypact_value_is_default(Yyn) \
  ((Yyn) == YYPACT_NINF)

#define YYTABLE_NINF (-103)

#define yytable_value_is_error(Yyn) \
  ((Yyn) == YYTABLE_NINF)

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
static const yytype_int16 yypact[] =
{
     357,   -41,  -223,  -223,  -223,  -223,  -223,  -223,  -223,  -223,
    -223,    18,    18,    18,   -59,   -33,   -33,  -223,  -223,   485,
       2,     2,   485,  -223,   485,   485,   485,   485,   485,  -223,
     485,    36,   876,  -223,  -223,  -223,  -223,  -223,  -223,  -223,
    -223,  -223,   -28,   -24,   -25,  -223,  -223,  -223,   -34,   -23,
    -223,  -223,  -223,  -223,  -223,  -223,  -223,   -12,  -223,  -223,
    -223,  -223,    12,    48,    11,    18,   286,   485,   286,   286,
     876,   -12,   485,    -4,   -24,    -4,    29,   876,    35,    -5,
     -45,   -45,   724,   876,    49,  1096,  -223,   485,   485,   485,
     485,   485,   485,   485,   485,   485,   485,   485,   485,   485,
     485,   485,   485,   485,   485,   485,    18,   485,   485,   485,
     485,   485,   485,   485,  -223,  -223,    18,   485,    18,  -223,
    -223,  -223,  -223,   485,   377,    18,    21,    24,    19,    18,
      30,    24,    18,    18,    83,   465,  -223,  -223,   762,    86,
    -223,  -223,   800,    18,  -223,  -223,    18,  -223,  -223,    41,
    -223,   485,   876,   907,   938,   969,   999,  1028,  1053,  1077,
     288,  1096,    26,  1096,    26,  1115,  1115,  1115,  1115,  1115,
    1115,  -223,    44,    91,    91,    -7,    -7,   -45,   -45,   -45,
      45,   524,  -223,   876,    32,    43,  -223,   -35,  -223,   564,
    -223,   396,    18,    18,    18,    18,    50,    52,  -223,    18,
      -2,  -223,    57,    46,    59,    61,    62,   485,   -49,  -223,
    -223,   286,    41,    45,    77,   485,   876,   485,    64,  -223,
     485,   485,  -223,  -223,   604,    72,    72,    72,    72,   286,
      24,   101,   485,  -223,    24,  -223,   286,    24,   200,   838,
     465,   465,  -223,  -223,  -223,  -223,   644,   684,  -223,  -223,
     876,  -223,  -223,  -223,  -223,  -223,  -223,  -223,   485,   876,
      73,  -223,   104,    24,  -223,  -223,  -223,   -46,  -223,  -223,
    -223,  -223,    64,    64,   876,   286,   286,  -223,   200,   200,
    -223,  -223,  -223,  -223,  -223,  -223,  -223
};

/* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
   Performed when YYTABLE does not specify something else to do.  Zero
   means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       0,     0,   158,   156,   147,   157,   153,   154,   151,   152,
     155,     0,     0,     0,     0,     0,     0,   143,   144,   145,
       0,     0,     0,     2,     0,     0,     0,     0,     0,     3,
       0,     0,   121,    44,    46,    45,    47,    48,    49,    50,
      51,    52,    54,    75,    53,   123,   124,   120,     0,     0,
     122,   114,   115,   116,   117,   118,   119,    68,    17,    16,
      14,    15,   102,   102,     0,     0,     0,     0,     0,     0,
     146,    68,     0,    23,     0,    25,    22,     0,    53,    68,
      20,    19,     0,   150,     0,    21,     1,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,    24,    26,     0,     0,     0,     4,
       5,     6,     7,    79,     0,    64,    69,    89,     0,     0,
     103,    89,     0,     0,     0,   108,   140,   107,     0,   137,
     142,    64,     0,     0,    12,    13,     0,     8,     9,    18,
     148,     0,    55,    56,    57,    59,    60,    61,    62,    63,
      58,    40,    41,    42,    43,    39,    38,    34,    35,    36,
      37,    78,     0,    28,    29,    30,    31,    27,    33,    32,
      76,     0,    77,    83,     0,    80,    81,    68,    71,     0,
     125,     0,     0,     0,     0,     0,     0,    90,    91,     0,
      93,   104,     0,     0,     0,     0,     0,     0,     0,   109,
     139,     0,     0,     0,     0,     0,   149,     0,    68,    74,
       0,     0,    70,    73,     0,   102,   102,   102,   102,     0,
       0,    95,     0,    97,    89,   105,     0,    89,   126,     0,
     113,   111,   106,   138,    10,    11,     0,     0,    65,    82,
      84,    72,    98,    99,   100,   101,    85,    92,     0,    94,
       0,    87,     0,     0,   134,   135,   132,     0,   127,   141,
     112,   110,    68,    68,    96,     0,     0,   133,   131,   129,
     136,    66,    67,    86,    88,   130,   128
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
    -223,  -223,    67,  -223,  -223,  -223,  -223,  -223,  -223,  -223,
    -223,  -223,   131,   -37,  -223,    -6,   163,  -223,  -223,   -32,
     189,   196,  -129,  -223,  -222,  -223,   -17,  -223,   -29,  -223,
    -223,     1,  -223,  -223,   -93,  -223,  -223,   181,  -223,  -223,
    -223,  -223,  -223,  -223,     0,  -223,  -223,  -223
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
       0,    31,    32,    33,    34,    35,    36,    37,    38,    39,
      40,    41,    42,   141,   126,    43,    44,   184,   185,   186,
     264,   265,   196,   197,   198,   199,   129,   130,   136,   208,
      47,   137,    49,   267,   268,    50,    51,    68,    52,    66,
      53,    54,    55,    56,    71,    58,    84,    59
};

/* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule whose
   number is the opposite.  If YYTABLE_NINF, syntax error.  */
static const yytype_int16 yytable[] =
{
      57,    48,   204,  -102,   221,   106,    65,     4,   257,   114,
     115,    62,    63,    64,    74,    74,   266,   113,   240,   241,
     125,   278,   279,     4,    79,   242,    60,    61,   280,     4,
     123,   124,    67,   119,   120,   116,    86,   232,   118,   139,
     140,   277,   117,   106,   121,   122,   132,   192,   193,   194,
     195,   111,   112,   123,   124,   113,   266,   266,   146,   143,
     123,   124,   147,   148,   128,   134,    96,    72,    98,   133,
     100,   101,   102,   103,   104,   105,   106,   127,   128,   106,
     107,   108,   109,   110,   111,   112,    70,   191,   113,    76,
     201,    77,    80,    81,    82,    83,   203,    85,   118,   207,
     171,   219,   144,   145,   211,   260,   172,   215,   262,   123,
     180,   217,   182,   131,   128,   220,   172,   235,   172,   229,
     150,   151,   234,   187,   230,   190,   237,   200,   236,   202,
     124,   200,   205,   206,   138,   238,   209,   213,   128,   142,
     258,   106,   275,   172,   244,   245,   214,   109,   110,   111,
     112,    73,    75,   113,   152,   153,   154,   155,   156,   157,
     158,   159,   160,   161,   162,   163,   164,   165,   166,   167,
     168,   169,   170,   276,   173,   174,   175,   176,   177,   178,
     179,   248,   243,   233,   181,   285,   286,    78,   249,    45,
     183,   189,   225,   226,   227,   228,    46,    69,     0,   231,
     256,     0,     0,     0,     0,     4,     0,   261,   252,   253,
     254,   255,    11,    12,     0,     0,     0,     0,   216,     0,
     187,     0,     0,   192,   193,   194,   195,   263,     0,     0,
     200,     0,     0,     0,   200,   281,   282,   200,   200,     0,
       0,   270,   271,     0,     0,     0,   283,   284,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   224,     0,
       0,     0,     0,   200,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   239,     0,     0,     0,   200,   200,
       0,     0,   246,     0,   247,     0,     0,   183,   250,     2,
       3,     4,     5,     6,     7,     8,     9,    10,     0,   259,
       0,    14,     0,    15,     0,    16,    17,    18,    19,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,    20,    21,   274,    22,    95,    96,    97,
      98,    99,   100,   101,   102,   103,   104,   105,   106,     0,
      25,    26,   107,   108,   109,   110,   111,   112,     0,     0,
     113,    27,    28,     0,     0,     0,    30,     0,     1,   135,
       2,     3,     4,     5,     6,     7,     8,     9,    10,    11,
      12,    13,    14,     0,    15,     0,    16,    17,    18,    19,
       2,     3,     4,     5,     6,     7,     8,     9,    10,     0,
       0,     0,     0,     0,    20,    21,     0,    22,     0,     2,
       3,     4,     5,     6,     7,     8,     9,    10,    23,    24,
       0,    25,    26,     0,    20,    21,     0,    22,     0,     0,
       0,     0,    27,    28,    29,     0,     0,    30,     0,     0,
       0,    25,    26,    20,    21,     0,    22,     0,     0,     0,
       0,     0,    27,    28,     0,     0,     0,    30,   188,     0,
      25,    26,     0,     0,     0,     0,     0,     0,     0,     0,
       0,    27,    28,     0,     0,     0,    30,   223,     2,     3,
       4,     5,     6,     7,     8,     9,    10,     0,     0,     0,
      14,     0,    15,     0,    16,    17,    18,    19,     2,     3,
       4,     5,     6,     7,     8,     9,    10,     0,     0,     0,
       0,     0,    20,    21,     0,    22,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,    25,
      26,     0,    20,    21,     0,    22,     0,     0,     0,     0,
      27,    28,     0,     0,     0,    30,     0,     0,     0,    25,
      26,     0,     0,     0,     0,     0,     0,     0,     0,     0,
      27,    28,    87,    88,    89,    30,    90,    91,    92,    93,
      94,     0,     0,    95,    96,    97,    98,    99,   100,   101,
     102,   103,   104,   105,   106,     0,     0,     0,   107,   108,
     109,   110,   111,   112,     0,     0,   113,     0,     0,     0,
       0,     0,    87,    88,    89,   218,    90,    91,    92,    93,
      94,     0,     0,    95,    96,    97,    98,    99,   100,   101,
     102,   103,   104,   105,   106,     0,     0,     0,   107,   108,
     109,   110,   111,   112,     0,     0,   113,     0,     0,     0,
       0,     0,    87,    88,    89,   222,    90,    91,    92,    93,
      94,     0,     0,    95,    96,    97,    98,    99,   100,   101,
     102,   103,   104,   105,   106,     0,     0,     0,   107,   108,
     109,   110,   111,   112,     0,     0,   113,     0,     0,     0,
       0,     0,    87,    88,    89,   251,    90,    91,    92,    93,
      94,     0,     0,    95,    96,    97,    98,    99,   100,   101,
     102,   103,   104,   105,   106,     0,     0,     0,   107,   108,
     109,   110,   111,   112,     0,     0,   113,     0,     0,     0,
       0,     0,    87,    88,    89,   272,    90,    91,    92,    93,
      94,     0,     0,    95,    96,    97,    98,    99,   100,   101,
     102,   103,   104,   105,   106,     0,     0,     0,   107,   108,
     109,   110,   111,   112,     0,     0,   113,     0,     0,     0,
       0,     0,    87,    88,    89,   273,    90,    91,    92,    93,
      94,     0,     0,    95,    96,    97,    98,    99,   100,   101,
     102,   103,   104,   105,   106,     0,     0,     0,   107,   108,
     109,   110,   111,   112,     0,     0,   113,     0,     0,     0,
      87,    88,    89,   149,    90,    91,    92,    93,    94,     0,
       0,    95,    96,    97,    98,    99,   100,   101,   102,   103,
     104,   105,   106,     0,     0,     0,   107,   108,   109,   110,
     111,   112,     0,     0,   113,     0,     0,     0,    87,    88,
      89,   210,    90,    91,    92,    93,    94,     0,     0,    95,
      96,    97,    98,    99,   100,   101,   102,   103,   104,   105,
     106,     0,     0,     0,   107,   108,   109,   110,   111,   112,
       0,     0,   113,     0,     0,     0,    87,    88,    89,   212,
      90,    91,    92,    93,    94,     0,     0,    95,    96,    97,
      98,    99,   100,   101,   102,   103,   104,   105,   106,     0,
       0,     0,   107,   108,   109,   110,   111,   112,     0,     0,
     113,     0,     0,     0,    87,    88,    89,   269,    90,    91,
      92,    93,    94,     0,     0,    95,    96,    97,    98,    99,
     100,   101,   102,   103,   104,   105,   106,     0,     0,     0,
     107,   108,   109,   110,   111,   112,    88,    89,   113,    90,
      91,    92,    93,    94,     0,     0,    95,    96,    97,    98,
      99,   100,   101,   102,   103,   104,   105,   106,     0,     0,
       0,   107,   108,   109,   110,   111,   112,     0,    89,   113,
      90,    91,    92,    93,    94,     0,     0,    95,    96,    97,
      98,    99,   100,   101,   102,   103,   104,   105,   106,     0,
       0,     0,   107,   108,   109,   110,   111,   112,     0,     0,
     113,    90,    91,    92,    93,    94,     0,     0,    95,    96,
      97,    98,    99,   100,   101,   102,   103,   104,   105,   106,
       0,     0,     0,   107,   108,   109,   110,   111,   112,     0,
       0,   113,    91,    92,    93,    94,     0,     0,    95,    96,
      97,    98,    99,   100,   101,   102,   103,   104,   105,   106,
       0,     0,     0,   107,   108,   109,   110,   111,   112,     0,
       0,   113,    92,    93,    94,     0,     0,    95,    96,    97,
      98,    99,   100,   101,   102,   103,   104,   105,   106,     0,
       0,     0,   107,   108,   109,   110,   111,   112,    93,    94,
     113,     0,    95,    96,    97,    98,    99,   100,   101,   102,
     103,   104,   105,   106,     0,     0,     0,   107,   108,   109,
     110,   111,   112,    94,     0,   113,    95,    96,    97,    98,
      99,   100,   101,   102,   103,   104,   105,   106,     0,     0,
       0,   107,   108,   109,   110,   111,   112,     0,     0,   113,
     100,   101,   102,   103,   104,   105,   106,     0,     0,     0,
     107,   108,   109,   110,   111,   112,     0,     0,   113,  -103,
    -103,  -103,  -103,  -103,  -103,   106,     0,     0,     0,   107,
     108,   109,   110,   111,   112,     0,     0,   113
};

static const yytype_int16 yycheck[] =
{
       0,     0,   131,     5,    39,    50,    65,     5,   230,    37,
      38,    11,    12,    13,    20,    21,   238,    62,    67,    68,
      57,    67,    68,     5,    24,    74,    67,    68,    74,     5,
      65,    66,    65,    67,    68,    63,     0,    39,    63,    68,
      69,   263,    66,    50,    67,    68,    63,    23,    24,    25,
      26,    58,    59,    65,    66,    62,   278,   279,    63,    63,
      65,    66,    67,    68,    66,    65,    40,    65,    42,    58,
      44,    45,    46,    47,    48,    49,    50,    65,    66,    50,
      54,    55,    56,    57,    58,    59,    19,    66,    62,    22,
      71,    24,    25,    26,    27,    28,    66,    30,    63,    16,
     106,    69,    67,    68,    18,   234,   106,    66,   237,    65,
     116,    66,   118,    65,    66,    72,   116,    71,   118,    69,
      71,    72,    65,   123,    72,   125,    65,   127,    69,   129,
      66,   131,   132,   133,    67,    73,   135,   143,    66,    72,
      39,    50,    69,   143,    67,    68,   146,    56,    57,    58,
      59,    20,    21,    62,    87,    88,    89,    90,    91,    92,
      93,    94,    95,    96,    97,    98,    99,   100,   101,   102,
     103,   104,   105,    69,   107,   108,   109,   110,   111,   112,
     113,   218,   211,   200,   117,   278,   279,    24,   220,     0,
     123,   124,   192,   193,   194,   195,     0,    16,    -1,   199,
     229,    -1,    -1,    -1,    -1,     5,    -1,   236,   225,   226,
     227,   228,    12,    13,    -1,    -1,    -1,    -1,   151,    -1,
     220,    -1,    -1,    23,    24,    25,    26,    27,    -1,    -1,
     230,    -1,    -1,    -1,   234,   272,   273,   237,   238,    -1,
      -1,   240,   241,    -1,    -1,    -1,   275,   276,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   191,    -1,
      -1,    -1,    -1,   263,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   207,    -1,    -1,    -1,   278,   279,
      -1,    -1,   215,    -1,   217,    -1,    -1,   220,   221,     3,
       4,     5,     6,     7,     8,     9,    10,    11,    -1,   232,
      -1,    15,    -1,    17,    -1,    19,    20,    21,    22,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    37,    38,   258,    40,    39,    40,    41,
      42,    43,    44,    45,    46,    47,    48,    49,    50,    -1,
      54,    55,    54,    55,    56,    57,    58,    59,    -1,    -1,
      62,    65,    66,    -1,    -1,    -1,    70,    -1,     1,    73,
       3,     4,     5,     6,     7,     8,     9,    10,    11,    12,
      13,    14,    15,    -1,    17,    -1,    19,    20,    21,    22,
       3,     4,     5,     6,     7,     8,     9,    10,    11,    -1,
      -1,    -1,    -1,    -1,    37,    38,    -1,    40,    -1,     3,
       4,     5,     6,     7,     8,     9,    10,    11,    51,    52,
      -1,    54,    55,    -1,    37,    38,    -1,    40,    -1,    -1,
      -1,    -1,    65,    66,    67,    -1,    -1,    70,    -1,    -1,
      -1,    54,    55,    37,    38,    -1,    40,    -1,    -1,    -1,
      -1,    -1,    65,    66,    -1,    -1,    -1,    70,    71,    -1,
      54,    55,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    65,    66,    -1,    -1,    -1,    70,    71,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    -1,    -1,    -1,
      15,    -1,    17,    -1,    19,    20,    21,    22,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    -1,    -1,    -1,
      -1,    -1,    37,    38,    -1,    40,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    54,
      55,    -1,    37,    38,    -1,    40,    -1,    -1,    -1,    -1,
      65,    66,    -1,    -1,    -1,    70,    -1,    -1,    -1,    54,
      55,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      65,    66,    28,    29,    30,    70,    32,    33,    34,    35,
      36,    -1,    -1,    39,    40,    41,    42,    43,    44,    45,
      46,    47,    48,    49,    50,    -1,    -1,    -1,    54,    55,
      56,    57,    58,    59,    -1,    -1,    62,    -1,    -1,    -1,
      -1,    -1,    28,    29,    30,    71,    32,    33,    34,    35,
      36,    -1,    -1,    39,    40,    41,    42,    43,    44,    45,
      46,    47,    48,    49,    50,    -1,    -1,    -1,    54,    55,
      56,    57,    58,    59,    -1,    -1,    62,    -1,    -1,    -1,
      -1,    -1,    28,    29,    30,    71,    32,    33,    34,    35,
      36,    -1,    -1,    39,    40,    41,    42,    43,    44,    45,
      46,    47,    48,    49,    50,    -1,    -1,    -1,    54,    55,
      56,    57,    58,    59,    -1,    -1,    62,    -1,    -1,    -1,
      -1,    -1,    28,    29,    30,    71,    32,    33,    34,    35,
      36,    -1,    -1,    39,    40,    41,    42,    43,    44,    45,
      46,    47,    48,    49,    50,    -1,    -1,    -1,    54,    55,
      56,    57,    58,    59,    -1,    -1,    62,    -1,    -1,    -1,
      -1,    -1,    28,    29,    30,    71,    32,    33,    34,    35,
      36,    -1,    -1,    39,    40,    41,    42,    43,    44,    45,
      46,    47,    48,    49,    50,    -1,    -1,    -1,    54,    55,
      56,    57,    58,    59,    -1,    -1,    62,    -1,    -1,    -1,
      -1,    -1,    28,    29,    30,    71,    32,    33,    34,    35,
      36,    -1,    -1,    39,    40,    41,    42,    43,    44,    45,
      46,    47,    48,    49,    50,    -1,    -1,    -1,    54,    55,
      56,    57,    58,    59,    -1,    -1,    62,    -1,    -1,    -1,
      28,    29,    30,    69,    32,    33,    34,    35,    36,    -1,
      -1,    39,    40,    41,    42,    43,    44,    45,    46,    47,
      48,    49,    50,    -1,    -1,    -1,    54,    55,    56,    57,
      58,    59,    -1,    -1,    62,    -1,    -1,    -1,    28,    29,
      30,    69,    32,    33,    34,    35,    36,    -1,    -1,    39,
      40,    41,    42,    43,    44,    45,    46,    47,    48,    49,
      50,    -1,    -1,    -1,    54,    55,    56,    57,    58,    59,
      -1,    -1,    62,    -1,    -1,    -1,    28,    29,    30,    69,
      32,    33,    34,    35,    36,    -1,    -1,    39,    40,    41,
      42,    43,    44,    45,    46,    47,    48,    49,    50,    -1,
      -1,    -1,    54,    55,    56,    57,    58,    59,    -1,    -1,
      62,    -1,    -1,    -1,    28,    29,    30,    69,    32,    33,
      34,    35,    36,    -1,    -1,    39,    40,    41,    42,    43,
      44,    45,    46,    47,    48,    49,    50,    -1,    -1,    -1,
      54,    55,    56,    57,    58,    59,    29,    30,    62,    32,
      33,    34,    35,    36,    -1,    -1,    39,    40,    41,    42,
      43,    44,    45,    46,    47,    48,    49,    50,    -1,    -1,
      -1,    54,    55,    56,    57,    58,    59,    -1,    30,    62,
      32,    33,    34,    35,    36,    -1,    -1,    39,    40,    41,
      42,    43,    44,    45,    46,    47,    48,    49,    50,    -1,
      -1,    -1,    54,    55,    56,    57,    58,    59,    -1,    -1,
      62,    32,    33,    34,    35,    36,    -1,    -1,    39,    40,
      41,    42,    43,    44,    45,    46,    47,    48,    49,    50,
      -1,    -1,    -1,    54,    55,    56,    57,    58,    59,    -1,
      -1,    62,    33,    34,    35,    36,    -1,    -1,    39,    40,
      41,    42,    43,    44,    45,    46,    47,    48,    49,    50,
      -1,    -1,    -1,    54,    55,    56,    57,    58,    59,    -1,
      -1,    62,    34,    35,    36,    -1,    -1,    39,    40,    41,
      42,    43,    44,    45,    46,    47,    48,    49,    50,    -1,
      -1,    -1,    54,    55,    56,    57,    58,    59,    35,    36,
      62,    -1,    39,    40,    41,    42,    43,    44,    45,    46,
      47,    48,    49,    50,    -1,    -1,    -1,    54,    55,    56,
      57,    58,    59,    36,    -1,    62,    39,    40,    41,    42,
      43,    44,    45,    46,    47,    48,    49,    50,    -1,    -1,
      -1,    54,    55,    56,    57,    58,    59,    -1,    -1,    62,
      44,    45,    46,    47,    48,    49,    50,    -1,    -1,    -1,
      54,    55,    56,    57,    58,    59,    -1,    -1,    62,    44,
      45,    46,    47,    48,    49,    50,    -1,    -1,    -1,    54,
      55,    56,    57,    58,    59,    -1,    -1,    62
};

/* YYSTOS[STATE-NUM] -- The symbol kind of the accessing symbol of
   state STATE-NUM.  */
static const yytype_int8 yystos[] =
{
       0,     1,     3,     4,     5,     6,     7,     8,     9,    10,
      11,    12,    13,    14,    15,    17,    19,    20,    21,    22,
      37,    38,    40,    51,    52,    54,    55,    65,    66,    67,
      70,    76,    77,    78,    79,    80,    81,    82,    83,    84,
      85,    86,    87,    90,    91,    95,    96,   105,   106,   107,
     110,   111,   113,   115,   116,   117,   118,   119,   120,   122,
      67,    68,   119,   119,   119,    65,   114,    65,   112,   112,
      77,   119,    65,    87,    90,    87,    77,    77,    91,   119,
      77,    77,    77,    77,   121,    77,     0,    28,    29,    30,
      32,    33,    34,    35,    36,    39,    40,    41,    42,    43,
      44,    45,    46,    47,    48,    49,    50,    54,    55,    56,
      57,    58,    59,    62,    37,    38,    63,    66,    63,    67,
      68,    67,    68,    65,    66,    88,    89,    65,    66,   101,
     102,    65,   101,    58,   119,    73,   103,   106,    77,   103,
     103,    88,    77,    63,    67,    68,    63,    67,    68,    69,
      71,    72,    77,    77,    77,    77,    77,    77,    77,    77,
      77,    77,    77,    77,    77,    77,    77,    77,    77,    77,
      77,    90,   119,    77,    77,    77,    77,    77,    77,    77,
      90,    77,    90,    77,    92,    93,    94,   119,    71,    77,
     119,    66,    23,    24,    25,    26,    97,    98,    99,   100,
     119,    71,   119,    66,    97,   119,   119,    16,   104,   106,
      69,    18,    69,    90,   119,    66,    77,    66,    71,    69,
      72,    39,    71,    71,    77,   119,   119,   119,   119,    69,
      72,   119,    39,   101,    65,    71,    69,    65,    73,    77,
      67,    68,    74,   103,    67,    68,    77,    77,    88,    94,
      77,    71,   101,   101,   101,   101,   103,    99,    39,    77,
      97,   103,    97,    27,    95,    96,    99,   108,   109,    69,
     106,   106,    71,    71,    77,    69,    69,    99,    67,    68,
      74,    88,    88,   103,   103,   109,   109
};

/* YYR1[RULE-NUM] -- Symbol kind of the left-hand side of rule RULE-NUM.  */
static const yytype_int8 yyr1[] =
{
       0,    75,    76,    76,    76,    76,    76,    76,    76,    76,
      76,    76,    76,    76,    76,    76,    77,    77,    77,    77,
      77,    77,    77,    77,    77,    77,    77,    77,    77,    77,
      77,    77,    77,    77,    77,    77,    77,    77,    77,    77,
      77,    77,    77,    77,    77,    77,    77,    77,    77,    77,
      77,    77,    77,    77,    77,    78,    79,    80,    81,    82,
      83,    84,    85,    86,    87,    87,    87,    87,    88,    88,
      89,    89,    89,    89,    90,    91,    91,    91,    91,    92,
      92,    93,    93,    94,    94,    95,    95,    96,    96,    97,
      97,    98,    98,    99,    99,    99,    99,   100,   100,   100,
     100,   100,   101,   101,   102,   102,   103,   103,   104,   104,
     104,   104,   104,   104,   105,   105,   105,   105,   105,   105,
     106,   106,   107,   107,   107,   107,   108,   108,   108,   108,
     108,   108,   109,   109,   109,   109,   110,   111,   111,   112,
     113,   114,   115,   116,   117,   118,   118,   119,   120,   121,
     121,   122,   122,   122,   122,   122,   122,   122,   122
};

/* YYR2[RULE-NUM] -- Number of symbols on the right-hand side of rule RULE-NUM.  */
static const yytype_int8 yyr2[] =
{
       0,     2,     1,     1,     2,     2,     2,     2,     3,     3,
       5,     5,     3,     3,     2,     2,     1,     1,     3,     2,
       2,     2,     2,     2,     2,     2,     2,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     2,     5,     7,     7,     0,     1,
       3,     2,     4,     3,     4,     1,     3,     3,     3,     0,
       1,     1,     3,     1,     3,     6,     8,     6,     8,     0,
       1,     1,     3,     1,     3,     2,     4,     2,     3,     3,
       3,     3,     0,     1,     2,     3,     3,     1,     0,     1,
       3,     2,     3,     2,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     3,     0,     1,     3,     2,
       3,     2,     1,     2,     1,     1,     7,     3,     5,     3,
       3,     5,     3,     1,     1,     1,     2,     1,     3,     3,
       1,     1,     1,     1,     1,     1,     1,     1,     1
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
#line 142 "./grammar.y"
            { delete (((*yyvaluep).syntaxElement)); }
#line 1473 "./grammar.tab.c"
        break;

    case YYSYMBOL_arrowAssign: /* arrowAssign  */
#line 143 "./grammar.y"
            { delete (((*yyvaluep).syntaxElement)); }
#line 1479 "./grammar.tab.c"
        break;

    case YYSYMBOL_tildeAssign: /* tildeAssign  */
#line 143 "./grammar.y"
            { delete (((*yyvaluep).syntaxElement)); }
#line 1485 "./grammar.tab.c"
        break;

    case YYSYMBOL_equationAssign: /* equationAssign  */
#line 143 "./grammar.y"
            { delete (((*yyvaluep).syntaxElement)); }
#line 1491 "./grammar.tab.c"
        break;

    case YYSYMBOL_workspaceAssign: /* workspaceAssign  */
#line 143 "./grammar.y"
            { delete (((*yyvaluep).syntaxElement)); }
#line 1497 "./grammar.tab.c"
        break;

    case YYSYMBOL_referenceAssign: /* referenceAssign  */
#line 144 "./grammar.y"
            { delete (((*yyvaluep).syntaxElement)); }
#line 1503 "./grammar.tab.c"
        break;

    case YYSYMBOL_additionAssign: /* additionAssign  */
#line 145 "./grammar.y"
            { delete (((*yyvaluep).syntaxElement)); }
#line 1509 "./grammar.tab.c"
        break;

    case YYSYMBOL_subtractionAssign: /* subtractionAssign  */
#line 145 "./grammar.y"
            { delete (((*yyvaluep).syntaxElement)); }
#line 1515 "./grammar.tab.c"
        break;

    case YYSYMBOL_multiplicationAssign: /* multiplicationAssign  */
#line 145 "./grammar.y"
            { delete (((*yyvaluep).syntaxElement)); }
#line 1521 "./grammar.tab.c"
        break;

    case YYSYMBOL_divisionAssign: /* divisionAssign  */
#line 145 "./grammar.y"
            { delete (((*yyvaluep).syntaxElement)); }
#line 1527 "./grammar.tab.c"
        break;

    case YYSYMBOL_variable: /* variable  */
#line 141 "./grammar.y"
            { delete (((*yyvaluep).syntaxElement)); }
#line 1533 "./grammar.tab.c"
        break;

    case YYSYMBOL_optElements: /* optElements  */
#line 137 "./grammar.y"
            { for (std::list<SyntaxElement*>::iterator it=((*yyvaluep).syntaxElementList)->begin(); it != ((*yyvaluep).syntaxElementList)->end(); it++) { SyntaxElement* theElement = *it; delete theElement; }; delete (((*yyvaluep).syntaxElementList)); }
#line 1539 "./grammar.tab.c"
        break;

    case YYSYMBOL_elementList: /* elementList  */
#line 137 "./grammar.y"
            { for (std::list<SyntaxElement*>::iterator it=((*yyvaluep).syntaxElementList)->begin(); it != ((*yyvaluep).syntaxElementList)->end(); it++) { SyntaxElement* theElement = *it; delete theElement; }; delete (((*yyvaluep).syntaxElementList)); }
#line 1545 "./grammar.tab.c"
        break;

    case YYSYMBOL_fxnCall: /* fxnCall  */
#line 141 "./grammar.y"
            { delete (((*yyvaluep).syntaxFunctionCall)); }
#line 1551 "./grammar.tab.c"
        break;

    case YYSYMBOL_functionCall: /* functionCall  */
#line 141 "./grammar.y"
            { delete (((*yyvaluep).syntaxFunctionCall)); }
#line 1557 "./grammar.tab.c"
        break;

    case YYSYMBOL_optArguments: /* optArguments  */
#line 138 "./grammar.y"
            { for (std::list<SyntaxLabeledExpr*>::iterator it=((*yyvaluep).argumentList)->begin(); it != ((*yyvaluep).argumentList)->end(); it++) { SyntaxLabeledExpr* theElement = *it; delete theElement; }; delete (((*yyvaluep).argumentList)); }
#line 1563 "./grammar.tab.c"
        break;

    case YYSYMBOL_argumentList: /* argumentList  */
#line 138 "./grammar.y"
            { for (std::list<SyntaxLabeledExpr*>::iterator it=((*yyvaluep).argumentList)->begin(); it != ((*yyvaluep).argumentList)->end(); it++) { SyntaxLabeledExpr* theElement = *it; delete theElement; }; delete (((*yyvaluep).argumentList)); }
#line 1569 "./grammar.tab.c"
        break;

    case YYSYMBOL_argument: /* argument  */
#line 141 "./grammar.y"
            { delete (((*yyvaluep).syntaxLabeledExpr)); }
#line 1575 "./grammar.tab.c"
        break;

    case YYSYMBOL_functionDef: /* functionDef  */
#line 147 "./grammar.y"
            { delete (((*yyvaluep).syntaxElement)); }
#line 1581 "./grammar.tab.c"
        break;

    case YYSYMBOL_procedureDef: /* procedureDef  */
#line 147 "./grammar.y"
            { delete (((*yyvaluep).syntaxElement)); }
#line 1587 "./grammar.tab.c"
        break;

    case YYSYMBOL_optFormals: /* optFormals  */
#line 139 "./grammar.y"
            { for (std::list<SyntaxFormal*>::iterator it=((*yyvaluep).formalList)->begin(); it != ((*yyvaluep).formalList)->end(); it++) { SyntaxFormal* theElement = *it; delete theElement; }; delete (((*yyvaluep).formalList)); }
#line 1593 "./grammar.tab.c"
        break;

    case YYSYMBOL_formalList: /* formalList  */
#line 139 "./grammar.y"
            { for (std::list<SyntaxFormal*>::iterator it=((*yyvaluep).formalList)->begin(); it != ((*yyvaluep).formalList)->end(); it++) { SyntaxFormal* theElement = *it; delete theElement; }; delete (((*yyvaluep).formalList)); }
#line 1599 "./grammar.tab.c"
        break;

    case YYSYMBOL_formal: /* formal  */
#line 141 "./grammar.y"
            { delete (((*yyvaluep).syntaxFormal)); }
#line 1605 "./grammar.tab.c"
        break;

    case YYSYMBOL_typeSpec: /* typeSpec  */
#line 140 "./grammar.y"
            { delete (((*yyvaluep).string)); }
#line 1611 "./grammar.tab.c"
        break;

    case YYSYMBOL_optDims: /* optDims  */
#line 140 "./grammar.y"
            { delete (((*yyvaluep).string)); }
#line 1617 "./grammar.tab.c"
        break;

    case YYSYMBOL_dimList: /* dimList  */
#line 140 "./grammar.y"
            { delete (((*yyvaluep).string)); }
#line 1623 "./grammar.tab.c"
        break;

    case YYSYMBOL_stmts: /* stmts  */
#line 137 "./grammar.y"
            { for (std::list<SyntaxElement*>::iterator it=((*yyvaluep).syntaxElementList)->begin(); it != ((*yyvaluep).syntaxElementList)->end(); it++) { SyntaxElement* theElement = *it; delete theElement; }; delete (((*yyvaluep).syntaxElementList)); }
#line 1629 "./grammar.tab.c"
        break;

    case YYSYMBOL_stmtList: /* stmtList  */
#line 137 "./grammar.y"
            { for (std::list<SyntaxElement*>::iterator it=((*yyvaluep).syntaxElementList)->begin(); it != ((*yyvaluep).syntaxElementList)->end(); it++) { SyntaxElement* theElement = *it; delete theElement; }; delete (((*yyvaluep).syntaxElementList)); }
#line 1635 "./grammar.tab.c"
        break;

    case YYSYMBOL_statement: /* statement  */
#line 142 "./grammar.y"
            { delete (((*yyvaluep).syntaxElement)); }
#line 1641 "./grammar.tab.c"
        break;

    case YYSYMBOL_stmt_or_expr: /* stmt_or_expr  */
#line 142 "./grammar.y"
            { delete (((*yyvaluep).syntaxElement)); }
#line 1647 "./grammar.tab.c"
        break;

    case YYSYMBOL_declaration: /* declaration  */
#line 146 "./grammar.y"
            { delete (((*yyvaluep).syntaxElement)); }
#line 1653 "./grammar.tab.c"
        break;

    case YYSYMBOL_memberDefs: /* memberDefs  */
#line 137 "./grammar.y"
            { for (std::list<SyntaxElement*>::iterator it=((*yyvaluep).syntaxElementList)->begin(); it != ((*yyvaluep).syntaxElementList)->end(); it++) { SyntaxElement* theElement = *it; delete theElement; }; delete (((*yyvaluep).syntaxElementList)); }
#line 1659 "./grammar.tab.c"
        break;

    case YYSYMBOL_memberDef: /* memberDef  */
#line 146 "./grammar.y"
            { delete (((*yyvaluep).syntaxElement)); }
#line 1665 "./grammar.tab.c"
        break;

    case YYSYMBOL_classDef: /* classDef  */
#line 146 "./grammar.y"
            { delete (((*yyvaluep).syntaxElement)); }
#line 1671 "./grammar.tab.c"
        break;

    case YYSYMBOL_ifStatement: /* ifStatement  */
#line 148 "./grammar.y"
            { delete (((*yyvaluep).syntaxElement)); }
#line 1677 "./grammar.tab.c"
        break;

    case YYSYMBOL_cond: /* cond  */
#line 149 "./grammar.y"
            { delete (((*yyvaluep).syntaxElement)); }
#line 1683 "./grammar.tab.c"
        break;

    case YYSYMBOL_forStatement: /* forStatement  */
#line 148 "./grammar.y"
            { delete (((*yyvaluep).syntaxElement)); }
#line 1689 "./grammar.tab.c"
        break;

    case YYSYMBOL_forCond: /* forCond  */
#line 149 "./grammar.y"
            { delete (((*yyvaluep).syntaxElement)); }
#line 1695 "./grammar.tab.c"
        break;

    case YYSYMBOL_whileStatement: /* whileStatement  */
#line 148 "./grammar.y"
            { delete (((*yyvaluep).syntaxElement)); }
#line 1701 "./grammar.tab.c"
        break;

    case YYSYMBOL_nextStatement: /* nextStatement  */
#line 150 "./grammar.y"
            { delete (((*yyvaluep).syntaxElement)); }
#line 1707 "./grammar.tab.c"
        break;

    case YYSYMBOL_breakStatement: /* breakStatement  */
#line 150 "./grammar.y"
            { delete (((*yyvaluep).syntaxElement)); }
#line 1713 "./grammar.tab.c"
        break;

    case YYSYMBOL_returnStatement: /* returnStatement  */
#line 149 "./grammar.y"
            { delete (((*yyvaluep).syntaxElement)); }
#line 1719 "./grammar.tab.c"
        break;

    case YYSYMBOL_identifier: /* identifier  */
#line 140 "./grammar.y"
            { delete (((*yyvaluep).string)); }
#line 1725 "./grammar.tab.c"
        break;

    case YYSYMBOL_vector: /* vector  */
#line 138 "./grammar.y"
            { for (std::list<SyntaxLabeledExpr*>::iterator it=((*yyvaluep).argumentList)->begin(); it != ((*yyvaluep).argumentList)->end(); it++) { SyntaxLabeledExpr* theElement = *it; delete theElement; }; delete (((*yyvaluep).argumentList)); }
#line 1731 "./grammar.tab.c"
        break;

    case YYSYMBOL_vectorList: /* vectorList  */
#line 138 "./grammar.y"
            { for (std::list<SyntaxLabeledExpr*>::iterator it=((*yyvaluep).argumentList)->begin(); it != ((*yyvaluep).argumentList)->end(); it++) { SyntaxLabeledExpr* theElement = *it; delete theElement; }; delete (((*yyvaluep).argumentList)); }
#line 1737 "./grammar.tab.c"
        break;

    case YYSYMBOL_constant: /* constant  */
#line 141 "./grammar.y"
            { delete (((*yyvaluep).syntaxElement)); }
#line 1743 "./grammar.tab.c"
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
#line 236 "./grammar.y"
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Bison encountered end_of_input; ignored\n");
#endif
                    return 0;
                }
#line 2043 "./grammar.tab.c"
    break;

  case 3: /* prog: '\n'  */
#line 243 "./grammar.y"
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Bison encountered newline; ignored\n");
#endif
                    return 0;
                }
#line 2054 "./grammar.tab.c"
    break;

  case 4: /* prog: stmt_or_expr '\n'  */
#line 250 "./grammar.y"
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Bison trying to execute statement or expression\n");
#endif
                    int rv = parser.execute((yyvsp[-1].syntaxElement), *executionEnvironment);
                    delete (yyvsp[-1].syntaxElement);
                    return rv;
                }
#line 2067 "./grammar.tab.c"
    break;

  case 5: /* prog: stmt_or_expr ';'  */
#line 259 "./grammar.y"
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Bison trying to execute statement or expression\n");
#endif
                    int rv =  parser.execute((yyvsp[-1].syntaxElement), *executionEnvironment);
                    delete (yyvsp[-1].syntaxElement);
                    return rv;
                }
#line 2080 "./grammar.tab.c"
    break;

  case 6: /* prog: declaration '\n'  */
#line 268 "./grammar.y"
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Bison trying to execute declaration\n");
#endif
                    int rv =  parser.execute((yyvsp[-1].syntaxElement), *executionEnvironment);
                    delete (yyvsp[-1].syntaxElement);
                    return rv;
                }
#line 2093 "./grammar.tab.c"
    break;

  case 7: /* prog: declaration ';'  */
#line 277 "./grammar.y"
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Bison trying to execute declaration\n");
#endif
                    int rv =  parser.execute((yyvsp[-1].syntaxElement), *executionEnvironment);
                    delete (yyvsp[-1].syntaxElement);
                    return rv;
                }
#line 2106 "./grammar.tab.c"
    break;

  case 8: /* prog: '?' identifier '\n'  */
#line 286 "./grammar.y"
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Bison trying to get help for symbol\n");
#endif
                    int rv =  parser.help(*(yyvsp[-1].string));
                    delete (yyvsp[-1].string);
                    return rv;
                }
#line 2119 "./grammar.tab.c"
    break;

  case 9: /* prog: '?' identifier ';'  */
#line 295 "./grammar.y"
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Bison trying to get help for symbol\n");
#endif
                    int rv =  parser.help(*(yyvsp[-1].string));
                    delete (yyvsp[-1].string);
                    return rv;
                }
#line 2132 "./grammar.tab.c"
    break;

  case 10: /* prog: '?' identifier '.' identifier '\n'  */
#line 304 "./grammar.y"
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Bison trying to get help for symbol\n");
#endif
                    int rv =  parser.help(*(yyvsp[-3].string), *(yyvsp[-1].string));
                    delete (yyvsp[-3].string);
                    delete (yyvsp[-1].string);
                    return rv;
                }
#line 2146 "./grammar.tab.c"
    break;

  case 11: /* prog: '?' identifier '.' identifier ';'  */
#line 314 "./grammar.y"
                {
#ifdef DEBUG_BISON_FLEX
                printf("Bison trying to get help for symbol\n");
#endif
                int rv =  parser.help(*(yyvsp[-3].string), *(yyvsp[-1].string));
                delete (yyvsp[-3].string);
                delete (yyvsp[-1].string);
                return rv;
                }
#line 2160 "./grammar.tab.c"
    break;

  case 12: /* prog: '?' functionCall '\n'  */
#line 324 "./grammar.y"
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Bison trying to get help for function call\n");
#endif
                    int rv =  parser.help((yyvsp[-1].syntaxFunctionCall));
                    delete (yyvsp[-1].syntaxFunctionCall);
                    return rv;
                }
#line 2173 "./grammar.tab.c"
    break;

  case 13: /* prog: '?' functionCall ';'  */
#line 333 "./grammar.y"
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Bison trying to get help for function call\n");
#endif
                    int rv =  parser.help((yyvsp[-1].syntaxFunctionCall));
                    delete (yyvsp[-1].syntaxFunctionCall);
                    return rv;
                }
#line 2186 "./grammar.tab.c"
    break;

  case 14: /* prog: error '\n'  */
#line 342 "./grammar.y"
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Bison error when reading line %d position %d to line %d position %d\n", (yyloc).first_line, (yyloc).first_column, (yyloc).last_line, (yyloc).last_column);
#endif
                    YYABORT;
                }
#line 2197 "./grammar.tab.c"
    break;

  case 15: /* prog: error ';'  */
#line 349 "./grammar.y"
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Bison error when reading line %d position %d to line %d position %d\n", (yyloc).first_line, (yyloc).first_column, (yyloc).last_line, (yyloc).last_column);
#endif
                    YYABORT;
                }
#line 2208 "./grammar.tab.c"
    break;

  case 16: /* expression: constant  */
#line 357 "./grammar.y"
                                            { (yyval.syntaxElement) = (yyvsp[0].syntaxElement); }
#line 2214 "./grammar.tab.c"
    break;

  case 17: /* expression: vector  */
#line 359 "./grammar.y"
                                            { (yyval.syntaxElement) = new SyntaxFunctionCall("v", (yyvsp[0].argumentList)); }
#line 2220 "./grammar.tab.c"
    break;

  case 18: /* expression: '(' expression ')'  */
#line 361 "./grammar.y"
                                            { (yyval.syntaxElement) = (yyvsp[-1].syntaxElement); }
#line 2226 "./grammar.tab.c"
    break;

  case 19: /* expression: '-' expression  */
#line 363 "./grammar.y"
                                            { (yyval.syntaxElement) = new SyntaxUnaryExpr(SyntaxUnaryExpr::UMinus, (yyvsp[0].syntaxElement)); }
#line 2232 "./grammar.tab.c"
    break;

  case 20: /* expression: '+' expression  */
#line 364 "./grammar.y"
                                            { (yyval.syntaxElement) = new SyntaxUnaryExpr(SyntaxUnaryExpr::UPlus, (yyvsp[0].syntaxElement)); }
#line 2238 "./grammar.tab.c"
    break;

  case 21: /* expression: '!' expression  */
#line 365 "./grammar.y"
                                            { (yyval.syntaxElement) = new SyntaxUnaryExpr(SyntaxUnaryExpr::UNot, (yyvsp[0].syntaxElement)); }
#line 2244 "./grammar.tab.c"
    break;

  case 22: /* expression: AND expression  */
#line 366 "./grammar.y"
                                            { (yyval.syntaxElement) = new SyntaxUnaryExpr(SyntaxUnaryExpr::UAnd, (yyvsp[0].syntaxElement)); }
#line 2250 "./grammar.tab.c"
    break;

  case 23: /* expression: DECREMENT variable  */
#line 368 "./grammar.y"
                                            { (yyval.syntaxElement) = new SyntaxDecrement( (yyvsp[0].syntaxElement), false ); }
#line 2256 "./grammar.tab.c"
    break;

  case 24: /* expression: variable DECREMENT  */
#line 369 "./grammar.y"
                                            { (yyval.syntaxElement) = new SyntaxDecrement( (yyvsp[-1].syntaxElement), true ); }
#line 2262 "./grammar.tab.c"
    break;

  case 25: /* expression: INCREMENT variable  */
#line 370 "./grammar.y"
                                            { (yyval.syntaxElement) = new SyntaxIncrement( (yyvsp[0].syntaxElement), false ); }
#line 2268 "./grammar.tab.c"
    break;

  case 26: /* expression: variable INCREMENT  */
#line 371 "./grammar.y"
                                            { (yyval.syntaxElement) = new SyntaxIncrement( (yyvsp[-1].syntaxElement), true ); }
#line 2274 "./grammar.tab.c"
    break;

  case 27: /* expression: expression ':' expression  */
#line 373 "./grammar.y"
                                            { (yyval.syntaxElement) = new SyntaxBinaryExpr(SyntaxBinaryExpr::Range, (yyvsp[-2].syntaxElement), (yyvsp[0].syntaxElement)); }
#line 2280 "./grammar.tab.c"
    break;

  case 28: /* expression: expression '+' expression  */
#line 375 "./grammar.y"
                                            { (yyval.syntaxElement) = new SyntaxBinaryExpr(SyntaxBinaryExpr::Add, (yyvsp[-2].syntaxElement), (yyvsp[0].syntaxElement)); }
#line 2286 "./grammar.tab.c"
    break;

  case 29: /* expression: expression '-' expression  */
#line 376 "./grammar.y"
                                            { (yyval.syntaxElement) = new SyntaxBinaryExpr(SyntaxBinaryExpr::Sub, (yyvsp[-2].syntaxElement), (yyvsp[0].syntaxElement)); }
#line 2292 "./grammar.tab.c"
    break;

  case 30: /* expression: expression '*' expression  */
#line 377 "./grammar.y"
                                            { (yyval.syntaxElement) = new SyntaxBinaryExpr(SyntaxBinaryExpr::Mul, (yyvsp[-2].syntaxElement), (yyvsp[0].syntaxElement)); }
#line 2298 "./grammar.tab.c"
    break;

  case 31: /* expression: expression '/' expression  */
#line 378 "./grammar.y"
                                            { (yyval.syntaxElement) = new SyntaxBinaryExpr(SyntaxBinaryExpr::Div, (yyvsp[-2].syntaxElement), (yyvsp[0].syntaxElement)); }
#line 2304 "./grammar.tab.c"
    break;

  case 32: /* expression: expression '^' expression  */
#line 379 "./grammar.y"
                                            { (yyval.syntaxElement) = new SyntaxBinaryExpr(SyntaxBinaryExpr::Exp, (yyvsp[-2].syntaxElement), (yyvsp[0].syntaxElement)); }
#line 2310 "./grammar.tab.c"
    break;

  case 33: /* expression: expression '%' expression  */
#line 380 "./grammar.y"
                                            { (yyval.syntaxElement) = new SyntaxBinaryExpr(SyntaxBinaryExpr::Mod, (yyvsp[-2].syntaxElement), (yyvsp[0].syntaxElement)); }
#line 2316 "./grammar.tab.c"
    break;

  case 34: /* expression: expression LT expression  */
#line 382 "./grammar.y"
                                            { (yyval.syntaxElement) = new SyntaxBinaryExpr(SyntaxBinaryExpr::Lt, (yyvsp[-2].syntaxElement), (yyvsp[0].syntaxElement)); }
#line 2322 "./grammar.tab.c"
    break;

  case 35: /* expression: expression LE expression  */
#line 383 "./grammar.y"
                                            { (yyval.syntaxElement) = new SyntaxBinaryExpr(SyntaxBinaryExpr::Le, (yyvsp[-2].syntaxElement), (yyvsp[0].syntaxElement)); }
#line 2328 "./grammar.tab.c"
    break;

  case 36: /* expression: expression EQ expression  */
#line 384 "./grammar.y"
                                            { (yyval.syntaxElement) = new SyntaxBinaryExpr(SyntaxBinaryExpr::Eq, (yyvsp[-2].syntaxElement), (yyvsp[0].syntaxElement)); }
#line 2334 "./grammar.tab.c"
    break;

  case 37: /* expression: expression NE expression  */
#line 385 "./grammar.y"
                                            { (yyval.syntaxElement) = new SyntaxBinaryExpr(SyntaxBinaryExpr::Ne, (yyvsp[-2].syntaxElement), (yyvsp[0].syntaxElement)); }
#line 2340 "./grammar.tab.c"
    break;

  case 38: /* expression: expression GE expression  */
#line 386 "./grammar.y"
                                            { (yyval.syntaxElement) = new SyntaxBinaryExpr(SyntaxBinaryExpr::Ge, (yyvsp[-2].syntaxElement), (yyvsp[0].syntaxElement)); }
#line 2346 "./grammar.tab.c"
    break;

  case 39: /* expression: expression GT expression  */
#line 387 "./grammar.y"
                                            { (yyval.syntaxElement) = new SyntaxBinaryExpr(SyntaxBinaryExpr::Gt, (yyvsp[-2].syntaxElement), (yyvsp[0].syntaxElement)); }
#line 2352 "./grammar.tab.c"
    break;

  case 40: /* expression: expression AND expression  */
#line 389 "./grammar.y"
                                            { (yyval.syntaxElement) = new SyntaxBinaryExpr(SyntaxBinaryExpr::And, (yyvsp[-2].syntaxElement), (yyvsp[0].syntaxElement)); }
#line 2358 "./grammar.tab.c"
    break;

  case 41: /* expression: expression OR expression  */
#line 390 "./grammar.y"
                                            { (yyval.syntaxElement) = new SyntaxBinaryExpr(SyntaxBinaryExpr::Or, (yyvsp[-2].syntaxElement), (yyvsp[0].syntaxElement)); }
#line 2364 "./grammar.tab.c"
    break;

  case 42: /* expression: expression AND2 expression  */
#line 391 "./grammar.y"
                                            { (yyval.syntaxElement) = new SyntaxBinaryExpr(SyntaxBinaryExpr::And2, (yyvsp[-2].syntaxElement), (yyvsp[0].syntaxElement)); }
#line 2370 "./grammar.tab.c"
    break;

  case 43: /* expression: expression OR2 expression  */
#line 392 "./grammar.y"
                                            { (yyval.syntaxElement) = new SyntaxBinaryExpr(SyntaxBinaryExpr::Or2, (yyvsp[-2].syntaxElement), (yyvsp[0].syntaxElement)); }
#line 2376 "./grammar.tab.c"
    break;

  case 44: /* expression: arrowAssign  */
#line 394 "./grammar.y"
                                            { (yyval.syntaxElement) = (yyvsp[0].syntaxElement); }
#line 2382 "./grammar.tab.c"
    break;

  case 45: /* expression: equationAssign  */
#line 395 "./grammar.y"
                                            { (yyval.syntaxElement) = (yyvsp[0].syntaxElement); }
#line 2388 "./grammar.tab.c"
    break;

  case 46: /* expression: tildeAssign  */
#line 396 "./grammar.y"
                                            { (yyval.syntaxElement) = (yyvsp[0].syntaxElement); }
#line 2394 "./grammar.tab.c"
    break;

  case 47: /* expression: workspaceAssign  */
#line 397 "./grammar.y"
                                            { (yyval.syntaxElement) = (yyvsp[0].syntaxElement); }
#line 2400 "./grammar.tab.c"
    break;

  case 48: /* expression: referenceAssign  */
#line 398 "./grammar.y"
                                            { (yyval.syntaxElement) = (yyvsp[0].syntaxElement); }
#line 2406 "./grammar.tab.c"
    break;

  case 49: /* expression: additionAssign  */
#line 400 "./grammar.y"
                                            { (yyval.syntaxElement) = (yyvsp[0].syntaxElement); }
#line 2412 "./grammar.tab.c"
    break;

  case 50: /* expression: subtractionAssign  */
#line 401 "./grammar.y"
                                            { (yyval.syntaxElement) = (yyvsp[0].syntaxElement); }
#line 2418 "./grammar.tab.c"
    break;

  case 51: /* expression: multiplicationAssign  */
#line 402 "./grammar.y"
                                            { (yyval.syntaxElement) = (yyvsp[0].syntaxElement); }
#line 2424 "./grammar.tab.c"
    break;

  case 52: /* expression: divisionAssign  */
#line 403 "./grammar.y"
                                            { (yyval.syntaxElement) = (yyvsp[0].syntaxElement); }
#line 2430 "./grammar.tab.c"
    break;

  case 53: /* expression: functionCall  */
#line 405 "./grammar.y"
                                            { (yyval.syntaxElement) = (yyvsp[0].syntaxFunctionCall); }
#line 2436 "./grammar.tab.c"
    break;

  case 54: /* expression: variable  */
#line 407 "./grammar.y"
                                            { (yyval.syntaxElement) = (yyvsp[0].syntaxElement); }
#line 2442 "./grammar.tab.c"
    break;

  case 55: /* arrowAssign: expression ARROW_ASSIGN expression  */
#line 411 "./grammar.y"
                    { 
#ifdef DEBUG_BISON_FLEX
                        printf("Parser inserting arrow assignment (ARROW_ASSIGN) in syntax tree\n");
#endif
                        (yyval.syntaxElement) = new SyntaxConstantAssignment((yyvsp[-2].syntaxElement), (yyvsp[0].syntaxElement));
                    }
#line 2453 "./grammar.tab.c"
    break;

  case 56: /* tildeAssign: expression TILDE_ASSIGN expression  */
#line 420 "./grammar.y"
                    {
#ifdef DEBUG_BISON_FLEX
                        printf("Parser inserting tilde assignment (TILDE_ASSIGN) in syntax tree\n");
#endif
                        (yyval.syntaxElement) = new SyntaxStochasticAssignment((yyvsp[-2].syntaxElement), (yyvsp[0].syntaxElement));
                    }
#line 2464 "./grammar.tab.c"
    break;

  case 57: /* equationAssign: expression EQUATION_ASSIGN expression  */
#line 429 "./grammar.y"
                    {
#ifdef DEBUG_BISON_FLEX
                        printf("Parser inserting equation assignment (EQUATION_ASSIGN) in syntax tree\n");
#endif
                        (yyval.syntaxElement) = new SyntaxDeterministicAssignment((yyvsp[-2].syntaxElement), (yyvsp[0].syntaxElement)); 
                    }
#line 2475 "./grammar.tab.c"
    break;

  case 58: /* workspaceAssign: expression EQUAL expression  */
#line 438 "./grammar.y"
                    {
#ifdef DEBUG_BISON_FLEX
                        printf("Parser inserting workspace assignment (WORKSPACE_ASSIGN) in syntax tree\n");
#endif
                        (yyval.syntaxElement) = new SyntaxWorkspaceVariableAssignment((yyvsp[-2].syntaxElement), (yyvsp[0].syntaxElement));
                    }
#line 2486 "./grammar.tab.c"
    break;

  case 59: /* referenceAssign: expression REFERENCE_ASSIGN expression  */
#line 447 "./grammar.y"
                    {
#ifdef DEBUG_BISON_FLEX
                        printf("Parser inserting reference assignment (REFERENCE_ASSIGN) in syntax tree\n");
#endif
                        (yyval.syntaxElement) = new SyntaxReferenceAssignment((yyvsp[-2].syntaxElement), (yyvsp[0].syntaxElement));
                    }
#line 2497 "./grammar.tab.c"
    break;

  case 60: /* additionAssign: expression ADDITION_ASSIGN expression  */
#line 456 "./grammar.y"
                    {
#ifdef DEBUG_BISON_FLEX
                        printf("Parser inserting addition assignment (ADDITION_ASSIGN) in syntax tree\n");
#endif
                        (yyval.syntaxElement) = new SyntaxAdditionAssignment((yyvsp[-2].syntaxElement), (yyvsp[0].syntaxElement)); 
                    }
#line 2508 "./grammar.tab.c"
    break;

  case 61: /* subtractionAssign: expression SUBTRACTION_ASSIGN expression  */
#line 465 "./grammar.y"
                        {
#ifdef DEBUG_BISON_FLEX
                            printf("Parser inserting subtraction assignment (SUBTRACTION_ASSIGN) in syntax tree\n");
#endif
                            (yyval.syntaxElement) = new SyntaxSubtractionAssignment((yyvsp[-2].syntaxElement), (yyvsp[0].syntaxElement)); 
                        }
#line 2519 "./grammar.tab.c"
    break;

  case 62: /* multiplicationAssign: expression MULTIPLICATION_ASSIGN expression  */
#line 474 "./grammar.y"
                            {
#ifdef DEBUG_BISON_FLEX
                                printf("Parser inserting multiplication assignment (MULTIPLICATION_ASSIGN) in syntax tree\n");
#endif
                                (yyval.syntaxElement) = new SyntaxMultiplicationAssignment((yyvsp[-2].syntaxElement), (yyvsp[0].syntaxElement)); 
                            }
#line 2530 "./grammar.tab.c"
    break;

  case 63: /* divisionAssign: expression DIVISION_ASSIGN expression  */
#line 483 "./grammar.y"
                    {
#ifdef DEBUG_BISON_FLEX
                        printf("Parser inserting division assignment (DIVISION_ASSIGN) in syntax tree\n");
#endif
                        (yyval.syntaxElement) = new SyntaxDivisionAssignment((yyvsp[-2].syntaxElement), (yyvsp[0].syntaxElement)); 
                    }
#line 2541 "./grammar.tab.c"
    break;

  case 64: /* variable: identifier optElements  */
#line 492 "./grammar.y"
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
#line 2558 "./grammar.tab.c"
    break;

  case 65: /* variable: fxnCall '[' expression ']' optElements  */
#line 505 "./grammar.y"
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Parser inserting variable (FUNCTION_VAR) in syntax tree\n");
#endif
                    (yyval.syntaxElement) = new SyntaxIndexOperation((yyvsp[-4].syntaxFunctionCall),(yyvsp[-2].syntaxElement));
                    for (std::list<SyntaxElement*>::iterator it=(yyvsp[0].syntaxElementList)->begin(); it!=(yyvsp[0].syntaxElementList)->end(); ++it)
                    {
                        (yyval.syntaxElement) = new SyntaxIndexOperation((yyval.syntaxElement),*it);
                    }
                    delete (yyvsp[0].syntaxElementList);
                }
#line 2574 "./grammar.tab.c"
    break;

  case 66: /* variable: '(' expression ')' '[' expression ']' optElements  */
#line 517 "./grammar.y"
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Parser inserting variable (EXPRESSION_VAR) in syntax tree\n");
#endif
                    (yyval.syntaxElement) = new SyntaxIndexOperation((yyvsp[-5].syntaxElement),(yyvsp[-2].syntaxElement));
                    for (std::list<SyntaxElement*>::iterator it=(yyvsp[0].syntaxElementList)->begin(); it!=(yyvsp[0].syntaxElementList)->end(); ++it)
                    {
                        (yyval.syntaxElement) = new SyntaxIndexOperation((yyval.syntaxElement),*it);
                    }
                    delete (yyvsp[0].syntaxElementList);
                }
#line 2590 "./grammar.tab.c"
    break;

  case 67: /* variable: variable '.' fxnCall '[' expression ']' optElements  */
#line 529 "./grammar.y"
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Parser inserting member variable (FUNCTION_VAR) in syntax tree\n");
#endif
                    (yyvsp[-4].syntaxFunctionCall)->setBaseVariable((yyvsp[-6].syntaxElement));
                    (yyval.syntaxElement) = new SyntaxIndexOperation((yyvsp[-4].syntaxFunctionCall),(yyvsp[-2].syntaxElement));
                    for (std::list<SyntaxElement*>::iterator it=(yyvsp[0].syntaxElementList)->begin(); it!=(yyvsp[0].syntaxElementList)->end(); ++it)
                    {
                        (yyval.syntaxElement) = new SyntaxIndexOperation((yyval.syntaxElement),*it);
                    }
                    delete (yyvsp[0].syntaxElementList);
                }
#line 2607 "./grammar.tab.c"
    break;

  case 68: /* optElements: %empty  */
#line 543 "./grammar.y"
                                                { (yyval.syntaxElementList) = new std::list<SyntaxElement*>(); }
#line 2613 "./grammar.tab.c"
    break;

  case 69: /* optElements: elementList  */
#line 544 "./grammar.y"
                                                { (yyval.syntaxElementList) = (yyvsp[0].syntaxElementList); }
#line 2619 "./grammar.tab.c"
    break;

  case 70: /* elementList: '[' expression ']'  */
#line 547 "./grammar.y"
                                                { (yyval.syntaxElementList) = new std::list<SyntaxElement*>(1, (yyvsp[-1].syntaxElement)); }
#line 2625 "./grammar.tab.c"
    break;

  case 71: /* elementList: '[' ']'  */
#line 548 "./grammar.y"
                                                { (yyval.syntaxElementList) = new std::list<SyntaxElement*>(); }
#line 2631 "./grammar.tab.c"
    break;

  case 72: /* elementList: elementList '[' expression ']'  */
#line 549 "./grammar.y"
                                                { (yyvsp[-3].syntaxElementList)->push_back((yyvsp[-1].syntaxElement)); (yyval.syntaxElementList) = (yyvsp[-3].syntaxElementList); }
#line 2637 "./grammar.tab.c"
    break;

  case 73: /* elementList: elementList '[' ']'  */
#line 550 "./grammar.y"
                                                { (yyvsp[-2].syntaxElementList)->push_back( NULL ); (yyval.syntaxElementList) = (yyvsp[-2].syntaxElementList); }
#line 2643 "./grammar.tab.c"
    break;

  case 74: /* fxnCall: identifier '(' optArguments ')'  */
#line 554 "./grammar.y"
                {
                    (yyval.syntaxFunctionCall) = new SyntaxFunctionCall(*(yyvsp[-3].string), (yyvsp[-1].argumentList));
                    delete (yyvsp[-3].string);
                }
#line 2652 "./grammar.tab.c"
    break;

  case 75: /* functionCall: fxnCall  */
#line 561 "./grammar.y"
                    {
#ifdef DEBUG_BISON_FLEX
                        printf("Parser inserting function call in syntax tree\n");
#endif
                        (yyval.syntaxFunctionCall) = (yyvsp[0].syntaxFunctionCall);
                    }
#line 2663 "./grammar.tab.c"
    break;

  case 76: /* functionCall: variable '.' fxnCall  */
#line 568 "./grammar.y"
                    {
#ifdef DEBUG_BISON_FLEX
                        printf("Parser inserting member call in syntax tree\n");
#endif
                        (yyvsp[0].syntaxFunctionCall)->setBaseVariable((yyvsp[-2].syntaxElement));
                        (yyval.syntaxFunctionCall) = (yyvsp[0].syntaxFunctionCall);
                    }
#line 2675 "./grammar.tab.c"
    break;

  case 77: /* functionCall: functionCall '.' fxnCall  */
#line 576 "./grammar.y"
                    {
#ifdef DEBUG_BISON_FLEX
                        printf("Parser inserting member call in syntax tree\n");
#endif
                        (yyvsp[0].syntaxFunctionCall)->setBaseVariable((yyvsp[-2].syntaxFunctionCall));
                        (yyval.syntaxFunctionCall) = (yyvsp[0].syntaxFunctionCall);
                    }
#line 2687 "./grammar.tab.c"
    break;

  case 78: /* functionCall: expression PIPE fxnCall  */
#line 584 "./grammar.y"
                    {
                        (yyval.syntaxFunctionCall) = (yyvsp[0].syntaxFunctionCall);
                        (yyval.syntaxFunctionCall)->pipeAddArg(new SyntaxLabeledExpr ("" , (yyvsp[-2].syntaxElement)));
                    }
#line 2696 "./grammar.tab.c"
    break;

  case 79: /* optArguments: %empty  */
#line 591 "./grammar.y"
                                            { (yyval.argumentList) = new std::list<SyntaxLabeledExpr*>(); }
#line 2702 "./grammar.tab.c"
    break;

  case 80: /* optArguments: argumentList  */
#line 592 "./grammar.y"
                                            { (yyval.argumentList) = (yyvsp[0].argumentList); }
#line 2708 "./grammar.tab.c"
    break;

  case 81: /* argumentList: argument  */
#line 595 "./grammar.y"
                                                { (yyval.argumentList) = new std::list<SyntaxLabeledExpr*>(1,(yyvsp[0].syntaxLabeledExpr)); }
#line 2714 "./grammar.tab.c"
    break;

  case 82: /* argumentList: argumentList ',' argument  */
#line 596 "./grammar.y"
                                                { (yyvsp[-2].argumentList)->push_back((yyvsp[0].syntaxLabeledExpr)); (yyval.argumentList) = (yyvsp[-2].argumentList); }
#line 2720 "./grammar.tab.c"
    break;

  case 83: /* argument: expression  */
#line 600 "./grammar.y"
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Parser inserting unlabeled argument in syntax tree\n");
#endif
                    (yyval.syntaxLabeledExpr) = new SyntaxLabeledExpr( "", (yyvsp[0].syntaxElement));
                }
#line 2731 "./grammar.tab.c"
    break;

  case 84: /* argument: identifier EQUAL expression  */
#line 607 "./grammar.y"
                { 
#ifdef DEBUG_BISON_FLEX
                    printf("Parser inserting labeled argument in syntax tree\n");
#endif
                    (yyval.syntaxLabeledExpr) = new SyntaxLabeledExpr(*(yyvsp[-2].string), (yyvsp[0].syntaxElement));
                    delete (yyvsp[-2].string);
                }
#line 2743 "./grammar.tab.c"
    break;

  case 85: /* functionDef: FUNCTION identifier '(' optFormals ')' stmts  */
#line 617 "./grammar.y"
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Parser inserting function definition in syntax tree\n");
#endif
                    (yyval.syntaxElement) = new SyntaxFunctionDef("", *(yyvsp[-4].string), (yyvsp[-2].formalList), (yyvsp[0].syntaxElementList));
                    delete (yyvsp[-4].string);
                }
#line 2755 "./grammar.tab.c"
    break;

  case 86: /* functionDef: FUNCTION identifier optDims identifier '(' optFormals ')' stmts  */
#line 626 "./grammar.y"
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
#line 2770 "./grammar.tab.c"
    break;

  case 87: /* procedureDef: PROCEDURE identifier '(' optFormals ')' stmts  */
#line 639 "./grammar.y"
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Parser inserting procedure definition in syntax tree\n");
#endif
                    (yyval.syntaxElement) = new SyntaxFunctionDef("", *(yyvsp[-4].string), (yyvsp[-2].formalList), (yyvsp[0].syntaxElementList), true);
                    delete (yyvsp[-4].string);
                }
#line 2782 "./grammar.tab.c"
    break;

  case 88: /* procedureDef: PROCEDURE identifier optDims identifier '(' optFormals ')' stmts  */
#line 648 "./grammar.y"
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
#line 2797 "./grammar.tab.c"
    break;

  case 89: /* optFormals: %empty  */
#line 660 "./grammar.y"
                                        { (yyval.formalList) = new std::list<SyntaxFormal*>(); }
#line 2803 "./grammar.tab.c"
    break;

  case 90: /* optFormals: formalList  */
#line 661 "./grammar.y"
                                        { (yyval.formalList) = (yyvsp[0].formalList); }
#line 2809 "./grammar.tab.c"
    break;

  case 91: /* formalList: formal  */
#line 664 "./grammar.y"
                                        { (yyval.formalList) = new std::list<SyntaxFormal*>(1, (yyvsp[0].syntaxFormal)); }
#line 2815 "./grammar.tab.c"
    break;

  case 92: /* formalList: formalList ',' formal  */
#line 665 "./grammar.y"
                                        { (yyvsp[-2].formalList)->push_back((yyvsp[0].syntaxFormal)); (yyval.formalList) = (yyvsp[-2].formalList); }
#line 2821 "./grammar.tab.c"
    break;

  case 93: /* formal: identifier  */
#line 669 "./grammar.y"
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Inserting labeled formal argument without default in syntax tree\n");
#endif
                    (yyval.syntaxFormal) = new SyntaxFormal(*(yyvsp[0].string), NULL );
                    delete (yyvsp[0].string);
                }
#line 2833 "./grammar.tab.c"
    break;

  case 94: /* formal: identifier EQUAL expression  */
#line 677 "./grammar.y"
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Inserting labeled formal argument with default in syntax tree\n");
#endif
                    (yyval.syntaxFormal) = new SyntaxFormal(*(yyvsp[-2].string), (yyvsp[0].syntaxElement));
                    delete (yyvsp[-2].string);
                }
#line 2845 "./grammar.tab.c"
    break;

  case 95: /* formal: typeSpec identifier  */
#line 685 "./grammar.y"
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Inserting typed labeled formal argument without default in syntax tree\n");
#endif
                    (yyval.syntaxFormal) = new SyntaxFormal(*(yyvsp[-1].string), *(yyvsp[0].string), NULL);
                    delete (yyvsp[-1].string);
                    delete (yyvsp[0].string);
                }
#line 2858 "./grammar.tab.c"
    break;

  case 96: /* formal: typeSpec identifier EQUAL expression  */
#line 694 "./grammar.y"
                {
#ifdef DEBUG_BISON_FLEX
                    printf("Inserting typed labeled formal argument with default in syntax tree\n");
#endif
                    (yyval.syntaxFormal) = new SyntaxFormal(*(yyvsp[-3].string), *(yyvsp[-2].string), (yyvsp[0].syntaxElement));
                    delete (yyvsp[-3].string);
                    delete (yyvsp[-2].string);
                }
#line 2871 "./grammar.tab.c"
    break;

  case 97: /* typeSpec: identifier optDims  */
#line 704 "./grammar.y"
                                                        { (yyvsp[-1].string)->append(*((yyvsp[0].string))); delete (yyvsp[0].string); (yyval.string) = (yyvsp[-1].string); }
#line 2877 "./grammar.tab.c"
    break;

  case 98: /* typeSpec: MOD_CONST identifier optDims  */
#line 705 "./grammar.y"
                                                        { (yyvsp[-1].string)->append(*((yyvsp[0].string))); (yyvsp[-1].string)->insert(0, "const ");           delete (yyvsp[0].string); (yyval.string) = (yyvsp[-1].string); }
#line 2883 "./grammar.tab.c"
    break;

  case 99: /* typeSpec: MOD_DYNAMIC identifier optDims  */
#line 706 "./grammar.y"
                                                        { (yyvsp[-1].string)->append(*((yyvsp[0].string))); (yyvsp[-1].string)->insert(0, "dynamic ");         delete (yyvsp[0].string); (yyval.string) = (yyvsp[-1].string); }
#line 2889 "./grammar.tab.c"
    break;

  case 100: /* typeSpec: MOD_STOCHASTIC identifier optDims  */
#line 707 "./grammar.y"
                                                        { (yyvsp[-1].string)->append(*((yyvsp[0].string))); (yyvsp[-1].string)->insert(0, "stochastic ");      delete (yyvsp[0].string); (yyval.string) = (yyvsp[-1].string); }
#line 2895 "./grammar.tab.c"
    break;

  case 101: /* typeSpec: MOD_DETERMINISTIC identifier optDims  */
#line 708 "./grammar.y"
                                                        { (yyvsp[-1].string)->append(*((yyvsp[0].string))); (yyvsp[-1].string)->insert(0, "deterministic ");   delete (yyvsp[0].string); (yyval.string) = (yyvsp[-1].string); }
#line 2901 "./grammar.tab.c"
    break;

  case 102: /* optDims: %empty  */
#line 711 "./grammar.y"
                                            { (yyval.string) = new std::string(""); }
#line 2907 "./grammar.tab.c"
    break;

  case 103: /* optDims: dimList  */
#line 712 "./grammar.y"
                                            { (yyval.string) = (yyvsp[0].string); }
#line 2913 "./grammar.tab.c"
    break;

  case 104: /* dimList: '[' ']'  */
#line 715 "./grammar.y"
                                            { (yyval.string) = new std::string("[]"); }
#line 2919 "./grammar.tab.c"
    break;

  case 105: /* dimList: dimList '[' ']'  */
#line 716 "./grammar.y"
                                            { (yyvsp[-2].string)->append("[]"); (yyval.string) = (yyvsp[-2].string); }
#line 2925 "./grammar.tab.c"
    break;

  case 106: /* stmts: '{' stmtList '}'  */
#line 719 "./grammar.y"
                                                { (yyval.syntaxElementList) = (yyvsp[-1].syntaxElementList); }
#line 2931 "./grammar.tab.c"
    break;

  case 107: /* stmts: stmt_or_expr  */
#line 721 "./grammar.y"
                {
                    std::list<SyntaxElement*>* stmts = new std::list<SyntaxElement*>();
                    stmts->push_back((yyvsp[0].syntaxElement));
                    (yyval.syntaxElementList) = stmts;
                }
#line 2941 "./grammar.tab.c"
    break;

  case 108: /* stmtList: %empty  */
#line 728 "./grammar.y"
                                            { (yyval.syntaxElementList) = new std::list<SyntaxElement*>(); }
#line 2947 "./grammar.tab.c"
    break;

  case 109: /* stmtList: stmt_or_expr  */
#line 729 "./grammar.y"
                                            { (yyval.syntaxElementList) = new std::list<SyntaxElement*>(1, (yyvsp[0].syntaxElement)); }
#line 2953 "./grammar.tab.c"
    break;

  case 110: /* stmtList: stmtList ';' stmt_or_expr  */
#line 730 "./grammar.y"
                                            { (yyvsp[-2].syntaxElementList)->push_back((yyvsp[0].syntaxElement)); (yyval.syntaxElementList) = (yyvsp[-2].syntaxElementList); }
#line 2959 "./grammar.tab.c"
    break;

  case 111: /* stmtList: stmtList ';'  */
#line 731 "./grammar.y"
                                            { (yyval.syntaxElementList) = (yyvsp[-1].syntaxElementList); }
#line 2965 "./grammar.tab.c"
    break;

  case 112: /* stmtList: stmtList '\n' stmt_or_expr  */
#line 732 "./grammar.y"
                                            { (yyvsp[-2].syntaxElementList)->push_back((yyvsp[0].syntaxElement)); (yyval.syntaxElementList) = (yyvsp[-2].syntaxElementList); }
#line 2971 "./grammar.tab.c"
    break;

  case 113: /* stmtList: stmtList '\n'  */
#line 733 "./grammar.y"
                                            { (yyval.syntaxElementList) = (yyvsp[-1].syntaxElementList); }
#line 2977 "./grammar.tab.c"
    break;

  case 114: /* statement: ifStatement  */
#line 736 "./grammar.y"
                                        { (yyval.syntaxElement) = (yyvsp[0].syntaxElement); }
#line 2983 "./grammar.tab.c"
    break;

  case 115: /* statement: forStatement  */
#line 737 "./grammar.y"
                                        { (yyval.syntaxElement) = (yyvsp[0].syntaxElement); }
#line 2989 "./grammar.tab.c"
    break;

  case 116: /* statement: whileStatement  */
#line 738 "./grammar.y"
                                        { (yyval.syntaxElement) = (yyvsp[0].syntaxElement); }
#line 2995 "./grammar.tab.c"
    break;

  case 117: /* statement: nextStatement  */
#line 739 "./grammar.y"
                                        { (yyval.syntaxElement) = (yyvsp[0].syntaxElement); }
#line 3001 "./grammar.tab.c"
    break;

  case 118: /* statement: breakStatement  */
#line 740 "./grammar.y"
                                        { (yyval.syntaxElement) = (yyvsp[0].syntaxElement); }
#line 3007 "./grammar.tab.c"
    break;

  case 119: /* statement: returnStatement  */
#line 741 "./grammar.y"
                                        { (yyval.syntaxElement) = (yyvsp[0].syntaxElement); }
#line 3013 "./grammar.tab.c"
    break;

  case 120: /* stmt_or_expr: statement  */
#line 744 "./grammar.y"
                                        { (yyval.syntaxElement) = (yyvsp[0].syntaxElement); }
#line 3019 "./grammar.tab.c"
    break;

  case 121: /* stmt_or_expr: expression  */
#line 745 "./grammar.y"
                                        { (yyval.syntaxElement) = (yyvsp[0].syntaxElement); }
#line 3025 "./grammar.tab.c"
    break;

  case 122: /* declaration: classDef  */
#line 748 "./grammar.y"
                                        { (yyval.syntaxElement) = (yyvsp[0].syntaxElement); }
#line 3031 "./grammar.tab.c"
    break;

  case 123: /* declaration: functionDef  */
#line 749 "./grammar.y"
                                        { (yyval.syntaxElement) = (yyvsp[0].syntaxElement); }
#line 3037 "./grammar.tab.c"
    break;

  case 124: /* declaration: procedureDef  */
#line 750 "./grammar.y"
                                        { (yyval.syntaxElement) = (yyvsp[0].syntaxElement); }
#line 3043 "./grammar.tab.c"
    break;

  case 125: /* declaration: identifier optElements identifier  */
#line 752 "./grammar.y"
                    {
#ifdef DEBUG_BISON_FLEX
                        printf("Parser inserting variable declaration in syntax tree\n");
#endif
                        (yyval.syntaxElement) = new SyntaxVariableDecl(*(yyvsp[-2].string), (yyvsp[-1].syntaxElementList), *(yyvsp[0].string));
                        delete (yyvsp[-2].string);
                        delete (yyvsp[0].string);
                    }
#line 3056 "./grammar.tab.c"
    break;

  case 126: /* memberDefs: %empty  */
#line 762 "./grammar.y"
                                                { (yyval.syntaxElementList) = new std::list<SyntaxElement*>(); }
#line 3062 "./grammar.tab.c"
    break;

  case 127: /* memberDefs: memberDef  */
#line 763 "./grammar.y"
                                                { (yyval.syntaxElementList) = new std::list<SyntaxElement*>(1, (yyvsp[0].syntaxElement)); }
#line 3068 "./grammar.tab.c"
    break;

  case 128: /* memberDefs: memberDefs ';' memberDef  */
#line 764 "./grammar.y"
                                                { (yyvsp[-2].syntaxElementList)->push_back((yyvsp[0].syntaxElement)); (yyval.syntaxElementList) = (yyvsp[-2].syntaxElementList); }
#line 3074 "./grammar.tab.c"
    break;

  case 129: /* memberDefs: memberDefs ';'  */
#line 765 "./grammar.y"
                                                { (yyval.syntaxElementList) = (yyvsp[-1].syntaxElementList); }
#line 3080 "./grammar.tab.c"
    break;

  case 130: /* memberDefs: memberDefs '\n' memberDef  */
#line 766 "./grammar.y"
                                                { (yyvsp[-2].syntaxElementList)->push_back((yyvsp[0].syntaxElement)); (yyval.syntaxElementList) = (yyvsp[-2].syntaxElementList); }
#line 3086 "./grammar.tab.c"
    break;

  case 131: /* memberDefs: memberDefs '\n'  */
#line 767 "./grammar.y"
                                                { (yyval.syntaxElementList) = (yyvsp[-1].syntaxElementList); }
#line 3092 "./grammar.tab.c"
    break;

  case 132: /* memberDef: formal  */
#line 770 "./grammar.y"
                                    { (yyval.syntaxElement) = (yyvsp[0].syntaxFormal); }
#line 3098 "./grammar.tab.c"
    break;

  case 133: /* memberDef: PROTECTED formal  */
#line 771 "./grammar.y"
                                    { (yyvsp[0].syntaxFormal)->setIsProtected(); (yyval.syntaxElement) = (yyvsp[0].syntaxFormal); }
#line 3104 "./grammar.tab.c"
    break;

  case 134: /* memberDef: functionDef  */
#line 772 "./grammar.y"
                                    { (yyval.syntaxElement) = (yyvsp[0].syntaxElement); }
#line 3110 "./grammar.tab.c"
    break;

  case 135: /* memberDef: procedureDef  */
#line 773 "./grammar.y"
                                    { (yyval.syntaxElement) = (yyvsp[0].syntaxElement); }
#line 3116 "./grammar.tab.c"
    break;

  case 136: /* classDef: CLASS identifier ':' identifier '{' memberDefs '}'  */
#line 777 "./grammar.y"
                {
#ifdef DEBUG_BISON_FLEX
                printf("Parser inserting class definition (CLASS_DEF) in syntax tree\n");
#endif
                    (yyval.syntaxElement) = new SyntaxClassDef(*(yyvsp[-5].string), *(yyvsp[-3].string), (yyvsp[-1].syntaxElementList));
                    delete (yyvsp[-5].string);
                    delete (yyvsp[-3].string);
                }
#line 3129 "./grammar.tab.c"
    break;

  case 137: /* ifStatement: IF cond stmts  */
#line 787 "./grammar.y"
                                            { (yyval.syntaxElement) = new SyntaxStatement(SyntaxStatement::If, (yyvsp[-1].syntaxElement), (yyvsp[0].syntaxElementList)); }
#line 3135 "./grammar.tab.c"
    break;

  case 138: /* ifStatement: IF cond stmts ELSE stmts  */
#line 788 "./grammar.y"
                                            { (yyval.syntaxElement) = new SyntaxStatement(SyntaxStatement::IfElse, (yyvsp[-3].syntaxElement), (yyvsp[-2].syntaxElementList), (yyvsp[0].syntaxElementList)); }
#line 3141 "./grammar.tab.c"
    break;

  case 139: /* cond: '(' expression ')'  */
#line 790 "./grammar.y"
                                  { (yyval.syntaxElement) = (yyvsp[-1].syntaxElement); }
#line 3147 "./grammar.tab.c"
    break;

  case 140: /* forStatement: FOR forCond stmts  */
#line 793 "./grammar.y"
                                        { (yyval.syntaxElement) = new SyntaxStatement(SyntaxStatement::For, (yyvsp[-1].syntaxElement), (yyvsp[0].syntaxElementList)); }
#line 3153 "./grammar.tab.c"
    break;

  case 141: /* forCond: '(' identifier IN expression ')'  */
#line 796 "./grammar.y"
                                                    { (yyval.syntaxElement) = new SyntaxForLoop(*(yyvsp[-3].string), (yyvsp[-1].syntaxElement)); delete (yyvsp[-3].string); }
#line 3159 "./grammar.tab.c"
    break;

  case 142: /* whileStatement: WHILE cond stmts  */
#line 799 "./grammar.y"
                                        { (yyval.syntaxElement) = new SyntaxStatement(SyntaxStatement::While, (yyvsp[-1].syntaxElement), (yyvsp[0].syntaxElementList)); }
#line 3165 "./grammar.tab.c"
    break;

  case 143: /* nextStatement: NEXT  */
#line 802 "./grammar.y"
                            { (yyval.syntaxElement) = new SyntaxStatement(SyntaxStatement::Next); }
#line 3171 "./grammar.tab.c"
    break;

  case 144: /* breakStatement: BREAK  */
#line 805 "./grammar.y"
                            { (yyval.syntaxElement) = new SyntaxStatement(SyntaxStatement::Break); }
#line 3177 "./grammar.tab.c"
    break;

  case 145: /* returnStatement: RETURN  */
#line 808 "./grammar.y"
                                        { (yyval.syntaxElement) = new SyntaxStatement(SyntaxStatement::Return); }
#line 3183 "./grammar.tab.c"
    break;

  case 146: /* returnStatement: RETURN expression  */
#line 809 "./grammar.y"
                                        { (yyval.syntaxElement) = new SyntaxStatement(SyntaxStatement::Return, (yyvsp[0].syntaxElement)); }
#line 3189 "./grammar.tab.c"
    break;

  case 147: /* identifier: NAME  */
#line 812 "./grammar.y"
                        { (yyval.string) = new std::string((yyvsp[0].c_string)); }
#line 3195 "./grammar.tab.c"
    break;

  case 148: /* vector: '[' vectorList ']'  */
#line 816 "./grammar.y"
                                        { (yyval.argumentList) = (yyvsp[-1].argumentList); }
#line 3201 "./grammar.tab.c"
    break;

  case 149: /* vectorList: vectorList ',' expression  */
#line 819 "./grammar.y"
                                            { (yyvsp[-2].argumentList)->push_back(new SyntaxLabeledExpr( "", (yyvsp[0].syntaxElement)) ); (yyval.argumentList) = (yyvsp[-2].argumentList); }
#line 3207 "./grammar.tab.c"
    break;

  case 150: /* vectorList: expression  */
#line 821 "./grammar.y"
                {
                (yyval.argumentList) = new std::list<SyntaxLabeledExpr*>(1, new SyntaxLabeledExpr("", (yyvsp[0].syntaxElement)) );
                }
#line 3215 "./grammar.tab.c"
    break;

  case 151: /* constant: FALSE  */
#line 827 "./grammar.y"
                {
                    (yyval.syntaxElement) = new SyntaxConstant(new RlBoolean(false) );
                }
#line 3223 "./grammar.tab.c"
    break;

  case 152: /* constant: TRUE  */
#line 831 "./grammar.y"
                {
                    (yyval.syntaxElement) = new SyntaxConstant(new RlBoolean(true) );
                }
#line 3231 "./grammar.tab.c"
    break;

  case 153: /* constant: RBNULL  */
#line 835 "./grammar.y"
                {
                    (yyval.syntaxElement) = new SyntaxConstant( NULL );
                }
#line 3239 "./grammar.tab.c"
    break;

  case 154: /* constant: RBTAB  */
#line 839 "./grammar.y"
                {
                    (yyval.syntaxElement) = new SyntaxConstant( new RlString("\t") );
                }
#line 3247 "./grammar.tab.c"
    break;

  case 155: /* constant: RBINF  */
#line 843 "./grammar.y"
                {
                    (yyval.syntaxElement) = new SyntaxConstant( new RealPos( RbConstants::Double::inf ) );
                }
#line 3255 "./grammar.tab.c"
    break;

  case 156: /* constant: INT  */
#line 847 "./grammar.y"
                {
                    if ( (yyvsp[0].longIntValue) < 0 ) {
                        (yyval.syntaxElement) = new SyntaxConstant(new Integer((yyvsp[0].longIntValue)) );
                    }
                    else { 
                        (yyval.syntaxElement) = new SyntaxConstant(new Natural((yyvsp[0].longIntValue)) );
                    }
                }
#line 3268 "./grammar.tab.c"
    break;

  case 157: /* constant: STRING  */
#line 856 "./grammar.y"
                {
                    (yyval.syntaxElement) = new SyntaxConstant(new RlString((yyvsp[0].c_string)) );
                }
#line 3276 "./grammar.tab.c"
    break;

  case 158: /* constant: REAL  */
#line 860 "./grammar.y"
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
#line 3339 "./grammar.tab.c"
    break;


#line 3343 "./grammar.tab.c"

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

#line 920 "./grammar.y"



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


