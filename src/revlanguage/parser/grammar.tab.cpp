/* A Bison parser, made by GNU Bison 2.3.  */

/* Skeleton implementation for Bison's Yacc-like parsers in C

   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003, 2004, 2005, 2006
   Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor,
   Boston, MA 02110-1301, USA.  */

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

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "2.3"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Using locations.  */
#define YYLSP_NEEDED 1



/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     REAL = 258,
     INT = 259,
     NAME = 260,
     STRING = 261,
     RBNULL = 262,
     RBTAB = 263,
     FALSE = 264,
     TRUE = 265,
     RBINF = 266,
     FUNCTION = 267,
     PROCEDURE = 268,
     CLASS = 269,
     FOR = 270,
     IN = 271,
     IF = 272,
     ELSE = 273,
     WHILE = 274,
     NEXT = 275,
     BREAK = 276,
     RETURN = 277,
     MOD_CONST = 278,
     MOD_DYNAMIC = 279,
     MOD_STOCHASTIC = 280,
     MOD_DETERMINISTIC = 281,
     PROTECTED = 282,
     ARROW_ASSIGN = 283,
     TILDE_ASSIGN = 284,
     EQUATION_ASSIGN = 285,
     WORKSPACE_ASSIGN = 286,
     REFERENCE_ASSIGN = 287,
     ADDITION_ASSIGN = 288,
     SUBTRACTION_ASSIGN = 289,
     MULTIPLICATION_ASSIGN = 290,
     DIVISION_ASSIGN = 291,
     DECREMENT = 292,
     INCREMENT = 293,
     EQUAL = 294,
     AND = 295,
     OR = 296,
     AND2 = 297,
     OR2 = 298,
     GT = 299,
     GE = 300,
     LT = 301,
     LE = 302,
     EQ = 303,
     NE = 304,
     END_OF_INPUT = 305,
     UNOT = 306,
     UPLUS = 307,
     UMINUS = 308,
     UAND = 309
   };
#endif
/* Tokens.  */
#define REAL 258
#define INT 259
#define NAME 260
#define STRING 261
#define RBNULL 262
#define RBTAB 263
#define FALSE 264
#define TRUE 265
#define RBINF 266
#define FUNCTION 267
#define PROCEDURE 268
#define CLASS 269
#define FOR 270
#define IN 271
#define IF 272
#define ELSE 273
#define WHILE 274
#define NEXT 275
#define BREAK 276
#define RETURN 277
#define MOD_CONST 278
#define MOD_DYNAMIC 279
#define MOD_STOCHASTIC 280
#define MOD_DETERMINISTIC 281
#define PROTECTED 282
#define ARROW_ASSIGN 283
#define TILDE_ASSIGN 284
#define EQUATION_ASSIGN 285
#define WORKSPACE_ASSIGN 286
#define REFERENCE_ASSIGN 287
#define ADDITION_ASSIGN 288
#define SUBTRACTION_ASSIGN 289
#define MULTIPLICATION_ASSIGN 290
#define DIVISION_ASSIGN 291
#define DECREMENT 292
#define INCREMENT 293
#define EQUAL 294
#define AND 295
#define OR 296
#define AND2 297
#define OR2 298
#define GT 299
#define GE 300
#define LT 301
#define LE 302
#define EQ 303
#define NE 304
#define END_OF_INPUT 305
#define UNOT 306
#define UPLUS 307
#define UMINUS 308
#define UAND 309




/* Copy the first part of user declarations.  */
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


/* Enabling traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

/* Enabling the token table.  */
#ifndef YYTOKEN_TABLE
# define YYTOKEN_TABLE 0
#endif

#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
#line 80 "./grammar.y"
{
    char*                                           c_string;
    std::string*                                    string;
    double                                          realValue;
    int                                             intValue;
    bool                                            boolValue;
    RevLanguage::SyntaxElement*                     syntaxElement;
    RevLanguage::SyntaxVariable*                    syntaxVariable;
    RevLanguage::SyntaxFunctionCall*                syntaxFunctionCall;
    RevLanguage::SyntaxLabeledExpr*                 syntaxLabeledExpr;
    RevLanguage::SyntaxFormal*                      syntaxFormal;
    std::list<RevLanguage::SyntaxElement*>*         syntaxElementList;
    std::list<RevLanguage::SyntaxLabeledExpr*>*     argumentList;
    std::list<RevLanguage::SyntaxFormal*>*          formalList;
}
/* Line 193 of yacc.c.  */
#line 298 "./grammar.tab.c"
	YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif

#if ! defined YYLTYPE && ! defined YYLTYPE_IS_DECLARED
typedef struct YYLTYPE
{
  int first_line;
  int first_column;
  int last_line;
  int last_column;
} YYLTYPE;
# define yyltype YYLTYPE /* obsolescent; will be withdrawn */
# define YYLTYPE_IS_DECLARED 1
# define YYLTYPE_IS_TRIVIAL 1
#endif


/* Copy the second part of user declarations.  */


/* Line 216 of yacc.c.  */
#line 323 "./grammar.tab.c"

#ifdef short
# undef short
#endif

#ifdef YYTYPE_UINT8
typedef YYTYPE_UINT8 yytype_uint8;
#else
typedef unsigned char yytype_uint8;
#endif

#ifdef YYTYPE_INT8
typedef YYTYPE_INT8 yytype_int8;
#elif (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
typedef signed char yytype_int8;
#else
typedef short int yytype_int8;
#endif

#ifdef YYTYPE_UINT16
typedef YYTYPE_UINT16 yytype_uint16;
#else
typedef unsigned short int yytype_uint16;
#endif

#ifdef YYTYPE_INT16
typedef YYTYPE_INT16 yytype_int16;
#else
typedef short int yytype_int16;
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif ! defined YYSIZE_T && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned int
# endif
#endif

#define YYSIZE_MAXIMUM ((YYSIZE_T) -1)

#ifndef YY_
# if defined YYENABLE_NLS && YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(msgid) dgettext ("bison-runtime", msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(msgid) msgid
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(e) ((void) (e))
#else
# define YYUSE(e) /* empty */
#endif

/* Identity function, used to suppress warnings about constant conditions.  */
#ifndef lint
# define YYID(n) (n)
#else
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static int
YYID (int i)
#else
static int
YYID (i)
    int i;
#endif
{
  return i;
}
#endif

#if ! defined yyoverflow || YYERROR_VERBOSE

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
#    if ! defined _ALLOCA_H && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#     ifndef _STDLIB_H
#      define _STDLIB_H 1
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's `empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (YYID (0))
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
#  if (defined __cplusplus && ! defined _STDLIB_H \
       && ! ((defined YYMALLOC || defined malloc) \
	     && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef _STDLIB_H
#    define _STDLIB_H 1
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* ! defined yyoverflow || YYERROR_VERBOSE */


#if (! defined yyoverflow \
     && (! defined __cplusplus \
	 || (defined YYLTYPE_IS_TRIVIAL && YYLTYPE_IS_TRIVIAL \
	     && defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yytype_int16 yyss;
  YYSTYPE yyvs;
    YYLTYPE yyls;
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (yytype_int16) + sizeof (YYSTYPE) + sizeof (YYLTYPE)) \
      + 2 * YYSTACK_GAP_MAXIMUM)

/* Copy COUNT objects from FROM to TO.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(To, From, Count) \
      __builtin_memcpy (To, From, (Count) * sizeof (*(From)))
#  else
#   define YYCOPY(To, From, Count)		\
      do					\
	{					\
	  YYSIZE_T yyi;				\
	  for (yyi = 0; yyi < (Count); yyi++)	\
	    (To)[yyi] = (From)[yyi];		\
	}					\
      while (YYID (0))
#  endif
# endif

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack)					\
    do									\
      {									\
	YYSIZE_T yynewbytes;						\
	YYCOPY (&yyptr->Stack, Stack, yysize);				\
	Stack = &yyptr->Stack;						\
	yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
	yyptr += yynewbytes / sizeof (*yyptr);				\
      }									\
    while (YYID (0))

#endif

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  86
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   1135

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  74
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  48
/* YYNRULES -- Number of rules.  */
#define YYNRULES  157
/* YYNRULES -- Number of states.  */
#define YYNSTATES  285

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   309

#define YYTRANSLATE(YYX)						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const yytype_uint8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
      66,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    69,     2,     2,     2,    58,     2,     2,
      64,    68,    55,    53,    71,    54,    62,    56,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,    57,    67,
       2,     2,     2,    51,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,    65,     2,    70,    61,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    72,     2,    73,     2,     2,     2,     2,
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
      45,    46,    47,    48,    49,    50,    52,    59,    60,    63
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const yytype_uint16 yyprhs[] =
{
       0,     0,     3,     5,     7,    10,    13,    16,    19,    23,
      27,    33,    39,    43,    47,    50,    53,    55,    57,    61,
      64,    67,    70,    73,    76,    79,    82,    85,    89,    93,
      97,   101,   105,   109,   113,   117,   121,   125,   129,   133,
     137,   141,   145,   149,   153,   155,   157,   159,   161,   163,
     165,   167,   169,   171,   173,   175,   179,   183,   187,   191,
     195,   199,   203,   207,   211,   214,   220,   228,   236,   237,
     239,   243,   246,   251,   255,   260,   262,   266,   270,   271,
     273,   275,   279,   281,   285,   292,   301,   308,   317,   318,
     320,   322,   326,   328,   332,   335,   340,   343,   347,   351,
     355,   359,   360,   362,   365,   369,   373,   375,   376,   378,
     382,   385,   389,   392,   394,   396,   398,   400,   402,   404,
     406,   408,   410,   412,   414,   418,   419,   421,   425,   428,
     432,   435,   437,   440,   442,   444,   452,   456,   462,   466,
     470,   476,   480,   482,   484,   486,   489,   491,   495,   499,
     501,   503,   505,   507,   509,   511,   513,   515
};

/* YYRHS -- A `-1'-separated list of the rules' RHS.  */
static const yytype_int8 yyrhs[] =
{
      75,     0,    -1,    50,    -1,    66,    -1,   105,    66,    -1,
     105,    67,    -1,   106,    66,    -1,   106,    67,    -1,    51,
     118,    66,    -1,    51,   118,    67,    -1,    51,   118,    62,
     118,    66,    -1,    51,   118,    62,   118,    67,    -1,    51,
      90,    66,    -1,    51,    90,    67,    -1,     1,    66,    -1,
       1,    67,    -1,   121,    -1,   119,    -1,    64,    76,    68,
      -1,    54,    76,    -1,    53,    76,    -1,    69,    76,    -1,
      40,    76,    -1,    37,    86,    -1,    86,    37,    -1,    38,
      86,    -1,    86,    38,    -1,    76,    57,    76,    -1,    76,
      53,    76,    -1,    76,    54,    76,    -1,    76,    55,    76,
      -1,    76,    56,    76,    -1,    76,    61,    76,    -1,    76,
      58,    76,    -1,    76,    46,    76,    -1,    76,    47,    76,
      -1,    76,    48,    76,    -1,    76,    49,    76,    -1,    76,
      45,    76,    -1,    76,    44,    76,    -1,    76,    40,    76,
      -1,    76,    41,    76,    -1,    76,    42,    76,    -1,    76,
      43,    76,    -1,    77,    -1,    79,    -1,    78,    -1,    80,
      -1,    81,    -1,    82,    -1,    83,    -1,    84,    -1,    85,
      -1,    90,    -1,    86,    -1,    76,    28,    76,    -1,    76,
      29,    76,    -1,    76,    30,    76,    -1,    76,    39,    76,
      -1,    76,    32,    76,    -1,    76,    33,    76,    -1,    76,
      34,    76,    -1,    76,    35,    76,    -1,    76,    36,    76,
      -1,   118,    87,    -1,    89,    65,    76,    70,    87,    -1,
      64,    76,    68,    65,    76,    70,    87,    -1,    86,    62,
      89,    65,    76,    70,    87,    -1,    -1,    88,    -1,    65,
      76,    70,    -1,    65,    70,    -1,    88,    65,    76,    70,
      -1,    88,    65,    70,    -1,   118,    64,    91,    68,    -1,
      89,    -1,    86,    62,    89,    -1,    90,    62,    89,    -1,
      -1,    92,    -1,    93,    -1,    92,    71,    93,    -1,    76,
      -1,   118,    39,    76,    -1,    12,   118,    64,    96,    68,
     102,    -1,    12,   118,   100,   118,    64,    96,    68,   102,
      -1,    13,   118,    64,    96,    68,   102,    -1,    13,   118,
     100,   118,    64,    96,    68,   102,    -1,    -1,    97,    -1,
      98,    -1,    97,    71,    98,    -1,   118,    -1,   118,    39,
      76,    -1,    99,   118,    -1,    99,   118,    39,    76,    -1,
     118,   100,    -1,    23,   118,   100,    -1,    24,   118,   100,
      -1,    25,   118,   100,    -1,    26,   118,   100,    -1,    -1,
     101,    -1,    65,    70,    -1,   101,    65,    70,    -1,    72,
     103,    73,    -1,   105,    -1,    -1,   105,    -1,   103,    67,
     105,    -1,   103,    67,    -1,   103,    66,   105,    -1,   103,
      66,    -1,   110,    -1,   112,    -1,   114,    -1,   115,    -1,
     116,    -1,   117,    -1,   104,    -1,    76,    -1,   109,    -1,
      94,    -1,    95,    -1,   118,    87,   118,    -1,    -1,   108,
      -1,   107,    67,   108,    -1,   107,    67,    -1,   107,    66,
     108,    -1,   107,    66,    -1,    98,    -1,    27,    98,    -1,
      94,    -1,    95,    -1,    14,   118,    57,   118,    72,   107,
      73,    -1,    17,   111,   102,    -1,    17,   111,   102,    18,
     102,    -1,    64,    76,    68,    -1,    15,   113,   102,    -1,
      64,   118,    16,    76,    68,    -1,    19,   111,   102,    -1,
      20,    -1,    21,    -1,    22,    -1,    22,    76,    -1,     5,
      -1,    65,   120,    70,    -1,   120,    71,    76,    -1,    76,
      -1,     9,    -1,    10,    -1,     7,    -1,     8,    -1,    11,
      -1,     4,    -1,     6,    -1,     3,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,   234,   234,   241,   248,   257,   266,   275,   284,   293,
     302,   312,   322,   331,   340,   347,   356,   358,   360,   362,
     363,   364,   365,   367,   368,   369,   370,   372,   374,   375,
     376,   377,   378,   379,   381,   382,   383,   384,   385,   386,
     388,   389,   390,   391,   393,   394,   395,   396,   397,   399,
     400,   401,   402,   404,   406,   409,   418,   427,   436,   445,
     454,   463,   472,   481,   490,   503,   515,   527,   542,   543,
     546,   547,   548,   549,   552,   560,   567,   575,   586,   587,
     590,   591,   594,   601,   611,   620,   633,   642,   655,   656,
     659,   660,   663,   671,   679,   688,   699,   700,   701,   702,
     703,   706,   707,   710,   711,   714,   715,   723,   724,   725,
     726,   727,   728,   731,   732,   733,   734,   735,   736,   739,
     740,   743,   744,   745,   746,   757,   758,   759,   760,   761,
     762,   765,   766,   767,   768,   771,   782,   783,   785,   788,
     791,   794,   797,   800,   803,   804,   807,   811,   814,   815,
     821,   825,   829,   833,   837,   841,   850,   854
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || YYTOKEN_TABLE
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "REAL", "INT", "NAME", "STRING",
  "RBNULL", "RBTAB", "FALSE", "TRUE", "RBINF", "FUNCTION", "PROCEDURE",
  "CLASS", "FOR", "IN", "IF", "ELSE", "WHILE", "NEXT", "BREAK", "RETURN",
  "MOD_CONST", "MOD_DYNAMIC", "MOD_STOCHASTIC", "MOD_DETERMINISTIC",
  "PROTECTED", "ARROW_ASSIGN", "TILDE_ASSIGN", "EQUATION_ASSIGN",
  "WORKSPACE_ASSIGN", "REFERENCE_ASSIGN", "ADDITION_ASSIGN",
  "SUBTRACTION_ASSIGN", "MULTIPLICATION_ASSIGN", "DIVISION_ASSIGN",
  "DECREMENT", "INCREMENT", "EQUAL", "AND", "OR", "AND2", "OR2", "GT",
  "GE", "LT", "LE", "EQ", "NE", "END_OF_INPUT", "'?'", "UNOT", "'+'",
  "'-'", "'*'", "'/'", "':'", "'%'", "UPLUS", "UMINUS", "'^'", "'.'",
  "UAND", "'('", "'['", "'\\n'", "';'", "')'", "'!'", "']'", "','", "'{'",
  "'}'", "$accept", "prog", "expression", "arrowAssign", "tildeAssign",
  "equationAssign", "workspaceAssign", "referenceAssign", "additionAssign",
  "subtractionAssign", "multiplicationAssign", "divisionAssign",
  "variable", "optElements", "elementList", "fxnCall", "functionCall",
  "optArguments", "argumentList", "argument", "functionDef",
  "procedureDef", "optFormals", "formalList", "formal", "typeSpec",
  "optDims", "dimList", "stmts", "stmtList", "statement", "stmt_or_expr",
  "declaration", "memberDefs", "memberDef", "classDef", "ifStatement",
  "cond", "forStatement", "forCond", "whileStatement", "nextStatement",
  "breakStatement", "returnStatement", "identifier", "vector",
  "vectorList", "constant", 0
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[YYLEX-NUM] -- Internal token number corresponding to
   token YYLEX-NUM.  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   283,   284,
     285,   286,   287,   288,   289,   290,   291,   292,   293,   294,
     295,   296,   297,   298,   299,   300,   301,   302,   303,   304,
     305,    63,   306,    43,    45,    42,    47,    58,    37,   307,
     308,    94,    46,   309,    40,    91,    10,    59,    41,    33,
      93,    44,   123,   125
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    74,    75,    75,    75,    75,    75,    75,    75,    75,
      75,    75,    75,    75,    75,    75,    76,    76,    76,    76,
      76,    76,    76,    76,    76,    76,    76,    76,    76,    76,
      76,    76,    76,    76,    76,    76,    76,    76,    76,    76,
      76,    76,    76,    76,    76,    76,    76,    76,    76,    76,
      76,    76,    76,    76,    76,    77,    78,    79,    80,    81,
      82,    83,    84,    85,    86,    86,    86,    86,    87,    87,
      88,    88,    88,    88,    89,    90,    90,    90,    91,    91,
      92,    92,    93,    93,    94,    94,    95,    95,    96,    96,
      97,    97,    98,    98,    98,    98,    99,    99,    99,    99,
      99,   100,   100,   101,   101,   102,   102,   103,   103,   103,
     103,   103,   103,   104,   104,   104,   104,   104,   104,   105,
     105,   106,   106,   106,   106,   107,   107,   107,   107,   107,
     107,   108,   108,   108,   108,   109,   110,   110,   111,   112,
     113,   114,   115,   116,   117,   117,   118,   119,   120,   120,
     121,   121,   121,   121,   121,   121,   121,   121
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     1,     1,     2,     2,     2,     2,     3,     3,
       5,     5,     3,     3,     2,     2,     1,     1,     3,     2,
       2,     2,     2,     2,     2,     2,     2,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     2,     5,     7,     7,     0,     1,
       3,     2,     4,     3,     4,     1,     3,     3,     0,     1,
       1,     3,     1,     3,     6,     8,     6,     8,     0,     1,
       1,     3,     1,     3,     2,     4,     2,     3,     3,     3,
       3,     0,     1,     2,     3,     3,     1,     0,     1,     3,
       2,     3,     2,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     3,     0,     1,     3,     2,     3,
       2,     1,     2,     1,     1,     7,     3,     5,     3,     3,
       5,     3,     1,     1,     1,     2,     1,     3,     3,     1,
       1,     1,     1,     1,     1,     1,     1,     1
};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       0,     0,   157,   155,   146,   156,   152,   153,   150,   151,
     154,     0,     0,     0,     0,     0,     0,   142,   143,   144,
       0,     0,     0,     2,     0,     0,     0,     0,     0,     3,
       0,     0,   120,    44,    46,    45,    47,    48,    49,    50,
      51,    52,    54,    75,    53,   122,   123,   119,     0,     0,
     121,   113,   114,   115,   116,   117,   118,    68,    17,    16,
      14,    15,   101,   101,     0,     0,     0,     0,     0,     0,
     145,    68,     0,    23,     0,    25,    22,     0,     0,     0,
      20,    19,     0,   149,     0,    21,     1,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,    24,    26,     0,     0,     0,     4,     5,
       6,     7,    78,     0,    64,    69,    88,     0,     0,   102,
      88,     0,     0,     0,   107,   139,   106,     0,   136,   141,
      64,     0,     0,    12,    13,     0,     8,     9,    18,   147,
       0,    55,    56,    57,    59,    60,    61,    62,    63,    58,
      40,    41,    42,    43,    39,    38,    34,    35,    36,    37,
      28,    29,    30,    31,    27,    33,    32,    76,     0,     0,
      77,    82,     0,    79,    80,    68,    71,     0,   124,     0,
       0,     0,     0,     0,     0,    89,    90,     0,    92,   103,
       0,     0,     0,     0,     0,     0,     0,   108,   138,     0,
       0,     0,     0,     0,   148,     0,    68,    74,     0,     0,
      70,    73,     0,   101,   101,   101,   101,     0,     0,    94,
       0,    96,    88,   104,     0,    88,   125,     0,   112,   110,
     105,   137,    10,    11,     0,     0,    65,    81,    83,    72,
      97,    98,    99,   100,    84,    91,     0,    93,     0,    86,
       0,     0,   133,   134,   131,     0,   126,   140,   111,   109,
      68,    68,    95,     0,     0,   132,   130,   128,   135,    66,
      67,    85,    87,   129,   127
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
      -1,    31,    32,    33,    34,    35,    36,    37,    38,    39,
      40,    41,    42,   140,   125,    43,    44,   182,   183,   184,
     262,   263,   194,   195,   196,   197,   128,   129,   135,   206,
      47,   136,    49,   265,   266,    50,    51,    68,    52,    66,
      53,    54,    55,    56,    71,    58,    84,    59
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -221
static const yytype_int16 yypact[] =
{
     353,   -24,  -221,  -221,  -221,  -221,  -221,  -221,  -221,  -221,
    -221,    18,    18,    18,   -47,   -36,   -36,  -221,  -221,   461,
       2,     2,   461,  -221,     2,   461,   461,   461,   461,  -221,
     461,    33,   842,  -221,  -221,  -221,  -221,  -221,  -221,  -221,
    -221,  -221,   -28,   -27,    -3,  -221,  -221,  -221,   -22,     5,
    -221,  -221,  -221,  -221,  -221,  -221,  -221,    10,  -221,  -221,
    -221,  -221,    22,    25,   -21,    18,   283,   461,   283,   283,
     842,    10,   461,    21,   -27,    21,  -221,    37,   -40,    54,
      19,    19,   694,   842,    35,   166,  -221,   461,   461,   461,
     461,   461,   461,   461,   461,   461,   461,   461,   461,   461,
     461,   461,   461,   461,   461,   461,   461,   461,   461,   461,
     461,   461,   461,  -221,  -221,    18,   461,    18,  -221,  -221,
    -221,  -221,   461,    44,    18,    30,    53,    -8,    18,    36,
      53,    18,    18,    87,   441,  -221,  -221,   731,    92,  -221,
    -221,   768,    18,  -221,  -221,    18,  -221,  -221,    58,  -221,
     461,   842,   872,   902,   932,   961,   989,  1013,  1036,   285,
     166,  1056,   166,  1056,  1074,  1074,  1074,  1074,  1074,  1074,
      12,    12,   -26,   -26,    19,    19,    19,    60,    70,   499,
    -221,   842,    73,    72,  -221,   -35,  -221,   538,  -221,   373,
      18,    18,    18,    18,    76,    75,  -221,    18,    -2,  -221,
      88,   113,   116,   121,   114,   461,   -52,  -221,  -221,   283,
      58,    60,    45,   461,   842,   461,   122,  -221,   461,   461,
    -221,  -221,   577,   129,   129,   129,   129,   283,    53,   156,
     461,  -221,    53,  -221,   283,    53,   124,   805,   441,   441,
    -221,  -221,  -221,  -221,   616,   655,  -221,  -221,   842,  -221,
    -221,  -221,  -221,  -221,  -221,  -221,   461,   842,   128,  -221,
     131,    53,  -221,  -221,  -221,   -48,  -221,  -221,  -221,  -221,
     122,   122,   842,   283,   283,  -221,   124,   124,  -221,  -221,
    -221,  -221,  -221,  -221,  -221
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
    -221,  -221,    66,  -221,  -221,  -221,  -221,  -221,  -221,  -221,
    -221,  -221,    40,   -37,  -221,   -15,   176,  -221,  -221,   -16,
     201,   203,  -128,  -221,  -220,  -221,   -17,  -221,   -29,  -221,
    -221,     1,  -221,  -221,  -137,  -221,  -221,   188,  -221,  -221,
    -221,  -221,  -221,  -221,     0,  -221,  -221,  -221
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -102
static const yytype_int16 yytable[] =
{
      57,    48,   202,  -101,   219,    74,    74,     4,   255,   113,
     114,    62,    63,    64,   238,   239,   264,    65,   276,   277,
     124,   240,   117,     4,    79,   278,   143,   144,    67,   122,
     123,   110,   111,    86,   115,   112,   132,   230,   116,   138,
     139,   275,    60,    61,   118,   119,   131,     2,     3,     4,
       5,     6,     7,     8,     9,    10,   264,   264,     4,   117,
      73,    75,   199,   127,    77,   133,    72,   108,   109,   110,
     111,   120,   121,   112,   122,   123,   190,   191,   192,   193,
     112,    20,    21,   142,    22,    70,   126,   127,    76,   130,
     127,    80,    81,    82,    83,   189,    85,    25,    26,   115,
     177,   201,   180,   205,   258,   149,   150,   260,    27,    28,
     209,   242,   243,    30,   186,   178,   145,   178,   122,   123,
     146,   147,   185,   213,   188,   215,   198,   211,   200,     4,
     198,   203,   204,   137,   122,   207,    11,    12,   141,   283,
     284,   217,   178,   218,   227,   212,   228,   190,   191,   192,
     193,   261,   232,   151,   152,   153,   154,   155,   156,   157,
     158,   159,   160,   161,   162,   163,   164,   165,   166,   167,
     168,   169,   170,   171,   172,   173,   174,   175,   176,   246,
     241,   231,   179,   233,   234,   235,   236,   123,   181,   187,
     223,   224,   225,   226,   127,   256,   273,   229,   254,   274,
      78,    45,   247,    46,    69,   259,   250,   251,   252,   253,
     100,   101,   102,   103,   104,   105,   214,     0,   185,   106,
     107,   108,   109,   110,   111,     0,     0,   112,   198,     0,
       0,     0,   198,   279,   280,   198,   198,     0,     0,   268,
     269,     0,     0,     0,   281,   282,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   222,     0,     0,     0,     0,
       0,   198,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   237,     0,     0,     0,     0,   198,   198,     0,   244,
       0,   245,     0,     0,   181,   248,     2,     3,     4,     5,
       6,     7,     8,     9,    10,     0,   257,     0,    14,     0,
      15,     0,    16,    17,    18,    19,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
      20,    21,   272,    22,    95,    96,    97,    98,    99,   100,
     101,   102,   103,   104,   105,     0,    25,    26,   106,   107,
     108,   109,   110,   111,     0,     0,   112,    27,    28,     0,
       0,     0,    30,     0,     1,   134,     2,     3,     4,     5,
       6,     7,     8,     9,    10,    11,    12,    13,    14,     0,
      15,     0,    16,    17,    18,    19,     2,     3,     4,     5,
       6,     7,     8,     9,    10,     0,     0,     0,     0,     0,
      20,    21,     0,    22,     0,     0,     0,     0,     0,     0,
       0,     0,     0,    23,    24,     0,    25,    26,     0,     0,
      20,    21,     0,    22,     0,     0,     0,    27,    28,    29,
       0,     0,    30,     0,     0,     0,    25,    26,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    27,    28,     0,
       0,     0,    30,   221,     2,     3,     4,     5,     6,     7,
       8,     9,    10,     0,     0,     0,    14,     0,    15,     0,
      16,    17,    18,    19,     2,     3,     4,     5,     6,     7,
       8,     9,    10,     0,     0,     0,     0,     0,    20,    21,
       0,    22,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,    25,    26,     0,     0,    20,    21,
       0,    22,     0,     0,     0,    27,    28,     0,     0,     0,
      30,     0,     0,     0,    25,    26,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    27,    28,    87,    88,    89,
      30,    90,    91,    92,    93,    94,     0,     0,    95,    96,
      97,    98,    99,   100,   101,   102,   103,   104,   105,     0,
       0,     0,   106,   107,   108,   109,   110,   111,     0,     0,
     112,     0,     0,     0,     0,     0,    87,    88,    89,   216,
      90,    91,    92,    93,    94,     0,     0,    95,    96,    97,
      98,    99,   100,   101,   102,   103,   104,   105,     0,     0,
       0,   106,   107,   108,   109,   110,   111,     0,     0,   112,
       0,     0,     0,     0,     0,    87,    88,    89,   220,    90,
      91,    92,    93,    94,     0,     0,    95,    96,    97,    98,
      99,   100,   101,   102,   103,   104,   105,     0,     0,     0,
     106,   107,   108,   109,   110,   111,     0,     0,   112,     0,
       0,     0,     0,     0,    87,    88,    89,   249,    90,    91,
      92,    93,    94,     0,     0,    95,    96,    97,    98,    99,
     100,   101,   102,   103,   104,   105,     0,     0,     0,   106,
     107,   108,   109,   110,   111,     0,     0,   112,     0,     0,
       0,     0,     0,    87,    88,    89,   270,    90,    91,    92,
      93,    94,     0,     0,    95,    96,    97,    98,    99,   100,
     101,   102,   103,   104,   105,     0,     0,     0,   106,   107,
     108,   109,   110,   111,     0,     0,   112,     0,     0,     0,
       0,     0,    87,    88,    89,   271,    90,    91,    92,    93,
      94,     0,     0,    95,    96,    97,    98,    99,   100,   101,
     102,   103,   104,   105,     0,     0,     0,   106,   107,   108,
     109,   110,   111,     0,     0,   112,     0,     0,     0,    87,
      88,    89,   148,    90,    91,    92,    93,    94,     0,     0,
      95,    96,    97,    98,    99,   100,   101,   102,   103,   104,
     105,     0,     0,     0,   106,   107,   108,   109,   110,   111,
       0,     0,   112,     0,     0,     0,    87,    88,    89,   208,
      90,    91,    92,    93,    94,     0,     0,    95,    96,    97,
      98,    99,   100,   101,   102,   103,   104,   105,     0,     0,
       0,   106,   107,   108,   109,   110,   111,     0,     0,   112,
       0,     0,     0,    87,    88,    89,   210,    90,    91,    92,
      93,    94,     0,     0,    95,    96,    97,    98,    99,   100,
     101,   102,   103,   104,   105,     0,     0,     0,   106,   107,
     108,   109,   110,   111,     0,     0,   112,     0,     0,     0,
      87,    88,    89,   267,    90,    91,    92,    93,    94,     0,
       0,    95,    96,    97,    98,    99,   100,   101,   102,   103,
     104,   105,     0,     0,     0,   106,   107,   108,   109,   110,
     111,    88,    89,   112,    90,    91,    92,    93,    94,     0,
       0,    95,    96,    97,    98,    99,   100,   101,   102,   103,
     104,   105,     0,     0,     0,   106,   107,   108,   109,   110,
     111,     0,    89,   112,    90,    91,    92,    93,    94,     0,
       0,    95,    96,    97,    98,    99,   100,   101,   102,   103,
     104,   105,     0,     0,     0,   106,   107,   108,   109,   110,
     111,     0,     0,   112,    90,    91,    92,    93,    94,     0,
       0,    95,    96,    97,    98,    99,   100,   101,   102,   103,
     104,   105,     0,     0,     0,   106,   107,   108,   109,   110,
     111,     0,     0,   112,    91,    92,    93,    94,     0,     0,
      95,    96,    97,    98,    99,   100,   101,   102,   103,   104,
     105,     0,     0,     0,   106,   107,   108,   109,   110,   111,
       0,     0,   112,    92,    93,    94,     0,     0,    95,    96,
      97,    98,    99,   100,   101,   102,   103,   104,   105,     0,
       0,     0,   106,   107,   108,   109,   110,   111,    93,    94,
     112,     0,    95,    96,    97,    98,    99,   100,   101,   102,
     103,   104,   105,     0,     0,     0,   106,   107,   108,   109,
     110,   111,    94,     0,   112,    95,    96,    97,    98,    99,
     100,   101,   102,   103,   104,   105,     0,     0,     0,   106,
     107,   108,   109,   110,   111,     0,    96,   112,    98,     0,
     100,   101,   102,   103,   104,   105,     0,     0,     0,   106,
     107,   108,   109,   110,   111,     0,     0,   112,  -102,  -102,
    -102,  -102,  -102,  -102,     0,     0,     0,   106,   107,   108,
     109,   110,   111,     0,     0,   112
};

static const yytype_int16 yycheck[] =
{
       0,     0,   130,     5,    39,    20,    21,     5,   228,    37,
      38,    11,    12,    13,    66,    67,   236,    64,    66,    67,
      57,    73,    62,     5,    24,    73,    66,    67,    64,    64,
      65,    57,    58,     0,    62,    61,    57,    39,    65,    68,
      69,   261,    66,    67,    66,    67,    63,     3,     4,     5,
       6,     7,     8,     9,    10,    11,   276,   277,     5,    62,
      20,    21,    70,    65,    24,    65,    64,    55,    56,    57,
      58,    66,    67,    61,    64,    65,    23,    24,    25,    26,
      61,    37,    38,    62,    40,    19,    64,    65,    22,    64,
      65,    25,    26,    27,    28,    65,    30,    53,    54,    62,
     115,    65,   117,    16,   232,    70,    71,   235,    64,    65,
      18,    66,    67,    69,    70,   115,    62,   117,    64,    65,
      66,    67,   122,    65,   124,    65,   126,   142,   128,     5,
     130,   131,   132,    67,    64,   134,    12,    13,    72,   276,
     277,    68,   142,    71,    68,   145,    71,    23,    24,    25,
      26,    27,    64,    87,    88,    89,    90,    91,    92,    93,
      94,    95,    96,    97,    98,    99,   100,   101,   102,   103,
     104,   105,   106,   107,   108,   109,   110,   111,   112,   216,
     209,   198,   116,    70,    68,    64,    72,    65,   122,   123,
     190,   191,   192,   193,    65,    39,    68,   197,   227,    68,
      24,     0,   218,     0,    16,   234,   223,   224,   225,   226,
      44,    45,    46,    47,    48,    49,   150,    -1,   218,    53,
      54,    55,    56,    57,    58,    -1,    -1,    61,   228,    -1,
      -1,    -1,   232,   270,   271,   235,   236,    -1,    -1,   238,
     239,    -1,    -1,    -1,   273,   274,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   189,    -1,    -1,    -1,    -1,
      -1,   261,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   205,    -1,    -1,    -1,    -1,   276,   277,    -1,   213,
      -1,   215,    -1,    -1,   218,   219,     3,     4,     5,     6,
       7,     8,     9,    10,    11,    -1,   230,    -1,    15,    -1,
      17,    -1,    19,    20,    21,    22,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      37,    38,   256,    40,    39,    40,    41,    42,    43,    44,
      45,    46,    47,    48,    49,    -1,    53,    54,    53,    54,
      55,    56,    57,    58,    -1,    -1,    61,    64,    65,    -1,
      -1,    -1,    69,    -1,     1,    72,     3,     4,     5,     6,
       7,     8,     9,    10,    11,    12,    13,    14,    15,    -1,
      17,    -1,    19,    20,    21,    22,     3,     4,     5,     6,
       7,     8,     9,    10,    11,    -1,    -1,    -1,    -1,    -1,
      37,    38,    -1,    40,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    50,    51,    -1,    53,    54,    -1,    -1,
      37,    38,    -1,    40,    -1,    -1,    -1,    64,    65,    66,
      -1,    -1,    69,    -1,    -1,    -1,    53,    54,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    64,    65,    -1,
      -1,    -1,    69,    70,     3,     4,     5,     6,     7,     8,
       9,    10,    11,    -1,    -1,    -1,    15,    -1,    17,    -1,
      19,    20,    21,    22,     3,     4,     5,     6,     7,     8,
       9,    10,    11,    -1,    -1,    -1,    -1,    -1,    37,    38,
      -1,    40,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    53,    54,    -1,    -1,    37,    38,
      -1,    40,    -1,    -1,    -1,    64,    65,    -1,    -1,    -1,
      69,    -1,    -1,    -1,    53,    54,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    64,    65,    28,    29,    30,
      69,    32,    33,    34,    35,    36,    -1,    -1,    39,    40,
      41,    42,    43,    44,    45,    46,    47,    48,    49,    -1,
      -1,    -1,    53,    54,    55,    56,    57,    58,    -1,    -1,
      61,    -1,    -1,    -1,    -1,    -1,    28,    29,    30,    70,
      32,    33,    34,    35,    36,    -1,    -1,    39,    40,    41,
      42,    43,    44,    45,    46,    47,    48,    49,    -1,    -1,
      -1,    53,    54,    55,    56,    57,    58,    -1,    -1,    61,
      -1,    -1,    -1,    -1,    -1,    28,    29,    30,    70,    32,
      33,    34,    35,    36,    -1,    -1,    39,    40,    41,    42,
      43,    44,    45,    46,    47,    48,    49,    -1,    -1,    -1,
      53,    54,    55,    56,    57,    58,    -1,    -1,    61,    -1,
      -1,    -1,    -1,    -1,    28,    29,    30,    70,    32,    33,
      34,    35,    36,    -1,    -1,    39,    40,    41,    42,    43,
      44,    45,    46,    47,    48,    49,    -1,    -1,    -1,    53,
      54,    55,    56,    57,    58,    -1,    -1,    61,    -1,    -1,
      -1,    -1,    -1,    28,    29,    30,    70,    32,    33,    34,
      35,    36,    -1,    -1,    39,    40,    41,    42,    43,    44,
      45,    46,    47,    48,    49,    -1,    -1,    -1,    53,    54,
      55,    56,    57,    58,    -1,    -1,    61,    -1,    -1,    -1,
      -1,    -1,    28,    29,    30,    70,    32,    33,    34,    35,
      36,    -1,    -1,    39,    40,    41,    42,    43,    44,    45,
      46,    47,    48,    49,    -1,    -1,    -1,    53,    54,    55,
      56,    57,    58,    -1,    -1,    61,    -1,    -1,    -1,    28,
      29,    30,    68,    32,    33,    34,    35,    36,    -1,    -1,
      39,    40,    41,    42,    43,    44,    45,    46,    47,    48,
      49,    -1,    -1,    -1,    53,    54,    55,    56,    57,    58,
      -1,    -1,    61,    -1,    -1,    -1,    28,    29,    30,    68,
      32,    33,    34,    35,    36,    -1,    -1,    39,    40,    41,
      42,    43,    44,    45,    46,    47,    48,    49,    -1,    -1,
      -1,    53,    54,    55,    56,    57,    58,    -1,    -1,    61,
      -1,    -1,    -1,    28,    29,    30,    68,    32,    33,    34,
      35,    36,    -1,    -1,    39,    40,    41,    42,    43,    44,
      45,    46,    47,    48,    49,    -1,    -1,    -1,    53,    54,
      55,    56,    57,    58,    -1,    -1,    61,    -1,    -1,    -1,
      28,    29,    30,    68,    32,    33,    34,    35,    36,    -1,
      -1,    39,    40,    41,    42,    43,    44,    45,    46,    47,
      48,    49,    -1,    -1,    -1,    53,    54,    55,    56,    57,
      58,    29,    30,    61,    32,    33,    34,    35,    36,    -1,
      -1,    39,    40,    41,    42,    43,    44,    45,    46,    47,
      48,    49,    -1,    -1,    -1,    53,    54,    55,    56,    57,
      58,    -1,    30,    61,    32,    33,    34,    35,    36,    -1,
      -1,    39,    40,    41,    42,    43,    44,    45,    46,    47,
      48,    49,    -1,    -1,    -1,    53,    54,    55,    56,    57,
      58,    -1,    -1,    61,    32,    33,    34,    35,    36,    -1,
      -1,    39,    40,    41,    42,    43,    44,    45,    46,    47,
      48,    49,    -1,    -1,    -1,    53,    54,    55,    56,    57,
      58,    -1,    -1,    61,    33,    34,    35,    36,    -1,    -1,
      39,    40,    41,    42,    43,    44,    45,    46,    47,    48,
      49,    -1,    -1,    -1,    53,    54,    55,    56,    57,    58,
      -1,    -1,    61,    34,    35,    36,    -1,    -1,    39,    40,
      41,    42,    43,    44,    45,    46,    47,    48,    49,    -1,
      -1,    -1,    53,    54,    55,    56,    57,    58,    35,    36,
      61,    -1,    39,    40,    41,    42,    43,    44,    45,    46,
      47,    48,    49,    -1,    -1,    -1,    53,    54,    55,    56,
      57,    58,    36,    -1,    61,    39,    40,    41,    42,    43,
      44,    45,    46,    47,    48,    49,    -1,    -1,    -1,    53,
      54,    55,    56,    57,    58,    -1,    40,    61,    42,    -1,
      44,    45,    46,    47,    48,    49,    -1,    -1,    -1,    53,
      54,    55,    56,    57,    58,    -1,    -1,    61,    44,    45,
      46,    47,    48,    49,    -1,    -1,    -1,    53,    54,    55,
      56,    57,    58,    -1,    -1,    61
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,     1,     3,     4,     5,     6,     7,     8,     9,    10,
      11,    12,    13,    14,    15,    17,    19,    20,    21,    22,
      37,    38,    40,    50,    51,    53,    54,    64,    65,    66,
      69,    75,    76,    77,    78,    79,    80,    81,    82,    83,
      84,    85,    86,    89,    90,    94,    95,   104,   105,   106,
     109,   110,   112,   114,   115,   116,   117,   118,   119,   121,
      66,    67,   118,   118,   118,    64,   113,    64,   111,   111,
      76,   118,    64,    86,    89,    86,    76,    86,    90,   118,
      76,    76,    76,    76,   120,    76,     0,    28,    29,    30,
      32,    33,    34,    35,    36,    39,    40,    41,    42,    43,
      44,    45,    46,    47,    48,    49,    53,    54,    55,    56,
      57,    58,    61,    37,    38,    62,    65,    62,    66,    67,
      66,    67,    64,    65,    87,    88,    64,    65,   100,   101,
      64,   100,    57,   118,    72,   102,   105,    76,   102,   102,
      87,    76,    62,    66,    67,    62,    66,    67,    68,    70,
      71,    76,    76,    76,    76,    76,    76,    76,    76,    76,
      76,    76,    76,    76,    76,    76,    76,    76,    76,    76,
      76,    76,    76,    76,    76,    76,    76,    89,   118,    76,
      89,    76,    91,    92,    93,   118,    70,    76,   118,    65,
      23,    24,    25,    26,    96,    97,    98,    99,   118,    70,
     118,    65,    96,   118,   118,    16,   103,   105,    68,    18,
      68,    89,   118,    65,    76,    65,    70,    68,    71,    39,
      70,    70,    76,   118,   118,   118,   118,    68,    71,   118,
      39,   100,    64,    70,    68,    64,    72,    76,    66,    67,
      73,   102,    66,    67,    76,    76,    87,    93,    76,    70,
     100,   100,   100,   100,   102,    98,    39,    76,    96,   102,
      96,    27,    94,    95,    98,   107,   108,    68,   105,   105,
      70,    70,    76,    68,    68,    98,    66,    67,    73,    87,
      87,   102,   102,   108,   108
};

#define yyerrok		(yyerrstatus = 0)
#define yyclearin	(yychar = YYEMPTY)
#define YYEMPTY		(-2)
#define YYEOF		0

#define YYACCEPT	goto yyacceptlab
#define YYABORT		goto yyabortlab
#define YYERROR		goto yyerrorlab


/* Like YYERROR except do call yyerror.  This remains here temporarily
   to ease the transition to the new meaning of YYERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  */

#define YYFAIL		goto yyerrlab

#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)					\
do								\
  if (yychar == YYEMPTY && yylen == 1)				\
    {								\
      yychar = (Token);						\
      yylval = (Value);						\
      yytoken = YYTRANSLATE (yychar);				\
      YYPOPSTACK (1);						\
      goto yybackup;						\
    }								\
  else								\
    {								\
      yyerror (YY_("syntax error: cannot back up")); \
      YYERROR;							\
    }								\
while (YYID (0))


#define YYTERROR	1
#define YYERRCODE	256


/* YYLLOC_DEFAULT -- Set CURRENT to span from RHS[1] to RHS[N].
   If N is 0, then set CURRENT to the empty location which ends
   the previous symbol: RHS[0] (always defined).  */

#define YYRHSLOC(Rhs, K) ((Rhs)[K])
#ifndef YYLLOC_DEFAULT
# define YYLLOC_DEFAULT(Current, Rhs, N)				\
    do									\
      if (YYID (N))                                                    \
	{								\
	  (Current).first_line   = YYRHSLOC (Rhs, 1).first_line;	\
	  (Current).first_column = YYRHSLOC (Rhs, 1).first_column;	\
	  (Current).last_line    = YYRHSLOC (Rhs, N).last_line;		\
	  (Current).last_column  = YYRHSLOC (Rhs, N).last_column;	\
	}								\
      else								\
	{								\
	  (Current).first_line   = (Current).last_line   =		\
	    YYRHSLOC (Rhs, 0).last_line;				\
	  (Current).first_column = (Current).last_column =		\
	    YYRHSLOC (Rhs, 0).last_column;				\
	}								\
    while (YYID (0))
#endif


/* YY_LOCATION_PRINT -- Print the location on the stream.
   This macro was not mandated originally: define only if we know
   we won't break user code: when these are the locations we know.  */

#ifndef YY_LOCATION_PRINT
# if defined YYLTYPE_IS_TRIVIAL && YYLTYPE_IS_TRIVIAL
#  define YY_LOCATION_PRINT(File, Loc)			\
     fprintf (File, "%d.%d-%d.%d",			\
	      (Loc).first_line, (Loc).first_column,	\
	      (Loc).last_line,  (Loc).last_column)
# else
#  define YY_LOCATION_PRINT(File, Loc) ((void) 0)
# endif
#endif


/* YYLEX -- calling `yylex' with the right arguments.  */

#ifdef YYLEX_PARAM
# define YYLEX yylex (YYLEX_PARAM)
#else
# define YYLEX yylex ()
#endif

/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)			\
do {						\
  if (yydebug)					\
    YYFPRINTF Args;				\
} while (YYID (0))

# define YY_SYMBOL_PRINT(Title, Type, Value, Location)			  \
do {									  \
  if (yydebug)								  \
    {									  \
      YYFPRINTF (stderr, "%s ", Title);					  \
      yy_symbol_print (stderr,						  \
		  Type, Value, Location); \
      YYFPRINTF (stderr, "\n");						  \
    }									  \
} while (YYID (0))


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep, YYLTYPE const * const yylocationp)
#else
static void
yy_symbol_value_print (yyoutput, yytype, yyvaluep, yylocationp)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
    YYLTYPE const * const yylocationp;
#endif
{
  if (!yyvaluep)
    return;
  YYUSE (yylocationp);
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# else
  YYUSE (yyoutput);
# endif
  switch (yytype)
    {
      default:
	break;
    }
}


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep, YYLTYPE const * const yylocationp)
#else
static void
yy_symbol_print (yyoutput, yytype, yyvaluep, yylocationp)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
    YYLTYPE const * const yylocationp;
#endif
{
  if (yytype < YYNTOKENS)
    YYFPRINTF (yyoutput, "token %s (", yytname[yytype]);
  else
    YYFPRINTF (yyoutput, "nterm %s (", yytname[yytype]);

  YY_LOCATION_PRINT (yyoutput, *yylocationp);
  YYFPRINTF (yyoutput, ": ");
  yy_symbol_value_print (yyoutput, yytype, yyvaluep, yylocationp);
  YYFPRINTF (yyoutput, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_stack_print (yytype_int16 *bottom, yytype_int16 *top)
#else
static void
yy_stack_print (bottom, top)
    yytype_int16 *bottom;
    yytype_int16 *top;
#endif
{
  YYFPRINTF (stderr, "Stack now");
  for (; bottom <= top; ++bottom)
    YYFPRINTF (stderr, " %d", *bottom);
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)				\
do {								\
  if (yydebug)							\
    yy_stack_print ((Bottom), (Top));				\
} while (YYID (0))


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_reduce_print (YYSTYPE *yyvsp, YYLTYPE *yylsp, int yyrule)
#else
static void
yy_reduce_print (yyvsp, yylsp, yyrule)
    YYSTYPE *yyvsp;
    YYLTYPE *yylsp;
    int yyrule;
#endif
{
  int yynrhs = yyr2[yyrule];
  int yyi;
  unsigned long int yylno = yyrline[yyrule];
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
	     yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      fprintf (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr, yyrhs[yyprhs[yyrule] + yyi],
		       &(yyvsp[(yyi + 1) - (yynrhs)])
		       , &(yylsp[(yyi + 1) - (yynrhs)])		       );
      fprintf (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)		\
do {					\
  if (yydebug)				\
    yy_reduce_print (yyvsp, yylsp, Rule); \
} while (YYID (0))

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef	YYINITDEPTH
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



#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined __GLIBC__ && defined _STRING_H
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static YYSIZE_T
yystrlen (const char *yystr)
#else
static YYSIZE_T
yystrlen (yystr)
    const char *yystr;
#endif
{
  YYSIZE_T yylen;
  for (yylen = 0; yystr[yylen]; yylen++)
    continue;
  return yylen;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined __GLIBC__ && defined _STRING_H && defined _GNU_SOURCE
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static char *
yystpcpy (char *yydest, const char *yysrc)
#else
static char *
yystpcpy (yydest, yysrc)
    char *yydest;
    const char *yysrc;
#endif
{
  char *yyd = yydest;
  const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

# ifndef yytnamerr
/* Copy to YYRES the contents of YYSTR after stripping away unnecessary
   quotes and backslashes, so that it's suitable for yyerror.  The
   heuristic is that double-quoting is unnecessary unless the string
   contains an apostrophe, a comma, or backslash (other than
   backslash-backslash).  YYSTR is taken from yytname.  If YYRES is
   null, do not copy; instead, return the length of what the result
   would have been.  */
static YYSIZE_T
yytnamerr (char *yyres, const char *yystr)
{
  if (*yystr == '"')
    {
      YYSIZE_T yyn = 0;
      char const *yyp = yystr;

      for (;;)
	switch (*++yyp)
	  {
	  case '\'':
	  case ',':
	    goto do_not_strip_quotes;

	  case '\\':
	    if (*++yyp != '\\')
	      goto do_not_strip_quotes;
	    /* Fall through.  */
	  default:
	    if (yyres)
	      yyres[yyn] = *yyp;
	    yyn++;
	    break;

	  case '"':
	    if (yyres)
	      yyres[yyn] = '\0';
	    return yyn;
	  }
    do_not_strip_quotes: ;
    }

  if (! yyres)
    return yystrlen (yystr);

  return yystpcpy (yyres, yystr) - yyres;
}
# endif

/* Copy into YYRESULT an error message about the unexpected token
   YYCHAR while in state YYSTATE.  Return the number of bytes copied,
   including the terminating null byte.  If YYRESULT is null, do not
   copy anything; just return the number of bytes that would be
   copied.  As a special case, return 0 if an ordinary "syntax error"
   message will do.  Return YYSIZE_MAXIMUM if overflow occurs during
   size calculation.  */
static YYSIZE_T
yysyntax_error (char *yyresult, int yystate, int yychar)
{
  int yyn = yypact[yystate];

  if (! (YYPACT_NINF < yyn && yyn <= YYLAST))
    return 0;
  else
    {
      int yytype = YYTRANSLATE (yychar);
      YYSIZE_T yysize0 = yytnamerr (0, yytname[yytype]);
      YYSIZE_T yysize = yysize0;
      YYSIZE_T yysize1;
      int yysize_overflow = 0;
      enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
      char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
      int yyx;

# if 0
      /* This is so xgettext sees the translatable formats that are
	 constructed on the fly.  */
      YY_("syntax error, unexpected %s");
      YY_("syntax error, unexpected %s, expecting %s");
      YY_("syntax error, unexpected %s, expecting %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s");
# endif
      char *yyfmt;
      char const *yyf;
      static char const yyunexpected[] = "syntax error, unexpected %s";
      static char const yyexpecting[] = ", expecting %s";
      static char const yyor[] = " or %s";
      char yyformat[sizeof yyunexpected
		    + sizeof yyexpecting - 1
		    + ((YYERROR_VERBOSE_ARGS_MAXIMUM - 2)
		       * (sizeof yyor - 1))];
      char const *yyprefix = yyexpecting;

      /* Start YYX at -YYN if negative to avoid negative indexes in
	 YYCHECK.  */
      int yyxbegin = yyn < 0 ? -yyn : 0;

      /* Stay within bounds of both yycheck and yytname.  */
      int yychecklim = YYLAST - yyn + 1;
      int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
      int yycount = 1;

      yyarg[0] = yytname[yytype];
      yyfmt = yystpcpy (yyformat, yyunexpected);

      for (yyx = yyxbegin; yyx < yyxend; ++yyx)
	if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR)
	  {
	    if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
	      {
		yycount = 1;
		yysize = yysize0;
		yyformat[sizeof yyunexpected - 1] = '\0';
		break;
	      }
	    yyarg[yycount++] = yytname[yyx];
	    yysize1 = yysize + yytnamerr (0, yytname[yyx]);
	    yysize_overflow |= (yysize1 < yysize);
	    yysize = yysize1;
	    yyfmt = yystpcpy (yyfmt, yyprefix);
	    yyprefix = yyor;
	  }

      yyf = YY_(yyformat);
      yysize1 = yysize + yystrlen (yyf);
      yysize_overflow |= (yysize1 < yysize);
      yysize = yysize1;

      if (yysize_overflow)
	return YYSIZE_MAXIMUM;

      if (yyresult)
	{
	  /* Avoid sprintf, as that infringes on the user's name space.
	     Don't have undefined behavior even if the translation
	     produced a string with the wrong number of "%s"s.  */
	  char *yyp = yyresult;
	  int yyi = 0;
	  while ((*yyp = *yyf) != '\0')
	    {
	      if (*yyp == '%' && yyf[1] == 's' && yyi < yycount)
		{
		  yyp += yytnamerr (yyp, yyarg[yyi++]);
		  yyf += 2;
		}
	      else
		{
		  yyp++;
		  yyf++;
		}
	    }
	}
      return yysize;
    }
}
#endif /* YYERROR_VERBOSE */


/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep, YYLTYPE *yylocationp)
#else
static void
yydestruct (yymsg, yytype, yyvaluep, yylocationp)
    const char *yymsg;
    int yytype;
    YYSTYPE *yyvaluep;
    YYLTYPE *yylocationp;
#endif
{
  YYUSE (yyvaluep);
  YYUSE (yylocationp);

  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  switch (yytype)
    {
      case 76: /* "expression" */
#line 141 "./grammar.y"
	{ delete ((yyvaluep->syntaxElement)); };
#line 1668 "./grammar.tab.c"
	break;
      case 77: /* "arrowAssign" */
#line 142 "./grammar.y"
	{ delete ((yyvaluep->syntaxElement)); };
#line 1673 "./grammar.tab.c"
	break;
      case 78: /* "tildeAssign" */
#line 142 "./grammar.y"
	{ delete ((yyvaluep->syntaxElement)); };
#line 1678 "./grammar.tab.c"
	break;
      case 79: /* "equationAssign" */
#line 142 "./grammar.y"
	{ delete ((yyvaluep->syntaxElement)); };
#line 1683 "./grammar.tab.c"
	break;
      case 80: /* "workspaceAssign" */
#line 142 "./grammar.y"
	{ delete ((yyvaluep->syntaxElement)); };
#line 1688 "./grammar.tab.c"
	break;
      case 81: /* "referenceAssign" */
#line 143 "./grammar.y"
	{ delete ((yyvaluep->syntaxElement)); };
#line 1693 "./grammar.tab.c"
	break;
      case 82: /* "additionAssign" */
#line 144 "./grammar.y"
	{ delete ((yyvaluep->syntaxElement)); };
#line 1698 "./grammar.tab.c"
	break;
      case 83: /* "subtractionAssign" */
#line 144 "./grammar.y"
	{ delete ((yyvaluep->syntaxElement)); };
#line 1703 "./grammar.tab.c"
	break;
      case 84: /* "multiplicationAssign" */
#line 144 "./grammar.y"
	{ delete ((yyvaluep->syntaxElement)); };
#line 1708 "./grammar.tab.c"
	break;
      case 85: /* "divisionAssign" */
#line 144 "./grammar.y"
	{ delete ((yyvaluep->syntaxElement)); };
#line 1713 "./grammar.tab.c"
	break;
      case 86: /* "variable" */
#line 140 "./grammar.y"
	{ delete ((yyvaluep->syntaxElement)); };
#line 1718 "./grammar.tab.c"
	break;
      case 87: /* "optElements" */
#line 136 "./grammar.y"
	{ for (std::list<SyntaxElement*>::iterator it=(yyvaluep->syntaxElementList)->begin(); it != (yyvaluep->syntaxElementList)->end(); it++) { SyntaxElement* theElement = *it; delete theElement; }; delete ((yyvaluep->syntaxElementList)); };
#line 1723 "./grammar.tab.c"
	break;
      case 88: /* "elementList" */
#line 136 "./grammar.y"
	{ for (std::list<SyntaxElement*>::iterator it=(yyvaluep->syntaxElementList)->begin(); it != (yyvaluep->syntaxElementList)->end(); it++) { SyntaxElement* theElement = *it; delete theElement; }; delete ((yyvaluep->syntaxElementList)); };
#line 1728 "./grammar.tab.c"
	break;
      case 89: /* "fxnCall" */
#line 140 "./grammar.y"
	{ delete ((yyvaluep->syntaxFunctionCall)); };
#line 1733 "./grammar.tab.c"
	break;
      case 90: /* "functionCall" */
#line 140 "./grammar.y"
	{ delete ((yyvaluep->syntaxFunctionCall)); };
#line 1738 "./grammar.tab.c"
	break;
      case 91: /* "optArguments" */
#line 137 "./grammar.y"
	{ for (std::list<SyntaxLabeledExpr*>::iterator it=(yyvaluep->argumentList)->begin(); it != (yyvaluep->argumentList)->end(); it++) { SyntaxLabeledExpr* theElement = *it; delete theElement; }; delete ((yyvaluep->argumentList)); };
#line 1743 "./grammar.tab.c"
	break;
      case 92: /* "argumentList" */
#line 137 "./grammar.y"
	{ for (std::list<SyntaxLabeledExpr*>::iterator it=(yyvaluep->argumentList)->begin(); it != (yyvaluep->argumentList)->end(); it++) { SyntaxLabeledExpr* theElement = *it; delete theElement; }; delete ((yyvaluep->argumentList)); };
#line 1748 "./grammar.tab.c"
	break;
      case 93: /* "argument" */
#line 140 "./grammar.y"
	{ delete ((yyvaluep->syntaxLabeledExpr)); };
#line 1753 "./grammar.tab.c"
	break;
      case 94: /* "functionDef" */
#line 146 "./grammar.y"
	{ delete ((yyvaluep->syntaxElement)); };
#line 1758 "./grammar.tab.c"
	break;
      case 95: /* "procedureDef" */
#line 146 "./grammar.y"
	{ delete ((yyvaluep->syntaxElement)); };
#line 1763 "./grammar.tab.c"
	break;
      case 96: /* "optFormals" */
#line 138 "./grammar.y"
	{ for (std::list<SyntaxFormal*>::iterator it=(yyvaluep->formalList)->begin(); it != (yyvaluep->formalList)->end(); it++) { SyntaxFormal* theElement = *it; delete theElement; }; delete ((yyvaluep->formalList)); };
#line 1768 "./grammar.tab.c"
	break;
      case 97: /* "formalList" */
#line 138 "./grammar.y"
	{ for (std::list<SyntaxFormal*>::iterator it=(yyvaluep->formalList)->begin(); it != (yyvaluep->formalList)->end(); it++) { SyntaxFormal* theElement = *it; delete theElement; }; delete ((yyvaluep->formalList)); };
#line 1773 "./grammar.tab.c"
	break;
      case 98: /* "formal" */
#line 140 "./grammar.y"
	{ delete ((yyvaluep->syntaxFormal)); };
#line 1778 "./grammar.tab.c"
	break;
      case 99: /* "typeSpec" */
#line 139 "./grammar.y"
	{ delete ((yyvaluep->string)); };
#line 1783 "./grammar.tab.c"
	break;
      case 100: /* "optDims" */
#line 139 "./grammar.y"
	{ delete ((yyvaluep->string)); };
#line 1788 "./grammar.tab.c"
	break;
      case 101: /* "dimList" */
#line 139 "./grammar.y"
	{ delete ((yyvaluep->string)); };
#line 1793 "./grammar.tab.c"
	break;
      case 102: /* "stmts" */
#line 136 "./grammar.y"
	{ for (std::list<SyntaxElement*>::iterator it=(yyvaluep->syntaxElementList)->begin(); it != (yyvaluep->syntaxElementList)->end(); it++) { SyntaxElement* theElement = *it; delete theElement; }; delete ((yyvaluep->syntaxElementList)); };
#line 1798 "./grammar.tab.c"
	break;
      case 103: /* "stmtList" */
#line 136 "./grammar.y"
	{ for (std::list<SyntaxElement*>::iterator it=(yyvaluep->syntaxElementList)->begin(); it != (yyvaluep->syntaxElementList)->end(); it++) { SyntaxElement* theElement = *it; delete theElement; }; delete ((yyvaluep->syntaxElementList)); };
#line 1803 "./grammar.tab.c"
	break;
      case 104: /* "statement" */
#line 141 "./grammar.y"
	{ delete ((yyvaluep->syntaxElement)); };
#line 1808 "./grammar.tab.c"
	break;
      case 105: /* "stmt_or_expr" */
#line 141 "./grammar.y"
	{ delete ((yyvaluep->syntaxElement)); };
#line 1813 "./grammar.tab.c"
	break;
      case 106: /* "declaration" */
#line 145 "./grammar.y"
	{ delete ((yyvaluep->syntaxElement)); };
#line 1818 "./grammar.tab.c"
	break;
      case 107: /* "memberDefs" */
#line 136 "./grammar.y"
	{ for (std::list<SyntaxElement*>::iterator it=(yyvaluep->syntaxElementList)->begin(); it != (yyvaluep->syntaxElementList)->end(); it++) { SyntaxElement* theElement = *it; delete theElement; }; delete ((yyvaluep->syntaxElementList)); };
#line 1823 "./grammar.tab.c"
	break;
      case 108: /* "memberDef" */
#line 145 "./grammar.y"
	{ delete ((yyvaluep->syntaxElement)); };
#line 1828 "./grammar.tab.c"
	break;
      case 109: /* "classDef" */
#line 145 "./grammar.y"
	{ delete ((yyvaluep->syntaxElement)); };
#line 1833 "./grammar.tab.c"
	break;
      case 110: /* "ifStatement" */
#line 147 "./grammar.y"
	{ delete ((yyvaluep->syntaxElement)); };
#line 1838 "./grammar.tab.c"
	break;
      case 111: /* "cond" */
#line 148 "./grammar.y"
	{ delete ((yyvaluep->syntaxElement)); };
#line 1843 "./grammar.tab.c"
	break;
      case 112: /* "forStatement" */
#line 147 "./grammar.y"
	{ delete ((yyvaluep->syntaxElement)); };
#line 1848 "./grammar.tab.c"
	break;
      case 113: /* "forCond" */
#line 148 "./grammar.y"
	{ delete ((yyvaluep->syntaxElement)); };
#line 1853 "./grammar.tab.c"
	break;
      case 114: /* "whileStatement" */
#line 147 "./grammar.y"
	{ delete ((yyvaluep->syntaxElement)); };
#line 1858 "./grammar.tab.c"
	break;
      case 115: /* "nextStatement" */
#line 149 "./grammar.y"
	{ delete ((yyvaluep->syntaxElement)); };
#line 1863 "./grammar.tab.c"
	break;
      case 116: /* "breakStatement" */
#line 149 "./grammar.y"
	{ delete ((yyvaluep->syntaxElement)); };
#line 1868 "./grammar.tab.c"
	break;
      case 117: /* "returnStatement" */
#line 148 "./grammar.y"
	{ delete ((yyvaluep->syntaxElement)); };
#line 1873 "./grammar.tab.c"
	break;
      case 118: /* "identifier" */
#line 139 "./grammar.y"
	{ delete ((yyvaluep->string)); };
#line 1878 "./grammar.tab.c"
	break;
      case 119: /* "vector" */
#line 137 "./grammar.y"
	{ for (std::list<SyntaxLabeledExpr*>::iterator it=(yyvaluep->argumentList)->begin(); it != (yyvaluep->argumentList)->end(); it++) { SyntaxLabeledExpr* theElement = *it; delete theElement; }; delete ((yyvaluep->argumentList)); };
#line 1883 "./grammar.tab.c"
	break;
      case 120: /* "vectorList" */
#line 137 "./grammar.y"
	{ for (std::list<SyntaxLabeledExpr*>::iterator it=(yyvaluep->argumentList)->begin(); it != (yyvaluep->argumentList)->end(); it++) { SyntaxLabeledExpr* theElement = *it; delete theElement; }; delete ((yyvaluep->argumentList)); };
#line 1888 "./grammar.tab.c"
	break;
      case 121: /* "constant" */
#line 140 "./grammar.y"
	{ delete ((yyvaluep->syntaxElement)); };
#line 1893 "./grammar.tab.c"
	break;

      default:
	break;
    }
}


/* Prevent warnings from -Wmissing-prototypes.  */

#ifdef YYPARSE_PARAM
#if defined __STDC__ || defined __cplusplus
int yyparse (void *YYPARSE_PARAM);
#else
int yyparse ();
#endif
#else /* ! YYPARSE_PARAM */
#if defined __STDC__ || defined __cplusplus
int yyparse (void);
#else
int yyparse ();
#endif
#endif /* ! YYPARSE_PARAM */



/* The look-ahead symbol.  */
int yychar;

/* The semantic value of the look-ahead symbol.  */
YYSTYPE yylval;

/* Number of syntax errors so far.  */
int yynerrs;
/* Location data for the look-ahead symbol.  */
YYLTYPE yylloc;



/*----------.
| yyparse.  |
`----------*/

#ifdef YYPARSE_PARAM
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void *YYPARSE_PARAM)
#else
int
yyparse (YYPARSE_PARAM)
    void *YYPARSE_PARAM;
#endif
#else /* ! YYPARSE_PARAM */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void)
#else
int
yyparse ()

#endif
#endif
{
  
  int yystate;
  int yyn;
  int yyresult;
  /* Number of tokens to shift before error messages enabled.  */
  int yyerrstatus;
  /* Look-ahead token as an internal (translated) token number.  */
  int yytoken = 0;
#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char yymsgbuf[128];
  char *yymsg = yymsgbuf;
  YYSIZE_T yymsg_alloc = sizeof yymsgbuf;
#endif

  /* Three stacks and their tools:
     `yyss': related to states,
     `yyvs': related to semantic values,
     `yyls': related to locations.

     Refer to the stacks thru separate pointers, to allow yyoverflow
     to reallocate them elsewhere.  */

  /* The state stack.  */
  yytype_int16 yyssa[YYINITDEPTH];
  yytype_int16 *yyss = yyssa;
  yytype_int16 *yyssp;

  /* The semantic value stack.  */
  YYSTYPE yyvsa[YYINITDEPTH];
  YYSTYPE *yyvs = yyvsa;
  YYSTYPE *yyvsp;

  /* The location stack.  */
  YYLTYPE yylsa[YYINITDEPTH];
  YYLTYPE *yyls = yylsa;
  YYLTYPE *yylsp;
  /* The locations where the error started and ended.  */
  YYLTYPE yyerror_range[2];

#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N), yylsp -= (N))

  YYSIZE_T yystacksize = YYINITDEPTH;

  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;
  YYLTYPE yyloc;

  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY;		/* Cause a token to be read.  */

  /* Initialize stack pointers.
     Waste one element of value and location stack
     so that they stay on the same level as the state stack.
     The wasted elements are never initialized.  */

  yyssp = yyss;
  yyvsp = yyvs;
  yylsp = yyls;
#if defined YYLTYPE_IS_TRIVIAL && YYLTYPE_IS_TRIVIAL
  /* Initialize the default location before parsing starts.  */
  yylloc.first_line   = yylloc.last_line   = 1;
  yylloc.first_column = yylloc.last_column = 0;
#endif

  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
	/* Give user a chance to reallocate the stack.  Use copies of
	   these so that the &'s don't force the real ones into
	   memory.  */
	YYSTYPE *yyvs1 = yyvs;
	yytype_int16 *yyss1 = yyss;
	YYLTYPE *yyls1 = yyls;

	/* Each stack pointer address is followed by the size of the
	   data in use in that stack, in bytes.  This used to be a
	   conditional around just the two extra args, but that might
	   be undefined if yyoverflow is a macro.  */
	yyoverflow (YY_("memory exhausted"),
		    &yyss1, yysize * sizeof (*yyssp),
		    &yyvs1, yysize * sizeof (*yyvsp),
		    &yyls1, yysize * sizeof (*yylsp),
		    &yystacksize);
	yyls = yyls1;
	yyss = yyss1;
	yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyexhaustedlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
	goto yyexhaustedlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
	yystacksize = YYMAXDEPTH;

      {
	yytype_int16 *yyss1 = yyss;
	union yyalloc *yyptr =
	  (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
	if (! yyptr)
	  goto yyexhaustedlab;
	YYSTACK_RELOCATE (yyss);
	YYSTACK_RELOCATE (yyvs);
	YYSTACK_RELOCATE (yyls);
#  undef YYSTACK_RELOCATE
	if (yyss1 != yyssa)
	  YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;
      yylsp = yyls + yysize - 1;

      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
		  (unsigned long int) yystacksize));

      if (yyss + yystacksize - 1 <= yyssp)
	YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

  /* Do appropriate processing given the current state.  Read a
     look-ahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to look-ahead token.  */
  yyn = yypact[yystate];
  if (yyn == YYPACT_NINF)
    goto yydefault;

  /* Not known => get a look-ahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid look-ahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = YYLEX;
    }

  if (yychar <= YYEOF)
    {
      yychar = yytoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
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
      if (yyn == 0 || yyn == YYTABLE_NINF)
	goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  if (yyn == YYFINAL)
    YYACCEPT;

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the look-ahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);

  /* Discard the shifted token unless it is eof.  */
  if (yychar != YYEOF)
    yychar = YYEMPTY;

  yystate = yyn;
  *++yyvsp = yylval;
  *++yylsp = yylloc;
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
| yyreduce -- Do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     `$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];

  /* Default location.  */
  YYLLOC_DEFAULT (yyloc, (yylsp - yylen), yylen);
  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
        case 2:
#line 235 "./grammar.y"
    {
#ifdef DEBUG_BISON_FLEX
                    printf("Bison encountered end_of_input; ignored\n");
#endif
                    return 0;
                }
    break;

  case 3:
#line 242 "./grammar.y"
    {
#ifdef DEBUG_BISON_FLEX
                    printf("Bison encountered newline; ignored\n");
#endif
                    return 0;
                }
    break;

  case 4:
#line 249 "./grammar.y"
    {
#ifdef DEBUG_BISON_FLEX
                    printf("Bison trying to execute statement or expression\n");
#endif
                    int rv = parser.execute((yyvsp[(1) - (2)].syntaxElement), *executionEnvironment);
                    delete (yyvsp[(1) - (2)].syntaxElement);
                    return rv;
                }
    break;

  case 5:
#line 258 "./grammar.y"
    {
#ifdef DEBUG_BISON_FLEX
                    printf("Bison trying to execute statement or expression\n");
#endif
                    int rv =  parser.execute((yyvsp[(1) - (2)].syntaxElement), *executionEnvironment);
                    delete (yyvsp[(1) - (2)].syntaxElement);
                    return rv;
                }
    break;

  case 6:
#line 267 "./grammar.y"
    {
#ifdef DEBUG_BISON_FLEX
                    printf("Bison trying to execute declaration\n");
#endif
                    int rv =  parser.execute((yyvsp[(1) - (2)].syntaxElement), *executionEnvironment);
                    delete (yyvsp[(1) - (2)].syntaxElement);
                    return rv;
                }
    break;

  case 7:
#line 276 "./grammar.y"
    {
#ifdef DEBUG_BISON_FLEX
                    printf("Bison trying to execute declaration\n");
#endif
                    int rv =  parser.execute((yyvsp[(1) - (2)].syntaxElement), *executionEnvironment);
                    delete (yyvsp[(1) - (2)].syntaxElement);
                    return rv;
                }
    break;

  case 8:
#line 285 "./grammar.y"
    {
#ifdef DEBUG_BISON_FLEX
                    printf("Bison trying to get help for symbol\n");
#endif
                    int rv =  parser.help(*(yyvsp[(2) - (3)].string));
                    delete (yyvsp[(2) - (3)].string);
                    return rv;
                }
    break;

  case 9:
#line 294 "./grammar.y"
    {
#ifdef DEBUG_BISON_FLEX
                    printf("Bison trying to get help for symbol\n");
#endif
                    int rv =  parser.help(*(yyvsp[(2) - (3)].string));
                    delete (yyvsp[(2) - (3)].string);
                    return rv;
                }
    break;

  case 10:
#line 303 "./grammar.y"
    {
#ifdef DEBUG_BISON_FLEX
                    printf("Bison trying to get help for symbol\n");
#endif
                    int rv =  parser.help(*(yyvsp[(2) - (5)].string), *(yyvsp[(4) - (5)].string));
                    delete (yyvsp[(2) - (5)].string);
                    delete (yyvsp[(4) - (5)].string);
                    return rv;
                }
    break;

  case 11:
#line 313 "./grammar.y"
    {
#ifdef DEBUG_BISON_FLEX
                printf("Bison trying to get help for symbol\n");
#endif
                int rv =  parser.help(*(yyvsp[(2) - (5)].string), *(yyvsp[(4) - (5)].string));
                delete (yyvsp[(2) - (5)].string);
                delete (yyvsp[(4) - (5)].string);
                return rv;
                }
    break;

  case 12:
#line 323 "./grammar.y"
    {
#ifdef DEBUG_BISON_FLEX
                    printf("Bison trying to get help for function call\n");
#endif
                    int rv =  parser.help((yyvsp[(2) - (3)].syntaxFunctionCall));
                    delete (yyvsp[(2) - (3)].syntaxFunctionCall);
                    return rv;
                }
    break;

  case 13:
#line 332 "./grammar.y"
    {
#ifdef DEBUG_BISON_FLEX
                    printf("Bison trying to get help for function call\n");
#endif
                    int rv =  parser.help((yyvsp[(2) - (3)].syntaxFunctionCall));
                    delete (yyvsp[(2) - (3)].syntaxFunctionCall);
                    return rv;
                }
    break;

  case 14:
#line 341 "./grammar.y"
    {
#ifdef DEBUG_BISON_FLEX
                    printf("Bison error when reading line %d position %d to line %d position %d\n", (yyloc).first_line, (yyloc).first_column, (yyloc).last_line, (yyloc).last_column);
#endif
                    YYABORT;
                }
    break;

  case 15:
#line 348 "./grammar.y"
    {
#ifdef DEBUG_BISON_FLEX
                    printf("Bison error when reading line %d position %d to line %d position %d\n", (yyloc).first_line, (yyloc).first_column, (yyloc).last_line, (yyloc).last_column);
#endif
                    YYABORT;
                }
    break;

  case 16:
#line 356 "./grammar.y"
    { (yyval.syntaxElement) = (yyvsp[(1) - (1)].syntaxElement); }
    break;

  case 17:
#line 358 "./grammar.y"
    { (yyval.syntaxElement) = new SyntaxFunctionCall("v", (yyvsp[(1) - (1)].argumentList)); }
    break;

  case 18:
#line 360 "./grammar.y"
    { (yyval.syntaxElement) = (yyvsp[(2) - (3)].syntaxElement); }
    break;

  case 19:
#line 362 "./grammar.y"
    { (yyval.syntaxElement) = new SyntaxUnaryExpr(SyntaxUnaryExpr::UMinus, (yyvsp[(2) - (2)].syntaxElement)); }
    break;

  case 20:
#line 363 "./grammar.y"
    { (yyval.syntaxElement) = new SyntaxUnaryExpr(SyntaxUnaryExpr::UPlus, (yyvsp[(2) - (2)].syntaxElement)); }
    break;

  case 21:
#line 364 "./grammar.y"
    { (yyval.syntaxElement) = new SyntaxUnaryExpr(SyntaxUnaryExpr::UNot, (yyvsp[(2) - (2)].syntaxElement)); }
    break;

  case 22:
#line 365 "./grammar.y"
    { (yyval.syntaxElement) = new SyntaxUnaryExpr(SyntaxUnaryExpr::UAnd, (yyvsp[(2) - (2)].syntaxElement)); }
    break;

  case 23:
#line 367 "./grammar.y"
    { (yyval.syntaxElement) = new SyntaxDecrement( (yyvsp[(2) - (2)].syntaxElement), false ); }
    break;

  case 24:
#line 368 "./grammar.y"
    { (yyval.syntaxElement) = new SyntaxDecrement( (yyvsp[(1) - (2)].syntaxElement), true ); }
    break;

  case 25:
#line 369 "./grammar.y"
    { (yyval.syntaxElement) = new SyntaxIncrement( (yyvsp[(2) - (2)].syntaxElement), false ); }
    break;

  case 26:
#line 370 "./grammar.y"
    { (yyval.syntaxElement) = new SyntaxIncrement( (yyvsp[(1) - (2)].syntaxElement), true ); }
    break;

  case 27:
#line 372 "./grammar.y"
    { (yyval.syntaxElement) = new SyntaxBinaryExpr(SyntaxBinaryExpr::Range, (yyvsp[(1) - (3)].syntaxElement), (yyvsp[(3) - (3)].syntaxElement)); }
    break;

  case 28:
#line 374 "./grammar.y"
    { (yyval.syntaxElement) = new SyntaxBinaryExpr(SyntaxBinaryExpr::Add, (yyvsp[(1) - (3)].syntaxElement), (yyvsp[(3) - (3)].syntaxElement)); }
    break;

  case 29:
#line 375 "./grammar.y"
    { (yyval.syntaxElement) = new SyntaxBinaryExpr(SyntaxBinaryExpr::Sub, (yyvsp[(1) - (3)].syntaxElement), (yyvsp[(3) - (3)].syntaxElement)); }
    break;

  case 30:
#line 376 "./grammar.y"
    { (yyval.syntaxElement) = new SyntaxBinaryExpr(SyntaxBinaryExpr::Mul, (yyvsp[(1) - (3)].syntaxElement), (yyvsp[(3) - (3)].syntaxElement)); }
    break;

  case 31:
#line 377 "./grammar.y"
    { (yyval.syntaxElement) = new SyntaxBinaryExpr(SyntaxBinaryExpr::Div, (yyvsp[(1) - (3)].syntaxElement), (yyvsp[(3) - (3)].syntaxElement)); }
    break;

  case 32:
#line 378 "./grammar.y"
    { (yyval.syntaxElement) = new SyntaxBinaryExpr(SyntaxBinaryExpr::Exp, (yyvsp[(1) - (3)].syntaxElement), (yyvsp[(3) - (3)].syntaxElement)); }
    break;

  case 33:
#line 379 "./grammar.y"
    { (yyval.syntaxElement) = new SyntaxBinaryExpr(SyntaxBinaryExpr::Mod, (yyvsp[(1) - (3)].syntaxElement), (yyvsp[(3) - (3)].syntaxElement)); }
    break;

  case 34:
#line 381 "./grammar.y"
    { (yyval.syntaxElement) = new SyntaxBinaryExpr(SyntaxBinaryExpr::Lt, (yyvsp[(1) - (3)].syntaxElement), (yyvsp[(3) - (3)].syntaxElement)); }
    break;

  case 35:
#line 382 "./grammar.y"
    { (yyval.syntaxElement) = new SyntaxBinaryExpr(SyntaxBinaryExpr::Le, (yyvsp[(1) - (3)].syntaxElement), (yyvsp[(3) - (3)].syntaxElement)); }
    break;

  case 36:
#line 383 "./grammar.y"
    { (yyval.syntaxElement) = new SyntaxBinaryExpr(SyntaxBinaryExpr::Eq, (yyvsp[(1) - (3)].syntaxElement), (yyvsp[(3) - (3)].syntaxElement)); }
    break;

  case 37:
#line 384 "./grammar.y"
    { (yyval.syntaxElement) = new SyntaxBinaryExpr(SyntaxBinaryExpr::Ne, (yyvsp[(1) - (3)].syntaxElement), (yyvsp[(3) - (3)].syntaxElement)); }
    break;

  case 38:
#line 385 "./grammar.y"
    { (yyval.syntaxElement) = new SyntaxBinaryExpr(SyntaxBinaryExpr::Ge, (yyvsp[(1) - (3)].syntaxElement), (yyvsp[(3) - (3)].syntaxElement)); }
    break;

  case 39:
#line 386 "./grammar.y"
    { (yyval.syntaxElement) = new SyntaxBinaryExpr(SyntaxBinaryExpr::Gt, (yyvsp[(1) - (3)].syntaxElement), (yyvsp[(3) - (3)].syntaxElement)); }
    break;

  case 40:
#line 388 "./grammar.y"
    { (yyval.syntaxElement) = new SyntaxBinaryExpr(SyntaxBinaryExpr::And, (yyvsp[(1) - (3)].syntaxElement), (yyvsp[(3) - (3)].syntaxElement)); }
    break;

  case 41:
#line 389 "./grammar.y"
    { (yyval.syntaxElement) = new SyntaxBinaryExpr(SyntaxBinaryExpr::Or, (yyvsp[(1) - (3)].syntaxElement), (yyvsp[(3) - (3)].syntaxElement)); }
    break;

  case 42:
#line 390 "./grammar.y"
    { (yyval.syntaxElement) = new SyntaxBinaryExpr(SyntaxBinaryExpr::And2, (yyvsp[(1) - (3)].syntaxElement), (yyvsp[(3) - (3)].syntaxElement)); }
    break;

  case 43:
#line 391 "./grammar.y"
    { (yyval.syntaxElement) = new SyntaxBinaryExpr(SyntaxBinaryExpr::Or2, (yyvsp[(1) - (3)].syntaxElement), (yyvsp[(3) - (3)].syntaxElement)); }
    break;

  case 44:
#line 393 "./grammar.y"
    { (yyval.syntaxElement) = (yyvsp[(1) - (1)].syntaxElement); }
    break;

  case 45:
#line 394 "./grammar.y"
    { (yyval.syntaxElement) = (yyvsp[(1) - (1)].syntaxElement); }
    break;

  case 46:
#line 395 "./grammar.y"
    { (yyval.syntaxElement) = (yyvsp[(1) - (1)].syntaxElement); }
    break;

  case 47:
#line 396 "./grammar.y"
    { (yyval.syntaxElement) = (yyvsp[(1) - (1)].syntaxElement); }
    break;

  case 48:
#line 397 "./grammar.y"
    { (yyval.syntaxElement) = (yyvsp[(1) - (1)].syntaxElement); }
    break;

  case 49:
#line 399 "./grammar.y"
    { (yyval.syntaxElement) = (yyvsp[(1) - (1)].syntaxElement); }
    break;

  case 50:
#line 400 "./grammar.y"
    { (yyval.syntaxElement) = (yyvsp[(1) - (1)].syntaxElement); }
    break;

  case 51:
#line 401 "./grammar.y"
    { (yyval.syntaxElement) = (yyvsp[(1) - (1)].syntaxElement); }
    break;

  case 52:
#line 402 "./grammar.y"
    { (yyval.syntaxElement) = (yyvsp[(1) - (1)].syntaxElement); }
    break;

  case 53:
#line 404 "./grammar.y"
    { (yyval.syntaxElement) = (yyvsp[(1) - (1)].syntaxFunctionCall); }
    break;

  case 54:
#line 406 "./grammar.y"
    { (yyval.syntaxElement) = (yyvsp[(1) - (1)].syntaxElement); }
    break;

  case 55:
#line 410 "./grammar.y"
    { 
#ifdef DEBUG_BISON_FLEX
                        printf("Parser inserting arrow assignment (ARROW_ASSIGN) in syntax tree\n");
#endif
                        (yyval.syntaxElement) = new SyntaxConstantAssignment((yyvsp[(1) - (3)].syntaxElement), (yyvsp[(3) - (3)].syntaxElement));
                    }
    break;

  case 56:
#line 419 "./grammar.y"
    {
#ifdef DEBUG_BISON_FLEX
                        printf("Parser inserting tilde assignment (TILDE_ASSIGN) in syntax tree\n");
#endif
                        (yyval.syntaxElement) = new SyntaxStochasticAssignment((yyvsp[(1) - (3)].syntaxElement), (yyvsp[(3) - (3)].syntaxElement));
                    }
    break;

  case 57:
#line 428 "./grammar.y"
    {
#ifdef DEBUG_BISON_FLEX
                        printf("Parser inserting equation assignment (EQUATION_ASSIGN) in syntax tree\n");
#endif
                        (yyval.syntaxElement) = new SyntaxDeterministicAssignment((yyvsp[(1) - (3)].syntaxElement), (yyvsp[(3) - (3)].syntaxElement)); 
                    }
    break;

  case 58:
#line 437 "./grammar.y"
    {
#ifdef DEBUG_BISON_FLEX
                        printf("Parser inserting workspace assignment (WORKSPACE_ASSIGN) in syntax tree\n");
#endif
                        (yyval.syntaxElement) = new SyntaxWorkspaceVariableAssignment((yyvsp[(1) - (3)].syntaxElement), (yyvsp[(3) - (3)].syntaxElement));
                    }
    break;

  case 59:
#line 446 "./grammar.y"
    {
#ifdef DEBUG_BISON_FLEX
                        printf("Parser inserting reference assignment (REFERENCE_ASSIGN) in syntax tree\n");
#endif
                        (yyval.syntaxElement) = new SyntaxReferenceAssignment((yyvsp[(1) - (3)].syntaxElement), (yyvsp[(3) - (3)].syntaxElement));
                    }
    break;

  case 60:
#line 455 "./grammar.y"
    {
#ifdef DEBUG_BISON_FLEX
                        printf("Parser inserting addition assignment (ADDITION_ASSIGN) in syntax tree\n");
#endif
                        (yyval.syntaxElement) = new SyntaxAdditionAssignment((yyvsp[(1) - (3)].syntaxElement), (yyvsp[(3) - (3)].syntaxElement)); 
                    }
    break;

  case 61:
#line 464 "./grammar.y"
    {
#ifdef DEBUG_BISON_FLEX
                            printf("Parser inserting subtraction assignment (SUBTRACTION_ASSIGN) in syntax tree\n");
#endif
                            (yyval.syntaxElement) = new SyntaxSubtractionAssignment((yyvsp[(1) - (3)].syntaxElement), (yyvsp[(3) - (3)].syntaxElement)); 
                        }
    break;

  case 62:
#line 473 "./grammar.y"
    {
#ifdef DEBUG_BISON_FLEX
                                printf("Parser inserting multiplication assignment (MULTIPLICATION_ASSIGN) in syntax tree\n");
#endif
                                (yyval.syntaxElement) = new SyntaxMultiplicationAssignment((yyvsp[(1) - (3)].syntaxElement), (yyvsp[(3) - (3)].syntaxElement)); 
                            }
    break;

  case 63:
#line 482 "./grammar.y"
    {
#ifdef DEBUG_BISON_FLEX
                        printf("Parser inserting division assignment (DIVISION_ASSIGN) in syntax tree\n");
#endif
                        (yyval.syntaxElement) = new SyntaxDivisionAssignment((yyvsp[(1) - (3)].syntaxElement), (yyvsp[(3) - (3)].syntaxElement)); 
                    }
    break;

  case 64:
#line 491 "./grammar.y"
    {
#ifdef DEBUG_BISON_FLEX
                    printf("Parser inserting variable (NAMED_VAR) in syntax tree\n");
#endif
                    (yyval.syntaxElement) = new SyntaxVariable(*(yyvsp[(1) - (2)].string));
                    for (std::list<SyntaxElement*>::iterator it=(yyvsp[(2) - (2)].syntaxElementList)->begin(); it!=(yyvsp[(2) - (2)].syntaxElementList)->end(); ++it)
                    {
                        (yyval.syntaxElement) = new SyntaxIndexOperation((yyval.syntaxElement),*it);
                    }
                    delete (yyvsp[(1) - (2)].string);
                    delete (yyvsp[(2) - (2)].syntaxElementList);
                }
    break;

  case 65:
#line 504 "./grammar.y"
    {
#ifdef DEBUG_BISON_FLEX
                    printf("Parser inserting variable (FUNCTION_VAR) in syntax tree\n");
#endif
                    (yyval.syntaxElement) = new SyntaxIndexOperation((yyvsp[(1) - (5)].syntaxFunctionCall),(yyvsp[(3) - (5)].syntaxElement));
                    for (std::list<SyntaxElement*>::iterator it=(yyvsp[(5) - (5)].syntaxElementList)->begin(); it!=(yyvsp[(5) - (5)].syntaxElementList)->end(); ++it)
                    {
                        (yyval.syntaxElement) = new SyntaxIndexOperation((yyval.syntaxElement),*it);
                    }
                    delete (yyvsp[(5) - (5)].syntaxElementList);
                }
    break;

  case 66:
#line 516 "./grammar.y"
    {
#ifdef DEBUG_BISON_FLEX
                    printf("Parser inserting variable (EXPRESSION_VAR) in syntax tree\n");
#endif
                    (yyval.syntaxElement) = new SyntaxIndexOperation((yyvsp[(2) - (7)].syntaxElement),(yyvsp[(5) - (7)].syntaxElement));
                    for (std::list<SyntaxElement*>::iterator it=(yyvsp[(7) - (7)].syntaxElementList)->begin(); it!=(yyvsp[(7) - (7)].syntaxElementList)->end(); ++it)
                    {
                        (yyval.syntaxElement) = new SyntaxIndexOperation((yyval.syntaxElement),*it);
                    }
                    delete (yyvsp[(7) - (7)].syntaxElementList);
                }
    break;

  case 67:
#line 528 "./grammar.y"
    {
#ifdef DEBUG_BISON_FLEX
                    printf("Parser inserting member variable (FUNCTION_VAR) in syntax tree\n");
#endif
                    (yyvsp[(3) - (7)].syntaxFunctionCall)->setBaseVariable((yyvsp[(1) - (7)].syntaxElement));
                    (yyval.syntaxElement) = new SyntaxIndexOperation((yyvsp[(3) - (7)].syntaxFunctionCall),(yyvsp[(5) - (7)].syntaxElement));
                    for (std::list<SyntaxElement*>::iterator it=(yyvsp[(7) - (7)].syntaxElementList)->begin(); it!=(yyvsp[(7) - (7)].syntaxElementList)->end(); ++it)
                    {
                        (yyval.syntaxElement) = new SyntaxIndexOperation((yyval.syntaxElement),*it);
                    }
                    delete (yyvsp[(7) - (7)].syntaxElementList);
                }
    break;

  case 68:
#line 542 "./grammar.y"
    { (yyval.syntaxElementList) = new std::list<SyntaxElement*>(); }
    break;

  case 69:
#line 543 "./grammar.y"
    { (yyval.syntaxElementList) = (yyvsp[(1) - (1)].syntaxElementList); }
    break;

  case 70:
#line 546 "./grammar.y"
    { (yyval.syntaxElementList) = new std::list<SyntaxElement*>(1, (yyvsp[(2) - (3)].syntaxElement)); }
    break;

  case 71:
#line 547 "./grammar.y"
    { (yyval.syntaxElementList) = new std::list<SyntaxElement*>(); }
    break;

  case 72:
#line 548 "./grammar.y"
    { (yyvsp[(1) - (4)].syntaxElementList)->push_back((yyvsp[(3) - (4)].syntaxElement)); (yyval.syntaxElementList) = (yyvsp[(1) - (4)].syntaxElementList); }
    break;

  case 73:
#line 549 "./grammar.y"
    { (yyvsp[(1) - (3)].syntaxElementList)->push_back( NULL ); (yyval.syntaxElementList) = (yyvsp[(1) - (3)].syntaxElementList); }
    break;

  case 74:
#line 553 "./grammar.y"
    {
                    (yyval.syntaxFunctionCall) = new SyntaxFunctionCall(*(yyvsp[(1) - (4)].string), (yyvsp[(3) - (4)].argumentList));
                    delete (yyvsp[(1) - (4)].string);
                }
    break;

  case 75:
#line 561 "./grammar.y"
    {
#ifdef DEBUG_BISON_FLEX
                        printf("Parser inserting function call in syntax tree\n");
#endif
                        (yyval.syntaxFunctionCall) = (yyvsp[(1) - (1)].syntaxFunctionCall);
                    }
    break;

  case 76:
#line 568 "./grammar.y"
    {
#ifdef DEBUG_BISON_FLEX
                        printf("Parser inserting member call in syntax tree\n");
#endif
                        (yyvsp[(3) - (3)].syntaxFunctionCall)->setBaseVariable((yyvsp[(1) - (3)].syntaxElement));
                        (yyval.syntaxFunctionCall) = (yyvsp[(3) - (3)].syntaxFunctionCall);
                    }
    break;

  case 77:
#line 576 "./grammar.y"
    {
#ifdef DEBUG_BISON_FLEX
                        printf("Parser inserting member call in syntax tree\n");
#endif
                        (yyvsp[(3) - (3)].syntaxFunctionCall)->setBaseVariable((yyvsp[(1) - (3)].syntaxFunctionCall));
                        (yyval.syntaxFunctionCall) = (yyvsp[(3) - (3)].syntaxFunctionCall);
                    }
    break;

  case 78:
#line 586 "./grammar.y"
    { (yyval.argumentList) = new std::list<SyntaxLabeledExpr*>(); }
    break;

  case 79:
#line 587 "./grammar.y"
    { (yyval.argumentList) = (yyvsp[(1) - (1)].argumentList); }
    break;

  case 80:
#line 590 "./grammar.y"
    { (yyval.argumentList) = new std::list<SyntaxLabeledExpr*>(1,(yyvsp[(1) - (1)].syntaxLabeledExpr)); }
    break;

  case 81:
#line 591 "./grammar.y"
    { (yyvsp[(1) - (3)].argumentList)->push_back((yyvsp[(3) - (3)].syntaxLabeledExpr)); (yyval.argumentList) = (yyvsp[(1) - (3)].argumentList); }
    break;

  case 82:
#line 595 "./grammar.y"
    {
#ifdef DEBUG_BISON_FLEX
                    printf("Parser inserting unlabeled argument in syntax tree\n");
#endif
                    (yyval.syntaxLabeledExpr) = new SyntaxLabeledExpr( "", (yyvsp[(1) - (1)].syntaxElement));
                }
    break;

  case 83:
#line 602 "./grammar.y"
    { 
#ifdef DEBUG_BISON_FLEX
                    printf("Parser inserting labeled argument in syntax tree\n");
#endif
                    (yyval.syntaxLabeledExpr) = new SyntaxLabeledExpr(*(yyvsp[(1) - (3)].string), (yyvsp[(3) - (3)].syntaxElement));
                    delete (yyvsp[(1) - (3)].string);
                }
    break;

  case 84:
#line 612 "./grammar.y"
    {
#ifdef DEBUG_BISON_FLEX
                    printf("Parser inserting function definition in syntax tree\n");
#endif
                    (yyval.syntaxElement) = new SyntaxFunctionDef("", *(yyvsp[(2) - (6)].string), (yyvsp[(4) - (6)].formalList), (yyvsp[(6) - (6)].syntaxElementList));
                    delete (yyvsp[(2) - (6)].string);
                }
    break;

  case 85:
#line 621 "./grammar.y"
    {
#ifdef DEBUG_BISON_FLEX
                    printf("Parser inserting typed function definition in syntax tree\n");
#endif
                    (yyvsp[(2) - (8)].string)->append(*(yyvsp[(3) - (8)].string));
                    (yyval.syntaxElement) = new SyntaxFunctionDef(*(yyvsp[(2) - (8)].string), *(yyvsp[(4) - (8)].string), (yyvsp[(6) - (8)].formalList), (yyvsp[(8) - (8)].syntaxElementList));
                    delete (yyvsp[(2) - (8)].string);
                    delete (yyvsp[(3) - (8)].string);
                    delete (yyvsp[(4) - (8)].string);
                }
    break;

  case 86:
#line 634 "./grammar.y"
    {
#ifdef DEBUG_BISON_FLEX
                    printf("Parser inserting procedure definition in syntax tree\n");
#endif
                    (yyval.syntaxElement) = new SyntaxFunctionDef("", *(yyvsp[(2) - (6)].string), (yyvsp[(4) - (6)].formalList), (yyvsp[(6) - (6)].syntaxElementList), true);
                    delete (yyvsp[(2) - (6)].string);
                }
    break;

  case 87:
#line 643 "./grammar.y"
    {
#ifdef DEBUG_BISON_FLEX
                    printf("Parser inserting typed procedure definition in syntax tree\n");
#endif
                    (yyvsp[(2) - (8)].string)->append(*(yyvsp[(3) - (8)].string));
                    (yyval.syntaxElement) = new SyntaxFunctionDef(*(yyvsp[(2) - (8)].string), *(yyvsp[(4) - (8)].string), (yyvsp[(6) - (8)].formalList), (yyvsp[(8) - (8)].syntaxElementList), true);
                    delete (yyvsp[(2) - (8)].string);
                    delete (yyvsp[(3) - (8)].string);
                    delete (yyvsp[(4) - (8)].string);
                }
    break;

  case 88:
#line 655 "./grammar.y"
    { (yyval.formalList) = new std::list<SyntaxFormal*>(); }
    break;

  case 89:
#line 656 "./grammar.y"
    { (yyval.formalList) = (yyvsp[(1) - (1)].formalList); }
    break;

  case 90:
#line 659 "./grammar.y"
    { (yyval.formalList) = new std::list<SyntaxFormal*>(1, (yyvsp[(1) - (1)].syntaxFormal)); }
    break;

  case 91:
#line 660 "./grammar.y"
    { (yyvsp[(1) - (3)].formalList)->push_back((yyvsp[(3) - (3)].syntaxFormal)); (yyval.formalList) = (yyvsp[(1) - (3)].formalList); }
    break;

  case 92:
#line 664 "./grammar.y"
    {
#ifdef DEBUG_BISON_FLEX
                    printf("Inserting labeled formal argument without default in syntax tree\n");
#endif
                    (yyval.syntaxFormal) = new SyntaxFormal(*(yyvsp[(1) - (1)].string), NULL );
                    delete (yyvsp[(1) - (1)].string);
                }
    break;

  case 93:
#line 672 "./grammar.y"
    {
#ifdef DEBUG_BISON_FLEX
                    printf("Inserting labeled formal argument with default in syntax tree\n");
#endif
                    (yyval.syntaxFormal) = new SyntaxFormal(*(yyvsp[(1) - (3)].string), (yyvsp[(3) - (3)].syntaxElement));
                    delete (yyvsp[(1) - (3)].string);
                }
    break;

  case 94:
#line 680 "./grammar.y"
    {
#ifdef DEBUG_BISON_FLEX
                    printf("Inserting typed labeled formal argument without default in syntax tree\n");
#endif
                    (yyval.syntaxFormal) = new SyntaxFormal(*(yyvsp[(1) - (2)].string), *(yyvsp[(2) - (2)].string), NULL);
                    delete (yyvsp[(1) - (2)].string);
                    delete (yyvsp[(2) - (2)].string);
                }
    break;

  case 95:
#line 689 "./grammar.y"
    {
#ifdef DEBUG_BISON_FLEX
                    printf("Inserting typed labeled formal argument with default in syntax tree\n");
#endif
                    (yyval.syntaxFormal) = new SyntaxFormal(*(yyvsp[(1) - (4)].string), *(yyvsp[(2) - (4)].string), (yyvsp[(4) - (4)].syntaxElement));
                    delete (yyvsp[(1) - (4)].string);
                    delete (yyvsp[(2) - (4)].string);
                }
    break;

  case 96:
#line 699 "./grammar.y"
    { (yyvsp[(1) - (2)].string)->append(*((yyvsp[(2) - (2)].string))); delete (yyvsp[(2) - (2)].string); (yyval.string) = (yyvsp[(1) - (2)].string); }
    break;

  case 97:
#line 700 "./grammar.y"
    { (yyvsp[(2) - (3)].string)->append(*((yyvsp[(3) - (3)].string))); (yyvsp[(2) - (3)].string)->insert(0, "const ");           delete (yyvsp[(3) - (3)].string); (yyval.string) = (yyvsp[(2) - (3)].string); }
    break;

  case 98:
#line 701 "./grammar.y"
    { (yyvsp[(2) - (3)].string)->append(*((yyvsp[(3) - (3)].string))); (yyvsp[(2) - (3)].string)->insert(0, "dynamic ");         delete (yyvsp[(3) - (3)].string); (yyval.string) = (yyvsp[(2) - (3)].string); }
    break;

  case 99:
#line 702 "./grammar.y"
    { (yyvsp[(2) - (3)].string)->append(*((yyvsp[(3) - (3)].string))); (yyvsp[(2) - (3)].string)->insert(0, "stochastic ");      delete (yyvsp[(3) - (3)].string); (yyval.string) = (yyvsp[(2) - (3)].string); }
    break;

  case 100:
#line 703 "./grammar.y"
    { (yyvsp[(2) - (3)].string)->append(*((yyvsp[(3) - (3)].string))); (yyvsp[(2) - (3)].string)->insert(0, "deterministic ");   delete (yyvsp[(3) - (3)].string); (yyval.string) = (yyvsp[(2) - (3)].string); }
    break;

  case 101:
#line 706 "./grammar.y"
    { (yyval.string) = new std::string(""); }
    break;

  case 102:
#line 707 "./grammar.y"
    { (yyval.string) = (yyvsp[(1) - (1)].string); }
    break;

  case 103:
#line 710 "./grammar.y"
    { (yyval.string) = new std::string("[]"); }
    break;

  case 104:
#line 711 "./grammar.y"
    { (yyvsp[(1) - (3)].string)->append("[]"); (yyval.string) = (yyvsp[(1) - (3)].string); }
    break;

  case 105:
#line 714 "./grammar.y"
    { (yyval.syntaxElementList) = (yyvsp[(2) - (3)].syntaxElementList); }
    break;

  case 106:
#line 716 "./grammar.y"
    {
                    std::list<SyntaxElement*>* stmts = new std::list<SyntaxElement*>();
                    stmts->push_back((yyvsp[(1) - (1)].syntaxElement));
                    (yyval.syntaxElementList) = stmts;
                }
    break;

  case 107:
#line 723 "./grammar.y"
    { (yyval.syntaxElementList) = new std::list<SyntaxElement*>(); }
    break;

  case 108:
#line 724 "./grammar.y"
    { (yyval.syntaxElementList) = new std::list<SyntaxElement*>(1, (yyvsp[(1) - (1)].syntaxElement)); }
    break;

  case 109:
#line 725 "./grammar.y"
    { (yyvsp[(1) - (3)].syntaxElementList)->push_back((yyvsp[(3) - (3)].syntaxElement)); (yyval.syntaxElementList) = (yyvsp[(1) - (3)].syntaxElementList); }
    break;

  case 110:
#line 726 "./grammar.y"
    { (yyval.syntaxElementList) = (yyvsp[(1) - (2)].syntaxElementList); }
    break;

  case 111:
#line 727 "./grammar.y"
    { (yyvsp[(1) - (3)].syntaxElementList)->push_back((yyvsp[(3) - (3)].syntaxElement)); (yyval.syntaxElementList) = (yyvsp[(1) - (3)].syntaxElementList); }
    break;

  case 112:
#line 728 "./grammar.y"
    { (yyval.syntaxElementList) = (yyvsp[(1) - (2)].syntaxElementList); }
    break;

  case 113:
#line 731 "./grammar.y"
    { (yyval.syntaxElement) = (yyvsp[(1) - (1)].syntaxElement); }
    break;

  case 114:
#line 732 "./grammar.y"
    { (yyval.syntaxElement) = (yyvsp[(1) - (1)].syntaxElement); }
    break;

  case 115:
#line 733 "./grammar.y"
    { (yyval.syntaxElement) = (yyvsp[(1) - (1)].syntaxElement); }
    break;

  case 116:
#line 734 "./grammar.y"
    { (yyval.syntaxElement) = (yyvsp[(1) - (1)].syntaxElement); }
    break;

  case 117:
#line 735 "./grammar.y"
    { (yyval.syntaxElement) = (yyvsp[(1) - (1)].syntaxElement); }
    break;

  case 118:
#line 736 "./grammar.y"
    { (yyval.syntaxElement) = (yyvsp[(1) - (1)].syntaxElement); }
    break;

  case 119:
#line 739 "./grammar.y"
    { (yyval.syntaxElement) = (yyvsp[(1) - (1)].syntaxElement); }
    break;

  case 120:
#line 740 "./grammar.y"
    { (yyval.syntaxElement) = (yyvsp[(1) - (1)].syntaxElement); }
    break;

  case 121:
#line 743 "./grammar.y"
    { (yyval.syntaxElement) = (yyvsp[(1) - (1)].syntaxElement); }
    break;

  case 122:
#line 744 "./grammar.y"
    { (yyval.syntaxElement) = (yyvsp[(1) - (1)].syntaxElement); }
    break;

  case 123:
#line 745 "./grammar.y"
    { (yyval.syntaxElement) = (yyvsp[(1) - (1)].syntaxElement); }
    break;

  case 124:
#line 747 "./grammar.y"
    {
#ifdef DEBUG_BISON_FLEX
                        printf("Parser inserting variable declaration in syntax tree\n");
#endif
                        (yyval.syntaxElement) = new SyntaxVariableDecl(*(yyvsp[(1) - (3)].string), (yyvsp[(2) - (3)].syntaxElementList), *(yyvsp[(3) - (3)].string));
                        delete (yyvsp[(1) - (3)].string);
                        delete (yyvsp[(3) - (3)].string);
                    }
    break;

  case 125:
#line 757 "./grammar.y"
    { (yyval.syntaxElementList) = new std::list<SyntaxElement*>(); }
    break;

  case 126:
#line 758 "./grammar.y"
    { (yyval.syntaxElementList) = new std::list<SyntaxElement*>(1, (yyvsp[(1) - (1)].syntaxElement)); }
    break;

  case 127:
#line 759 "./grammar.y"
    { (yyvsp[(1) - (3)].syntaxElementList)->push_back((yyvsp[(3) - (3)].syntaxElement)); (yyval.syntaxElementList) = (yyvsp[(1) - (3)].syntaxElementList); }
    break;

  case 128:
#line 760 "./grammar.y"
    { (yyval.syntaxElementList) = (yyvsp[(1) - (2)].syntaxElementList); }
    break;

  case 129:
#line 761 "./grammar.y"
    { (yyvsp[(1) - (3)].syntaxElementList)->push_back((yyvsp[(3) - (3)].syntaxElement)); (yyval.syntaxElementList) = (yyvsp[(1) - (3)].syntaxElementList); }
    break;

  case 130:
#line 762 "./grammar.y"
    { (yyval.syntaxElementList) = (yyvsp[(1) - (2)].syntaxElementList); }
    break;

  case 131:
#line 765 "./grammar.y"
    { (yyval.syntaxElement) = (yyvsp[(1) - (1)].syntaxFormal); }
    break;

  case 132:
#line 766 "./grammar.y"
    { (yyvsp[(2) - (2)].syntaxFormal)->setIsProtected(); (yyval.syntaxElement) = (yyvsp[(2) - (2)].syntaxFormal); }
    break;

  case 133:
#line 767 "./grammar.y"
    { (yyval.syntaxElement) = (yyvsp[(1) - (1)].syntaxElement); }
    break;

  case 134:
#line 768 "./grammar.y"
    { (yyval.syntaxElement) = (yyvsp[(1) - (1)].syntaxElement); }
    break;

  case 135:
#line 772 "./grammar.y"
    {
#ifdef DEBUG_BISON_FLEX
                printf("Parser inserting class definition (CLASS_DEF) in syntax tree\n");
#endif
                    (yyval.syntaxElement) = new SyntaxClassDef(*(yyvsp[(2) - (7)].string), *(yyvsp[(4) - (7)].string), (yyvsp[(6) - (7)].syntaxElementList));
                    delete (yyvsp[(2) - (7)].string);
                    delete (yyvsp[(4) - (7)].string);
                }
    break;

  case 136:
#line 782 "./grammar.y"
    { (yyval.syntaxElement) = new SyntaxStatement(SyntaxStatement::If, (yyvsp[(2) - (3)].syntaxElement), (yyvsp[(3) - (3)].syntaxElementList)); }
    break;

  case 137:
#line 783 "./grammar.y"
    { (yyval.syntaxElement) = new SyntaxStatement(SyntaxStatement::IfElse, (yyvsp[(2) - (5)].syntaxElement), (yyvsp[(3) - (5)].syntaxElementList), (yyvsp[(5) - (5)].syntaxElementList)); }
    break;

  case 138:
#line 785 "./grammar.y"
    { (yyval.syntaxElement) = (yyvsp[(2) - (3)].syntaxElement); }
    break;

  case 139:
#line 788 "./grammar.y"
    { (yyval.syntaxElement) = new SyntaxStatement(SyntaxStatement::For, (yyvsp[(2) - (3)].syntaxElement), (yyvsp[(3) - (3)].syntaxElementList)); }
    break;

  case 140:
#line 791 "./grammar.y"
    { (yyval.syntaxElement) = new SyntaxForLoop(*(yyvsp[(2) - (5)].string), (yyvsp[(4) - (5)].syntaxElement)); delete (yyvsp[(2) - (5)].string); }
    break;

  case 141:
#line 794 "./grammar.y"
    { (yyval.syntaxElement) = new SyntaxStatement(SyntaxStatement::While, (yyvsp[(2) - (3)].syntaxElement), (yyvsp[(3) - (3)].syntaxElementList)); }
    break;

  case 142:
#line 797 "./grammar.y"
    { (yyval.syntaxElement) = new SyntaxStatement(SyntaxStatement::Next); }
    break;

  case 143:
#line 800 "./grammar.y"
    { (yyval.syntaxElement) = new SyntaxStatement(SyntaxStatement::Break); }
    break;

  case 144:
#line 803 "./grammar.y"
    { (yyval.syntaxElement) = new SyntaxStatement(SyntaxStatement::Return); }
    break;

  case 145:
#line 804 "./grammar.y"
    { (yyval.syntaxElement) = new SyntaxStatement(SyntaxStatement::Return, (yyvsp[(2) - (2)].syntaxElement)); }
    break;

  case 146:
#line 807 "./grammar.y"
    { (yyval.string) = new std::string((yyvsp[(1) - (1)].c_string)); }
    break;

  case 147:
#line 811 "./grammar.y"
    { (yyval.argumentList) = (yyvsp[(2) - (3)].argumentList); }
    break;

  case 148:
#line 814 "./grammar.y"
    { (yyvsp[(1) - (3)].argumentList)->push_back(new SyntaxLabeledExpr( "", (yyvsp[(3) - (3)].syntaxElement)) ); (yyval.argumentList) = (yyvsp[(1) - (3)].argumentList); }
    break;

  case 149:
#line 816 "./grammar.y"
    {
                (yyval.argumentList) = new std::list<SyntaxLabeledExpr*>(1, new SyntaxLabeledExpr("", (yyvsp[(1) - (1)].syntaxElement)) );
                }
    break;

  case 150:
#line 822 "./grammar.y"
    {
                    (yyval.syntaxElement) = new SyntaxConstant(new RlBoolean(false) );
                }
    break;

  case 151:
#line 826 "./grammar.y"
    {
                    (yyval.syntaxElement) = new SyntaxConstant(new RlBoolean(true) );
                }
    break;

  case 152:
#line 830 "./grammar.y"
    {
                    (yyval.syntaxElement) = new SyntaxConstant( NULL );
                }
    break;

  case 153:
#line 834 "./grammar.y"
    {
                    (yyval.syntaxElement) = new SyntaxConstant( new RlString("\t") );
                }
    break;

  case 154:
#line 838 "./grammar.y"
    {
                    (yyval.syntaxElement) = new SyntaxConstant( new RealPos( RbConstants::Double::inf ) );
                }
    break;

  case 155:
#line 842 "./grammar.y"
    {
                    if ( (yyvsp[(1) - (1)].intValue) < 0 ) {
                        (yyval.syntaxElement) = new SyntaxConstant(new Integer((yyvsp[(1) - (1)].intValue)) );
                    }
                    else { 
                        (yyval.syntaxElement) = new SyntaxConstant(new Natural((yyvsp[(1) - (1)].intValue)) );
                    }
                }
    break;

  case 156:
#line 851 "./grammar.y"
    {
                    (yyval.syntaxElement) = new SyntaxConstant(new RlString((yyvsp[(1) - (1)].c_string)) );
                }
    break;

  case 157:
#line 855 "./grammar.y"
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
                    
                    if ((yyvsp[(1) - (1)].realValue) >= 0.0 && (yyvsp[(1) - (1)].realValue) <= 1.0) {
                        #ifdef DEBUG_BISON_FLEX
                        printf("Parser inserting RealPos constant in syntax tree \n");
                        #endif
                        (yyval.syntaxElement) = new SyntaxConstant(new Probability((yyvsp[(1) - (1)].realValue)) );
                    }
                    else if ((yyvsp[(1) - (1)].realValue) >= 0.0) {
#ifdef DEBUG_BISON_FLEX
                        printf("Parser inserting RealPos constant in syntax tree \n");
#endif
                        (yyval.syntaxElement) = new SyntaxConstant(new RealPos((yyvsp[(1) - (1)].realValue)) );
                    }
                    else {
#ifdef DEBUG_BISON_FLEX
                        printf("Parser inserting Real constant in syntax tree \n");
#endif
                        (yyval.syntaxElement) = new SyntaxConstant(new Real((yyvsp[(1) - (1)].realValue)) );
                    }
                }
    break;


/* Line 1267 of yacc.c.  */
#line 3359 "./grammar.tab.c"
      default: break;
    }
  YY_SYMBOL_PRINT ("-> $$ =", yyr1[yyn], &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;
  *++yylsp = yyloc;

  /* Now `shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*------------------------------------.
| yyerrlab -- here on detecting error |
`------------------------------------*/
yyerrlab:
  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if ! YYERROR_VERBOSE
      yyerror (YY_("syntax error"));
#else
      {
	YYSIZE_T yysize = yysyntax_error (0, yystate, yychar);
	if (yymsg_alloc < yysize && yymsg_alloc < YYSTACK_ALLOC_MAXIMUM)
	  {
	    YYSIZE_T yyalloc = 2 * yysize;
	    if (! (yysize <= yyalloc && yyalloc <= YYSTACK_ALLOC_MAXIMUM))
	      yyalloc = YYSTACK_ALLOC_MAXIMUM;
	    if (yymsg != yymsgbuf)
	      YYSTACK_FREE (yymsg);
	    yymsg = (char *) YYSTACK_ALLOC (yyalloc);
	    if (yymsg)
	      yymsg_alloc = yyalloc;
	    else
	      {
		yymsg = yymsgbuf;
		yymsg_alloc = sizeof yymsgbuf;
	      }
	  }

	if (0 < yysize && yysize <= yymsg_alloc)
	  {
	    (void) yysyntax_error (yymsg, yystate, yychar);
	    yyerror (yymsg);
	  }
	else
	  {
	    yyerror (YY_("syntax error"));
	    if (yysize != 0)
	      goto yyexhaustedlab;
	  }
      }
#endif
    }

  yyerror_range[0] = yylloc;

  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse look-ahead token after an
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

  /* Else will try to reuse look-ahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:

  /* Pacify compilers like GCC when the user code never invokes
     YYERROR and the label yyerrorlab therefore never appears in user
     code.  */
  if (/*CONSTCOND*/ 0)
     goto yyerrorlab;

  yyerror_range[0] = yylsp[1-yylen];
  /* Do not reclaim the symbols of the rule which action triggered
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
  yyerrstatus = 3;	/* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (yyn != YYPACT_NINF)
	{
	  yyn += YYTERROR;
	  if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
	    {
	      yyn = yytable[yyn];
	      if (0 < yyn)
		break;
	    }
	}

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
	YYABORT;

      yyerror_range[0] = *yylsp;
      yydestruct ("Error: popping",
		  yystos[yystate], yyvsp, yylsp);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  if (yyn == YYFINAL)
    YYACCEPT;

  *++yyvsp = yylval;

  yyerror_range[1] = yylloc;
  /* Using YYLLOC is tempting, but would change the location of
     the look-ahead.  YYLOC is available though.  */
  YYLLOC_DEFAULT (yyloc, (yyerror_range - 1), 2);
  *++yylsp = yyloc;

  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", yystos[yyn], yyvsp, yylsp);

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;

/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;

#ifndef yyoverflow
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
  if (yychar != YYEOF && yychar != YYEMPTY)
     yydestruct ("Cleanup: discarding lookahead",
		 yytoken, &yylval, &yylloc);
  /* Do not reclaim the symbols of the rule which action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
		  yystos[*yyssp], yyvsp, yylsp);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
#if YYERROR_VERBOSE
  if (yymsg != yymsgbuf)
    YYSTACK_FREE (yymsg);
#endif
  /* Make sure YYID is used.  */
  return YYID (yyresult);
}


#line 915 "./grammar.y"



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



