/* A Bison parser, made by GNU Bison 3.8.2.  */

/* Bison interface for Yacc-like parsers in C

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

/* DO NOT RELY ON FEATURES THAT ARE NOT DOCUMENTED in the manual,
   especially those whose name start with YY_ or yy_.  They are
   private implementation details that can be changed or removed.  */

#ifndef YY_YY_GRAMMAR_TAB_H_INCLUDED
# define YY_YY_GRAMMAR_TAB_H_INCLUDED
/* Debug traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif
#if YYDEBUG
extern int yydebug;
#endif

/* Token kinds.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
  enum yytokentype
  {
    YYEMPTY = -2,
    YYEOF = 0,                     /* "end of file"  */
    YYerror = 256,                 /* error  */
    YYUNDEF = 257,                 /* "invalid token"  */
    REAL = 258,                    /* REAL  */
    INT = 259,                     /* INT  */
    NAME = 260,                    /* NAME  */
    STRING = 261,                  /* STRING  */
    RBNULL = 262,                  /* RBNULL  */
    RBTAB = 263,                   /* RBTAB  */
    FALSE = 264,                   /* FALSE  */
    TRUE = 265,                    /* TRUE  */
    RBINF = 266,                   /* RBINF  */
    FUNCTION = 267,                /* FUNCTION  */
    PROCEDURE = 268,               /* PROCEDURE  */
    CLASS = 269,                   /* CLASS  */
    FOR = 270,                     /* FOR  */
    IN = 271,                      /* IN  */
    IF = 272,                      /* IF  */
    ELSE = 273,                    /* ELSE  */
    WHILE = 274,                   /* WHILE  */
    NEXT = 275,                    /* NEXT  */
    BREAK = 276,                   /* BREAK  */
    RETURN = 277,                  /* RETURN  */
    MOD_CONST = 278,               /* MOD_CONST  */
    MOD_DYNAMIC = 279,             /* MOD_DYNAMIC  */
    MOD_STOCHASTIC = 280,          /* MOD_STOCHASTIC  */
    MOD_DETERMINISTIC = 281,       /* MOD_DETERMINISTIC  */
    PROTECTED = 282,               /* PROTECTED  */
    ARROW_ASSIGN = 283,            /* ARROW_ASSIGN  */
    TILDE_ASSIGN = 284,            /* TILDE_ASSIGN  */
    EQUATION_ASSIGN = 285,         /* EQUATION_ASSIGN  */
    WORKSPACE_ASSIGN = 286,        /* WORKSPACE_ASSIGN  */
    REFERENCE_ASSIGN = 287,        /* REFERENCE_ASSIGN  */
    ADDITION_ASSIGN = 288,         /* ADDITION_ASSIGN  */
    SUBTRACTION_ASSIGN = 289,      /* SUBTRACTION_ASSIGN  */
    MULTIPLICATION_ASSIGN = 290,   /* MULTIPLICATION_ASSIGN  */
    DIVISION_ASSIGN = 291,         /* DIVISION_ASSIGN  */
    DECREMENT = 292,               /* DECREMENT  */
    INCREMENT = 293,               /* INCREMENT  */
    EQUAL = 294,                   /* EQUAL  */
    AND = 295,                     /* AND  */
    OR = 296,                      /* OR  */
    AND2 = 297,                    /* AND2  */
    OR2 = 298,                     /* OR2  */
    GT = 299,                      /* GT  */
    GE = 300,                      /* GE  */
    LT = 301,                      /* LT  */
    LE = 302,                      /* LE  */
    EQ = 303,                      /* EQ  */
    NE = 304,                      /* NE  */
    PIPE = 305,                    /* PIPE  */
    PIPE_PLACEHOLDER = 306,        /* PIPE_PLACEHOLDER  */
    END_OF_INPUT = 307,            /* END_OF_INPUT  */
    UNOT = 308,                    /* UNOT  */
    UMINUS = 309,                  /* UMINUS  */
    UPLUS = 310,                   /* UPLUS  */
    UAND = 311                     /* UAND  */
  };
  typedef enum yytokentype yytoken_kind_t;
#endif
/* Token kinds.  */
#define YYEMPTY -2
#define YYEOF 0
#define YYerror 256
#define YYUNDEF 257
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
#define PIPE 305
#define PIPE_PLACEHOLDER 306
#define END_OF_INPUT 307
#define UNOT 308
#define UMINUS 309
#define UPLUS 310
#define UAND 311

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
union YYSTYPE
{
#line 83 "./grammar.y"

    char*                                           c_string;
    std::string*                                    string;
    double                                          realValue;
    std::int64_t                                    longIntValue;
    bool                                            boolValue;
    RevLanguage::SyntaxElement*                     syntaxElement;
    RevLanguage::SyntaxVariable*                    syntaxVariable;
    RevLanguage::SyntaxFunctionCall*                syntaxFunctionCall;
    RevLanguage::SyntaxLabeledExpr*                 syntaxLabeledExpr;
    RevLanguage::SyntaxFormal*                      syntaxFormal;
    std::list<RevLanguage::SyntaxElement*>*         syntaxElementList;
    std::list<RevLanguage::SyntaxLabeledExpr*>*     argumentList;
    std::list<RevLanguage::SyntaxFormal*>*          formalList;

#line 195 "./grammar.tab.h"

};
typedef union YYSTYPE YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif

/* Location type.  */
#if ! defined YYLTYPE && ! defined YYLTYPE_IS_DECLARED
typedef struct YYLTYPE YYLTYPE;
struct YYLTYPE
{
  int first_line;
  int first_column;
  int last_line;
  int last_column;
};
# define YYLTYPE_IS_DECLARED 1
# define YYLTYPE_IS_TRIVIAL 1
#endif


extern YYSTYPE yylval;
extern YYLTYPE yylloc;

int yyparse (void);


#endif /* !YY_YY_GRAMMAR_TAB_H_INCLUDED  */
