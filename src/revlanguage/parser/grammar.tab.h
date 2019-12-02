/* A Bison parser, made by GNU Bison 2.3.  */

/* Skeleton interface for Bison's Yacc-like parsers in C

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
/* Line 1529 of yacc.c.  */
#line 173 "./grammar.tab.h"
	YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif

extern YYSTYPE yylval;

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

extern YYLTYPE yylloc;
