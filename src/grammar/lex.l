%{
/**
 * @file
 * Lexicographical analysis for the Rev language.
 *
 * @brief Lexicographical analysis using flex
 *
 * @author Fredrik Ronquist and Sebastian Hoehna
 */

/* The following statements go into the resulting cpp code */

/* Define the call-back function for receiving input */
#define YY_INPUT rrinput

#include <cstring>

/* Hack for Windows */
#ifdef WIN32
#define strdup(a) _strdup(a)
#define fileno _fileno
#define YY_SKIP_YYWRAP
#define yywrap() 1
int isatty(int) { return 0; }
#else
extern "C" { 
   char *strdup(const char *s);
}
#define YY_SKIP_YYWRAP
#define yywrap() 1
#endif
    
    
#ifndef WIN32
#include <unistd.h> // for isatty()
#endif

/* Include objects needed to understand grammar.tab.h */
#include "RbOptions.h"     // To catch any parser debug flags
#include "Parser.h"
#include "RlString.h"
#include "SyntaxElement.h"
#include "SyntaxFormal.h"
#include "SyntaxFunctionCall.h"
#include "SyntaxLabeledExpr.h"

#include "grammar.tab.h"

#define YY_NEVER_INTERACTIVE 1
extern void rrinput(char *,size_t &,size_t);
#define YY_NO_UNISTD_H

size_t yycolumn;

#define YY_USER_ACTION yycolumn += yyleng;

%}

DIGIT    [0-9] 
EXP      [Ee][+-]?[0-9]+
ID       [_]*[a-zA-Z][a-zA-Z0-9_]*
STRING   \x22(\\.|[^\"])*\x22
COMMENT  \x23.*

%%
function                                {
                                            return FUNCTION;
                                        }
procedure                               {
                                            return PROCEDURE;
                                        }
class                                   {
                                            return CLASS;
                                        }
for                                     {
                                            return FOR;
                                        }
in                                      {
                                            return IN;
                                        }
if                                      {
                                            return IF;
                                        }
else                                    {
                                            return ELSE;
                                        }
while                                   {
                                            return WHILE;
                                        }
next                                    {
                                            return NEXT;
                                        }
break                                   {
                                            return BREAK;
                                        }
return                                  {
                                            return RETURN;
                                        }
null                                    {
                                            return RBNULL;
                                        }
tab                                     {
                                            return RBTAB;
                                        }
false                                   {
                                            return FALSE;
                                        }
true                                    {
                                            return TRUE;
                                        }
NULL                                    {
                                            return RBNULL;
                                        }
TAB                                     {
                                            return RBTAB;
                                        }
Inf                                     {
                                            return RBINF;
                                        }
FALSE                                   {
                                            return FALSE;
                                        }
TRUE                                    {
                                            return TRUE;
                                        }
const                                   {
                                            return MOD_CONST;
                                        }
dynamic                                 {
                                            return MOD_DYNAMIC;
                                        }
stochastic                              {
                                            return MOD_STOCHASTIC;
                                        }
deterministic                           {
                                            return MOD_DETERMINISTIC;
                                        }
protected                               {
                                            return PROTECTED;
                                        }
{DIGIT}+"."{DIGIT}*({EXP})?             {
#ifdef DEBUG_BISON_FLEX
                                            printf("Flex found real %s -- ", yytext);
#endif
                                            yylval.realValue = strtod(yytext,NULL); 
#ifdef DEBUG_BISON_FLEX
                                            printf("double value = %f\n",yylval.realValue);
#endif
                                            return REAL;
                                        }
{DIGIT}*"."{DIGIT}+({EXP})?             {
                                            yylval.realValue = strtod(yytext,NULL); 
                                            return REAL;
                                        }
{DIGIT}+{EXP}                           { 
                                            /* Treat ints with exponents as floating point values */
                                            yylval.realValue = strtod(yytext,NULL); 
                                            return REAL;
                                        }
{DIGIT}+                                {
                                            yylval.longIntValue = strtoll(yytext,NULL,10);
                                            if ((yylval.longIntValue == LLONG_MIN || yylval.longIntValue == LLONG_MAX) && errno == ERANGE) {
                                                yy_fatal_error("Integer value out of range.");
                                            }
                                            return INT;
                                        }
{ID}                                    {
                                            yylval.c_string = yytext;
                                            return NAME;
                                        }
{STRING}                                {
                                            /* Strip quotes from string */
                                            yytext[strlen(yytext)-1] = '\0';
                                            yylval.c_string = yytext + 1;
                                            return STRING;
                                        }
{COMMENT}                               {
                                            /* Ignore comment */
                                        }
"<-"                                    {
                                            return ARROW_ASSIGN;
                                        }
":="                                    {
                                            return EQUATION_ASSIGN;
                                        }
"~"                                     {
                                            return TILDE_ASSIGN;
                                        }
"<-&"                                   {
                                            return REFERENCE_ASSIGN;
                                        }
"+="                                    {
                                            return ADDITION_ASSIGN;
                                        }
"-="                                    {
                                            return SUBTRACTION_ASSIGN;
                                        }
"*="                                    {
                                            return MULTIPLICATION_ASSIGN;
                                        }
"/="                                    {
                                            return DIVISION_ASSIGN;
                                        }
"++"                                    {
                                            return INCREMENT;
                                        }
"--"                                    {
                                            return DECREMENT;
                                        }
"="                                     {
                                            return EQUAL;
                                        }
"&"                                     {
                                            return AND;
                                        }
"|"                                     {
                                            return OR;
                                        }
"&&"                                    {
                                            return AND;
                                        }
"||"                                    {
                                            return OR;
                                        }
">"                                     {
                                            return GT;
                                        }
">="                                    {
                                            return GE;
                                        }
"<"                                     {
                                            return LT;
                                        }
"<="                                    {
                                            return LE;
                                        }
"=="                                    {
                                            return EQ;
                                        }
"!="                                    {
                                            return NE;
                                        }
"|>"                                    {
                                            return PIPE;
                                        }
"_"                                     {
                                            return PIPE_PLACEHOLDER;
                                        }
[ \t]+                                  {
                                            /* Ignore white space */
                                        }
"\n"                                    {
                                            foundNewline = true;
                                            return yytext[0];
                                        }
";"                                     {
                                            return yytext[0];
                                        }
.                                       {
                                            /* any character not matching existing rules */
                                            /* pass it onto the parser and let the parser handle it */
                                            return yytext[0];
                                        }
%%
