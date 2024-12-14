%{
#include <iostream>
#include <cstring>
#include "sites.hpp"
#include "functions.hpp"
#include "radial.hpp"
#include "basis.hpp"
#include "locproj.tab.h"
using namespace std;

#define YY_DECL extern "C" int yylex()
int line_num = 1;
%}
%option noyywrap
%x COMMENT
%%
<INITIAL>{
   ([-+]+)?[0-9]+\.([0-9]+)? { yylval.fval = atof(yytext); return FLOAT; };
   ([-+]+)?\.[0-9]+          { yylval.fval = atof(yytext); return FLOAT; };
   ([-+]+)?[0-9]+            { yylval.ival = atoi(yytext); return INT; };
   (?i:PR)                   { return PAW_PR_MOD; };
   (?i:PS)                   { return PAW_PS_MOD; };
   (?i:HY)                   { return RHYDRO_MOD; };
   [a-zA-Z]                  { yylval.sval = strdup(yytext); return STRING; };
   [a-zA-Z][a-zA-Z0-9\-)(]+  { yylval.sval = strdup(yytext); return STRING; };
   ":"            { return DELIM1; };
   "-"            { return DASH; };
   "("            { return OPENKL; };
   ")"            { return CLOSEKL; };
   ","            { return KOMMA; };
   ";"            { return DELIM2; };
   "\""           { return DELIM2; };
   \n             { ++line_num; return DELIM2; };
}
<COMMENT>{
      \n             { ++line_num; BEGIN(INITIAL); };
}
<*>{
      \!             { BEGIN(COMMENT); };
      \#             { BEGIN(COMMENT); };
      [ \t]           ;
      .               ;
}
%%
