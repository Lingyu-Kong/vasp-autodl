%{
#include <cstdio>
#include <iostream>
#include <cstdlib>
#include <string.h>
using namespace std;

#include "sites.hpp"
#include "functions.hpp"
#include "radial.hpp"
#include "basis.hpp"

char* strtmp;

// stuff from flex that bison needs to know about:
extern "C" int yylex();
extern "C" FILE *yyin;
extern int line_num;
// for support to read from string
typedef struct yy_buffer_state * YY_BUFFER_STATE;
extern YY_BUFFER_STATE yy_scan_bytes(const char *bytes, int len);
extern YY_BUFFER_STATE yy_scan_string(const char *string);
extern void yy_delete_buffer(YY_BUFFER_STATE buffer);

int yyparse(t_basis_vector* tbv);

void yyerror(t_basis_vector*, const char *s);
%}

// Bison fundamentally works by asking flex to get the next token, which it
// returns as an object of type "yystype".  But tokens could be of any
// arbitrary data type!  So we deal with that in Bison by defining a C union
// holding each of the types of tokens that Flex could return, and have Bison
// use that union instead of "int" for the definition of "yystype":
%union {
    int ival;
    float fval;
    char *sval;

    t_modifier* pradial;
    t_func_vector* pfuncs;
    t_site_vector* psitevector;

    t_basis_vector* pbasisvector;

    site_info* sinfo;
}

%parse-param {t_basis_vector* tbv}

%define parse.error verbose

// define the constant-string tokens:
%token PAW_PR_MOD
%token PAW_PS_MOD
%token RHYDRO_MOD
%token DELIM1
%token DELIM2
%token DASH
%token OPENKL
%token CLOSEKL
%token KOMMA

// define the "terminal symbol" token types I'm going to use (in CAPS
// by convention), and associate each with a field of the union:
%token <ival> INT
%token <fval> FLOAT
%token <sval> STRING

%type <pradial> modifiers
%type <pfuncs> functions
%type <sinfo> vector
%type <psitevector> sites
%type <pbasisvector> compound

%%
// the first rule defined is the highest-level rule

match:
        match command       { }
|                           { }
;

command:
    compound DELIM2 | DELIM2 
;

compound:
        sites DELIM1 functions DELIM1 modifiers {
      $$ = tbv; assemble_all(tbv,$1,$3,*$5); }
      | sites DELIM1 functions {
      $$ = tbv; assemble_all(tbv,$1,$3,t_modifier(1,-1,-1)); };
        
modifiers:
      PAW_PR_MOD             { $$ = new t_modifier(1,-1,-1); }
    | PAW_PR_MOD INT         { $$ = new t_modifier(1,-1,$2); }
    | PAW_PR_MOD INT INT     { $$ = new t_modifier(1,$3,$2); }
    | PAW_PS_MOD             { $$ = new t_modifier(2,-1,-1); }
    | PAW_PS_MOD INT         { $$ = new t_modifier(2,-1,$2); }
    | PAW_PS_MOD INT INT     { $$ = new t_modifier(2,$3,$2); }
    | RHYDRO_MOD             { $$ = new t_modifier(3,-1, 1); }
    | RHYDRO_MOD INT         { $$ = new t_modifier(3,-1,$2); }
    | RHYDRO_MOD INT FLOAT   { $$ = new t_modifier(3,-1,$2,$3); }
;

functions:
functions STRING     { $$ = $1;  interpret_function_string($$, $2);}
        | STRING     { $$ = new t_func_vector; interpret_function_string($$, $1);}
;

vector:
         OPENKL FLOAT KOMMA FLOAT KOMMA FLOAT CLOSEKL  {$$ = new site_info($2,$4,$6);}
       | OPENKL FLOAT KOMMA FLOAT KOMMA  INT  CLOSEKL  {$$ = new site_info($2,$4,$6);}
       | OPENKL FLOAT KOMMA  INT  KOMMA FLOAT CLOSEKL  {$$ = new site_info($2,$4,$6);}
       | OPENKL FLOAT KOMMA  INT  KOMMA  INT  CLOSEKL  {$$ = new site_info($2,$4,$6);}
       | OPENKL  INT  KOMMA FLOAT KOMMA FLOAT CLOSEKL  {$$ = new site_info($2,$4,$6);}
       | OPENKL  INT  KOMMA FLOAT KOMMA  INT  CLOSEKL  {$$ = new site_info($2,$4,$6);}
       | OPENKL  INT  KOMMA  INT  KOMMA FLOAT CLOSEKL  {$$ = new site_info($2,$4,$6);}
       | OPENKL  INT  KOMMA  INT  KOMMA  INT  CLOSEKL  {$$ = new site_info($2,$4,$6);}
       ;
        

sites:
 sites INT              { $1->push_back(site_info($2)); $$ = $1;}
 | sites INT DASH INT   { $$ = $1; push_back_fromto($$,$2,$4);}
 | sites vector         { $$ = $1; $$->push_back(*$2); delete $2;}
 | INT                  { $$ = new t_site_vector; $$->push_back(site_info($1)); }
 | INT DASH INT         { if ($1 > $3) {
                            yyerror(tbv,"first index larger than second in range expression");
                          }
                          $$ = new t_site_vector;
                          push_back_fromto($$,$1,$3);
                        }
 | vector               { $$ = new t_site_vector; $$->push_back(*$1); delete $1;}
;
%%

t_basis_vector tbvec;

extern "C" void parse_file_C(int& n_basis, int& strlen, char* filename){

    tbvec.clear();

    std::string s;

    strtmp = new char[strlen+1];

    for (int i=0;i<strlen;i++)
    strtmp[i] = filename[i];
    strtmp[strlen] = '\0';

    // open a file handle to a particular file:
    FILE *myfile = fopen(strtmp, "r");

    // make sure it's valid:
    if (!myfile) {
    cout << "I can't open input file " <<  strtmp << endl;
    delete[] strtmp;
    exit(0);
    }

//  std::cerr << "Opened file \"" << strtmp << "\".\n"; 

//  delete[] strtmp;

    // set lex to read from it instead of defaulting to STDIN:
    yyin = myfile;

    // parse through the input until there is no more:

    do {
    yyparse(&tbvec);
    } while (!feof(yyin));

//    output_basisfunctions(&tbvec);

    n_basis = tbvec.size();

    delete[] strtmp;

    fclose(myfile);
}

// refer to https://stackoverflow.com/questions/780676/string-input-to-flex-lexer (retrieved 8/05/2020)
extern "C" void parse_string_C(int& n_basis, int& strlen, char* string){

    //yydebug = 1;
    tbvec.clear();
    YY_BUFFER_STATE buffer = yy_scan_bytes(string,strlen);
    yyparse(&tbvec);
    yy_delete_buffer(buffer);
    n_basis = tbvec.size();
}

extern "C" void fill_basis_info_C(t_basis* result, const int& i){
    (*result) = tbvec[i-1];
}

extern "C" void free_parser_C(){
    t_basis_vector().swap(tbvec);
}

void yyerror(t_basis_vector*, const char *s) {
    cout << "Parse error on line " << line_num << ": " << s << endl;
    // might as well halt now:
    exit(-1);
}
