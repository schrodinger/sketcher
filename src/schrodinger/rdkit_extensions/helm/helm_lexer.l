/*
 * This is a flex file which is used to generate a C++ lexer for parsing a HELM
 * string.
 */

%{

#include <string_view>
#include "schrodinger/rdkit_extensions/helm/token_scanner.h"

#undef  YY_DECL
#define YY_DECL int helm::TokenScanner::lex( helm::TokenParser::semantic_type * const lval, helm::TokenParser::location_type *loc )

using token = helm::TokenParser::token;

#define yyterminate() return( token::END )

#define YY_NO_UNISTD_H

/* update location on matching */
#define YY_USER_ACTION loc->step(); loc->columns(yyleng);

%}

%option debug
%option nodefault
%option yyclass="helm::TokenScanner"
%option noyywrap
%option c++
%option outfile="lex.yy.cpp"

%s POLYMER

positive_int                [1-9][0-9]*
peptide_polymer_id          PEPTIDE{positive_int}
nucleotide_polymer_id       RNA{positive_int}
chem_polymer_id             CHEM{positive_int}
blob_polymer_id             BLOB{positive_int}

%%
%{
            yylval = lval;
%}

    /* Helm polymers */
{peptide_polymer_id}|{nucleotide_polymer_id}|{chem_polymer_id}|{blob_polymer_id} {
    BEGIN POLYMER;
    yylval->build<std::string_view>(yytext);
    if (yytext[0] == 'B') { return token::BLOB_ID; }
    else if (yytext[0] == 'R') { return token::RNA_ID; }
    else if (yytext[0] == 'C') { return token::CHEM_ID; }
    else { return token::PEPTIDE_ID; }
}

<POLYMER>BEAD { return token::BEAD;}

<POLYMER>[a-zA-Z][a-zA-Z]* {
    yylval->build<std::string_view>(yytext);
    return yyleng > 1 ? token::MULTI_CHARACTER_MONOMER : token::SINGLE_CHARACTER_MONOMER;
}

<POLYMER>.|\n { return yytext[0]; }

    /* General tokens */
V2\.0 { return token::VERSION_TOKEN; }
\$ |
\n |
. { BEGIN 0; return yytext[0]; }

%%
