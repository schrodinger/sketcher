/*
 * This is a flex file which is used to generate a C++ lexer for parsing a HELM
 * string.
 *
 * NOTE:
 *      To generate the C++ lexer (lex.yy.cpp), run the following command:
 *
 *          flex helm_lexer.l
 */

/* We must include stdint.h before the generated file tries to define its own
 * integral limits.  Otherwise, MSVC will complain about the macros being
 * redefined when it gets transitively included later
 */
%top{
#include <stdint.h>
}

%{

#include <string>
#include <string_view>
#include "schrodinger/rdkit_extensions/helm/helm_parser.h"
#include "schrodinger/rdkit_extensions/helm/token_scanner.h"

#undef  YY_DECL
#define YY_DECL int helm::TokenScanner::lex( helm::TokenParser::semantic_type * const lval, helm::TokenParser::location_type *loc )

using token = helm::TokenParser::token;

#define yyterminate() return( token::END )

#define YY_NO_UNISTD_H

#define YY_USER_ACTION \
        loc->step(); loc->columns(yyleng); \
        yylval->build<std::string_view>(ref_string_view.substr(loc->begin.column - 1, yyleng));

%}

%option debug
%option nodefault
%option yyclass="helm::TokenScanner"
%option noline
%option noyywrap
%option c++
%option outfile="lex.yy.cpp"
%option full

/* We're defining different token recognition states because token identities differ
 * based on what type of section we're parsing. This is mainly to allow the
 * following:
 *
 *     * Capture unknown monomers since the definition is X for peptides and
 *       N for nucleotides
 *     * Capture unknown sequences for BLOB polymers since they can be anything
 *       including whitespace characters.
 *
 */
%x PEPTIDE NUCLEOTIDE BLOB CHEM CONNECTIONS POLYMER_GROUPS

float_token                 [0-9]*\.[0-9]+
positive_int                [1-9][0-9]*
polymer_id                  (PEPTIDE|RNA|CHEM|BLOB){positive_int}

multicharacter_monomer      (\[[^\[\]\n]+\]|[^\[\]\n]+)+
inline_annotation           \"[^\"]+\"
num_repetitions             {positive_int}|{positive_int}\-{positive_int}
ratio                       ({float_token}|{positive_int})


%%
%{
            yylval = lval;
%}

    /* Helm polymers */
<INITIAL>{polymer_id} {
    if (yytext[0] == 'B') { BEGIN BLOB; return token::BLOB_ID; }
    else if (yytext[0] == 'R') { BEGIN NUCLEOTIDE; return token::RNA_ID; }
    else if (yytext[0] == 'C') { BEGIN CHEM; return token::CHEM_ID; }
    else { BEGIN PEPTIDE; return token::PEPTIDE_ID; }
}
<BLOB>[a-zA-Z ]+ { return token::UNKNOWN_SEQUENCE; }
    /* Special monomer recognition for peptides and nucleotide */
<PEPTIDE>X |
<NUCLEOTIDE>N { return token::UNKNOWN_MONOMER; }

<PEPTIDE,NUCLEOTIDE,CHEM>[a-zA-Z] { return token::SINGLE_CHARACTER_MONOMER; }
<PEPTIDE,NUCLEOTIDE,CHEM>_ { return token::MISSING_MONOMER; }
<PEPTIDE,NUCLEOTIDE,CHEM>\?|{ratio} { return token::MONOMER_RATIO; }
<PEPTIDE,NUCLEOTIDE,CHEM>\* { return token::MONOMER_WILDCARD; }
<PEPTIDE,NUCLEOTIDE,CHEM>\[{multicharacter_monomer}\] {
    return helm::is_smiles_monomer(yytext) ? token::INLINE_SMILES_MONOMER
                                           : token::MULTI_CHARACTER_MONOMER;
}
<PEPTIDE,NUCLEOTIDE,BLOB,CHEM>\| { BEGIN INITIAL; return yytext[0]; }
<PEPTIDE,NUCLEOTIDE,BLOB,CHEM>{inline_annotation} { return token::ANNOTATION; }
<PEPTIDE,NUCLEOTIDE>\'{num_repetitions}\' { return token::REPETITIONS; }
<PEPTIDE,NUCLEOTIDE,BLOB,CHEM>\$ { BEGIN CONNECTIONS; return yytext[0]; }


<CONNECTIONS>{polymer_id} { return token::POLYMER_ID; }
<CONNECTIONS>pair { return token::HYDROGEN_PAIRING; }
<CONNECTIONS>\? { return token::UNDEFINED_RESIDUE_NUMBER_OR_RGROUP; }
<CONNECTIONS>R{positive_int} { return token::RGROUP; }
<CONNECTIONS>{positive_int}|[a-zA-Z]+ { return token::CONNECTION_RESIDUE; }
<CONNECTIONS>{inline_annotation} { return token::ANNOTATION; }
<CONNECTIONS>G{positive_int} { return token::POLYMER_GROUP_ID; }
<CONNECTIONS>\$ { BEGIN POLYMER_GROUPS; return yytext[0]; }


<POLYMER_GROUPS>G{positive_int} { return token::POLYMER_GROUP_ID; }
<POLYMER_GROUPS>{polymer_id} { return token::POLYMER_ID;
}
<POLYMER_GROUPS>({ratio}|{ratio}\-{ratio}) { return token::POLYMER_GROUP_RATIO; }
<POLYMER_GROUPS>\$ { return yytext[0]; }

    /* General tokens */
<*>.|\n { return yytext[0]; }

%%
