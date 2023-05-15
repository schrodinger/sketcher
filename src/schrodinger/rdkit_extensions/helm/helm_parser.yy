%{
/*
 * This is a bison file which is used to generate a C++ parser for parsing a HELM
 * string.
 *
 * The grammar defines the specification for a valid HELMV2.0 string and uses
 * a helm::HelmParser instance to add polymer and monomer information to an
 * ROMol object.
 */
%}

%skeleton "lalr1.cc"
%require  "3.0"
%debug
%defines
%define api.namespace {helm}
%define api.parser.class {TokenParser}
%defines "helm_parser.tab.hh"
%output "helm_parser.tab.cpp"

%code requires{
   namespace helm {
      class HelmParser;
      class TokenScanner;
   }
}

%parse-param { TokenScanner  &scanner  }
%parse-param { HelmParser  &helm_parser}

%code{
   #include <string>
   #include <string_view>

   #include "schrodinger/rdkit_extensions/helm/helm_parser.h"

#undef yylex
#define yylex scanner.lex
}

%define api.value.type variant
%define parse.assert


%token BEAD VERSION_TOKEN
%token END 0
%token <std::string_view> BLOB_ID CHEM_ID PEPTIDE_ID RNA_ID SINGLE_CHARACTER_MONOMER MULTI_CHARACTER_MONOMER

%locations

%start helm


%%

helm: polymers '$' connections '$' polymer_groups '$' extended_annotations '$' VERSION_TOKEN;

polymers: polymer
       ;

polymer: PEPTIDE_ID '{' monomers '}' { helm_parser.add_polymer($1); }
       | RNA_ID '{' monomers '}'{ helm_parser.add_polymer($1); }
       | CHEM_ID '{' monomer '}'{ helm_parser.add_polymer($1); }
       | BLOB_ID '{' bead '}'{ helm_parser.add_polymer($1); }
       ;

monomers: monomer
       | monomers '.' monomer
       ;

monomer: SINGLE_CHARACTER_MONOMER { helm_parser.add_monomer($1); }
       | '[' MULTI_CHARACTER_MONOMER ']' { helm_parser.add_monomer($2);}
       ;

bead: BEAD { helm_parser.add_monomer("BEAD"); }
    ;

connections: /* currently unsupported */;
polymer_groups: /* currently unsupported */;
extended_annotations: /* currently unsupported */;


%%


void
helm::TokenParser::error(const location_type &l, const std::string& err_message )
{
    // save the pointer to the beginning of the current token
    helm_parser.saveErrorInformation(l.begin.column);
}
