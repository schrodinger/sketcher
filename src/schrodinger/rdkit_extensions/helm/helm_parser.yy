%{
/*
 * This is a bison file which is used to generate a C++ parser for parsing a HELM
 * string.
 *
 * The grammar defines the specification for a valid HELMV2.0 string and uses
 * a helm::HelmParser instance to retrieve the parsed helm information.
 *
 *
 * NOTE:
 *      To generate the C++ parser, run the following command:
 *
 *          bison helm_parser.yy
 *
 *      In addition to the helm_parser.tab.hh and helm_parser.tab.cpp files, we
 *      should also generate the following files:
 *          * location.hh
 *          * position.hh
 *          * stack.hh
 */
%}

%skeleton "lalr1.cc"
%require  "3.8.2"
%debug
%defines
%define api.namespace {helm}
%define api.parser.class {TokenParser}
%define api.value.automove
%defines "helm_parser.tab.hh"
%output "helm_parser.tab.cpp"
%no-lines

%code requires{
   #include <string>
   #include <string_view>

   namespace helm {
      class HelmParser;
      class TokenScanner;
   };
}

%parse-param { TokenScanner  &scanner  }
%parse-param { HelmParser  &helm_parser}

%code{
   #include <algorithm>
   #include <boost/algorithm/string.hpp>
   #include <string>
   #include <string_view>

   #include "schrodinger/rdkit_extensions/helm/helm_parser.h"

#undef yylex
#define yylex scanner.lex

    static const std::string branch_monomer_group_err_message{
      "Only one branch monomer is allowed. If you intended to create a branched"
      " group, create a second polymer with the monomer sequence and define "
      "a custom connection with the desired linkage. E.g.: PEPTIDE1{A(C.G.K)P} "
      "could be either "
      "PEPTIDE1{A(C)P}|PEPTIDE2{G.K}$PEPTIDE1,PEPTIDE2,2:R2-1:R1$$$$V2.0 or "
      "PEPTIDE1{A.P}|PEPTIDE2{C.G.K}$PEPTIDE1,PEPTIDE2,1:R3-1:R1$$$$V2.0"};
}

%define api.value.type variant
%locations


/* token to indicate end of parsing */
%token END 0

/* Polymer id tokens */
%token <std::string_view> BLOB_ID CHEM_ID PEPTIDE_ID RNA_ID POLYMER_ID
/* Monomer id tokens */
%token <std::string_view> SINGLE_CHARACTER_MONOMER MULTI_CHARACTER_MONOMER
%token <std::string_view> UNKNOWN_MONOMER MISSING_MONOMER MONOMER_WILDCARD
%token <std::string_view> MONOMER_RATIO INLINE_SMILES_MONOMER ANNOTATION REPETITIONS
/* top-level tokens */
%token <std::string_view> UNKNOWN_SEQUENCE

/* Connection tokens */
%token <std::string_view> HYDROGEN_PAIRING RGROUP UNDEFINED_RESIDUE_NUMBER_OR_RGROUP
%token <std::string_view> CONNECTION_RESIDUE

/* Polymer group tokens */
%token <std::string_view> POLYMER_GROUP_RATIO POLYMER_GROUP_ID

%type <std::string_view> annotation

%type <std::string_view> smiles_monomer monomer_id
%type <std::string_view> monomer_list_item monomer_and_list monomer_or_list monomer_list
%type <size_t> monomer_sequence branch_monomer_group
%type <std::string_view> polymer blob repetitions

%type <std::string_view> connection_polymer attachment_point residue_and_list
%type <std::string_view> residue_or_list connection_monomer
%type <std::string_view> polymer_group_item polymer_or_list polymer_and_list

// defining these types to check bad syntax
%type <size_t> connections polymer_groups

%start helm

%%


/**
  The overall structure of a HELMV2 string

  POLYMERS$CONNECTIONS$POLYMER_GROUPS$EXTEND_ANNOTATIONS$VERSION

  NOTE: The EXTENDED_ANNOTATIONS and VERSION sections are parsed outside of Bison
  because unusually large JSON inputs can significantly impact the parser
  performance.
*/
helm: polymers '$' connections '$' polymer_groups '$' extended_annotations_and_version;

/**
 * Defining a recursive rule for sections. The base case should for most is the
 * empty case since all but the polymers section is optional.
 *
 * NOTE: Condensing the non-empty and empty rules will allow parsing a
 * valid section that starts with a '|' token, e.g., |BLOB1,BLOB2,*:*-*:*.
 * We should handle this
 */
polymers: polymer_unit
        | polymers '|' polymer_unit
        ;
connections: /* connections are optional */ { $$ = 0; }
           | connection { $$ = 1; }
           | connections '|' connection {
                if ($1 == 0) { TokenParser::error(@1, "syntax error"); YYABORT; }
                ++$$;
           };
polymer_groups:  /* polymer groups are optional */ { $$ = 0; }
              | polymer_group { $$ = 1; }
              | polymer_groups '|' polymer_group {
                    if ($1 == 0) { TokenParser::error(@1, "syntax error"); YYABORT; }
                    ++$$;
              };
extended_annotations_and_version:  /* This is parsed outside of Bison */;


/* A structure of a polymer is as follows:
 *     * PEPTIDE and RNA polymers can consist of multiple monomer groups
 *       separated by a dot.
 *     * CHEM polymers can only have a single monomer. This monomer must conform
 *       to the monomer specifications for the PEPTIDE and RNA monomers.
 *     * BLOB polymers can only have a single monomer. This monomer can be
 *       anything. Examples are "BEAD" and "Gold particle".
 *     * Polymers can have an inline annotation
 *
 *
 * For PEPTIDE and RNA polymers, we define a monomer group as a substructure
 * that can for a backbone linkage to another subsstructure within the same polymer.
 * These are:
 *     * A single monomer, e.g., X, [dX], (X.X.X.X), *, _
 *     * A branch monomer group, e.g., X(X), X(X)X
 *     * A repeated monomer sequence, e.g., X'N', (X.X.X)'N'
 *
 */
polymer_unit: polymer { helm_parser.add_polymer($1, {}); }
            | polymer annotation { helm_parser.add_polymer($1, $2); }
            ;
polymer: PEPTIDE_ID '{' monomers '}' { $$ = $1; }
       | RNA_ID '{' monomers '}'{ $$ = $1; }
       | CHEM_ID '{' monomer_unit '}' { $$ = $1; }
       | BLOB_ID '{' blob '}' { $$ = $1; }
       ;
monomers: monomer_group
        | monomers '.' monomer_group
        ;
monomer_group: monomer_unit
             | repeated_monomers
             | branch_monomer_group { helm_parser.mark_branch_monomer($1); }
             /* Verbosely handling this since it's a likely case */
             | monomer_unit '(' monomer_sequence ')' {
                   helm_parser.saveError(@2.begin.column,
                                         branch_monomer_group_err_message);
                   YYABORT;
             }
             | monomer_unit '(' monomer_sequence ')' monomer_unit {
                   helm_parser.saveError(@2.begin.column,
                                         branch_monomer_group_err_message);
                   YYABORT;
             }
             ;
repeated_monomers: monomer_unit repetitions {
                    helm_parser.mark_last_n_monomers_as_repeated(1, $2, {});
                 }
                 | monomer_unit repetitions annotation {
                    helm_parser.mark_last_n_monomers_as_repeated(1, $2, $3);
                 }
                 | '(' monomer_sequence ')' repetitions {
                    helm_parser.mark_last_n_monomers_as_repeated($2, $4, {});
                 }
                 | '(' monomer_sequence ')' repetitions annotation {
                    helm_parser.mark_last_n_monomers_as_repeated($2, $4, $5);
                 }
                 ;
repetitions: REPETITIONS { $$ = $1; $$.remove_prefix(1); $$.remove_suffix(1); }
monomer_sequence: monomer_unit '.' monomer_unit { $$ = 2; }
                | branch_monomer_group { $$ = $1; helm_parser.mark_branch_monomer($$); }
                | monomer_sequence '.' monomer_unit { $$ = $1 + 1; }
                | monomer_sequence '.' branch_monomer_group {
                    helm_parser.mark_branch_monomer($3);
                    $$ = $1 + $3;
                }
                ;
branch_monomer_group: monomer_unit '(' monomer_unit ')' { $$ = 2; }
                    | monomer_unit '(' monomer_unit ')' monomer_unit { $$ = 3; };
monomer_unit: monomer_id { helm_parser.add_monomer_with_id($1, {}); };
       | monomer_id annotation { helm_parser.add_monomer_with_id($1,  $2); }
       | smiles_monomer { helm_parser.add_smiles_monomer($1,  {}); }
       | smiles_monomer annotation { helm_parser.add_smiles_monomer($1, $2); }
       | '(' monomer_list ')'  { helm_parser.add_monomer_list($2, {}); }
       | '(' monomer_list ')' annotation { helm_parser.add_monomer_list($2, $4); }
       ;
monomer_list: monomer_and_list | monomer_or_list;
monomer_and_list: monomer_list_item '+' monomer_list_item {
                    $$ = {$1.data(), $1.size() + 1 + $3.size()};
               }
               | monomer_and_list '+' monomer_list_item {
                    $$ = {$1.data(), $1.size() + 1 + $3.size()};
               }
               ;
monomer_or_list: monomer_list_item ',' monomer_list_item {
                    $$ = {$1.data(), $1.size() + 1 + $3.size()};
               }
               | monomer_or_list ',' monomer_list_item {
                    $$ = {$1.data(), $1.size() + 1 + $3.size()};
               }
               ;

monomer_list_item: monomer_id { $$ = $1; }
                 | monomer_id ':' MONOMER_RATIO { $$ = {$1.data(), $1.size() + $3.size() + 1};}
                 ;

smiles_monomer: INLINE_SMILES_MONOMER { $$ = $1; $$.remove_prefix(1); $$.remove_suffix(1); }

monomer_id: SINGLE_CHARACTER_MONOMER { helm_parser.add_residue_name($1); $$ = $1; }
       | MULTI_CHARACTER_MONOMER {
            auto tmp = std::string_view{$1};
            tmp.remove_prefix(1);
            tmp.remove_suffix(1);
            helm_parser.add_residue_name(tmp);
            $$ = $1;
       }
       | UNKNOWN_MONOMER {
            helm_parser.add_wildcard_or_unknown_residue();
            $$ = $1;
       }
       | MISSING_MONOMER { $$ = $1; }
       | MONOMER_WILDCARD {
            helm_parser.add_wildcard_or_unknown_residue();
            $$ = $1;
       }
       ;
blob: UNKNOWN_SEQUENCE { helm_parser.add_monomer($1, false, false, false, {}); }


/*
 * The structure of a connection is as follows
 *     * A custom connection can be formed between two monomers
 *     * Rgroups that partake in a custom connection must be specified except in
 *       the case of BLOBs or ambiguous connections.
 *     * If multiple polymers form the same type of custom bond, you can alias
 *       the polymer ids through a polymer group and use the polymer group id
 *       as the connection polymer.
 *     * You can define multiple connection monomers is they for the same type
 *       of bond, e.g., POLYMER,POLYMER,(MONOMER,MONOMER):R3-MONOMER:R1
 *
 */
connection: connection_polymer ',' connection_polymer ','
            connection_monomer ':' attachment_point '-'
            connection_monomer ':' attachment_point {
                helm_parser.add_connection({
                                            $1,  // from_id
                                            $5,  // from_res
                                            $7,  // from_rgroup
                                            $3,  // to_id
                                            $9,  // to_res
                                            $11, // to_rgroup
                                            {}   // annotation
                                          });
        }
        | connection_polymer ',' connection_polymer ','
            connection_monomer ':' attachment_point '-'
            connection_monomer ':' attachment_point annotation {
                helm_parser.add_connection({
                                            $1,  // from_id
                                            $5,  // from_res
                                            $7,  // from_rgroup
                                            $3,  // to_id
                                            $9,  // to_res
                                            $11, // to_rgroup
                                            $12  // annotation
                                          });
        }
        ;
connection_polymer: POLYMER_ID | POLYMER_GROUP_ID;
attachment_point: HYDROGEN_PAIRING | RGROUP | UNDEFINED_RESIDUE_NUMBER_OR_RGROUP;
connection_monomer: CONNECTION_RESIDUE { $$ = $1; }
                  | UNDEFINED_RESIDUE_NUMBER_OR_RGROUP { $$ = $1; }
                  | '(' residue_and_list ')' { $$ = $2; }
                  | '(' residue_or_list ')' { $$ = $2; }
                  ;
residue_and_list: CONNECTION_RESIDUE '+' CONNECTION_RESIDUE {
                    $$ = {$1.data(), $1.size() + 1 + $3.size()};
                }
                | residue_and_list '+' CONNECTION_RESIDUE {
                    $$ = {$1.data(), $1.size() + 1 + $3.size()};
                }
                ;
residue_or_list: CONNECTION_RESIDUE ',' CONNECTION_RESIDUE {
                    $$ = {$1.data(), $1.size() + 1 + $3.size()};
                }
                | residue_or_list ',' CONNECTION_RESIDUE {
                    $$ = {$1.data(), $1.size() + 1 + $3.size()};
                }
                ;

/* Polymer group */
polymer_group: POLYMER_GROUP_ID '(' polymer_and_list ')' { helm_parser.add_polymer_group($1, $3, true); }
             | POLYMER_GROUP_ID '(' polymer_or_list ')' { helm_parser.add_polymer_group($1, $3, false); }
             ;
polymer_or_list: polymer_group_item ',' polymer_group_item {
                $$ = {$1.data(), $1.size() + 1 + $3.size()};
            }
            | polymer_or_list ',' polymer_group_item {
                $$ = {$1.data(), $1.size() + 1 + $3.size()};
            }
polymer_and_list: polymer_group_item '+' polymer_group_item {
                $$ = {$1.data(), $1.size() + 1 + $3.size()};
            }
            | polymer_and_list '+' polymer_group_item {
                $$ = {$1.data(), $1.size() + 1 + $3.size()};
            }
            ;
polymer_group_item: POLYMER_ID { $$ = $1; }
                  | POLYMER_GROUP_ID { $$  = $1; }
                  | POLYMER_ID ':' POLYMER_GROUP_RATIO {
                    $$ = {$1.data(), $1.size() + 1 + $3.size()};
                  }
                  | POLYMER_GROUP_ID ':' POLYMER_GROUP_RATIO {
                    $$ = {$1.data(), $1.size() + 1 + $3.size()};
                  }
                  ;

/* General tokens */
annotation: ANNOTATION { $$ = $1; $$.remove_prefix(1); $$.remove_suffix(1); };
%%


void
helm::TokenParser::error(const location_type &l, const std::string& err_message )
{
    // save the pointer to the beginning of the current token
    helm_parser.saveError(l.begin.column, err_message);
}
