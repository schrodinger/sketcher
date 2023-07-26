// A Bison parser, made by GNU Bison 3.8.2.

// Skeleton implementation for Bison LALR(1) parsers in C++

// Copyright (C) 2002-2015, 2018-2021 Free Software Foundation, Inc.

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

// As a special exception, you may create a larger work that contains
// part or all of the Bison parser skeleton and distribute that work
// under terms of your choice, so long as that work isn't itself a
// parser generator using the skeleton or a modified version thereof
// as a parser skeleton.  Alternatively, if you modify or redistribute
// the parser skeleton itself, you may (at your option) remove this
// special exception, which will cause the skeleton and the resulting
// Bison output files to be licensed under the GNU General Public
// License without this special exception.

// This special exception was added by the Free Software Foundation in
// version 2.2 of Bison.

// DO NOT RELY ON FEATURES THAT ARE NOT DOCUMENTED in the manual,
// especially those whose name start with YY_ or yy_.  They are
// private implementation details that can be changed or removed.

// First part of user prologue.
#line 1 "../helm_parser.yy"

/*
 * This is a bison file which is used to generate a C++ parser for parsing a
 * HELM string.
 *
 * The grammar defines the specification for a valid HELMV2.0 string and uses
 * a helm::HelmParser instance to retrieve the parsed helm information.
 */

#line 51 "helm_parser.tab.cpp"

#include "helm_parser.tab.hh"

// Unqualified %code blocks.
#line 34 "../helm_parser.yy"

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

#line 78 "helm_parser.tab.cpp"

#ifndef YY_
#if defined YYENABLE_NLS && YYENABLE_NLS
#if ENABLE_NLS
#include <libintl.h> // FIXME: INFRINGES ON USER NAME SPACE.
#define YY_(msgid) dgettext("bison-runtime", msgid)
#endif
#endif
#ifndef YY_
#define YY_(msgid) msgid
#endif
#endif

// Whether we are compiled with exception support.
#ifndef YY_EXCEPTIONS
#if defined __GNUC__ && !defined __EXCEPTIONS
#define YY_EXCEPTIONS 0
#else
#define YY_EXCEPTIONS 1
#endif
#endif

#define YYRHSLOC(Rhs, K) ((Rhs)[K].location)
/* YYLLOC_DEFAULT -- Set CURRENT to span from RHS[1] to RHS[N].
   If N is 0, then set CURRENT to the empty location which ends
   the previous symbol: RHS[0] (always defined).  */

#ifndef YYLLOC_DEFAULT
#define YYLLOC_DEFAULT(Current, Rhs, N)                             \
    do                                                              \
        if (N) {                                                    \
            (Current).begin = YYRHSLOC(Rhs, 1).begin;               \
            (Current).end = YYRHSLOC(Rhs, N).end;                   \
        } else {                                                    \
            (Current).begin = (Current).end = YYRHSLOC(Rhs, 0).end; \
        }                                                           \
    while (false)
#endif

// Enable debugging if requested.
#if YYDEBUG

// A pseudo ostream that takes yydebug_ into account.
#define YYCDEBUG  \
    if (yydebug_) \
    (*yycdebug_)

#define YY_SYMBOL_PRINT(Title, Symbol)     \
    do {                                   \
        if (yydebug_) {                    \
            *yycdebug_ << Title << ' ';    \
            yy_print_(*yycdebug_, Symbol); \
            *yycdebug_ << '\n';            \
        }                                  \
    } while (false)

#define YY_REDUCE_PRINT(Rule)       \
    do {                            \
        if (yydebug_)               \
            yy_reduce_print_(Rule); \
    } while (false)

#define YY_STACK_PRINT()       \
    do {                       \
        if (yydebug_)          \
            yy_stack_print_(); \
    } while (false)

#else // !YYDEBUG

#define YYCDEBUG \
    if (false)   \
    std::cerr
#define YY_SYMBOL_PRINT(Title, Symbol) YY_USE(Symbol)
#define YY_REDUCE_PRINT(Rule) static_cast<void>(0)
#define YY_STACK_PRINT() static_cast<void>(0)

#endif // !YYDEBUG

#define yyerrok (yyerrstatus_ = 0)
#define yyclearin (yyla.clear())

#define YYACCEPT goto yyacceptlab
#define YYABORT goto yyabortlab
#define YYERROR goto yyerrorlab
#define YYRECOVERING() (!!yyerrstatus_)

#line 15 "../helm_parser.yy"
namespace helm
{
#line 171 "helm_parser.tab.cpp"

/// Build a parser object.
TokenParser::TokenParser(TokenScanner& scanner_yyarg,
                         HelmParser& helm_parser_yyarg)
#if YYDEBUG
    :
    yydebug_(false),
    yycdebug_(&std::cerr),
#else
    :
#endif
    scanner(scanner_yyarg),
    helm_parser(helm_parser_yyarg)
{
}

TokenParser::~TokenParser()
{
}

TokenParser::syntax_error::~syntax_error() YY_NOEXCEPT YY_NOTHROW
{
}

/*---------.
| symbol.  |
`---------*/

// basic_symbol.
template <typename Base>
TokenParser::basic_symbol<Base>::basic_symbol(const basic_symbol& that) :
    Base(that),
    value(),
    location(that.location)
{
    switch (this->kind()) {
        case symbol_kind::S_connections:          // connections
        case symbol_kind::S_polymer_groups:       // polymer_groups
        case symbol_kind::S_monomer_sequence:     // monomer_sequence
        case symbol_kind::S_branch_monomer_group: // branch_monomer_group
            value.copy<size_t>(YY_MOVE(that.value));
            break;

        case symbol_kind::S_BLOB_ID:    // BLOB_ID
        case symbol_kind::S_CHEM_ID:    // CHEM_ID
        case symbol_kind::S_PEPTIDE_ID: // PEPTIDE_ID
        case symbol_kind::S_RNA_ID:     // RNA_ID
        case symbol_kind::S_POLYMER_ID: // POLYMER_ID
        case symbol_kind::
            S_SINGLE_CHARACTER_MONOMER:              // SINGLE_CHARACTER_MONOMER
        case symbol_kind::S_MULTI_CHARACTER_MONOMER: // MULTI_CHARACTER_MONOMER
        case symbol_kind::S_UNKNOWN_MONOMER:         // UNKNOWN_MONOMER
        case symbol_kind::S_MISSING_MONOMER:         // MISSING_MONOMER
        case symbol_kind::S_MONOMER_WILDCARD:        // MONOMER_WILDCARD
        case symbol_kind::S_MONOMER_RATIO:           // MONOMER_RATIO
        case symbol_kind::S_INLINE_SMILES_MONOMER:   // INLINE_SMILES_MONOMER
        case symbol_kind::S_ANNOTATION:              // ANNOTATION
        case symbol_kind::S_REPETITIONS:             // REPETITIONS
        case symbol_kind::S_UNKNOWN_SEQUENCE:        // UNKNOWN_SEQUENCE
        case symbol_kind::S_VERSION_TOKEN:           // VERSION_TOKEN
        case symbol_kind::S_HYDROGEN_PAIRING:        // HYDROGEN_PAIRING
        case symbol_kind::S_RGROUP:                  // RGROUP
        case symbol_kind::
            S_UNDEFINED_RESIDUE_NUMBER_OR_RGROUP: // UNDEFINED_RESIDUE_NUMBER_OR_RGROUP
        case symbol_kind::S_CONNECTION_RESIDUE:  // CONNECTION_RESIDUE
        case symbol_kind::S_POLYMER_GROUP_RATIO: // POLYMER_GROUP_RATIO
        case symbol_kind::S_POLYMER_GROUP_ID:    // POLYMER_GROUP_ID
        case symbol_kind::
            S_EXTENDED_ANNOTATIONS_TOKEN:         // EXTENDED_ANNOTATIONS_TOKEN
        case symbol_kind::S_extended_annotations: // extended_annotations
        case symbol_kind::S_version:              // version
        case symbol_kind::S_polymer:              // polymer
        case symbol_kind::S_repetitions:          // repetitions
        case symbol_kind::S_monomer_list:         // monomer_list
        case symbol_kind::S_monomer_and_list:     // monomer_and_list
        case symbol_kind::S_monomer_or_list:      // monomer_or_list
        case symbol_kind::S_monomer_list_item:    // monomer_list_item
        case symbol_kind::S_smiles_monomer:       // smiles_monomer
        case symbol_kind::S_monomer_id:           // monomer_id
        case symbol_kind::S_blob:                 // blob
        case symbol_kind::S_connection_polymer:   // connection_polymer
        case symbol_kind::S_attachment_point:     // attachment_point
        case symbol_kind::S_connection_monomer:   // connection_monomer
        case symbol_kind::S_residue_and_list:     // residue_and_list
        case symbol_kind::S_residue_or_list:      // residue_or_list
        case symbol_kind::S_polymer_or_list:      // polymer_or_list
        case symbol_kind::S_polymer_and_list:     // polymer_and_list
        case symbol_kind::S_polymer_group_item:   // polymer_group_item
        case symbol_kind::S_annotation:           // annotation
            value.copy<std::string_view>(YY_MOVE(that.value));
            break;

        default:
            break;
    }
}

template <typename Base> TokenParser::symbol_kind_type
TokenParser::basic_symbol<Base>::type_get() const YY_NOEXCEPT
{
    return this->kind();
}

template <typename Base>
bool TokenParser::basic_symbol<Base>::empty() const YY_NOEXCEPT
{
    return this->kind() == symbol_kind::S_YYEMPTY;
}

template <typename Base>
void TokenParser::basic_symbol<Base>::move(basic_symbol& s)
{
    super_type::move(s);
    switch (this->kind()) {
        case symbol_kind::S_connections:          // connections
        case symbol_kind::S_polymer_groups:       // polymer_groups
        case symbol_kind::S_monomer_sequence:     // monomer_sequence
        case symbol_kind::S_branch_monomer_group: // branch_monomer_group
            value.move<size_t>(YY_MOVE(s.value));
            break;

        case symbol_kind::S_BLOB_ID:    // BLOB_ID
        case symbol_kind::S_CHEM_ID:    // CHEM_ID
        case symbol_kind::S_PEPTIDE_ID: // PEPTIDE_ID
        case symbol_kind::S_RNA_ID:     // RNA_ID
        case symbol_kind::S_POLYMER_ID: // POLYMER_ID
        case symbol_kind::
            S_SINGLE_CHARACTER_MONOMER:              // SINGLE_CHARACTER_MONOMER
        case symbol_kind::S_MULTI_CHARACTER_MONOMER: // MULTI_CHARACTER_MONOMER
        case symbol_kind::S_UNKNOWN_MONOMER:         // UNKNOWN_MONOMER
        case symbol_kind::S_MISSING_MONOMER:         // MISSING_MONOMER
        case symbol_kind::S_MONOMER_WILDCARD:        // MONOMER_WILDCARD
        case symbol_kind::S_MONOMER_RATIO:           // MONOMER_RATIO
        case symbol_kind::S_INLINE_SMILES_MONOMER:   // INLINE_SMILES_MONOMER
        case symbol_kind::S_ANNOTATION:              // ANNOTATION
        case symbol_kind::S_REPETITIONS:             // REPETITIONS
        case symbol_kind::S_UNKNOWN_SEQUENCE:        // UNKNOWN_SEQUENCE
        case symbol_kind::S_VERSION_TOKEN:           // VERSION_TOKEN
        case symbol_kind::S_HYDROGEN_PAIRING:        // HYDROGEN_PAIRING
        case symbol_kind::S_RGROUP:                  // RGROUP
        case symbol_kind::
            S_UNDEFINED_RESIDUE_NUMBER_OR_RGROUP: // UNDEFINED_RESIDUE_NUMBER_OR_RGROUP
        case symbol_kind::S_CONNECTION_RESIDUE:  // CONNECTION_RESIDUE
        case symbol_kind::S_POLYMER_GROUP_RATIO: // POLYMER_GROUP_RATIO
        case symbol_kind::S_POLYMER_GROUP_ID:    // POLYMER_GROUP_ID
        case symbol_kind::
            S_EXTENDED_ANNOTATIONS_TOKEN:         // EXTENDED_ANNOTATIONS_TOKEN
        case symbol_kind::S_extended_annotations: // extended_annotations
        case symbol_kind::S_version:              // version
        case symbol_kind::S_polymer:              // polymer
        case symbol_kind::S_repetitions:          // repetitions
        case symbol_kind::S_monomer_list:         // monomer_list
        case symbol_kind::S_monomer_and_list:     // monomer_and_list
        case symbol_kind::S_monomer_or_list:      // monomer_or_list
        case symbol_kind::S_monomer_list_item:    // monomer_list_item
        case symbol_kind::S_smiles_monomer:       // smiles_monomer
        case symbol_kind::S_monomer_id:           // monomer_id
        case symbol_kind::S_blob:                 // blob
        case symbol_kind::S_connection_polymer:   // connection_polymer
        case symbol_kind::S_attachment_point:     // attachment_point
        case symbol_kind::S_connection_monomer:   // connection_monomer
        case symbol_kind::S_residue_and_list:     // residue_and_list
        case symbol_kind::S_residue_or_list:      // residue_or_list
        case symbol_kind::S_polymer_or_list:      // polymer_or_list
        case symbol_kind::S_polymer_and_list:     // polymer_and_list
        case symbol_kind::S_polymer_group_item:   // polymer_group_item
        case symbol_kind::S_annotation:           // annotation
            value.move<std::string_view>(YY_MOVE(s.value));
            break;

        default:
            break;
    }

    location = YY_MOVE(s.location);
}

// by_kind.
TokenParser::by_kind::by_kind() YY_NOEXCEPT : kind_(symbol_kind::S_YYEMPTY)
{
}

#if 201103L <= YY_CPLUSPLUS
TokenParser::by_kind::by_kind(by_kind&& that) YY_NOEXCEPT : kind_(that.kind_)
{
    that.clear();
}
#endif

TokenParser::by_kind::by_kind(const by_kind& that) YY_NOEXCEPT
    : kind_(that.kind_)
{
}

TokenParser::by_kind::by_kind(token_kind_type t) YY_NOEXCEPT
    : kind_(yytranslate_(t))
{
}

void TokenParser::by_kind::clear() YY_NOEXCEPT
{
    kind_ = symbol_kind::S_YYEMPTY;
}

void TokenParser::by_kind::move(by_kind& that)
{
    kind_ = that.kind_;
    that.clear();
}

TokenParser::symbol_kind_type TokenParser::by_kind::kind() const YY_NOEXCEPT
{
    return kind_;
}

TokenParser::symbol_kind_type TokenParser::by_kind::type_get() const YY_NOEXCEPT
{
    return this->kind();
}

// by_state.
TokenParser::by_state::by_state() YY_NOEXCEPT : state(empty_state)
{
}

TokenParser::by_state::by_state(const by_state& that) YY_NOEXCEPT
    : state(that.state)
{
}

void TokenParser::by_state::clear() YY_NOEXCEPT
{
    state = empty_state;
}

void TokenParser::by_state::move(by_state& that)
{
    state = that.state;
    that.clear();
}

TokenParser::by_state::by_state(state_type s) YY_NOEXCEPT : state(s)
{
}

TokenParser::symbol_kind_type TokenParser::by_state::kind() const YY_NOEXCEPT
{
    if (state == empty_state)
        return symbol_kind::S_YYEMPTY;
    else
        return YY_CAST(symbol_kind_type, yystos_[+state]);
}

TokenParser::stack_symbol_type::stack_symbol_type()
{
}

TokenParser::stack_symbol_type::stack_symbol_type(YY_RVREF(stack_symbol_type)
                                                      that) :
    super_type(YY_MOVE(that.state), YY_MOVE(that.location))
{
    switch (that.kind()) {
        case symbol_kind::S_connections:          // connections
        case symbol_kind::S_polymer_groups:       // polymer_groups
        case symbol_kind::S_monomer_sequence:     // monomer_sequence
        case symbol_kind::S_branch_monomer_group: // branch_monomer_group
            value.YY_MOVE_OR_COPY<size_t>(YY_MOVE(that.value));
            break;

        case symbol_kind::S_BLOB_ID:    // BLOB_ID
        case symbol_kind::S_CHEM_ID:    // CHEM_ID
        case symbol_kind::S_PEPTIDE_ID: // PEPTIDE_ID
        case symbol_kind::S_RNA_ID:     // RNA_ID
        case symbol_kind::S_POLYMER_ID: // POLYMER_ID
        case symbol_kind::
            S_SINGLE_CHARACTER_MONOMER:              // SINGLE_CHARACTER_MONOMER
        case symbol_kind::S_MULTI_CHARACTER_MONOMER: // MULTI_CHARACTER_MONOMER
        case symbol_kind::S_UNKNOWN_MONOMER:         // UNKNOWN_MONOMER
        case symbol_kind::S_MISSING_MONOMER:         // MISSING_MONOMER
        case symbol_kind::S_MONOMER_WILDCARD:        // MONOMER_WILDCARD
        case symbol_kind::S_MONOMER_RATIO:           // MONOMER_RATIO
        case symbol_kind::S_INLINE_SMILES_MONOMER:   // INLINE_SMILES_MONOMER
        case symbol_kind::S_ANNOTATION:              // ANNOTATION
        case symbol_kind::S_REPETITIONS:             // REPETITIONS
        case symbol_kind::S_UNKNOWN_SEQUENCE:        // UNKNOWN_SEQUENCE
        case symbol_kind::S_VERSION_TOKEN:           // VERSION_TOKEN
        case symbol_kind::S_HYDROGEN_PAIRING:        // HYDROGEN_PAIRING
        case symbol_kind::S_RGROUP:                  // RGROUP
        case symbol_kind::
            S_UNDEFINED_RESIDUE_NUMBER_OR_RGROUP: // UNDEFINED_RESIDUE_NUMBER_OR_RGROUP
        case symbol_kind::S_CONNECTION_RESIDUE:  // CONNECTION_RESIDUE
        case symbol_kind::S_POLYMER_GROUP_RATIO: // POLYMER_GROUP_RATIO
        case symbol_kind::S_POLYMER_GROUP_ID:    // POLYMER_GROUP_ID
        case symbol_kind::
            S_EXTENDED_ANNOTATIONS_TOKEN:         // EXTENDED_ANNOTATIONS_TOKEN
        case symbol_kind::S_extended_annotations: // extended_annotations
        case symbol_kind::S_version:              // version
        case symbol_kind::S_polymer:              // polymer
        case symbol_kind::S_repetitions:          // repetitions
        case symbol_kind::S_monomer_list:         // monomer_list
        case symbol_kind::S_monomer_and_list:     // monomer_and_list
        case symbol_kind::S_monomer_or_list:      // monomer_or_list
        case symbol_kind::S_monomer_list_item:    // monomer_list_item
        case symbol_kind::S_smiles_monomer:       // smiles_monomer
        case symbol_kind::S_monomer_id:           // monomer_id
        case symbol_kind::S_blob:                 // blob
        case symbol_kind::S_connection_polymer:   // connection_polymer
        case symbol_kind::S_attachment_point:     // attachment_point
        case symbol_kind::S_connection_monomer:   // connection_monomer
        case symbol_kind::S_residue_and_list:     // residue_and_list
        case symbol_kind::S_residue_or_list:      // residue_or_list
        case symbol_kind::S_polymer_or_list:      // polymer_or_list
        case symbol_kind::S_polymer_and_list:     // polymer_and_list
        case symbol_kind::S_polymer_group_item:   // polymer_group_item
        case symbol_kind::S_annotation:           // annotation
            value.YY_MOVE_OR_COPY<std::string_view>(YY_MOVE(that.value));
            break;

        default:
            break;
    }

#if 201103L <= YY_CPLUSPLUS
    // that is emptied.
    that.state = empty_state;
#endif
}

TokenParser::stack_symbol_type::stack_symbol_type(state_type s,
                                                  YY_MOVE_REF(symbol_type)
                                                      that) :
    super_type(s, YY_MOVE(that.location))
{
    switch (that.kind()) {
        case symbol_kind::S_connections:          // connections
        case symbol_kind::S_polymer_groups:       // polymer_groups
        case symbol_kind::S_monomer_sequence:     // monomer_sequence
        case symbol_kind::S_branch_monomer_group: // branch_monomer_group
            value.move<size_t>(YY_MOVE(that.value));
            break;

        case symbol_kind::S_BLOB_ID:    // BLOB_ID
        case symbol_kind::S_CHEM_ID:    // CHEM_ID
        case symbol_kind::S_PEPTIDE_ID: // PEPTIDE_ID
        case symbol_kind::S_RNA_ID:     // RNA_ID
        case symbol_kind::S_POLYMER_ID: // POLYMER_ID
        case symbol_kind::
            S_SINGLE_CHARACTER_MONOMER:              // SINGLE_CHARACTER_MONOMER
        case symbol_kind::S_MULTI_CHARACTER_MONOMER: // MULTI_CHARACTER_MONOMER
        case symbol_kind::S_UNKNOWN_MONOMER:         // UNKNOWN_MONOMER
        case symbol_kind::S_MISSING_MONOMER:         // MISSING_MONOMER
        case symbol_kind::S_MONOMER_WILDCARD:        // MONOMER_WILDCARD
        case symbol_kind::S_MONOMER_RATIO:           // MONOMER_RATIO
        case symbol_kind::S_INLINE_SMILES_MONOMER:   // INLINE_SMILES_MONOMER
        case symbol_kind::S_ANNOTATION:              // ANNOTATION
        case symbol_kind::S_REPETITIONS:             // REPETITIONS
        case symbol_kind::S_UNKNOWN_SEQUENCE:        // UNKNOWN_SEQUENCE
        case symbol_kind::S_VERSION_TOKEN:           // VERSION_TOKEN
        case symbol_kind::S_HYDROGEN_PAIRING:        // HYDROGEN_PAIRING
        case symbol_kind::S_RGROUP:                  // RGROUP
        case symbol_kind::
            S_UNDEFINED_RESIDUE_NUMBER_OR_RGROUP: // UNDEFINED_RESIDUE_NUMBER_OR_RGROUP
        case symbol_kind::S_CONNECTION_RESIDUE:  // CONNECTION_RESIDUE
        case symbol_kind::S_POLYMER_GROUP_RATIO: // POLYMER_GROUP_RATIO
        case symbol_kind::S_POLYMER_GROUP_ID:    // POLYMER_GROUP_ID
        case symbol_kind::
            S_EXTENDED_ANNOTATIONS_TOKEN:         // EXTENDED_ANNOTATIONS_TOKEN
        case symbol_kind::S_extended_annotations: // extended_annotations
        case symbol_kind::S_version:              // version
        case symbol_kind::S_polymer:              // polymer
        case symbol_kind::S_repetitions:          // repetitions
        case symbol_kind::S_monomer_list:         // monomer_list
        case symbol_kind::S_monomer_and_list:     // monomer_and_list
        case symbol_kind::S_monomer_or_list:      // monomer_or_list
        case symbol_kind::S_monomer_list_item:    // monomer_list_item
        case symbol_kind::S_smiles_monomer:       // smiles_monomer
        case symbol_kind::S_monomer_id:           // monomer_id
        case symbol_kind::S_blob:                 // blob
        case symbol_kind::S_connection_polymer:   // connection_polymer
        case symbol_kind::S_attachment_point:     // attachment_point
        case symbol_kind::S_connection_monomer:   // connection_monomer
        case symbol_kind::S_residue_and_list:     // residue_and_list
        case symbol_kind::S_residue_or_list:      // residue_or_list
        case symbol_kind::S_polymer_or_list:      // polymer_or_list
        case symbol_kind::S_polymer_and_list:     // polymer_and_list
        case symbol_kind::S_polymer_group_item:   // polymer_group_item
        case symbol_kind::S_annotation:           // annotation
            value.move<std::string_view>(YY_MOVE(that.value));
            break;

        default:
            break;
    }

    // that is emptied.
    that.kind_ = symbol_kind::S_YYEMPTY;
}

#if YY_CPLUSPLUS < 201103L
TokenParser::stack_symbol_type&
TokenParser::stack_symbol_type::operator=(const stack_symbol_type& that)
{
    state = that.state;
    switch (that.kind()) {
        case symbol_kind::S_connections:          // connections
        case symbol_kind::S_polymer_groups:       // polymer_groups
        case symbol_kind::S_monomer_sequence:     // monomer_sequence
        case symbol_kind::S_branch_monomer_group: // branch_monomer_group
            value.copy<size_t>(that.value);
            break;

        case symbol_kind::S_BLOB_ID:    // BLOB_ID
        case symbol_kind::S_CHEM_ID:    // CHEM_ID
        case symbol_kind::S_PEPTIDE_ID: // PEPTIDE_ID
        case symbol_kind::S_RNA_ID:     // RNA_ID
        case symbol_kind::S_POLYMER_ID: // POLYMER_ID
        case symbol_kind::
            S_SINGLE_CHARACTER_MONOMER:              // SINGLE_CHARACTER_MONOMER
        case symbol_kind::S_MULTI_CHARACTER_MONOMER: // MULTI_CHARACTER_MONOMER
        case symbol_kind::S_UNKNOWN_MONOMER:         // UNKNOWN_MONOMER
        case symbol_kind::S_MISSING_MONOMER:         // MISSING_MONOMER
        case symbol_kind::S_MONOMER_WILDCARD:        // MONOMER_WILDCARD
        case symbol_kind::S_MONOMER_RATIO:           // MONOMER_RATIO
        case symbol_kind::S_INLINE_SMILES_MONOMER:   // INLINE_SMILES_MONOMER
        case symbol_kind::S_ANNOTATION:              // ANNOTATION
        case symbol_kind::S_REPETITIONS:             // REPETITIONS
        case symbol_kind::S_UNKNOWN_SEQUENCE:        // UNKNOWN_SEQUENCE
        case symbol_kind::S_VERSION_TOKEN:           // VERSION_TOKEN
        case symbol_kind::S_HYDROGEN_PAIRING:        // HYDROGEN_PAIRING
        case symbol_kind::S_RGROUP:                  // RGROUP
        case symbol_kind::
            S_UNDEFINED_RESIDUE_NUMBER_OR_RGROUP: // UNDEFINED_RESIDUE_NUMBER_OR_RGROUP
        case symbol_kind::S_CONNECTION_RESIDUE:  // CONNECTION_RESIDUE
        case symbol_kind::S_POLYMER_GROUP_RATIO: // POLYMER_GROUP_RATIO
        case symbol_kind::S_POLYMER_GROUP_ID:    // POLYMER_GROUP_ID
        case symbol_kind::
            S_EXTENDED_ANNOTATIONS_TOKEN:         // EXTENDED_ANNOTATIONS_TOKEN
        case symbol_kind::S_extended_annotations: // extended_annotations
        case symbol_kind::S_version:              // version
        case symbol_kind::S_polymer:              // polymer
        case symbol_kind::S_repetitions:          // repetitions
        case symbol_kind::S_monomer_list:         // monomer_list
        case symbol_kind::S_monomer_and_list:     // monomer_and_list
        case symbol_kind::S_monomer_or_list:      // monomer_or_list
        case symbol_kind::S_monomer_list_item:    // monomer_list_item
        case symbol_kind::S_smiles_monomer:       // smiles_monomer
        case symbol_kind::S_monomer_id:           // monomer_id
        case symbol_kind::S_blob:                 // blob
        case symbol_kind::S_connection_polymer:   // connection_polymer
        case symbol_kind::S_attachment_point:     // attachment_point
        case symbol_kind::S_connection_monomer:   // connection_monomer
        case symbol_kind::S_residue_and_list:     // residue_and_list
        case symbol_kind::S_residue_or_list:      // residue_or_list
        case symbol_kind::S_polymer_or_list:      // polymer_or_list
        case symbol_kind::S_polymer_and_list:     // polymer_and_list
        case symbol_kind::S_polymer_group_item:   // polymer_group_item
        case symbol_kind::S_annotation:           // annotation
            value.copy<std::string_view>(that.value);
            break;

        default:
            break;
    }

    location = that.location;
    return *this;
}

TokenParser::stack_symbol_type&
TokenParser::stack_symbol_type::operator=(stack_symbol_type& that)
{
    state = that.state;
    switch (that.kind()) {
        case symbol_kind::S_connections:          // connections
        case symbol_kind::S_polymer_groups:       // polymer_groups
        case symbol_kind::S_monomer_sequence:     // monomer_sequence
        case symbol_kind::S_branch_monomer_group: // branch_monomer_group
            value.move<size_t>(that.value);
            break;

        case symbol_kind::S_BLOB_ID:    // BLOB_ID
        case symbol_kind::S_CHEM_ID:    // CHEM_ID
        case symbol_kind::S_PEPTIDE_ID: // PEPTIDE_ID
        case symbol_kind::S_RNA_ID:     // RNA_ID
        case symbol_kind::S_POLYMER_ID: // POLYMER_ID
        case symbol_kind::
            S_SINGLE_CHARACTER_MONOMER:              // SINGLE_CHARACTER_MONOMER
        case symbol_kind::S_MULTI_CHARACTER_MONOMER: // MULTI_CHARACTER_MONOMER
        case symbol_kind::S_UNKNOWN_MONOMER:         // UNKNOWN_MONOMER
        case symbol_kind::S_MISSING_MONOMER:         // MISSING_MONOMER
        case symbol_kind::S_MONOMER_WILDCARD:        // MONOMER_WILDCARD
        case symbol_kind::S_MONOMER_RATIO:           // MONOMER_RATIO
        case symbol_kind::S_INLINE_SMILES_MONOMER:   // INLINE_SMILES_MONOMER
        case symbol_kind::S_ANNOTATION:              // ANNOTATION
        case symbol_kind::S_REPETITIONS:             // REPETITIONS
        case symbol_kind::S_UNKNOWN_SEQUENCE:        // UNKNOWN_SEQUENCE
        case symbol_kind::S_VERSION_TOKEN:           // VERSION_TOKEN
        case symbol_kind::S_HYDROGEN_PAIRING:        // HYDROGEN_PAIRING
        case symbol_kind::S_RGROUP:                  // RGROUP
        case symbol_kind::
            S_UNDEFINED_RESIDUE_NUMBER_OR_RGROUP: // UNDEFINED_RESIDUE_NUMBER_OR_RGROUP
        case symbol_kind::S_CONNECTION_RESIDUE:  // CONNECTION_RESIDUE
        case symbol_kind::S_POLYMER_GROUP_RATIO: // POLYMER_GROUP_RATIO
        case symbol_kind::S_POLYMER_GROUP_ID:    // POLYMER_GROUP_ID
        case symbol_kind::
            S_EXTENDED_ANNOTATIONS_TOKEN:         // EXTENDED_ANNOTATIONS_TOKEN
        case symbol_kind::S_extended_annotations: // extended_annotations
        case symbol_kind::S_version:              // version
        case symbol_kind::S_polymer:              // polymer
        case symbol_kind::S_repetitions:          // repetitions
        case symbol_kind::S_monomer_list:         // monomer_list
        case symbol_kind::S_monomer_and_list:     // monomer_and_list
        case symbol_kind::S_monomer_or_list:      // monomer_or_list
        case symbol_kind::S_monomer_list_item:    // monomer_list_item
        case symbol_kind::S_smiles_monomer:       // smiles_monomer
        case symbol_kind::S_monomer_id:           // monomer_id
        case symbol_kind::S_blob:                 // blob
        case symbol_kind::S_connection_polymer:   // connection_polymer
        case symbol_kind::S_attachment_point:     // attachment_point
        case symbol_kind::S_connection_monomer:   // connection_monomer
        case symbol_kind::S_residue_and_list:     // residue_and_list
        case symbol_kind::S_residue_or_list:      // residue_or_list
        case symbol_kind::S_polymer_or_list:      // polymer_or_list
        case symbol_kind::S_polymer_and_list:     // polymer_and_list
        case symbol_kind::S_polymer_group_item:   // polymer_group_item
        case symbol_kind::S_annotation:           // annotation
            value.move<std::string_view>(that.value);
            break;

        default:
            break;
    }

    location = that.location;
    // that is emptied.
    that.state = empty_state;
    return *this;
}
#endif

template <typename Base> void
TokenParser::yy_destroy_(const char* yymsg, basic_symbol<Base>& yysym) const
{
    if (yymsg)
        YY_SYMBOL_PRINT(yymsg, yysym);
}

#if YYDEBUG
template <typename Base> void
TokenParser::yy_print_(std::ostream& yyo, const basic_symbol<Base>& yysym) const
{
    std::ostream& yyoutput = yyo;
    YY_USE(yyoutput);
    if (yysym.empty())
        yyo << "empty symbol";
    else {
        symbol_kind_type yykind = yysym.kind();
        yyo << (yykind < YYNTOKENS ? "token" : "nterm") << ' ' << yysym.name()
            << " (" << yysym.location << ": ";
        YY_USE(yykind);
        yyo << ')';
    }
}
#endif

void TokenParser::yypush_(const char* m, YY_MOVE_REF(stack_symbol_type) sym)
{
    if (m)
        YY_SYMBOL_PRINT(m, sym);
    yystack_.push(YY_MOVE(sym));
}

void TokenParser::yypush_(const char* m, state_type s,
                          YY_MOVE_REF(symbol_type) sym)
{
#if 201103L <= YY_CPLUSPLUS
    yypush_(m, stack_symbol_type(s, std::move(sym)));
#else
    stack_symbol_type ss(s, sym);
    yypush_(m, ss);
#endif
}

void TokenParser::yypop_(int n) YY_NOEXCEPT
{
    yystack_.pop(n);
}

#if YYDEBUG
std::ostream& TokenParser::debug_stream() const
{
    return *yycdebug_;
}

void TokenParser::set_debug_stream(std::ostream& o)
{
    yycdebug_ = &o;
}

TokenParser::debug_level_type TokenParser::debug_level() const
{
    return yydebug_;
}

void TokenParser::set_debug_level(debug_level_type l)
{
    yydebug_ = l;
}
#endif // YYDEBUG

TokenParser::state_type TokenParser::yy_lr_goto_state_(state_type yystate,
                                                       int yysym)
{
    int yyr = yypgoto_[yysym - YYNTOKENS] + yystate;
    if (0 <= yyr && yyr <= yylast_ && yycheck_[yyr] == yystate)
        return yytable_[yyr];
    else
        return yydefgoto_[yysym - YYNTOKENS];
}

bool TokenParser::yy_pact_value_is_default_(int yyvalue) YY_NOEXCEPT
{
    return yyvalue == yypact_ninf_;
}

bool TokenParser::yy_table_value_is_error_(int yyvalue) YY_NOEXCEPT
{
    return yyvalue == yytable_ninf_;
}

int TokenParser::operator()()
{
    return parse();
}

int TokenParser::parse()
{
    int yyn;
    /// Length of the RHS of the rule being reduced.
    int yylen = 0;

    // Error handling.
    int yynerrs_ = 0;
    int yyerrstatus_ = 0;

    /// The lookahead symbol.
    symbol_type yyla;

    /// The locations where the error started and ended.
    stack_symbol_type yyerror_range[3];

    /// The return value of parse ().
    int yyresult;

#if YY_EXCEPTIONS
    try
#endif // YY_EXCEPTIONS
    {
        YYCDEBUG << "Starting parse\n";

        /* Initialize the stack.  The initial state will be set in
           yynewstate, since the latter expects the semantical and the
           location values to have been already stored, initialize these
           stacks with a primary value.  */
        yystack_.clear();
        yypush_(YY_NULLPTR, 0, YY_MOVE(yyla));

    /*-----------------------------------------------.
    | yynewstate -- push a new symbol on the stack.  |
    `-----------------------------------------------*/
    yynewstate:
        YYCDEBUG << "Entering state " << int(yystack_[0].state) << '\n';
        YY_STACK_PRINT();

        // Accept?
        if (yystack_[0].state == yyfinal_)
            YYACCEPT;

        goto yybackup;

    /*-----------.
    | yybackup.  |
    `-----------*/
    yybackup:
        // Try to take a decision without lookahead.
        yyn = yypact_[+yystack_[0].state];
        if (yy_pact_value_is_default_(yyn))
            goto yydefault;

        // Read a lookahead token.
        if (yyla.empty()) {
            YYCDEBUG << "Reading a token\n";
#if YY_EXCEPTIONS
            try
#endif // YY_EXCEPTIONS
            {
                yyla.kind_ = yytranslate_(yylex(&yyla.value, &yyla.location));
            }
#if YY_EXCEPTIONS
            catch (const syntax_error& yyexc) {
                YYCDEBUG << "Caught exception: " << yyexc.what() << '\n';
                error(yyexc);
                goto yyerrlab1;
            }
#endif // YY_EXCEPTIONS
        }
        YY_SYMBOL_PRINT("Next token is", yyla);

        if (yyla.kind() == symbol_kind::S_YYerror) {
            // The scanner already issued an error message, process directly
            // to error recovery.  But do not keep the error token as
            // lookahead, it is too special and may lead us to an endless
            // loop in error recovery. */
            yyla.kind_ = symbol_kind::S_YYUNDEF;
            goto yyerrlab1;
        }

        /* If the proper action on seeing token YYLA.TYPE is to reduce or
           to detect an error, take that action.  */
        yyn += yyla.kind();
        if (yyn < 0 || yylast_ < yyn || yycheck_[yyn] != yyla.kind()) {
            goto yydefault;
        }

        // Reduce or error.
        yyn = yytable_[yyn];
        if (yyn <= 0) {
            if (yy_table_value_is_error_(yyn))
                goto yyerrlab;
            yyn = -yyn;
            goto yyreduce;
        }

        // Count tokens shifted since error; after three, turn off error status.
        if (yyerrstatus_)
            --yyerrstatus_;

        // Shift the lookahead token.
        yypush_("Shifting", state_type(yyn), YY_MOVE(yyla));
        goto yynewstate;

    /*-----------------------------------------------------------.
    | yydefault -- do the default action for the current state.  |
    `-----------------------------------------------------------*/
    yydefault:
        yyn = yydefact_[+yystack_[0].state];
        if (yyn == 0)
            goto yyerrlab;
        goto yyreduce;

    /*-----------------------------.
    | yyreduce -- do a reduction.  |
    `-----------------------------*/
    yyreduce:
        yylen = yyr2_[yyn];
        {
            stack_symbol_type yylhs;
            yylhs.state = yy_lr_goto_state_(yystack_[yylen].state, yyr1_[yyn]);
            /* Variants are always initialized to an empty instance of the
               correct type. The default '$$ = $1' action is NOT applied
               when using variants.  */
            switch (yyr1_[yyn]) {
                case symbol_kind::S_connections:      // connections
                case symbol_kind::S_polymer_groups:   // polymer_groups
                case symbol_kind::S_monomer_sequence: // monomer_sequence
                case symbol_kind::
                    S_branch_monomer_group: // branch_monomer_group
                    yylhs.value.emplace<size_t>();
                    break;

                case symbol_kind::S_BLOB_ID:    // BLOB_ID
                case symbol_kind::S_CHEM_ID:    // CHEM_ID
                case symbol_kind::S_PEPTIDE_ID: // PEPTIDE_ID
                case symbol_kind::S_RNA_ID:     // RNA_ID
                case symbol_kind::S_POLYMER_ID: // POLYMER_ID
                case symbol_kind::
                    S_SINGLE_CHARACTER_MONOMER: // SINGLE_CHARACTER_MONOMER
                case symbol_kind::
                    S_MULTI_CHARACTER_MONOMER:        // MULTI_CHARACTER_MONOMER
                case symbol_kind::S_UNKNOWN_MONOMER:  // UNKNOWN_MONOMER
                case symbol_kind::S_MISSING_MONOMER:  // MISSING_MONOMER
                case symbol_kind::S_MONOMER_WILDCARD: // MONOMER_WILDCARD
                case symbol_kind::S_MONOMER_RATIO:    // MONOMER_RATIO
                case symbol_kind::
                    S_INLINE_SMILES_MONOMER:          // INLINE_SMILES_MONOMER
                case symbol_kind::S_ANNOTATION:       // ANNOTATION
                case symbol_kind::S_REPETITIONS:      // REPETITIONS
                case symbol_kind::S_UNKNOWN_SEQUENCE: // UNKNOWN_SEQUENCE
                case symbol_kind::S_VERSION_TOKEN:    // VERSION_TOKEN
                case symbol_kind::S_HYDROGEN_PAIRING: // HYDROGEN_PAIRING
                case symbol_kind::S_RGROUP:           // RGROUP
                case symbol_kind::
                    S_UNDEFINED_RESIDUE_NUMBER_OR_RGROUP: // UNDEFINED_RESIDUE_NUMBER_OR_RGROUP
                case symbol_kind::S_CONNECTION_RESIDUE:  // CONNECTION_RESIDUE
                case symbol_kind::S_POLYMER_GROUP_RATIO: // POLYMER_GROUP_RATIO
                case symbol_kind::S_POLYMER_GROUP_ID:    // POLYMER_GROUP_ID
                case symbol_kind::
                    S_EXTENDED_ANNOTATIONS_TOKEN: // EXTENDED_ANNOTATIONS_TOKEN
                case symbol_kind::
                    S_extended_annotations:             // extended_annotations
                case symbol_kind::S_version:            // version
                case symbol_kind::S_polymer:            // polymer
                case symbol_kind::S_repetitions:        // repetitions
                case symbol_kind::S_monomer_list:       // monomer_list
                case symbol_kind::S_monomer_and_list:   // monomer_and_list
                case symbol_kind::S_monomer_or_list:    // monomer_or_list
                case symbol_kind::S_monomer_list_item:  // monomer_list_item
                case symbol_kind::S_smiles_monomer:     // smiles_monomer
                case symbol_kind::S_monomer_id:         // monomer_id
                case symbol_kind::S_blob:               // blob
                case symbol_kind::S_connection_polymer: // connection_polymer
                case symbol_kind::S_attachment_point:   // attachment_point
                case symbol_kind::S_connection_monomer: // connection_monomer
                case symbol_kind::S_residue_and_list:   // residue_and_list
                case symbol_kind::S_residue_or_list:    // residue_or_list
                case symbol_kind::S_polymer_or_list:    // polymer_or_list
                case symbol_kind::S_polymer_and_list:   // polymer_and_list
                case symbol_kind::S_polymer_group_item: // polymer_group_item
                case symbol_kind::S_annotation:         // annotation
                    yylhs.value.emplace<std::string_view>();
                    break;

                default:
                    break;
            }

            // Default location.
            {
                stack_type::slice range(yystack_, yylen);
                YYLLOC_DEFAULT(yylhs.location, range, yylen);
                yyerror_range[1].location = yylhs.location;
            }

            // Perform the reduction.
            YY_REDUCE_PRINT(yyn);
#if YY_EXCEPTIONS
            try
#endif // YY_EXCEPTIONS
            {
                switch (yyn) {
                    case 2: // helm: polymers '$' connections '$' polymer_groups
                            // '$' extended_annotations '$' version
#line 107 "../helm_parser.yy"
                    {
                        if (YY_MOVE(yystack_[0].value.as<std::string_view>())
                                .empty() &&
                            !YY_MOVE(yystack_[2].value.as<std::string_view>())
                                 .empty()) {
                            helm_parser.saveError(
                                yystack_[2].location.begin.column,
                                "HELM annotations are currently unsupported");
                            YYABORT;
                        } else if (!YY_MOVE(yystack_[0]
                                                .value.as<std::string_view>())
                                        .empty() &&
                                   YY_MOVE(yystack_[0]
                                               .value.as<std::string_view>()) !=
                                       "V2.0") {
                            helm_parser.saveError(
                                yystack_[0].location.begin.column,
                                "Only HELM and HELMV2.0 versions are currently "
                                "supported");
                            YYABORT;
                        }
                    }
#line 1042 "helm_parser.tab.cpp"
                    break;

                    case 5: // connections: %empty
#line 130 "../helm_parser.yy"
                    {
                        yylhs.value.as<size_t>() = 0;
                    }
#line 1048 "helm_parser.tab.cpp"
                    break;

                    case 6: // connections: connection
#line 131 "../helm_parser.yy"
                    {
                        yylhs.value.as<size_t>() = 1;
                    }
#line 1054 "helm_parser.tab.cpp"
                    break;

                    case 7: // connections: connections '|' connection
#line 132 "../helm_parser.yy"
                    {
                        if (YY_MOVE(yystack_[2].value.as<size_t>()) == 0) {
                            TokenParser::error(yystack_[2].location,
                                               "syntax error");
                            YYABORT;
                        }
                        ++yylhs.value.as<size_t>();
                    }
#line 1063 "helm_parser.tab.cpp"
                    break;

                    case 8: // polymer_groups: %empty
#line 136 "../helm_parser.yy"
                    {
                        yylhs.value.as<size_t>() = 0;
                    }
#line 1069 "helm_parser.tab.cpp"
                    break;

                    case 9: // polymer_groups: polymer_group
#line 137 "../helm_parser.yy"
                    {
                        yylhs.value.as<size_t>() = 1;
                    }
#line 1075 "helm_parser.tab.cpp"
                    break;

                    case 10: // polymer_groups: polymer_groups '|' polymer_group
#line 138 "../helm_parser.yy"
                    {
                        if (YY_MOVE(yystack_[2].value.as<size_t>()) == 0) {
                            TokenParser::error(yystack_[2].location,
                                               "syntax error");
                            YYABORT;
                        }
                        ++yylhs.value.as<size_t>();
                    }
#line 1084 "helm_parser.tab.cpp"
                    break;

                    case 11: // extended_annotations: %empty
#line 142 "../helm_parser.yy"
                    {
                        yylhs.value.as<std::string_view>() = {};
                    }
#line 1090 "helm_parser.tab.cpp"
                    break;

                    case 12: // extended_annotations: EXTENDED_ANNOTATIONS_TOKEN
#line 143 "../helm_parser.yy"
                    {
                        helm_parser.add_extended_annotation(
                            YY_MOVE(yystack_[0].value.as<std::string_view>()));
                        yylhs.value.as<std::string_view>() =
                            YY_MOVE(yystack_[0].value.as<std::string_view>());
                    }
#line 1099 "helm_parser.tab.cpp"
                    break;

                    case 13: // version: %empty
#line 147 "../helm_parser.yy"
                    {
                        yylhs.value.as<std::string_view>() = {};
                    }
#line 1105 "helm_parser.tab.cpp"
                    break;

                    case 14: // version: VERSION_TOKEN
#line 147 "../helm_parser.yy"
                    {
                        yylhs.value.as<std::string_view>() =
                            YY_MOVE(yystack_[0].value.as<std::string_view>());
                    }
#line 1111 "helm_parser.tab.cpp"
                    break;

                    case 15: // polymer_unit: polymer
#line 170 "../helm_parser.yy"
                    {
                        helm_parser.add_polymer(
                            YY_MOVE(yystack_[0].value.as<std::string_view>()),
                            {});
                    }
#line 1117 "helm_parser.tab.cpp"
                    break;

                    case 16: // polymer_unit: polymer annotation
#line 171 "../helm_parser.yy"
                    {
                        helm_parser.add_polymer(
                            YY_MOVE(yystack_[1].value.as<std::string_view>()),
                            YY_MOVE(yystack_[0].value.as<std::string_view>()));
                    }
#line 1123 "helm_parser.tab.cpp"
                    break;

                    case 17: // polymer: PEPTIDE_ID '{' monomers '}'
#line 173 "../helm_parser.yy"
                    {
                        yylhs.value.as<std::string_view>() =
                            YY_MOVE(yystack_[3].value.as<std::string_view>());
                    }
#line 1129 "helm_parser.tab.cpp"
                    break;

                    case 18: // polymer: RNA_ID '{' monomers '}'
#line 174 "../helm_parser.yy"
                    {
                        yylhs.value.as<std::string_view>() =
                            YY_MOVE(yystack_[3].value.as<std::string_view>());
                    }
#line 1135 "helm_parser.tab.cpp"
                    break;

                    case 19: // polymer: CHEM_ID '{' monomer_unit '}'
#line 175 "../helm_parser.yy"
                    {
                        yylhs.value.as<std::string_view>() =
                            YY_MOVE(yystack_[3].value.as<std::string_view>());
                    }
#line 1141 "helm_parser.tab.cpp"
                    break;

                    case 20: // polymer: BLOB_ID '{' blob '}'
#line 176 "../helm_parser.yy"
                    {
                        yylhs.value.as<std::string_view>() =
                            YY_MOVE(yystack_[3].value.as<std::string_view>());
                    }
#line 1147 "helm_parser.tab.cpp"
                    break;

                    case 25: // monomer_group: branch_monomer_group
#line 183 "../helm_parser.yy"
                    {
                        helm_parser.mark_branch_monomer(
                            YY_MOVE(yystack_[0].value.as<size_t>()));
                    }
#line 1153 "helm_parser.tab.cpp"
                    break;

                    case 26: // monomer_group: monomer_unit '(' monomer_sequence
                             // ')'
#line 185 "../helm_parser.yy"
                    {
                        helm_parser.saveError(yystack_[2].location.begin.column,
                                              branch_monomer_group_err_message);
                        YYABORT;
                    }
#line 1163 "helm_parser.tab.cpp"
                    break;

                    case 27: // monomer_group: monomer_unit '(' monomer_sequence
                             // ')' monomer_unit
#line 190 "../helm_parser.yy"
                    {
                        helm_parser.saveError(yystack_[3].location.begin.column,
                                              branch_monomer_group_err_message);
                        YYABORT;
                    }
#line 1173 "helm_parser.tab.cpp"
                    break;

                    case 28: // repeated_monomers: monomer_unit repetitions
#line 196 "../helm_parser.yy"
                    {
                        helm_parser.mark_last_n_monomers_as_repeated(
                            1,
                            YY_MOVE(yystack_[0].value.as<std::string_view>()),
                            {});
                    }
#line 1181 "helm_parser.tab.cpp"
                    break;

                    case 29: // repeated_monomers: monomer_unit repetitions
                             // annotation
#line 199 "../helm_parser.yy"
                    {
                        helm_parser.mark_last_n_monomers_as_repeated(
                            1,
                            YY_MOVE(yystack_[1].value.as<std::string_view>()),
                            YY_MOVE(yystack_[0].value.as<std::string_view>()));
                    }
#line 1189 "helm_parser.tab.cpp"
                    break;

                    case 30: // repeated_monomers: '(' monomer_sequence ')'
                             // repetitions
#line 202 "../helm_parser.yy"
                    {
                        helm_parser.mark_last_n_monomers_as_repeated(
                            YY_MOVE(yystack_[2].value.as<size_t>()),
                            YY_MOVE(yystack_[0].value.as<std::string_view>()),
                            {});
                    }
#line 1197 "helm_parser.tab.cpp"
                    break;

                    case 31: // repeated_monomers: '(' monomer_sequence ')'
                             // repetitions annotation
#line 205 "../helm_parser.yy"
                    {
                        helm_parser.mark_last_n_monomers_as_repeated(
                            YY_MOVE(yystack_[3].value.as<size_t>()),
                            YY_MOVE(yystack_[1].value.as<std::string_view>()),
                            YY_MOVE(yystack_[0].value.as<std::string_view>()));
                    }
#line 1205 "helm_parser.tab.cpp"
                    break;

                    case 32: // repetitions: REPETITIONS
#line 209 "../helm_parser.yy"
                    {
                        yylhs.value.as<std::string_view>() =
                            YY_MOVE(yystack_[0].value.as<std::string_view>());
                        yylhs.value.as<std::string_view>().remove_prefix(1);
                        yylhs.value.as<std::string_view>().remove_suffix(1);
                    }
#line 1211 "helm_parser.tab.cpp"
                    break;

                    case 33: // monomer_sequence: monomer_unit '.' monomer_unit
#line 210 "../helm_parser.yy"
                    {
                        yylhs.value.as<size_t>() = 2;
                    }
#line 1217 "helm_parser.tab.cpp"
                    break;

                    case 34: // monomer_sequence: branch_monomer_group
#line 211 "../helm_parser.yy"
                    {
                        yylhs.value.as<size_t>() =
                            YY_MOVE(yystack_[0].value.as<size_t>());
                        helm_parser.mark_branch_monomer(
                            yylhs.value.as<size_t>());
                    }
#line 1223 "helm_parser.tab.cpp"
                    break;

                    case 35: // monomer_sequence: monomer_sequence '.'
                             // monomer_unit
#line 212 "../helm_parser.yy"
                    {
                        yylhs.value.as<size_t>() =
                            YY_MOVE(yystack_[2].value.as<size_t>()) + 1;
                    }
#line 1229 "helm_parser.tab.cpp"
                    break;

                    case 36: // monomer_sequence: monomer_sequence '.'
                             // branch_monomer_group
#line 213 "../helm_parser.yy"
                    {
                        helm_parser.mark_branch_monomer(
                            YY_MOVE(yystack_[0].value.as<size_t>()));
                        yylhs.value.as<size_t>() =
                            YY_MOVE(yystack_[2].value.as<size_t>()) +
                            YY_MOVE(yystack_[0].value.as<size_t>());
                    }
#line 1238 "helm_parser.tab.cpp"
                    break;

                    case 37: // branch_monomer_group: monomer_unit '('
                             // monomer_unit ')'
#line 218 "../helm_parser.yy"
                    {
                        yylhs.value.as<size_t>() = 2;
                    }
#line 1244 "helm_parser.tab.cpp"
                    break;

                    case 38: // branch_monomer_group: monomer_unit '('
                             // monomer_unit ')' monomer_unit
#line 219 "../helm_parser.yy"
                    {
                        yylhs.value.as<size_t>() = 3;
                    }
#line 1250 "helm_parser.tab.cpp"
                    break;

                    case 39: // monomer_unit: monomer_id
#line 220 "../helm_parser.yy"
                    {
                        helm_parser.add_monomer_with_id(
                            YY_MOVE(yystack_[0].value.as<std::string_view>()),
                            {});
                    }
#line 1256 "helm_parser.tab.cpp"
                    break;

                    case 40: // monomer_unit: monomer_id annotation
#line 221 "../helm_parser.yy"
                    {
                        helm_parser.add_monomer_with_id(
                            YY_MOVE(yystack_[1].value.as<std::string_view>()),
                            YY_MOVE(yystack_[0].value.as<std::string_view>()));
                    }
#line 1262 "helm_parser.tab.cpp"
                    break;

                    case 41: // monomer_unit: smiles_monomer
#line 222 "../helm_parser.yy"
                    {
                        helm_parser.add_smiles_monomer(
                            YY_MOVE(yystack_[0].value.as<std::string_view>()),
                            {});
                    }
#line 1268 "helm_parser.tab.cpp"
                    break;

                    case 42: // monomer_unit: smiles_monomer annotation
#line 223 "../helm_parser.yy"
                    {
                        helm_parser.add_smiles_monomer(
                            YY_MOVE(yystack_[1].value.as<std::string_view>()),
                            YY_MOVE(yystack_[0].value.as<std::string_view>()));
                    }
#line 1274 "helm_parser.tab.cpp"
                    break;

                    case 43: // monomer_unit: '(' monomer_list ')'
#line 224 "../helm_parser.yy"
                    {
                        helm_parser.add_monomer_list(
                            YY_MOVE(yystack_[1].value.as<std::string_view>()),
                            {});
                    }
#line 1280 "helm_parser.tab.cpp"
                    break;

                    case 44: // monomer_unit: '(' monomer_list ')' annotation
#line 225 "../helm_parser.yy"
                    {
                        helm_parser.add_monomer_list(
                            YY_MOVE(yystack_[2].value.as<std::string_view>()),
                            YY_MOVE(yystack_[0].value.as<std::string_view>()));
                    }
#line 1286 "helm_parser.tab.cpp"
                    break;

                    case 45: // monomer_list: monomer_and_list
#line 227 "../helm_parser.yy"
                    {
                        yylhs.value.as<std::string_view>() =
                            YY_MOVE(yystack_[0].value.as<std::string_view>());
                    }
#line 1292 "helm_parser.tab.cpp"
                    break;

                    case 46: // monomer_list: monomer_or_list
#line 227 "../helm_parser.yy"
                    {
                        yylhs.value.as<std::string_view>() =
                            YY_MOVE(yystack_[0].value.as<std::string_view>());
                    }
#line 1298 "helm_parser.tab.cpp"
                    break;

                    case 47: // monomer_and_list: monomer_list_item '+'
                             // monomer_list_item
#line 228 "../helm_parser.yy"
                    {
                        yylhs.value.as<std::string_view>() = {
                            YY_MOVE(yystack_[2].value.as<std::string_view>())
                                .data(),
                            YY_MOVE(yystack_[2].value.as<std::string_view>())
                                    .size() +
                                1 +
                                YY_MOVE(
                                    yystack_[0].value.as<std::string_view>())
                                    .size()};
                    }
#line 1306 "helm_parser.tab.cpp"
                    break;

                    case 48: // monomer_and_list: monomer_and_list '+'
                             // monomer_list_item
#line 231 "../helm_parser.yy"
                    {
                        yylhs.value.as<std::string_view>() = {
                            YY_MOVE(yystack_[2].value.as<std::string_view>())
                                .data(),
                            YY_MOVE(yystack_[2].value.as<std::string_view>())
                                    .size() +
                                1 +
                                YY_MOVE(
                                    yystack_[0].value.as<std::string_view>())
                                    .size()};
                    }
#line 1314 "helm_parser.tab.cpp"
                    break;

                    case 49: // monomer_or_list: monomer_list_item ','
                             // monomer_list_item
#line 235 "../helm_parser.yy"
                    {
                        yylhs.value.as<std::string_view>() = {
                            YY_MOVE(yystack_[2].value.as<std::string_view>())
                                .data(),
                            YY_MOVE(yystack_[2].value.as<std::string_view>())
                                    .size() +
                                1 +
                                YY_MOVE(
                                    yystack_[0].value.as<std::string_view>())
                                    .size()};
                    }
#line 1322 "helm_parser.tab.cpp"
                    break;

                    case 50: // monomer_or_list: monomer_or_list ','
                             // monomer_list_item
#line 238 "../helm_parser.yy"
                    {
                        yylhs.value.as<std::string_view>() = {
                            YY_MOVE(yystack_[2].value.as<std::string_view>())
                                .data(),
                            YY_MOVE(yystack_[2].value.as<std::string_view>())
                                    .size() +
                                1 +
                                YY_MOVE(
                                    yystack_[0].value.as<std::string_view>())
                                    .size()};
                    }
#line 1330 "helm_parser.tab.cpp"
                    break;

                    case 51: // monomer_list_item: monomer_id
#line 243 "../helm_parser.yy"
                    {
                        yylhs.value.as<std::string_view>() =
                            YY_MOVE(yystack_[0].value.as<std::string_view>());
                    }
#line 1336 "helm_parser.tab.cpp"
                    break;

                    case 52: // monomer_list_item: monomer_id ':' MONOMER_RATIO
#line 244 "../helm_parser.yy"
                    {
                        yylhs.value.as<std::string_view>() = {
                            YY_MOVE(yystack_[2].value.as<std::string_view>())
                                .data(),
                            YY_MOVE(yystack_[2].value.as<std::string_view>())
                                    .size() +
                                YY_MOVE(
                                    yystack_[0].value.as<std::string_view>())
                                    .size() +
                                1};
                    }
#line 1342 "helm_parser.tab.cpp"
                    break;

                    case 53: // smiles_monomer: INLINE_SMILES_MONOMER
#line 247 "../helm_parser.yy"
                    {
                        yylhs.value.as<std::string_view>() =
                            YY_MOVE(yystack_[0].value.as<std::string_view>());
                        yylhs.value.as<std::string_view>().remove_prefix(1);
                        yylhs.value.as<std::string_view>().remove_suffix(1);
                    }
#line 1348 "helm_parser.tab.cpp"
                    break;

                    case 54: // monomer_id: SINGLE_CHARACTER_MONOMER
#line 249 "../helm_parser.yy"
                    {
                        helm_parser.add_residue_name(
                            YY_MOVE(yystack_[0].value.as<std::string_view>()));
                        yylhs.value.as<std::string_view>() =
                            YY_MOVE(yystack_[0].value.as<std::string_view>());
                    }
#line 1354 "helm_parser.tab.cpp"
                    break;

                    case 55: // monomer_id: MULTI_CHARACTER_MONOMER
#line 250 "../helm_parser.yy"
                    {
                        auto tmp = std::string_view{
                            YY_MOVE(yystack_[0].value.as<std::string_view>())};
                        tmp.remove_prefix(1);
                        tmp.remove_suffix(1);
                        helm_parser.add_residue_name(tmp);
                        yylhs.value.as<std::string_view>() =
                            YY_MOVE(yystack_[0].value.as<std::string_view>());
                    }
#line 1366 "helm_parser.tab.cpp"
                    break;

                    case 56: // monomer_id: UNKNOWN_MONOMER
#line 257 "../helm_parser.yy"
                    {
                        helm_parser.add_wildcard_or_unknown_residue();
                        yylhs.value.as<std::string_view>() =
                            YY_MOVE(yystack_[0].value.as<std::string_view>());
                    }
#line 1375 "helm_parser.tab.cpp"
                    break;

                    case 57: // monomer_id: MISSING_MONOMER
#line 261 "../helm_parser.yy"
                    {
                        yylhs.value.as<std::string_view>() =
                            YY_MOVE(yystack_[0].value.as<std::string_view>());
                    }
#line 1381 "helm_parser.tab.cpp"
                    break;

                    case 58: // monomer_id: MONOMER_WILDCARD
#line 262 "../helm_parser.yy"
                    {
                        helm_parser.add_wildcard_or_unknown_residue();
                        yylhs.value.as<std::string_view>() =
                            YY_MOVE(yystack_[0].value.as<std::string_view>());
                    }
#line 1390 "helm_parser.tab.cpp"
                    break;

                    case 59: // blob: UNKNOWN_SEQUENCE
#line 267 "../helm_parser.yy"
                    {
                        helm_parser.add_monomer(
                            YY_MOVE(yystack_[0].value.as<std::string_view>()),
                            false, false, false, {});
                    }
#line 1396 "helm_parser.tab.cpp"
                    break;

                    case 60: // connection: connection_polymer ','
                             // connection_polymer ',' connection_monomer ':'
                             // attachment_point '-' connection_monomer ':'
                             // attachment_point
#line 284 "../helm_parser.yy"
                    {
                        helm_parser.add_connection({
                            YY_MOVE(
                                yystack_[10]
                                    .value.as<std::string_view>()), // from_id
                            YY_MOVE(
                                yystack_[6]
                                    .value.as<std::string_view>()), // from_res
                            YY_MOVE(yystack_[4]
                                        .value
                                        .as<std::string_view>()), // from_rgroup
                            YY_MOVE(yystack_[8]
                                        .value.as<std::string_view>()), // to_id
                            YY_MOVE(
                                yystack_[2]
                                    .value.as<std::string_view>()), // to_res
                            YY_MOVE(
                                yystack_[0]
                                    .value.as<std::string_view>()), // to_rgroup
                            {} // annotation
                        });
                    }
#line 1412 "helm_parser.tab.cpp"
                    break;

                    case 61: // connection: connection_polymer ','
                             // connection_polymer ',' connection_monomer ':'
                             // attachment_point '-' connection_monomer ':'
                             // attachment_point annotation
#line 297 "../helm_parser.yy"
                    {
                        helm_parser.add_connection({
                            YY_MOVE(
                                yystack_[11]
                                    .value.as<std::string_view>()), // from_id
                            YY_MOVE(
                                yystack_[7]
                                    .value.as<std::string_view>()), // from_res
                            YY_MOVE(yystack_[5]
                                        .value
                                        .as<std::string_view>()), // from_rgroup
                            YY_MOVE(yystack_[9]
                                        .value.as<std::string_view>()), // to_id
                            YY_MOVE(
                                yystack_[3]
                                    .value.as<std::string_view>()), // to_res
                            YY_MOVE(
                                yystack_[1]
                                    .value.as<std::string_view>()), // to_rgroup
                            YY_MOVE(
                                yystack_[0]
                                    .value.as<std::string_view>()) // annotation
                        });
                    }
#line 1428 "helm_parser.tab.cpp"
                    break;

                    case 62: // connection_polymer: POLYMER_ID
#line 309 "../helm_parser.yy"
                    {
                        yylhs.value.as<std::string_view>() =
                            YY_MOVE(yystack_[0].value.as<std::string_view>());
                    }
#line 1434 "helm_parser.tab.cpp"
                    break;

                    case 63: // connection_polymer: POLYMER_GROUP_ID
#line 309 "../helm_parser.yy"
                    {
                        yylhs.value.as<std::string_view>() =
                            YY_MOVE(yystack_[0].value.as<std::string_view>());
                    }
#line 1440 "helm_parser.tab.cpp"
                    break;

                    case 64: // attachment_point: HYDROGEN_PAIRING
#line 310 "../helm_parser.yy"
                    {
                        yylhs.value.as<std::string_view>() =
                            YY_MOVE(yystack_[0].value.as<std::string_view>());
                    }
#line 1446 "helm_parser.tab.cpp"
                    break;

                    case 65: // attachment_point: RGROUP
#line 310 "../helm_parser.yy"
                    {
                        yylhs.value.as<std::string_view>() =
                            YY_MOVE(yystack_[0].value.as<std::string_view>());
                    }
#line 1452 "helm_parser.tab.cpp"
                    break;

                    case 66: // attachment_point:
                             // UNDEFINED_RESIDUE_NUMBER_OR_RGROUP
#line 310 "../helm_parser.yy"
                    {
                        yylhs.value.as<std::string_view>() =
                            YY_MOVE(yystack_[0].value.as<std::string_view>());
                    }
#line 1458 "helm_parser.tab.cpp"
                    break;

                    case 67: // connection_monomer: CONNECTION_RESIDUE
#line 311 "../helm_parser.yy"
                    {
                        yylhs.value.as<std::string_view>() =
                            YY_MOVE(yystack_[0].value.as<std::string_view>());
                    }
#line 1464 "helm_parser.tab.cpp"
                    break;

                    case 68: // connection_monomer:
                             // UNDEFINED_RESIDUE_NUMBER_OR_RGROUP
#line 312 "../helm_parser.yy"
                    {
                        yylhs.value.as<std::string_view>() =
                            YY_MOVE(yystack_[0].value.as<std::string_view>());
                    }
#line 1470 "helm_parser.tab.cpp"
                    break;

                    case 69: // connection_monomer: '(' residue_and_list ')'
#line 313 "../helm_parser.yy"
                    {
                        yylhs.value.as<std::string_view>() =
                            YY_MOVE(yystack_[1].value.as<std::string_view>());
                    }
#line 1476 "helm_parser.tab.cpp"
                    break;

                    case 70: // connection_monomer: '(' residue_or_list ')'
#line 314 "../helm_parser.yy"
                    {
                        yylhs.value.as<std::string_view>() =
                            YY_MOVE(yystack_[1].value.as<std::string_view>());
                    }
#line 1482 "helm_parser.tab.cpp"
                    break;

                    case 71: // residue_and_list: CONNECTION_RESIDUE '+'
                             // CONNECTION_RESIDUE
#line 316 "../helm_parser.yy"
                    {
                        yylhs.value.as<std::string_view>() = {
                            YY_MOVE(yystack_[2].value.as<std::string_view>())
                                .data(),
                            YY_MOVE(yystack_[2].value.as<std::string_view>())
                                    .size() +
                                1 +
                                YY_MOVE(
                                    yystack_[0].value.as<std::string_view>())
                                    .size()};
                    }
#line 1490 "helm_parser.tab.cpp"
                    break;

                    case 72: // residue_and_list: residue_and_list '+'
                             // CONNECTION_RESIDUE
#line 319 "../helm_parser.yy"
                    {
                        yylhs.value.as<std::string_view>() = {
                            YY_MOVE(yystack_[2].value.as<std::string_view>())
                                .data(),
                            YY_MOVE(yystack_[2].value.as<std::string_view>())
                                    .size() +
                                1 +
                                YY_MOVE(
                                    yystack_[0].value.as<std::string_view>())
                                    .size()};
                    }
#line 1498 "helm_parser.tab.cpp"
                    break;

                    case 73: // residue_or_list: CONNECTION_RESIDUE ','
                             // CONNECTION_RESIDUE
#line 323 "../helm_parser.yy"
                    {
                        yylhs.value.as<std::string_view>() = {
                            YY_MOVE(yystack_[2].value.as<std::string_view>())
                                .data(),
                            YY_MOVE(yystack_[2].value.as<std::string_view>())
                                    .size() +
                                1 +
                                YY_MOVE(
                                    yystack_[0].value.as<std::string_view>())
                                    .size()};
                    }
#line 1506 "helm_parser.tab.cpp"
                    break;

                    case 74: // residue_or_list: residue_or_list ','
                             // CONNECTION_RESIDUE
#line 326 "../helm_parser.yy"
                    {
                        yylhs.value.as<std::string_view>() = {
                            YY_MOVE(yystack_[2].value.as<std::string_view>())
                                .data(),
                            YY_MOVE(yystack_[2].value.as<std::string_view>())
                                    .size() +
                                1 +
                                YY_MOVE(
                                    yystack_[0].value.as<std::string_view>())
                                    .size()};
                    }
#line 1514 "helm_parser.tab.cpp"
                    break;

                    case 75: // polymer_group: POLYMER_GROUP_ID '('
                             // polymer_and_list ')'
#line 332 "../helm_parser.yy"
                    {
                        helm_parser.add_polymer_group(
                            YY_MOVE(yystack_[3].value.as<std::string_view>()),
                            YY_MOVE(yystack_[1].value.as<std::string_view>()),
                            true);
                    }
#line 1520 "helm_parser.tab.cpp"
                    break;

                    case 76: // polymer_group: POLYMER_GROUP_ID '('
                             // polymer_or_list ')'
#line 333 "../helm_parser.yy"
                    {
                        helm_parser.add_polymer_group(
                            YY_MOVE(yystack_[3].value.as<std::string_view>()),
                            YY_MOVE(yystack_[1].value.as<std::string_view>()),
                            false);
                    }
#line 1526 "helm_parser.tab.cpp"
                    break;

                    case 77: // polymer_or_list: polymer_group_item ','
                             // polymer_group_item
#line 335 "../helm_parser.yy"
                    {
                        yylhs.value.as<std::string_view>() = {
                            YY_MOVE(yystack_[2].value.as<std::string_view>())
                                .data(),
                            YY_MOVE(yystack_[2].value.as<std::string_view>())
                                    .size() +
                                1 +
                                YY_MOVE(
                                    yystack_[0].value.as<std::string_view>())
                                    .size()};
                    }
#line 1534 "helm_parser.tab.cpp"
                    break;

                    case 78: // polymer_or_list: polymer_or_list ','
                             // polymer_group_item
#line 338 "../helm_parser.yy"
                    {
                        yylhs.value.as<std::string_view>() = {
                            YY_MOVE(yystack_[2].value.as<std::string_view>())
                                .data(),
                            YY_MOVE(yystack_[2].value.as<std::string_view>())
                                    .size() +
                                1 +
                                YY_MOVE(
                                    yystack_[0].value.as<std::string_view>())
                                    .size()};
                    }
#line 1542 "helm_parser.tab.cpp"
                    break;

                    case 79: // polymer_and_list: polymer_group_item '+'
                             // polymer_group_item
#line 341 "../helm_parser.yy"
                    {
                        yylhs.value.as<std::string_view>() = {
                            YY_MOVE(yystack_[2].value.as<std::string_view>())
                                .data(),
                            YY_MOVE(yystack_[2].value.as<std::string_view>())
                                    .size() +
                                1 +
                                YY_MOVE(
                                    yystack_[0].value.as<std::string_view>())
                                    .size()};
                    }
#line 1550 "helm_parser.tab.cpp"
                    break;

                    case 80: // polymer_and_list: polymer_and_list '+'
                             // polymer_group_item
#line 344 "../helm_parser.yy"
                    {
                        yylhs.value.as<std::string_view>() = {
                            YY_MOVE(yystack_[2].value.as<std::string_view>())
                                .data(),
                            YY_MOVE(yystack_[2].value.as<std::string_view>())
                                    .size() +
                                1 +
                                YY_MOVE(
                                    yystack_[0].value.as<std::string_view>())
                                    .size()};
                    }
#line 1558 "helm_parser.tab.cpp"
                    break;

                    case 81: // polymer_group_item: POLYMER_ID
#line 348 "../helm_parser.yy"
                    {
                        yylhs.value.as<std::string_view>() =
                            YY_MOVE(yystack_[0].value.as<std::string_view>());
                    }
#line 1564 "helm_parser.tab.cpp"
                    break;

                    case 82: // polymer_group_item: POLYMER_GROUP_ID
#line 349 "../helm_parser.yy"
                    {
                        yylhs.value.as<std::string_view>() =
                            YY_MOVE(yystack_[0].value.as<std::string_view>());
                    }
#line 1570 "helm_parser.tab.cpp"
                    break;

                    case 83: // polymer_group_item: POLYMER_ID ':'
                             // POLYMER_GROUP_RATIO
#line 350 "../helm_parser.yy"
                    {
                        yylhs.value.as<std::string_view>() = {
                            YY_MOVE(yystack_[2].value.as<std::string_view>())
                                .data(),
                            YY_MOVE(yystack_[2].value.as<std::string_view>())
                                    .size() +
                                1 +
                                YY_MOVE(
                                    yystack_[0].value.as<std::string_view>())
                                    .size()};
                    }
#line 1578 "helm_parser.tab.cpp"
                    break;

                    case 84: // polymer_group_item: POLYMER_GROUP_ID ':'
                             // POLYMER_GROUP_RATIO
#line 353 "../helm_parser.yy"
                    {
                        yylhs.value.as<std::string_view>() = {
                            YY_MOVE(yystack_[2].value.as<std::string_view>())
                                .data(),
                            YY_MOVE(yystack_[2].value.as<std::string_view>())
                                    .size() +
                                1 +
                                YY_MOVE(
                                    yystack_[0].value.as<std::string_view>())
                                    .size()};
                    }
#line 1586 "helm_parser.tab.cpp"
                    break;

                    case 85: // annotation: ANNOTATION
#line 359 "../helm_parser.yy"
                    {
                        yylhs.value.as<std::string_view>() =
                            YY_MOVE(yystack_[0].value.as<std::string_view>());
                        yylhs.value.as<std::string_view>().remove_prefix(1);
                        yylhs.value.as<std::string_view>().remove_suffix(1);
                    }
#line 1592 "helm_parser.tab.cpp"
                    break;

#line 1596 "helm_parser.tab.cpp"

                    default:
                        break;
                }
            }
#if YY_EXCEPTIONS
            catch (const syntax_error& yyexc) {
                YYCDEBUG << "Caught exception: " << yyexc.what() << '\n';
                error(yyexc);
                YYERROR;
            }
#endif // YY_EXCEPTIONS
            YY_SYMBOL_PRINT("-> $$ =", yylhs);
            yypop_(yylen);
            yylen = 0;

            // Shift the result of the reduction.
            yypush_(YY_NULLPTR, YY_MOVE(yylhs));
        }
        goto yynewstate;

    /*--------------------------------------.
    | yyerrlab -- here on detecting error.  |
    `--------------------------------------*/
    yyerrlab:
        // If not already recovering from an error, report this error.
        if (!yyerrstatus_) {
            ++yynerrs_;
            std::string msg = YY_("syntax error");
            error(yyla.location, YY_MOVE(msg));
        }

        yyerror_range[1].location = yyla.location;
        if (yyerrstatus_ == 3) {
            /* If just tried and failed to reuse lookahead token after an
               error, discard it.  */

            // Return failure if at end of input.
            if (yyla.kind() == symbol_kind::S_YYEOF)
                YYABORT;
            else if (!yyla.empty()) {
                yy_destroy_("Error: discarding", yyla);
                yyla.clear();
            }
        }

        // Else will try to reuse lookahead token after shifting the error
        // token.
        goto yyerrlab1;

    /*---------------------------------------------------.
    | yyerrorlab -- error raised explicitly by YYERROR.  |
    `---------------------------------------------------*/
    yyerrorlab:
        /* Pacify compilers when the user code never invokes YYERROR and
           the label yyerrorlab therefore never appears in user code.  */
        if (false)
            YYERROR;

        /* Do not reclaim the symbols of the rule whose action triggered
           this YYERROR.  */
        yypop_(yylen);
        yylen = 0;
        YY_STACK_PRINT();
        goto yyerrlab1;

    /*-------------------------------------------------------------.
    | yyerrlab1 -- common code for both syntax error and YYERROR.  |
    `-------------------------------------------------------------*/
    yyerrlab1:
        yyerrstatus_ = 3; // Each real token shifted decrements this.
        // Pop stack until we find a state that shifts the error token.
        for (;;) {
            yyn = yypact_[+yystack_[0].state];
            if (!yy_pact_value_is_default_(yyn)) {
                yyn += symbol_kind::S_YYerror;
                if (0 <= yyn && yyn <= yylast_ &&
                    yycheck_[yyn] == symbol_kind::S_YYerror) {
                    yyn = yytable_[yyn];
                    if (0 < yyn)
                        break;
                }
            }

            // Pop the current state because it cannot handle the error token.
            if (yystack_.size() == 1)
                YYABORT;

            yyerror_range[1].location = yystack_[0].location;
            yy_destroy_("Error: popping", yystack_[0]);
            yypop_();
            YY_STACK_PRINT();
        }
        {
            stack_symbol_type error_token;

            yyerror_range[2].location = yyla.location;
            YYLLOC_DEFAULT(error_token.location, yyerror_range, 2);

            // Shift the error token.
            error_token.state = state_type(yyn);
            yypush_("Shifting", YY_MOVE(error_token));
        }
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

    /*-----------------------------------------------------.
    | yyreturn -- parsing is finished, return the result.  |
    `-----------------------------------------------------*/
    yyreturn:
        if (!yyla.empty())
            yy_destroy_("Cleanup: discarding lookahead", yyla);

        /* Do not reclaim the symbols of the rule whose action triggered
           this YYABORT or YYACCEPT.  */
        yypop_(yylen);
        YY_STACK_PRINT();
        while (1 < yystack_.size()) {
            yy_destroy_("Cleanup: popping", yystack_[0]);
            yypop_();
        }

        return yyresult;
    }
#if YY_EXCEPTIONS
    catch (...) {
        YYCDEBUG << "Exception caught: cleaning lookahead and stack\n";
        // Do not try to display the values of the reclaimed symbols,
        // as their printers might throw an exception.
        if (!yyla.empty())
            yy_destroy_(YY_NULLPTR, yyla);

        while (1 < yystack_.size()) {
            yy_destroy_(YY_NULLPTR, yystack_[0]);
            yypop_();
        }
        throw;
    }
#endif // YY_EXCEPTIONS
}

void TokenParser::error(const syntax_error& yyexc)
{
    error(yyexc.location, yyexc.what());
}

#if YYDEBUG || 0
const char* TokenParser::symbol_name(symbol_kind_type yysymbol)
{
    return yytname_[yysymbol];
}
#endif // #if YYDEBUG || 0

const signed char TokenParser::yypact_ninf_ = -69;

const signed char TokenParser::yytable_ninf_ = -52;

const signed char TokenParser::yypact_[] = {
    55,  21,  28,  35,  57,  92,  56,  -69, 93,  90,  -1,  17,  17,  -69, 8,
    55,  -69, -69, -69, 80,  -69, -69, -69, -69, -69, -69, 58,  81,  93,  93,
    -1,  6,   -69, -69, -69, -13, 59,  -69, -69, 64,  -69, 77,  -69, -69, 82,
    79,  83,  60,  78,  -69, -69, -69, -18, -69, 65,  -11, -69, 17,  -69, -1,
    93,  -69, 91,  8,   8,   93,  58,  58,  58,  58,  103, -1,  102, -1,  -1,
    -69, 48,  45,  -69, 88,  71,  -69, -69, 86,  -69, -69, -69, -69, -69, -69,
    -69, 94,  93,  -69, 89,  -1,  -1,  9,   97,  91,  16,  -69, -69, -69, 95,
    96,  -15, 67,  68,  -69, 98,  -69, -69, -69, 101, 99,  104, 105, -69, 9,
    -69, 9,   9,   9,   108, 70,  73,  47,  20,  -69, -69, -69, -69, -69, -69,
    -69, -69, 107, 110, -69, 111, -69, 113, -69, -69, -69, 100, -69, -69, -69,
    -69, 16,  106, 20,  93,  -69};

const signed char TokenParser::yydefact_[] = {
    0,  0,  0,  0,  0,  0,  0,  3,  15, 0,  0,  0,  0,  1,  5,  0,  85, 16,
    59, 0,  54, 55, 56, 57, 58, 53, 0,  0,  41, 39, 0,  0,  21, 24, 25, 23,
    0,  62, 63, 0,  6,  0,  4,  20, 0,  45, 46, 0,  51, 19, 42, 40, 0,  34,
    0,  39, 17, 0,  32, 0,  28, 18, 8,  0,  0,  43, 0,  0,  0,  0,  0,  0,
    0,  0,  0,  22, 0,  0,  29, 0,  0,  9,  7,  0,  44, 48, 50, 47, 49, 52,
    36, 35, 30, 33, 0,  26, 37, 0,  11, 0,  0,  31, 27, 38, 81, 82, 0,  0,
    0,  12, 0,  10, 68, 67, 0,  0,  0,  0,  76, 0,  75, 0,  0,  0,  13, 0,
    0,  0,  0,  83, 84, 78, 80, 79, 77, 14, 2,  0,  0,  69, 0,  70, 0,  64,
    65, 66, 0,  71, 73, 72, 74, 0,  0,  0,  60, 61};

const signed char TokenParser::yypgoto_[] = {
    -69, -69, -69, -69, -69, -69, -69, 122, -69, 126, 85,  -69,
    72,  84,  -25, -9,  -69, -69, -69, 5,   -69, -24, -69, 76,
    87,  -6,  -3,  -69, -69, 41,  -69, -69, -68, -8};

const unsigned char TokenParser::yydefgoto_[] = {
    0,  5,   6,   39,  80,  110, 136, 7,   8,   31, 32, 33,
    60, 52,  34,  35,  44,  45,  46,  47,  28,  29, 19, 40,
    41, 146, 115, 126, 127, 81,  106, 107, 108, 51};

const short TokenParser::yytable_[] = {
    17,  27,  48,  58,  16,  53,  55,  20,  21,  22,  23,  24,  71,  25,
    72,  37,  104, 118, 59,  119, 50,  54,  -51, -51, 70,  20,  21,  22,
    23,  24,  26,  25,  38,  105, 53,  56,  57,  112, 113, 143, 144, 145,
    48,  48,  48,  48,  90,  114, 30,  9,   77,  131, 78,  132, 133, 134,
    10,  84,  1,   2,   3,   4,   91,  11,  93,  94,  20,  21,  22,  23,
    24,  85,  86,  87,  88,  73,  74,  96,  71,  141, 95,  142, 14,  15,
    101, 12,  102, 103, 61,  57,  62,  63,  13,  68,  69,  73,  74,  98,
    99,  120, 121, 122, 123, 137, 138, 139, 140, 18,  16,  43,  49,  64,
    66,  70,  65,  79,  89,  67,  58,  97,  100, 96,  109, 125, 124, 74,
    135, 129, 130, 147, 116, 117, 148, 149, 128, 150, 151, 42,  36,  82,
    111, 153, 75,  76,  92,  0,   155, 154, 152, 0,   0,   83};

const short TokenParser::yycheck_[] = {
    8,   10, 26,  16,  15,  30, 30, 8,  9,  10, 11,  12,  30,  14, 32, 7,  7,
    32,  31, 34,  28,  30,  33, 34, 35, 8,  9,  10,  11,  12,  31, 14, 24, 24,
    59,  29, 30,  21,  22,  19, 20, 21, 66, 67, 68,  69,  71,  31, 31, 28, 59,
    119, 60, 121, 122, 123, 28, 65, 3,  4,  5,  6,   71,  28,  73, 74, 8,  9,
    10,  11, 12,  66,  67,  68, 69, 30, 31, 32, 30,  32,  32,  34, 26, 27, 92,
    28,  95, 96,  29,  30,  26, 27, 0,  33, 34, 30,  31,  26,  27, 32, 33, 33,
    34,  33, 34,  32,  33,  17, 15, 29, 29, 34, 33,  35,  32,  24, 13, 34, 16,
    31,  34, 32,  25,  22,  26, 31, 18, 23, 23, 22,  35,  35,  22, 22, 35, 22,
    36,  15, 12,  63,  99,  35, 57, 59, 72, -1, 154, 153, 151, -1, -1, 64};

const signed char TokenParser::yystos_[] = {
    0,  3,  4,  5,  6,  38, 39, 44, 45, 28, 28, 28, 28, 0,  26, 27, 15, 70,
    17, 59, 8,  9,  10, 11, 12, 14, 31, 52, 57, 58, 31, 46, 47, 48, 51, 52,
    46, 7,  24, 40, 60, 61, 44, 29, 53, 54, 55, 56, 58, 29, 70, 70, 50, 51,
    52, 58, 29, 30, 16, 31, 49, 29, 26, 27, 34, 32, 33, 34, 33, 34, 35, 30,
    32, 30, 31, 47, 50, 52, 70, 24, 41, 66, 60, 61, 70, 56, 56, 56, 56, 13,
    51, 52, 49, 52, 52, 32, 32, 31, 26, 27, 34, 70, 52, 52, 7,  24, 67, 68,
    69, 25, 42, 66, 21, 22, 31, 63, 35, 35, 32, 34, 32, 33, 33, 34, 26, 22,
    64, 65, 35, 23, 23, 69, 69, 69, 69, 18, 43, 33, 34, 32, 33, 32, 34, 19,
    20, 21, 62, 22, 22, 22, 22, 36, 63, 35, 62, 70};

const signed char TokenParser::yyr1_[] = {
    0,  37, 38, 39, 39, 40, 40, 40, 41, 41, 41, 42, 42, 43, 43, 44, 44, 45,
    45, 45, 45, 46, 46, 47, 47, 47, 47, 47, 48, 48, 48, 48, 49, 50, 50, 50,
    50, 51, 51, 52, 52, 52, 52, 52, 52, 53, 53, 54, 54, 55, 55, 56, 56, 57,
    58, 58, 58, 58, 58, 59, 60, 60, 61, 61, 62, 62, 62, 63, 63, 63, 63, 64,
    64, 65, 65, 66, 66, 67, 67, 68, 68, 69, 69, 69, 69, 70};

const signed char TokenParser::yyr2_[] = {
    0, 2, 9, 1, 3, 0, 1, 3, 0, 1, 3, 0, 1, 0, 1, 1, 2,  4,  4, 4, 4, 1,
    3, 1, 1, 1, 4, 5, 2, 3, 4, 5, 1, 3, 1, 3, 3, 4, 5,  1,  2, 1, 2, 3,
    4, 1, 1, 3, 3, 3, 3, 1, 3, 1, 1, 1, 1, 1, 1, 1, 11, 12, 1, 1, 1, 1,
    1, 1, 1, 3, 3, 3, 3, 3, 3, 4, 4, 3, 3, 3, 3, 1, 1,  3,  3, 1};

#if YYDEBUG
// YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
// First, the terminals, then, starting at \a YYNTOKENS, nonterminals.
const char* const TokenParser::yytname_[] = {
    "END",
    "error",
    "\"invalid token\"",
    "BLOB_ID",
    "CHEM_ID",
    "PEPTIDE_ID",
    "RNA_ID",
    "POLYMER_ID",
    "SINGLE_CHARACTER_MONOMER",
    "MULTI_CHARACTER_MONOMER",
    "UNKNOWN_MONOMER",
    "MISSING_MONOMER",
    "MONOMER_WILDCARD",
    "MONOMER_RATIO",
    "INLINE_SMILES_MONOMER",
    "ANNOTATION",
    "REPETITIONS",
    "UNKNOWN_SEQUENCE",
    "VERSION_TOKEN",
    "HYDROGEN_PAIRING",
    "RGROUP",
    "UNDEFINED_RESIDUE_NUMBER_OR_RGROUP",
    "CONNECTION_RESIDUE",
    "POLYMER_GROUP_RATIO",
    "POLYMER_GROUP_ID",
    "EXTENDED_ANNOTATIONS_TOKEN",
    "'$'",
    "'|'",
    "'{'",
    "'}'",
    "'.'",
    "'('",
    "')'",
    "'+'",
    "','",
    "':'",
    "'-'",
    "$accept",
    "helm",
    "polymers",
    "connections",
    "polymer_groups",
    "extended_annotations",
    "version",
    "polymer_unit",
    "polymer",
    "monomers",
    "monomer_group",
    "repeated_monomers",
    "repetitions",
    "monomer_sequence",
    "branch_monomer_group",
    "monomer_unit",
    "monomer_list",
    "monomer_and_list",
    "monomer_or_list",
    "monomer_list_item",
    "smiles_monomer",
    "monomer_id",
    "blob",
    "connection",
    "connection_polymer",
    "attachment_point",
    "connection_monomer",
    "residue_and_list",
    "residue_or_list",
    "polymer_group",
    "polymer_or_list",
    "polymer_and_list",
    "polymer_group_item",
    "annotation",
    YY_NULLPTR};
#endif

#if YYDEBUG
const short TokenParser::yyrline_[] = {
    0,   107, 107, 127, 128, 130, 131, 132, 136, 137, 138, 142, 143, 147, 147,
    170, 171, 173, 174, 175, 176, 178, 179, 181, 182, 183, 185, 190, 196, 199,
    202, 205, 209, 210, 211, 212, 213, 218, 219, 220, 221, 222, 223, 224, 225,
    227, 227, 228, 231, 235, 238, 243, 244, 247, 249, 250, 257, 261, 262, 267,
    282, 295, 309, 309, 310, 310, 310, 311, 312, 313, 314, 316, 319, 323, 326,
    332, 333, 335, 338, 341, 344, 348, 349, 350, 353, 359};

void TokenParser::yy_stack_print_() const
{
    *yycdebug_ << "Stack now";
    for (stack_type::const_iterator i = yystack_.begin(),
                                    i_end = yystack_.end();
         i != i_end; ++i)
        *yycdebug_ << ' ' << int(i->state);
    *yycdebug_ << '\n';
}

void TokenParser::yy_reduce_print_(int yyrule) const
{
    int yylno = yyrline_[yyrule];
    int yynrhs = yyr2_[yyrule];
    // Print the symbols being reduced, and their result.
    *yycdebug_ << "Reducing stack by rule " << yyrule - 1 << " (line " << yylno
               << "):\n";
    // The symbols being reduced.
    for (int yyi = 0; yyi < yynrhs; yyi++)
        YY_SYMBOL_PRINT("   $" << yyi + 1 << " =",
                        yystack_[(yynrhs) - (yyi + 1)]);
}
#endif // YYDEBUG

TokenParser::symbol_kind_type TokenParser::yytranslate_(int t) YY_NOEXCEPT
{
    // YYTRANSLATE[TOKEN-NUM] -- Symbol number corresponding to
    // TOKEN-NUM as returned by yylex.
    static const signed char translate_table[] = {
        0,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 2, 2, 2,  2,
        2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 2, 2, 26, 2,
        2,  2,  31, 32, 2,  33, 34, 36, 30, 2,  2,  2,  2,  2,  2, 2, 2, 2,  2,
        2,  35, 2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 2, 2, 2,  2,
        2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 2, 2, 2,  2,
        2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 2, 2, 2,  2,
        2,  2,  2,  2,  2,  2,  2,  2,  2,  28, 27, 29, 2,  2,  2, 2, 2, 2,  2,
        2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 2, 2, 2,  2,
        2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 2, 2, 2,  2,
        2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 2, 2, 2,  2,
        2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 2, 2, 2,  2,
        2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 2, 2, 2,  2,
        2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 2, 2, 2,  2,
        2,  2,  2,  2,  2,  2,  2,  2,  2,  1,  2,  3,  4,  5,  6, 7, 8, 9,  10,
        11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25};
    // Last valid token kind.
    const int code_max = 280;

    if (t <= 0)
        return symbol_kind::S_YYEOF;
    else if (t <= code_max)
        return static_cast<symbol_kind_type>(translate_table[t]);
    else
        return symbol_kind::S_YYUNDEF;
}

#line 15 "../helm_parser.yy"
} // namespace helm
#line 2068 "helm_parser.tab.cpp"

#line 360 "../helm_parser.yy"

void helm::TokenParser::error(const location_type& l,
                              const std::string& err_message)
{
    // save the pointer to the beginning of the current token
    helm_parser.saveError(l.begin.column, err_message);
}
