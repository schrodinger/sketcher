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

/*
 * This is a bison file which is used to generate a C++ parser for parsing a
 * HELM string.
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

#include "helm_parser.tab.hh"

// Unqualified %code blocks.

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

namespace helm
{

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
        case symbol_kind::S_HYDROGEN_PAIRING:        // HYDROGEN_PAIRING
        case symbol_kind::S_RGROUP:                  // RGROUP
        case symbol_kind::
            S_UNDEFINED_RESIDUE_NUMBER_OR_RGROUP: // UNDEFINED_RESIDUE_NUMBER_OR_RGROUP
        case symbol_kind::S_CONNECTION_RESIDUE:  // CONNECTION_RESIDUE
        case symbol_kind::S_POLYMER_GROUP_RATIO: // POLYMER_GROUP_RATIO
        case symbol_kind::S_POLYMER_GROUP_ID:    // POLYMER_GROUP_ID
        case symbol_kind::S_polymer:             // polymer
        case symbol_kind::S_repetitions:         // repetitions
        case symbol_kind::S_monomer_list:        // monomer_list
        case symbol_kind::S_monomer_and_list:    // monomer_and_list
        case symbol_kind::S_monomer_or_list:     // monomer_or_list
        case symbol_kind::S_monomer_list_item:   // monomer_list_item
        case symbol_kind::S_smiles_monomer:      // smiles_monomer
        case symbol_kind::S_monomer_id:          // monomer_id
        case symbol_kind::S_blob:                // blob
        case symbol_kind::S_connection_polymer:  // connection_polymer
        case symbol_kind::S_attachment_point:    // attachment_point
        case symbol_kind::S_connection_monomer:  // connection_monomer
        case symbol_kind::S_residue_and_list:    // residue_and_list
        case symbol_kind::S_residue_or_list:     // residue_or_list
        case symbol_kind::S_polymer_or_list:     // polymer_or_list
        case symbol_kind::S_polymer_and_list:    // polymer_and_list
        case symbol_kind::S_polymer_group_item:  // polymer_group_item
        case symbol_kind::S_annotation:          // annotation
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
        case symbol_kind::S_HYDROGEN_PAIRING:        // HYDROGEN_PAIRING
        case symbol_kind::S_RGROUP:                  // RGROUP
        case symbol_kind::
            S_UNDEFINED_RESIDUE_NUMBER_OR_RGROUP: // UNDEFINED_RESIDUE_NUMBER_OR_RGROUP
        case symbol_kind::S_CONNECTION_RESIDUE:  // CONNECTION_RESIDUE
        case symbol_kind::S_POLYMER_GROUP_RATIO: // POLYMER_GROUP_RATIO
        case symbol_kind::S_POLYMER_GROUP_ID:    // POLYMER_GROUP_ID
        case symbol_kind::S_polymer:             // polymer
        case symbol_kind::S_repetitions:         // repetitions
        case symbol_kind::S_monomer_list:        // monomer_list
        case symbol_kind::S_monomer_and_list:    // monomer_and_list
        case symbol_kind::S_monomer_or_list:     // monomer_or_list
        case symbol_kind::S_monomer_list_item:   // monomer_list_item
        case symbol_kind::S_smiles_monomer:      // smiles_monomer
        case symbol_kind::S_monomer_id:          // monomer_id
        case symbol_kind::S_blob:                // blob
        case symbol_kind::S_connection_polymer:  // connection_polymer
        case symbol_kind::S_attachment_point:    // attachment_point
        case symbol_kind::S_connection_monomer:  // connection_monomer
        case symbol_kind::S_residue_and_list:    // residue_and_list
        case symbol_kind::S_residue_or_list:     // residue_or_list
        case symbol_kind::S_polymer_or_list:     // polymer_or_list
        case symbol_kind::S_polymer_and_list:    // polymer_and_list
        case symbol_kind::S_polymer_group_item:  // polymer_group_item
        case symbol_kind::S_annotation:          // annotation
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
        case symbol_kind::S_HYDROGEN_PAIRING:        // HYDROGEN_PAIRING
        case symbol_kind::S_RGROUP:                  // RGROUP
        case symbol_kind::
            S_UNDEFINED_RESIDUE_NUMBER_OR_RGROUP: // UNDEFINED_RESIDUE_NUMBER_OR_RGROUP
        case symbol_kind::S_CONNECTION_RESIDUE:  // CONNECTION_RESIDUE
        case symbol_kind::S_POLYMER_GROUP_RATIO: // POLYMER_GROUP_RATIO
        case symbol_kind::S_POLYMER_GROUP_ID:    // POLYMER_GROUP_ID
        case symbol_kind::S_polymer:             // polymer
        case symbol_kind::S_repetitions:         // repetitions
        case symbol_kind::S_monomer_list:        // monomer_list
        case symbol_kind::S_monomer_and_list:    // monomer_and_list
        case symbol_kind::S_monomer_or_list:     // monomer_or_list
        case symbol_kind::S_monomer_list_item:   // monomer_list_item
        case symbol_kind::S_smiles_monomer:      // smiles_monomer
        case symbol_kind::S_monomer_id:          // monomer_id
        case symbol_kind::S_blob:                // blob
        case symbol_kind::S_connection_polymer:  // connection_polymer
        case symbol_kind::S_attachment_point:    // attachment_point
        case symbol_kind::S_connection_monomer:  // connection_monomer
        case symbol_kind::S_residue_and_list:    // residue_and_list
        case symbol_kind::S_residue_or_list:     // residue_or_list
        case symbol_kind::S_polymer_or_list:     // polymer_or_list
        case symbol_kind::S_polymer_and_list:    // polymer_and_list
        case symbol_kind::S_polymer_group_item:  // polymer_group_item
        case symbol_kind::S_annotation:          // annotation
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
        case symbol_kind::S_HYDROGEN_PAIRING:        // HYDROGEN_PAIRING
        case symbol_kind::S_RGROUP:                  // RGROUP
        case symbol_kind::
            S_UNDEFINED_RESIDUE_NUMBER_OR_RGROUP: // UNDEFINED_RESIDUE_NUMBER_OR_RGROUP
        case symbol_kind::S_CONNECTION_RESIDUE:  // CONNECTION_RESIDUE
        case symbol_kind::S_POLYMER_GROUP_RATIO: // POLYMER_GROUP_RATIO
        case symbol_kind::S_POLYMER_GROUP_ID:    // POLYMER_GROUP_ID
        case symbol_kind::S_polymer:             // polymer
        case symbol_kind::S_repetitions:         // repetitions
        case symbol_kind::S_monomer_list:        // monomer_list
        case symbol_kind::S_monomer_and_list:    // monomer_and_list
        case symbol_kind::S_monomer_or_list:     // monomer_or_list
        case symbol_kind::S_monomer_list_item:   // monomer_list_item
        case symbol_kind::S_smiles_monomer:      // smiles_monomer
        case symbol_kind::S_monomer_id:          // monomer_id
        case symbol_kind::S_blob:                // blob
        case symbol_kind::S_connection_polymer:  // connection_polymer
        case symbol_kind::S_attachment_point:    // attachment_point
        case symbol_kind::S_connection_monomer:  // connection_monomer
        case symbol_kind::S_residue_and_list:    // residue_and_list
        case symbol_kind::S_residue_or_list:     // residue_or_list
        case symbol_kind::S_polymer_or_list:     // polymer_or_list
        case symbol_kind::S_polymer_and_list:    // polymer_and_list
        case symbol_kind::S_polymer_group_item:  // polymer_group_item
        case symbol_kind::S_annotation:          // annotation
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
        case symbol_kind::S_HYDROGEN_PAIRING:        // HYDROGEN_PAIRING
        case symbol_kind::S_RGROUP:                  // RGROUP
        case symbol_kind::
            S_UNDEFINED_RESIDUE_NUMBER_OR_RGROUP: // UNDEFINED_RESIDUE_NUMBER_OR_RGROUP
        case symbol_kind::S_CONNECTION_RESIDUE:  // CONNECTION_RESIDUE
        case symbol_kind::S_POLYMER_GROUP_RATIO: // POLYMER_GROUP_RATIO
        case symbol_kind::S_POLYMER_GROUP_ID:    // POLYMER_GROUP_ID
        case symbol_kind::S_polymer:             // polymer
        case symbol_kind::S_repetitions:         // repetitions
        case symbol_kind::S_monomer_list:        // monomer_list
        case symbol_kind::S_monomer_and_list:    // monomer_and_list
        case symbol_kind::S_monomer_or_list:     // monomer_or_list
        case symbol_kind::S_monomer_list_item:   // monomer_list_item
        case symbol_kind::S_smiles_monomer:      // smiles_monomer
        case symbol_kind::S_monomer_id:          // monomer_id
        case symbol_kind::S_blob:                // blob
        case symbol_kind::S_connection_polymer:  // connection_polymer
        case symbol_kind::S_attachment_point:    // attachment_point
        case symbol_kind::S_connection_monomer:  // connection_monomer
        case symbol_kind::S_residue_and_list:    // residue_and_list
        case symbol_kind::S_residue_or_list:     // residue_or_list
        case symbol_kind::S_polymer_or_list:     // polymer_or_list
        case symbol_kind::S_polymer_and_list:    // polymer_and_list
        case symbol_kind::S_polymer_group_item:  // polymer_group_item
        case symbol_kind::S_annotation:          // annotation
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
        case symbol_kind::S_HYDROGEN_PAIRING:        // HYDROGEN_PAIRING
        case symbol_kind::S_RGROUP:                  // RGROUP
        case symbol_kind::
            S_UNDEFINED_RESIDUE_NUMBER_OR_RGROUP: // UNDEFINED_RESIDUE_NUMBER_OR_RGROUP
        case symbol_kind::S_CONNECTION_RESIDUE:  // CONNECTION_RESIDUE
        case symbol_kind::S_POLYMER_GROUP_RATIO: // POLYMER_GROUP_RATIO
        case symbol_kind::S_POLYMER_GROUP_ID:    // POLYMER_GROUP_ID
        case symbol_kind::S_polymer:             // polymer
        case symbol_kind::S_repetitions:         // repetitions
        case symbol_kind::S_monomer_list:        // monomer_list
        case symbol_kind::S_monomer_and_list:    // monomer_and_list
        case symbol_kind::S_monomer_or_list:     // monomer_or_list
        case symbol_kind::S_monomer_list_item:   // monomer_list_item
        case symbol_kind::S_smiles_monomer:      // smiles_monomer
        case symbol_kind::S_monomer_id:          // monomer_id
        case symbol_kind::S_blob:                // blob
        case symbol_kind::S_connection_polymer:  // connection_polymer
        case symbol_kind::S_attachment_point:    // attachment_point
        case symbol_kind::S_connection_monomer:  // connection_monomer
        case symbol_kind::S_residue_and_list:    // residue_and_list
        case symbol_kind::S_residue_or_list:     // residue_or_list
        case symbol_kind::S_polymer_or_list:     // polymer_or_list
        case symbol_kind::S_polymer_and_list:    // polymer_and_list
        case symbol_kind::S_polymer_group_item:  // polymer_group_item
        case symbol_kind::S_annotation:          // annotation
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
                case symbol_kind::S_HYDROGEN_PAIRING: // HYDROGEN_PAIRING
                case symbol_kind::S_RGROUP:           // RGROUP
                case symbol_kind::
                    S_UNDEFINED_RESIDUE_NUMBER_OR_RGROUP: // UNDEFINED_RESIDUE_NUMBER_OR_RGROUP
                case symbol_kind::S_CONNECTION_RESIDUE:  // CONNECTION_RESIDUE
                case symbol_kind::S_POLYMER_GROUP_RATIO: // POLYMER_GROUP_RATIO
                case symbol_kind::S_POLYMER_GROUP_ID:    // POLYMER_GROUP_ID
                case symbol_kind::S_polymer:             // polymer
                case symbol_kind::S_repetitions:         // repetitions
                case symbol_kind::S_monomer_list:        // monomer_list
                case symbol_kind::S_monomer_and_list:    // monomer_and_list
                case symbol_kind::S_monomer_or_list:     // monomer_or_list
                case symbol_kind::S_monomer_list_item:   // monomer_list_item
                case symbol_kind::S_smiles_monomer:      // smiles_monomer
                case symbol_kind::S_monomer_id:          // monomer_id
                case symbol_kind::S_blob:                // blob
                case symbol_kind::S_connection_polymer:  // connection_polymer
                case symbol_kind::S_attachment_point:    // attachment_point
                case symbol_kind::S_connection_monomer:  // connection_monomer
                case symbol_kind::S_residue_and_list:    // residue_and_list
                case symbol_kind::S_residue_or_list:     // residue_or_list
                case symbol_kind::S_polymer_or_list:     // polymer_or_list
                case symbol_kind::S_polymer_and_list:    // polymer_and_list
                case symbol_kind::S_polymer_group_item:  // polymer_group_item
                case symbol_kind::S_annotation:          // annotation
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
                    case 5: // connections: %empty
                    {
                        yylhs.value.as<size_t>() = 0;
                    } break;

                    case 6: // connections: connection
                    {
                        yylhs.value.as<size_t>() = 1;
                    } break;

                    case 7: // connections: connections '|' connection
                    {
                        if (YY_MOVE(yystack_[2].value.as<size_t>()) == 0) {
                            TokenParser::error(yystack_[2].location,
                                               "syntax error");
                            YYABORT;
                        }
                        ++yylhs.value.as<size_t>();
                    } break;

                    case 8: // polymer_groups: %empty
                    {
                        yylhs.value.as<size_t>() = 0;
                    } break;

                    case 9: // polymer_groups: polymer_group
                    {
                        yylhs.value.as<size_t>() = 1;
                    } break;

                    case 10: // polymer_groups: polymer_groups '|' polymer_group
                    {
                        if (YY_MOVE(yystack_[2].value.as<size_t>()) == 0) {
                            TokenParser::error(yystack_[2].location,
                                               "syntax error");
                            YYABORT;
                        }
                        ++yylhs.value.as<size_t>();
                    } break;

                    case 12: // polymer_unit: polymer
                    {
                        helm_parser.addPolymer(
                            YY_MOVE(yystack_[0].value.as<std::string_view>()),
                            {});
                    } break;

                    case 13: // polymer_unit: polymer annotation
                    {
                        helm_parser.addPolymer(
                            YY_MOVE(yystack_[1].value.as<std::string_view>()),
                            YY_MOVE(yystack_[0].value.as<std::string_view>()));
                    } break;

                    case 14: // polymer: PEPTIDE_ID '{' monomers '}'
                    {
                        yylhs.value.as<std::string_view>() =
                            YY_MOVE(yystack_[3].value.as<std::string_view>());
                    } break;

                    case 15: // polymer: RNA_ID '{' monomers '}'
                    {
                        yylhs.value.as<std::string_view>() =
                            YY_MOVE(yystack_[3].value.as<std::string_view>());
                    } break;

                    case 16: // polymer: CHEM_ID '{' monomer_unit '}'
                    {
                        yylhs.value.as<std::string_view>() =
                            YY_MOVE(yystack_[3].value.as<std::string_view>());
                    } break;

                    case 17: // polymer: BLOB_ID '{' blob '}'
                    {
                        yylhs.value.as<std::string_view>() =
                            YY_MOVE(yystack_[3].value.as<std::string_view>());
                    } break;

                    case 22: // monomer_group: branch_monomer_group
                    {
                        helm_parser.markBranchMonomer(
                            YY_MOVE(yystack_[0].value.as<size_t>()));
                    } break;

                    case 23: // monomer_group: monomer_unit '(' monomer_sequence
                             // ')'
                    {
                        helm_parser.saveError(yystack_[2].location.begin.column,
                                              branch_monomer_group_err_message);
                        YYABORT;
                    } break;

                    case 24: // monomer_group: monomer_unit '(' monomer_sequence
                             // ')' monomer_unit
                    {
                        helm_parser.saveError(yystack_[3].location.begin.column,
                                              branch_monomer_group_err_message);
                        YYABORT;
                    } break;

                    case 25: // repeated_monomers: monomer_unit repetitions
                    {
                        helm_parser.markLastNMonomersAsRepeated(
                            1,
                            YY_MOVE(yystack_[0].value.as<std::string_view>()),
                            {});
                    } break;

                    case 26: // repeated_monomers: monomer_unit repetitions
                             // annotation
                    {
                        helm_parser.markLastNMonomersAsRepeated(
                            1,
                            YY_MOVE(yystack_[1].value.as<std::string_view>()),
                            YY_MOVE(yystack_[0].value.as<std::string_view>()));
                    } break;

                    case 27: // repeated_monomers: '(' monomer_sequence ')'
                             // repetitions
                    {
                        helm_parser.markLastNMonomersAsRepeated(
                            YY_MOVE(yystack_[2].value.as<size_t>()),
                            YY_MOVE(yystack_[0].value.as<std::string_view>()),
                            {});
                    } break;

                    case 28: // repeated_monomers: '(' monomer_sequence ')'
                             // repetitions annotation
                    {
                        helm_parser.markLastNMonomersAsRepeated(
                            YY_MOVE(yystack_[3].value.as<size_t>()),
                            YY_MOVE(yystack_[1].value.as<std::string_view>()),
                            YY_MOVE(yystack_[0].value.as<std::string_view>()));
                    } break;

                    case 29: // repetitions: REPETITIONS
                    {
                        yylhs.value.as<std::string_view>() =
                            YY_MOVE(yystack_[0].value.as<std::string_view>());
                        yylhs.value.as<std::string_view>().remove_prefix(1);
                        yylhs.value.as<std::string_view>().remove_suffix(1);
                    } break;

                    case 30: // monomer_sequence: monomer_unit '.' monomer_unit
                    {
                        yylhs.value.as<size_t>() = 2;
                    } break;

                    case 31: // monomer_sequence: branch_monomer_group
                    {
                        yylhs.value.as<size_t>() =
                            YY_MOVE(yystack_[0].value.as<size_t>());
                        helm_parser.markBranchMonomer(yylhs.value.as<size_t>());
                    } break;

                    case 32: // monomer_sequence: monomer_sequence '.'
                             // monomer_unit
                    {
                        yylhs.value.as<size_t>() =
                            YY_MOVE(yystack_[2].value.as<size_t>()) + 1;
                    } break;

                    case 33: // monomer_sequence: monomer_sequence '.'
                             // branch_monomer_group
                    {
                        helm_parser.markBranchMonomer(
                            YY_MOVE(yystack_[0].value.as<size_t>()));
                        yylhs.value.as<size_t>() =
                            YY_MOVE(yystack_[2].value.as<size_t>()) +
                            YY_MOVE(yystack_[0].value.as<size_t>());
                    } break;

                    case 34: // branch_monomer_group: monomer_unit '('
                             // monomer_unit ')'
                    {
                        yylhs.value.as<size_t>() = 2;
                    } break;

                    case 35: // branch_monomer_group: monomer_unit '('
                             // monomer_unit ')' monomer_unit
                    {
                        yylhs.value.as<size_t>() = 3;
                    } break;

                    case 36: // monomer_unit: monomer_id
                    {
                        helm_parser.addMonomerWithId(
                            YY_MOVE(yystack_[0].value.as<std::string_view>()),
                            {});
                    } break;

                    case 37: // monomer_unit: monomer_id annotation
                    {
                        helm_parser.addMonomerWithId(
                            YY_MOVE(yystack_[1].value.as<std::string_view>()),
                            YY_MOVE(yystack_[0].value.as<std::string_view>()));
                    } break;

                    case 38: // monomer_unit: smiles_monomer
                    {
                        helm_parser.addSmilesMonomer(
                            YY_MOVE(yystack_[0].value.as<std::string_view>()),
                            {});
                    } break;

                    case 39: // monomer_unit: smiles_monomer annotation
                    {
                        helm_parser.addSmilesMonomer(
                            YY_MOVE(yystack_[1].value.as<std::string_view>()),
                            YY_MOVE(yystack_[0].value.as<std::string_view>()));
                    } break;

                    case 40: // monomer_unit: '(' monomer_list ')'
                    {
                        helm_parser.addMonomerList(
                            YY_MOVE(yystack_[1].value.as<std::string_view>()),
                            {});
                    } break;

                    case 41: // monomer_unit: '(' monomer_list ')' annotation
                    {
                        helm_parser.addMonomerList(
                            YY_MOVE(yystack_[2].value.as<std::string_view>()),
                            YY_MOVE(yystack_[0].value.as<std::string_view>()));
                    } break;

                    case 42: // monomer_list: monomer_and_list
                    {
                        yylhs.value.as<std::string_view>() =
                            YY_MOVE(yystack_[0].value.as<std::string_view>());
                    } break;

                    case 43: // monomer_list: monomer_or_list
                    {
                        yylhs.value.as<std::string_view>() =
                            YY_MOVE(yystack_[0].value.as<std::string_view>());
                    } break;

                    case 44: // monomer_and_list: monomer_list_item '+'
                             // monomer_list_item
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
                    } break;

                    case 45: // monomer_and_list: monomer_and_list '+'
                             // monomer_list_item
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
                    } break;

                    case 46: // monomer_or_list: monomer_list_item ','
                             // monomer_list_item
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
                    } break;

                    case 47: // monomer_or_list: monomer_or_list ','
                             // monomer_list_item
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
                    } break;

                    case 48: // monomer_list_item: monomer_id
                    {
                        yylhs.value.as<std::string_view>() =
                            YY_MOVE(yystack_[0].value.as<std::string_view>());
                    } break;

                    case 49: // monomer_list_item: monomer_id ':' MONOMER_RATIO
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
                    } break;

                    case 50: // smiles_monomer: INLINE_SMILES_MONOMER
                    {
                        yylhs.value.as<std::string_view>() =
                            YY_MOVE(yystack_[0].value.as<std::string_view>());
                        yylhs.value.as<std::string_view>().remove_prefix(1);
                        yylhs.value.as<std::string_view>().remove_suffix(1);
                    } break;

                    case 51: // monomer_id: SINGLE_CHARACTER_MONOMER
                    {
                        helm_parser.addResidueName(
                            YY_MOVE(yystack_[0].value.as<std::string_view>()));
                        yylhs.value.as<std::string_view>() =
                            YY_MOVE(yystack_[0].value.as<std::string_view>());
                    } break;

                    case 52: // monomer_id: MULTI_CHARACTER_MONOMER
                    {
                        auto tmp = std::string_view{
                            YY_MOVE(yystack_[0].value.as<std::string_view>())};
                        tmp.remove_prefix(1);
                        tmp.remove_suffix(1);
                        helm_parser.addResidueName(tmp);
                        yylhs.value.as<std::string_view>() =
                            YY_MOVE(yystack_[0].value.as<std::string_view>());
                    } break;

                    case 53: // monomer_id: UNKNOWN_MONOMER
                    {
                        helm_parser.addWildcardOrUnknownResidue();
                        yylhs.value.as<std::string_view>() =
                            YY_MOVE(yystack_[0].value.as<std::string_view>());
                    } break;

                    case 54: // monomer_id: MISSING_MONOMER
                    {
                        yylhs.value.as<std::string_view>() =
                            YY_MOVE(yystack_[0].value.as<std::string_view>());
                    } break;

                    case 55: // monomer_id: MONOMER_WILDCARD
                    {
                        helm_parser.addWildcardOrUnknownResidue();
                        yylhs.value.as<std::string_view>() =
                            YY_MOVE(yystack_[0].value.as<std::string_view>());
                    } break;

                    case 56: // blob: UNKNOWN_SEQUENCE
                    {
                        helm_parser.addMonomer(
                            YY_MOVE(yystack_[0].value.as<std::string_view>()),
                            false, false, false, {});
                    } break;

                    case 57: // connection: connection_polymer ','
                             // connection_polymer ',' connection_monomer ':'
                             // attachment_point '-' connection_monomer ':'
                             // attachment_point
                    {
                        helm_parser.addConnection({
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
                    } break;

                    case 58: // connection: connection_polymer ','
                             // connection_polymer ',' connection_monomer ':'
                             // attachment_point '-' connection_monomer ':'
                             // attachment_point annotation
                    {
                        helm_parser.addConnection({
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
                    } break;

                    case 59: // connection_polymer: POLYMER_ID
                    {
                        yylhs.value.as<std::string_view>() =
                            YY_MOVE(yystack_[0].value.as<std::string_view>());
                    } break;

                    case 60: // connection_polymer: POLYMER_GROUP_ID
                    {
                        yylhs.value.as<std::string_view>() =
                            YY_MOVE(yystack_[0].value.as<std::string_view>());
                    } break;

                    case 61: // attachment_point: HYDROGEN_PAIRING
                    {
                        yylhs.value.as<std::string_view>() =
                            YY_MOVE(yystack_[0].value.as<std::string_view>());
                    } break;

                    case 62: // attachment_point: RGROUP
                    {
                        yylhs.value.as<std::string_view>() =
                            YY_MOVE(yystack_[0].value.as<std::string_view>());
                    } break;

                    case 63: // attachment_point:
                             // UNDEFINED_RESIDUE_NUMBER_OR_RGROUP
                    {
                        yylhs.value.as<std::string_view>() =
                            YY_MOVE(yystack_[0].value.as<std::string_view>());
                    } break;

                    case 64: // connection_monomer: CONNECTION_RESIDUE
                    {
                        yylhs.value.as<std::string_view>() =
                            YY_MOVE(yystack_[0].value.as<std::string_view>());
                    } break;

                    case 65: // connection_monomer:
                             // UNDEFINED_RESIDUE_NUMBER_OR_RGROUP
                    {
                        yylhs.value.as<std::string_view>() =
                            YY_MOVE(yystack_[0].value.as<std::string_view>());
                    } break;

                    case 66: // connection_monomer: '(' residue_and_list ')'
                    {
                        yylhs.value.as<std::string_view>() =
                            YY_MOVE(yystack_[1].value.as<std::string_view>());
                    } break;

                    case 67: // connection_monomer: '(' residue_or_list ')'
                    {
                        yylhs.value.as<std::string_view>() =
                            YY_MOVE(yystack_[1].value.as<std::string_view>());
                    } break;

                    case 68: // residue_and_list: CONNECTION_RESIDUE '+'
                             // CONNECTION_RESIDUE
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
                    } break;

                    case 69: // residue_and_list: residue_and_list '+'
                             // CONNECTION_RESIDUE
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
                    } break;

                    case 70: // residue_or_list: CONNECTION_RESIDUE ','
                             // CONNECTION_RESIDUE
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
                    } break;

                    case 71: // residue_or_list: residue_or_list ','
                             // CONNECTION_RESIDUE
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
                    } break;

                    case 72: // polymer_group: POLYMER_GROUP_ID '('
                             // polymer_and_list ')'
                    {
                        helm_parser.addPolymerGroup(
                            YY_MOVE(yystack_[3].value.as<std::string_view>()),
                            YY_MOVE(yystack_[1].value.as<std::string_view>()),
                            true);
                    } break;

                    case 73: // polymer_group: POLYMER_GROUP_ID '('
                             // polymer_or_list ')'
                    {
                        helm_parser.addPolymerGroup(
                            YY_MOVE(yystack_[3].value.as<std::string_view>()),
                            YY_MOVE(yystack_[1].value.as<std::string_view>()),
                            false);
                    } break;

                    case 74: // polymer_or_list: polymer_group_item ','
                             // polymer_group_item
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
                    } break;

                    case 75: // polymer_or_list: polymer_or_list ','
                             // polymer_group_item
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
                    } break;

                    case 76: // polymer_and_list: polymer_group_item '+'
                             // polymer_group_item
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
                    } break;

                    case 77: // polymer_and_list: polymer_and_list '+'
                             // polymer_group_item
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
                    } break;

                    case 78: // polymer_group_item: POLYMER_ID
                    {
                        yylhs.value.as<std::string_view>() =
                            YY_MOVE(yystack_[0].value.as<std::string_view>());
                    } break;

                    case 79: // polymer_group_item: POLYMER_GROUP_ID
                    {
                        yylhs.value.as<std::string_view>() =
                            YY_MOVE(yystack_[0].value.as<std::string_view>());
                    } break;

                    case 80: // polymer_group_item: POLYMER_ID ':'
                             // POLYMER_GROUP_RATIO
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
                    } break;

                    case 81: // polymer_group_item: POLYMER_GROUP_ID ':'
                             // POLYMER_GROUP_RATIO
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
                    } break;

                    case 82: // annotation: ANNOTATION
                    {
                        yylhs.value.as<std::string_view>() =
                            YY_MOVE(yystack_[0].value.as<std::string_view>());
                        yylhs.value.as<std::string_view>().remove_prefix(1);
                        yylhs.value.as<std::string_view>().remove_suffix(1);
                    } break;

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

const signed char TokenParser::yypact_ninf_ = -53;

const signed char TokenParser::yytable_ninf_ = -49;

const signed char TokenParser::yypact_[] = {
    55,  41,  53,  59,  83,  12,  58,  -53, 25,  93,  -1,  27,  27,  -53,
    -4,  55,  -53, -53, -53, 54,  -53, -53, -53, -53, -53, -53, 14,  67,
    25,  25,  -1,  61,  -53, -53, -53, -12, 63,  -53, -53, 68,  -53, 79,
    -53, -53, 82,  84,  81,  64,  85,  -53, -53, -53, 23,  -53, 69,  0,
    -53, 27,  -53, -1,  25,  -53, 91,  -4,  -4,  25,  14,  14,  14,  14,
    103, -1,  101, -1,  -1,  -53, 50,  19,  -53, 90,  75,  -53, -53, 88,
    -53, -53, -53, -53, -53, -53, -53, 92,  25,  -53, 94,  -1,  -1,  7,
    -53, 91,  34,  -53, -53, -53, 89,  95,  -14, 71,  72,  -53, -53, -53,
    -53, 102, 96,  104, 105, -53, 7,   -53, 7,   7,   7,   74,  77,  -3,
    57,  -53, -53, -53, -53, -53, -53, 109, 110, -53, 111, -53, 112, -53,
    -53, -53, 100, -53, -53, -53, -53, 34,  106, 57,  25,  -53};

const signed char TokenParser::yydefact_[] = {
    0,  0,  0,  0,  0,  0,  0,  3,  12, 0,  0,  0,  0,  1,  5,  0,  82, 13, 56,
    0,  51, 52, 53, 54, 55, 50, 0,  0,  38, 36, 0,  0,  18, 21, 22, 20, 0,  59,
    60, 0,  6,  0,  4,  17, 0,  42, 43, 0,  48, 16, 39, 37, 0,  31, 0,  36, 14,
    0,  29, 0,  25, 15, 8,  0,  0,  40, 0,  0,  0,  0,  0,  0,  0,  0,  0,  19,
    0,  0,  26, 0,  0,  9,  7,  0,  41, 45, 47, 44, 46, 49, 33, 32, 27, 30, 0,
    23, 34, 0,  11, 0,  0,  28, 24, 35, 78, 79, 0,  0,  0,  2,  10, 65, 64, 0,
    0,  0,  0,  73, 0,  72, 0,  0,  0,  0,  0,  0,  0,  80, 81, 75, 77, 76, 74,
    0,  0,  66, 0,  67, 0,  61, 62, 63, 0,  68, 70, 69, 71, 0,  0,  0,  57, 58};

const signed char TokenParser::yypgoto_[] = {
    -53, -53, -53, -53, -53, -53, 120, -53, 113, 80,  -53,
    66,  86,  -25, -9,  -53, -53, -53, 5,   -53, -24, -53,
    73,  76,  -6,  1,   -53, -53, 42,  -53, -53, -52, -8};

const unsigned char TokenParser::yydefgoto_[] = {
    0,  5,  6,  39, 80, 109, 7,  8,   31,  32,  33,  60, 52,  34,  35,  44, 45,
    46, 47, 28, 29, 19, 40,  41, 142, 114, 124, 125, 81, 106, 107, 108, 51};

const short TokenParser::yytable_[] = {
    17,  27,  48,  37,  58,  53, 55,  20,  21,  22,  23,  24,  13,  25,  104,
    16,  117, 59,  118, 38,  50, 54,  20,  21,  22,  23,  24,  137, 26,  138,
    105, -48, -48, 70,  53,  20, 21,  22,  23,  24,  16,  25,  48,  48,  48,
    48,  90,  73,  74,  96,  77, 71,  78,  72,  111, 112, 30,  84,  1,   2,
    3,   4,   91,  113, 93,  94, 129, 9,   130, 131, 132, 85,  86,  87,  88,
    139, 140, 141, 71,  10,  95, 43,  14,  15,  101, 11,  102, 103, 56,  57,
    61,  57,  62,  63,  49,  68, 69,  73,  74,  98,  99,  119, 120, 121, 122,
    133, 134, 135, 136, 12,  18, 64,  65,  67,  79,  66,  89,  58,  70,  97,
    100, 74,  115, 123, 96,  36, 127, 128, 116, 126, 143, 144, 145, 146, 147,
    42,  82,  75,  92,  149, 83, 110, 151, 150, 0,   76,  0,   0,   148};

const short TokenParser::yycheck_[] = {
    8,   10,  26,  7,  16, 30, 30,  8,   9,  10, 11, 12, 0,  14, 7,  15,  30,
    29,  32,  23,  28, 30, 8,  9,   10,  11, 12, 30, 29, 32, 23, 31, 32,  33,
    59,  8,   9,   10, 11, 12, 15,  14,  66, 67, 68, 69, 71, 28, 29, 30,  59,
    28,  60,  30,  20, 21, 29, 65,  3,   4,  5,  6,  71, 29, 73, 74, 118, 26,
    120, 121, 122, 66, 67, 68, 69,  18,  19, 20, 28, 26, 30, 27, 24, 25,  92,
    26,  95,  96,  27, 28, 27, 28,  24,  25, 27, 31, 32, 28, 29, 24, 25,  30,
    31,  31,  32,  31, 32, 30, 31,  26,  17, 32, 30, 32, 23, 31, 13, 16,  33,
    29,  32,  29,  33, 21, 30, 12,  22,  22, 33, 33, 21, 21, 21, 21, 34,  15,
    63,  57,  72,  33, 64, 99, 150, 149, -1, 59, -1, -1, 147};

const signed char TokenParser::yystos_[] = {
    0,  3,  4,  5,  6,  36, 37, 41, 42, 26, 26, 26, 26, 0,  24, 25, 15, 67, 17,
    56, 8,  9,  10, 11, 12, 14, 29, 49, 54, 55, 29, 43, 44, 45, 48, 49, 43, 7,
    23, 38, 57, 58, 41, 27, 50, 51, 52, 53, 55, 27, 67, 67, 47, 48, 49, 55, 27,
    28, 16, 29, 46, 27, 24, 25, 32, 30, 31, 32, 31, 32, 33, 28, 30, 28, 29, 44,
    47, 49, 67, 23, 39, 63, 57, 58, 67, 53, 53, 53, 53, 13, 48, 49, 46, 49, 49,
    30, 30, 29, 24, 25, 32, 67, 49, 49, 7,  23, 64, 65, 66, 40, 63, 20, 21, 29,
    60, 33, 33, 30, 32, 30, 31, 31, 32, 21, 61, 62, 33, 22, 22, 66, 66, 66, 66,
    31, 32, 30, 31, 30, 32, 18, 19, 20, 59, 21, 21, 21, 21, 34, 60, 33, 59, 67};

const signed char TokenParser::yyr1_[] = {
    0,  35, 36, 37, 37, 38, 38, 38, 39, 39, 39, 40, 41, 41, 42, 42, 42,
    42, 43, 43, 44, 44, 44, 44, 44, 45, 45, 45, 45, 46, 47, 47, 47, 47,
    48, 48, 49, 49, 49, 49, 49, 49, 50, 50, 51, 51, 52, 52, 53, 53, 54,
    55, 55, 55, 55, 55, 56, 57, 57, 58, 58, 59, 59, 59, 60, 60, 60, 60,
    61, 61, 62, 62, 63, 63, 64, 64, 65, 65, 66, 66, 66, 66, 67};

const signed char TokenParser::yyr2_[] = {
    0, 2, 7, 1, 3, 0, 1, 3, 0, 1, 3, 0, 1, 2, 4, 4,  4,  4, 1, 3, 1,
    1, 1, 4, 5, 2, 3, 4, 5, 1, 3, 1, 3, 3, 4, 5, 1,  2,  1, 2, 3, 4,
    1, 1, 3, 3, 3, 3, 1, 3, 1, 1, 1, 1, 1, 1, 1, 11, 12, 1, 1, 1, 1,
    1, 1, 1, 3, 3, 3, 3, 3, 3, 4, 4, 3, 3, 3, 3, 1,  1,  3, 3, 1};

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
    "HYDROGEN_PAIRING",
    "RGROUP",
    "UNDEFINED_RESIDUE_NUMBER_OR_RGROUP",
    "CONNECTION_RESIDUE",
    "POLYMER_GROUP_RATIO",
    "POLYMER_GROUP_ID",
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
    "extended_annotations_and_version",
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
    0,   118, 118, 128, 129, 131, 132, 133, 137, 138, 139, 143, 164, 165,
    167, 168, 169, 170, 172, 173, 175, 176, 177, 179, 184, 190, 193, 196,
    199, 203, 204, 205, 206, 207, 212, 213, 214, 215, 216, 217, 218, 219,
    221, 221, 222, 225, 229, 232, 237, 238, 241, 243, 244, 251, 255, 256,
    261, 276, 289, 303, 303, 304, 304, 304, 305, 306, 307, 308, 310, 313,
    317, 320, 326, 327, 329, 332, 335, 338, 342, 343, 344, 347, 353};

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
        0,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 2, 2, 2, 2, 2,  2,
        2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 2, 2, 2, 2, 24, 2,
        2,  2,  29, 30, 2,  31, 32, 34, 28, 2,  2,  2,  2, 2, 2, 2, 2, 2,  2,
        2,  33, 2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 2, 2, 2, 2, 2,  2,
        2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 2, 2, 2, 2, 2,  2,
        2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 2, 2, 2, 2, 2,  2,
        2,  2,  2,  2,  2,  2,  2,  2,  2,  26, 25, 27, 2, 2, 2, 2, 2, 2,  2,
        2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 2, 2, 2, 2, 2,  2,
        2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 2, 2, 2, 2, 2,  2,
        2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 2, 2, 2, 2, 2,  2,
        2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 2, 2, 2, 2, 2,  2,
        2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 2, 2, 2, 2, 2,  2,
        2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 2, 2, 2, 2, 2,  2,
        2,  2,  2,  2,  2,  2,  2,  2,  2,  1,  2,  3,  4, 5, 6, 7, 8, 9,  10,
        11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23};
    // Last valid token kind.
    const int code_max = 278;

    if (t <= 0)
        return symbol_kind::S_YYEOF;
    else if (t <= code_max)
        return static_cast<symbol_kind_type>(translate_table[t]);
    else
        return symbol_kind::S_YYUNDEF;
}

} // namespace helm

void helm::TokenParser::error(const location_type& l,
                              const std::string& err_message)
{
    // save the pointer to the beginning of the current token
    helm_parser.saveError(l.begin.column, err_message);
}
