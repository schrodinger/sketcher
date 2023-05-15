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
 * a helm::HelmParser instance to add polymer and monomer information to an
 * ROMol object.
 */

#line 52 "helm_parser.tab.cpp"

#include "helm_parser.tab.hh"

// Unqualified %code blocks.
#line 31 "../helm_parser.yy"

#include <string>
#include <string_view>

#include "schrodinger/rdkit_extensions/helm/helm_parser.h"

#undef yylex
#define yylex scanner.lex

#line 69 "helm_parser.tab.cpp"

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

#line 16 "../helm_parser.yy"
namespace helm
{
#line 162 "helm_parser.tab.cpp"

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
        case symbol_kind::S_BLOB_ID:    // BLOB_ID
        case symbol_kind::S_CHEM_ID:    // CHEM_ID
        case symbol_kind::S_PEPTIDE_ID: // PEPTIDE_ID
        case symbol_kind::S_RNA_ID:     // RNA_ID
        case symbol_kind::
            S_SINGLE_CHARACTER_MONOMER:              // SINGLE_CHARACTER_MONOMER
        case symbol_kind::S_MULTI_CHARACTER_MONOMER: // MULTI_CHARACTER_MONOMER
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
        case symbol_kind::S_BLOB_ID:    // BLOB_ID
        case symbol_kind::S_CHEM_ID:    // CHEM_ID
        case symbol_kind::S_PEPTIDE_ID: // PEPTIDE_ID
        case symbol_kind::S_RNA_ID:     // RNA_ID
        case symbol_kind::
            S_SINGLE_CHARACTER_MONOMER:              // SINGLE_CHARACTER_MONOMER
        case symbol_kind::S_MULTI_CHARACTER_MONOMER: // MULTI_CHARACTER_MONOMER
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
        case symbol_kind::S_BLOB_ID:    // BLOB_ID
        case symbol_kind::S_CHEM_ID:    // CHEM_ID
        case symbol_kind::S_PEPTIDE_ID: // PEPTIDE_ID
        case symbol_kind::S_RNA_ID:     // RNA_ID
        case symbol_kind::
            S_SINGLE_CHARACTER_MONOMER:              // SINGLE_CHARACTER_MONOMER
        case symbol_kind::S_MULTI_CHARACTER_MONOMER: // MULTI_CHARACTER_MONOMER
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
        case symbol_kind::S_BLOB_ID:    // BLOB_ID
        case symbol_kind::S_CHEM_ID:    // CHEM_ID
        case symbol_kind::S_PEPTIDE_ID: // PEPTIDE_ID
        case symbol_kind::S_RNA_ID:     // RNA_ID
        case symbol_kind::
            S_SINGLE_CHARACTER_MONOMER:              // SINGLE_CHARACTER_MONOMER
        case symbol_kind::S_MULTI_CHARACTER_MONOMER: // MULTI_CHARACTER_MONOMER
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
        case symbol_kind::S_BLOB_ID:    // BLOB_ID
        case symbol_kind::S_CHEM_ID:    // CHEM_ID
        case symbol_kind::S_PEPTIDE_ID: // PEPTIDE_ID
        case symbol_kind::S_RNA_ID:     // RNA_ID
        case symbol_kind::
            S_SINGLE_CHARACTER_MONOMER:              // SINGLE_CHARACTER_MONOMER
        case symbol_kind::S_MULTI_CHARACTER_MONOMER: // MULTI_CHARACTER_MONOMER
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
        case symbol_kind::S_BLOB_ID:    // BLOB_ID
        case symbol_kind::S_CHEM_ID:    // CHEM_ID
        case symbol_kind::S_PEPTIDE_ID: // PEPTIDE_ID
        case symbol_kind::S_RNA_ID:     // RNA_ID
        case symbol_kind::
            S_SINGLE_CHARACTER_MONOMER:              // SINGLE_CHARACTER_MONOMER
        case symbol_kind::S_MULTI_CHARACTER_MONOMER: // MULTI_CHARACTER_MONOMER
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
                case symbol_kind::S_BLOB_ID:    // BLOB_ID
                case symbol_kind::S_CHEM_ID:    // CHEM_ID
                case symbol_kind::S_PEPTIDE_ID: // PEPTIDE_ID
                case symbol_kind::S_RNA_ID:     // RNA_ID
                case symbol_kind::
                    S_SINGLE_CHARACTER_MONOMER: // SINGLE_CHARACTER_MONOMER
                case symbol_kind::
                    S_MULTI_CHARACTER_MONOMER: // MULTI_CHARACTER_MONOMER
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
                    case 4: // polymer: PEPTIDE_ID '{' monomers '}'
#line 61 "../helm_parser.yy"
                    {
                        helm_parser.add_polymer(
                            yystack_[3].value.as<std::string_view>());
                    }
#line 715 "helm_parser.tab.cpp"
                    break;

                    case 5: // polymer: RNA_ID '{' monomers '}'
#line 62 "../helm_parser.yy"
                    {
                        helm_parser.add_polymer(
                            yystack_[3].value.as<std::string_view>());
                    }
#line 721 "helm_parser.tab.cpp"
                    break;

                    case 6: // polymer: CHEM_ID '{' monomer '}'
#line 63 "../helm_parser.yy"
                    {
                        helm_parser.add_polymer(
                            yystack_[3].value.as<std::string_view>());
                    }
#line 727 "helm_parser.tab.cpp"
                    break;

                    case 7: // polymer: BLOB_ID '{' bead '}'
#line 64 "../helm_parser.yy"
                    {
                        helm_parser.add_polymer(
                            yystack_[3].value.as<std::string_view>());
                    }
#line 733 "helm_parser.tab.cpp"
                    break;

                    case 10: // monomer: SINGLE_CHARACTER_MONOMER
#line 71 "../helm_parser.yy"
                    {
                        helm_parser.add_monomer(
                            yystack_[0].value.as<std::string_view>());
                    }
#line 739 "helm_parser.tab.cpp"
                    break;

                    case 11: // monomer: '[' MULTI_CHARACTER_MONOMER ']'
#line 72 "../helm_parser.yy"
                    {
                        helm_parser.add_monomer(
                            yystack_[1].value.as<std::string_view>());
                    }
#line 745 "helm_parser.tab.cpp"
                    break;

                    case 12: // bead: BEAD
#line 75 "../helm_parser.yy"
                    {
                        helm_parser.add_monomer("BEAD");
                    }
#line 751 "helm_parser.tab.cpp"
                    break;

#line 755 "helm_parser.tab.cpp"

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

const signed char TokenParser::yypact_ninf_ = -10;

const signed char TokenParser::yytable_ninf_ = -1;

const signed char TokenParser::yypact_[] = {
    -3,  -6,  1,   2,   3,   10,  5,  -10, 14, -8, -8,  -8, -10,
    -10, -10, 6,   -10, 11,  7,   -5, -10, -2, 12, -10, 8,  -10,
    -10, -8,  -10, -10, -10, -10, 15, -10, 16, 18, -10};

const signed char TokenParser::yydefact_[] = {
    0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0,  0,  1, 13, 12, 0, 10, 0, 0,
    0, 8, 0, 0, 7, 0, 6, 4, 0, 5, 14, 11, 9, 0,  15, 0, 0,  2};

const signed char TokenParser::yypgoto_[] = {-10, -10, -10, -10, 17,
                                             -9,  -10, -10, -10, -10};

const signed char TokenParser::yydefgoto_[] = {0,  5,  6,  7,  19,
                                               20, 15, 22, 32, 34};

const signed char TokenParser::yytable_[] = {
    18, 16, 1,  2,  3,  4,  8,  17, 26, 27, 12, 28, 27, 9, 10,
    11, 13, 14, 31, 23, 25, 24, 36, 29, 30, 0,  33, 35, 21};

const signed char TokenParser::yycheck_[] = {
    9,  9,  5, 6,  7,  8,  12, 15, 13, 14, 0,  13, 14, 12, 12,
    12, 11, 3, 27, 13, 13, 10, 4,  11, 16, -1, 11, 11, 11};

const signed char TokenParser::yystos_[] = {
    0,  5,  6,  7,  8,  18, 19, 20, 12, 12, 12, 12, 0,  11, 3,  23, 9,  15, 22,
    21, 22, 21, 24, 13, 10, 13, 13, 14, 13, 11, 16, 22, 25, 11, 26, 11, 4};

const signed char TokenParser::yyr1_[] = {0,  17, 18, 19, 20, 20, 20, 20,
                                          21, 21, 22, 22, 23, 24, 25, 26};

const signed char TokenParser::yyr2_[] = {0, 2, 9, 1, 4, 4, 4, 4,
                                          1, 3, 1, 3, 1, 0, 0, 0};

#if YYDEBUG
// YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
// First, the terminals, then, starting at \a YYNTOKENS, nonterminals.
const char* const TokenParser::yytname_[] = {"END",
                                             "error",
                                             "\"invalid token\"",
                                             "BEAD",
                                             "VERSION_TOKEN",
                                             "BLOB_ID",
                                             "CHEM_ID",
                                             "PEPTIDE_ID",
                                             "RNA_ID",
                                             "SINGLE_CHARACTER_MONOMER",
                                             "MULTI_CHARACTER_MONOMER",
                                             "'$'",
                                             "'{'",
                                             "'}'",
                                             "'.'",
                                             "'['",
                                             "']'",
                                             "$accept",
                                             "helm",
                                             "polymers",
                                             "polymer",
                                             "monomers",
                                             "monomer",
                                             "bead",
                                             "connections",
                                             "polymer_groups",
                                             "extended_annotations",
                                             YY_NULLPTR};
#endif

#if YYDEBUG
const signed char TokenParser::yyrline_[] = {0,  56, 56, 58, 61, 62, 63, 64,
                                             67, 68, 71, 72, 75, 78, 79, 80};

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
        0, 2, 2, 2, 2,  2, 2, 2,  2, 2,  2, 2, 2, 2, 2, 2,  2, 2, 2,  2, 2,
        2, 2, 2, 2, 2,  2, 2, 2,  2, 2,  2, 2, 2, 2, 2, 11, 2, 2, 2,  2, 2,
        2, 2, 2, 2, 14, 2, 2, 2,  2, 2,  2, 2, 2, 2, 2, 2,  2, 2, 2,  2, 2,
        2, 2, 2, 2, 2,  2, 2, 2,  2, 2,  2, 2, 2, 2, 2, 2,  2, 2, 2,  2, 2,
        2, 2, 2, 2, 2,  2, 2, 15, 2, 16, 2, 2, 2, 2, 2, 2,  2, 2, 2,  2, 2,
        2, 2, 2, 2, 2,  2, 2, 2,  2, 2,  2, 2, 2, 2, 2, 2,  2, 2, 12, 2, 13,
        2, 2, 2, 2, 2,  2, 2, 2,  2, 2,  2, 2, 2, 2, 2, 2,  2, 2, 2,  2, 2,
        2, 2, 2, 2, 2,  2, 2, 2,  2, 2,  2, 2, 2, 2, 2, 2,  2, 2, 2,  2, 2,
        2, 2, 2, 2, 2,  2, 2, 2,  2, 2,  2, 2, 2, 2, 2, 2,  2, 2, 2,  2, 2,
        2, 2, 2, 2, 2,  2, 2, 2,  2, 2,  2, 2, 2, 2, 2, 2,  2, 2, 2,  2, 2,
        2, 2, 2, 2, 2,  2, 2, 2,  2, 2,  2, 2, 2, 2, 2, 2,  2, 2, 2,  2, 2,
        2, 2, 2, 2, 2,  2, 2, 2,  2, 2,  2, 2, 2, 2, 2, 2,  2, 2, 2,  2, 2,
        2, 2, 2, 2, 1,  2, 3, 4,  5, 6,  7, 8, 9, 10};
    // Last valid token kind.
    const int code_max = 265;

    if (t <= 0)
        return symbol_kind::S_YYEOF;
    else if (t <= code_max)
        return static_cast<symbol_kind_type>(translate_table[t]);
    else
        return symbol_kind::S_YYUNDEF;
}

#line 16 "../helm_parser.yy"
} // namespace helm
#line 1122 "helm_parser.tab.cpp"

#line 83 "../helm_parser.yy"

void helm::TokenParser::error(const location_type& l,
                              const std::string& err_message)
{
    // save the pointer to the beginning of the current token
    helm_parser.saveErrorInformation(l.begin.column);
}
