// A Bison parser, made by GNU Bison 3.8.2.

// Skeleton interface for Bison LALR(1) parsers in C++

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


/**
 ** \file helm_parser.tab.hh
 ** Define the helm::parser class.
 */

// C++ LALR(1) parser skeleton written by Akim Demaille.

// DO NOT RELY ON FEATURES THAT ARE NOT DOCUMENTED in the manual,
// especially those whose name start with YY_ or yy_.  They are
// private implementation details that can be changed or removed.

#ifndef YY_YY_HELM_PARSER_TAB_HH_INCLUDED
# define YY_YY_HELM_PARSER_TAB_HH_INCLUDED
// "%code requires" blocks.
#line 21 "../helm_parser.yy"

   #include <string>
   #include <string_view>

   namespace helm {
      class HelmParser;
      class TokenScanner;
   };

#line 59 "helm_parser.tab.hh"


# include <cstdlib> // std::abort
# include <iostream>
# include <stdexcept>
# include <string>
# include <vector>

#if defined __cplusplus
# define YY_CPLUSPLUS __cplusplus
#else
# define YY_CPLUSPLUS 199711L
#endif

// Support move semantics when possible.
#if 201103L <= YY_CPLUSPLUS
# define YY_MOVE           std::move
# define YY_MOVE_OR_COPY   move
# define YY_MOVE_REF(Type) Type&&
# define YY_RVREF(Type)    Type&&
# define YY_COPY(Type)     Type
#else
# define YY_MOVE
# define YY_MOVE_OR_COPY   copy
# define YY_MOVE_REF(Type) Type&
# define YY_RVREF(Type)    const Type&
# define YY_COPY(Type)     const Type&
#endif

// Support noexcept when possible.
#if 201103L <= YY_CPLUSPLUS
# define YY_NOEXCEPT noexcept
# define YY_NOTHROW
#else
# define YY_NOEXCEPT
# define YY_NOTHROW throw ()
#endif

// Support constexpr when possible.
#if 201703 <= YY_CPLUSPLUS
# define YY_CONSTEXPR constexpr
#else
# define YY_CONSTEXPR
#endif
# include "location.hh"


#ifndef YY_ATTRIBUTE_PURE
# if defined __GNUC__ && 2 < __GNUC__ + (96 <= __GNUC_MINOR__)
#  define YY_ATTRIBUTE_PURE __attribute__ ((__pure__))
# else
#  define YY_ATTRIBUTE_PURE
# endif
#endif

#ifndef YY_ATTRIBUTE_UNUSED
# if defined __GNUC__ && 2 < __GNUC__ + (7 <= __GNUC_MINOR__)
#  define YY_ATTRIBUTE_UNUSED __attribute__ ((__unused__))
# else
#  define YY_ATTRIBUTE_UNUSED
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YY_USE(E) ((void) (E))
#else
# define YY_USE(E) /* empty */
#endif

/* Suppress an incorrect diagnostic about yylval being uninitialized.  */
#if defined __GNUC__ && ! defined __ICC && 406 <= __GNUC__ * 100 + __GNUC_MINOR__
# if __GNUC__ * 100 + __GNUC_MINOR__ < 407
#  define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN                           \
    _Pragma ("GCC diagnostic push")                                     \
    _Pragma ("GCC diagnostic ignored \"-Wuninitialized\"")
# else
#  define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN                           \
    _Pragma ("GCC diagnostic push")                                     \
    _Pragma ("GCC diagnostic ignored \"-Wuninitialized\"")              \
    _Pragma ("GCC diagnostic ignored \"-Wmaybe-uninitialized\"")
# endif
# define YY_IGNORE_MAYBE_UNINITIALIZED_END      \
    _Pragma ("GCC diagnostic pop")
#else
# define YY_INITIAL_VALUE(Value) Value
#endif
#ifndef YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_END
#endif
#ifndef YY_INITIAL_VALUE
# define YY_INITIAL_VALUE(Value) /* Nothing. */
#endif

#if defined __cplusplus && defined __GNUC__ && ! defined __ICC && 6 <= __GNUC__
# define YY_IGNORE_USELESS_CAST_BEGIN                          \
    _Pragma ("GCC diagnostic push")                            \
    _Pragma ("GCC diagnostic ignored \"-Wuseless-cast\"")
# define YY_IGNORE_USELESS_CAST_END            \
    _Pragma ("GCC diagnostic pop")
#endif
#ifndef YY_IGNORE_USELESS_CAST_BEGIN
# define YY_IGNORE_USELESS_CAST_BEGIN
# define YY_IGNORE_USELESS_CAST_END
#endif

# ifndef YY_CAST
#  ifdef __cplusplus
#   define YY_CAST(Type, Val) static_cast<Type> (Val)
#   define YY_REINTERPRET_CAST(Type, Val) reinterpret_cast<Type> (Val)
#  else
#   define YY_CAST(Type, Val) ((Type) (Val))
#   define YY_REINTERPRET_CAST(Type, Val) ((Type) (Val))
#  endif
# endif
# ifndef YY_NULLPTR
#  if defined __cplusplus
#   if 201103L <= __cplusplus
#    define YY_NULLPTR nullptr
#   else
#    define YY_NULLPTR 0
#   endif
#  else
#   define YY_NULLPTR ((void*)0)
#  endif
# endif

/* Debug traces.  */
#ifndef YYDEBUG
# define YYDEBUG 1
#endif

#line 15 "../helm_parser.yy"
namespace helm {
#line 195 "helm_parser.tab.hh"




  /// A Bison parser.
  class TokenParser
  {
  public:
#ifdef YYSTYPE
# ifdef __GNUC__
#  pragma GCC message "bison: do not #define YYSTYPE in C++, use %define api.value.type"
# endif
    typedef YYSTYPE value_type;
#else
  /// A buffer to store and retrieve objects.
  ///
  /// Sort of a variant, but does not keep track of the nature
  /// of the stored data, since that knowledge is available
  /// via the current parser state.
  class value_type
  {
  public:
    /// Type of *this.
    typedef value_type self_type;

    /// Empty construction.
    value_type () YY_NOEXCEPT
      : yyraw_ ()
    {}

    /// Construct and fill.
    template <typename T>
    value_type (YY_RVREF (T) t)
    {
      new (yyas_<T> ()) T (YY_MOVE (t));
    }

#if 201103L <= YY_CPLUSPLUS
    /// Non copyable.
    value_type (const self_type&) = delete;
    /// Non copyable.
    self_type& operator= (const self_type&) = delete;
#endif

    /// Destruction, allowed only if empty.
    ~value_type () YY_NOEXCEPT
    {}

# if 201103L <= YY_CPLUSPLUS
    /// Instantiate a \a T in here from \a t.
    template <typename T, typename... U>
    T&
    emplace (U&&... u)
    {
      return *new (yyas_<T> ()) T (std::forward <U>(u)...);
    }
# else
    /// Instantiate an empty \a T in here.
    template <typename T>
    T&
    emplace ()
    {
      return *new (yyas_<T> ()) T ();
    }

    /// Instantiate a \a T in here from \a t.
    template <typename T>
    T&
    emplace (const T& t)
    {
      return *new (yyas_<T> ()) T (t);
    }
# endif

    /// Instantiate an empty \a T in here.
    /// Obsolete, use emplace.
    template <typename T>
    T&
    build ()
    {
      return emplace<T> ();
    }

    /// Instantiate a \a T in here from \a t.
    /// Obsolete, use emplace.
    template <typename T>
    T&
    build (const T& t)
    {
      return emplace<T> (t);
    }

    /// Accessor to a built \a T.
    template <typename T>
    T&
    as () YY_NOEXCEPT
    {
      return *yyas_<T> ();
    }

    /// Const accessor to a built \a T (for %printer).
    template <typename T>
    const T&
    as () const YY_NOEXCEPT
    {
      return *yyas_<T> ();
    }

    /// Swap the content with \a that, of same type.
    ///
    /// Both variants must be built beforehand, because swapping the actual
    /// data requires reading it (with as()), and this is not possible on
    /// unconstructed variants: it would require some dynamic testing, which
    /// should not be the variant's responsibility.
    /// Swapping between built and (possibly) non-built is done with
    /// self_type::move ().
    template <typename T>
    void
    swap (self_type& that) YY_NOEXCEPT
    {
      std::swap (as<T> (), that.as<T> ());
    }

    /// Move the content of \a that to this.
    ///
    /// Destroys \a that.
    template <typename T>
    void
    move (self_type& that)
    {
# if 201103L <= YY_CPLUSPLUS
      emplace<T> (std::move (that.as<T> ()));
# else
      emplace<T> ();
      swap<T> (that);
# endif
      that.destroy<T> ();
    }

# if 201103L <= YY_CPLUSPLUS
    /// Move the content of \a that to this.
    template <typename T>
    void
    move (self_type&& that)
    {
      emplace<T> (std::move (that.as<T> ()));
      that.destroy<T> ();
    }
#endif

    /// Copy the content of \a that to this.
    template <typename T>
    void
    copy (const self_type& that)
    {
      emplace<T> (that.as<T> ());
    }

    /// Destroy the stored \a T.
    template <typename T>
    void
    destroy ()
    {
      as<T> ().~T ();
    }

  private:
#if YY_CPLUSPLUS < 201103L
    /// Non copyable.
    value_type (const self_type&);
    /// Non copyable.
    self_type& operator= (const self_type&);
#endif

    /// Accessor to raw memory as \a T.
    template <typename T>
    T*
    yyas_ () YY_NOEXCEPT
    {
      void *yyp = yyraw_;
      return static_cast<T*> (yyp);
     }

    /// Const accessor to raw memory as \a T.
    template <typename T>
    const T*
    yyas_ () const YY_NOEXCEPT
    {
      const void *yyp = yyraw_;
      return static_cast<const T*> (yyp);
     }

    /// An auxiliary type to compute the largest semantic type.
    union union_type
    {
      // connections
      // polymer_groups
      // monomer_sequence
      // branch_monomer_group
      char dummy1[sizeof (size_t)];

      // BLOB_ID
      // CHEM_ID
      // PEPTIDE_ID
      // RNA_ID
      // POLYMER_ID
      // SINGLE_CHARACTER_MONOMER
      // MULTI_CHARACTER_MONOMER
      // UNKNOWN_MONOMER
      // MISSING_MONOMER
      // MONOMER_WILDCARD
      // MONOMER_RATIO
      // INLINE_SMILES_MONOMER
      // ANNOTATION
      // REPETITIONS
      // UNKNOWN_SEQUENCE
      // VERSION_TOKEN
      // HYDROGEN_PAIRING
      // RGROUP
      // UNDEFINED_RESIDUE_NUMBER_OR_RGROUP
      // CONNECTION_RESIDUE
      // POLYMER_GROUP_RATIO
      // POLYMER_GROUP_ID
      // EXTENDED_ANNOTATIONS_TOKEN
      // extended_annotations
      // version
      // polymer
      // repetitions
      // monomer_list
      // monomer_and_list
      // monomer_or_list
      // monomer_list_item
      // smiles_monomer
      // monomer_id
      // blob
      // connection_polymer
      // attachment_point
      // connection_monomer
      // residue_and_list
      // residue_or_list
      // polymer_or_list
      // polymer_and_list
      // polymer_group_item
      // annotation
      char dummy2[sizeof (std::string_view)];
    };

    /// The size of the largest semantic type.
    enum { size = sizeof (union_type) };

    /// A buffer to store semantic values.
    union
    {
      /// Strongest alignment constraints.
      long double yyalign_me_;
      /// A buffer large enough to store any of the semantic values.
      char yyraw_[size];
    };
  };

#endif
    /// Backward compatibility (Bison 3.8).
    typedef value_type semantic_type;

    /// Symbol locations.
    typedef location location_type;

    /// Syntax errors thrown from user actions.
    struct syntax_error : std::runtime_error
    {
      syntax_error (const location_type& l, const std::string& m)
        : std::runtime_error (m)
        , location (l)
      {}

      syntax_error (const syntax_error& s)
        : std::runtime_error (s.what ())
        , location (s.location)
      {}

      ~syntax_error () YY_NOEXCEPT YY_NOTHROW;

      location_type location;
    };

    /// Token kinds.
    struct token
    {
      enum token_kind_type
      {
        YYEMPTY = -2,
    END = 0,                       // END
    YYerror = 256,                 // error
    YYUNDEF = 257,                 // "invalid token"
    BLOB_ID = 258,                 // BLOB_ID
    CHEM_ID = 259,                 // CHEM_ID
    PEPTIDE_ID = 260,              // PEPTIDE_ID
    RNA_ID = 261,                  // RNA_ID
    POLYMER_ID = 262,              // POLYMER_ID
    SINGLE_CHARACTER_MONOMER = 263, // SINGLE_CHARACTER_MONOMER
    MULTI_CHARACTER_MONOMER = 264, // MULTI_CHARACTER_MONOMER
    UNKNOWN_MONOMER = 265,         // UNKNOWN_MONOMER
    MISSING_MONOMER = 266,         // MISSING_MONOMER
    MONOMER_WILDCARD = 267,        // MONOMER_WILDCARD
    MONOMER_RATIO = 268,           // MONOMER_RATIO
    INLINE_SMILES_MONOMER = 269,   // INLINE_SMILES_MONOMER
    ANNOTATION = 270,              // ANNOTATION
    REPETITIONS = 271,             // REPETITIONS
    UNKNOWN_SEQUENCE = 272,        // UNKNOWN_SEQUENCE
    VERSION_TOKEN = 273,           // VERSION_TOKEN
    HYDROGEN_PAIRING = 274,        // HYDROGEN_PAIRING
    RGROUP = 275,                  // RGROUP
    UNDEFINED_RESIDUE_NUMBER_OR_RGROUP = 276, // UNDEFINED_RESIDUE_NUMBER_OR_RGROUP
    CONNECTION_RESIDUE = 277,      // CONNECTION_RESIDUE
    POLYMER_GROUP_RATIO = 278,     // POLYMER_GROUP_RATIO
    POLYMER_GROUP_ID = 279,        // POLYMER_GROUP_ID
    EXTENDED_ANNOTATIONS_TOKEN = 280 // EXTENDED_ANNOTATIONS_TOKEN
      };
      /// Backward compatibility alias (Bison 3.6).
      typedef token_kind_type yytokentype;
    };

    /// Token kind, as returned by yylex.
    typedef token::token_kind_type token_kind_type;

    /// Backward compatibility alias (Bison 3.6).
    typedef token_kind_type token_type;

    /// Symbol kinds.
    struct symbol_kind
    {
      enum symbol_kind_type
      {
        YYNTOKENS = 37, ///< Number of tokens.
        S_YYEMPTY = -2,
        S_YYEOF = 0,                             // END
        S_YYerror = 1,                           // error
        S_YYUNDEF = 2,                           // "invalid token"
        S_BLOB_ID = 3,                           // BLOB_ID
        S_CHEM_ID = 4,                           // CHEM_ID
        S_PEPTIDE_ID = 5,                        // PEPTIDE_ID
        S_RNA_ID = 6,                            // RNA_ID
        S_POLYMER_ID = 7,                        // POLYMER_ID
        S_SINGLE_CHARACTER_MONOMER = 8,          // SINGLE_CHARACTER_MONOMER
        S_MULTI_CHARACTER_MONOMER = 9,           // MULTI_CHARACTER_MONOMER
        S_UNKNOWN_MONOMER = 10,                  // UNKNOWN_MONOMER
        S_MISSING_MONOMER = 11,                  // MISSING_MONOMER
        S_MONOMER_WILDCARD = 12,                 // MONOMER_WILDCARD
        S_MONOMER_RATIO = 13,                    // MONOMER_RATIO
        S_INLINE_SMILES_MONOMER = 14,            // INLINE_SMILES_MONOMER
        S_ANNOTATION = 15,                       // ANNOTATION
        S_REPETITIONS = 16,                      // REPETITIONS
        S_UNKNOWN_SEQUENCE = 17,                 // UNKNOWN_SEQUENCE
        S_VERSION_TOKEN = 18,                    // VERSION_TOKEN
        S_HYDROGEN_PAIRING = 19,                 // HYDROGEN_PAIRING
        S_RGROUP = 20,                           // RGROUP
        S_UNDEFINED_RESIDUE_NUMBER_OR_RGROUP = 21, // UNDEFINED_RESIDUE_NUMBER_OR_RGROUP
        S_CONNECTION_RESIDUE = 22,               // CONNECTION_RESIDUE
        S_POLYMER_GROUP_RATIO = 23,              // POLYMER_GROUP_RATIO
        S_POLYMER_GROUP_ID = 24,                 // POLYMER_GROUP_ID
        S_EXTENDED_ANNOTATIONS_TOKEN = 25,       // EXTENDED_ANNOTATIONS_TOKEN
        S_26_ = 26,                              // '$'
        S_27_ = 27,                              // '|'
        S_28_ = 28,                              // '{'
        S_29_ = 29,                              // '}'
        S_30_ = 30,                              // '.'
        S_31_ = 31,                              // '('
        S_32_ = 32,                              // ')'
        S_33_ = 33,                              // '+'
        S_34_ = 34,                              // ','
        S_35_ = 35,                              // ':'
        S_36_ = 36,                              // '-'
        S_YYACCEPT = 37,                         // $accept
        S_helm = 38,                             // helm
        S_polymers = 39,                         // polymers
        S_connections = 40,                      // connections
        S_polymer_groups = 41,                   // polymer_groups
        S_extended_annotations = 42,             // extended_annotations
        S_version = 43,                          // version
        S_polymer_unit = 44,                     // polymer_unit
        S_polymer = 45,                          // polymer
        S_monomers = 46,                         // monomers
        S_monomer_group = 47,                    // monomer_group
        S_repeated_monomers = 48,                // repeated_monomers
        S_repetitions = 49,                      // repetitions
        S_monomer_sequence = 50,                 // monomer_sequence
        S_branch_monomer_group = 51,             // branch_monomer_group
        S_monomer_unit = 52,                     // monomer_unit
        S_monomer_list = 53,                     // monomer_list
        S_monomer_and_list = 54,                 // monomer_and_list
        S_monomer_or_list = 55,                  // monomer_or_list
        S_monomer_list_item = 56,                // monomer_list_item
        S_smiles_monomer = 57,                   // smiles_monomer
        S_monomer_id = 58,                       // monomer_id
        S_blob = 59,                             // blob
        S_connection = 60,                       // connection
        S_connection_polymer = 61,               // connection_polymer
        S_attachment_point = 62,                 // attachment_point
        S_connection_monomer = 63,               // connection_monomer
        S_residue_and_list = 64,                 // residue_and_list
        S_residue_or_list = 65,                  // residue_or_list
        S_polymer_group = 66,                    // polymer_group
        S_polymer_or_list = 67,                  // polymer_or_list
        S_polymer_and_list = 68,                 // polymer_and_list
        S_polymer_group_item = 69,               // polymer_group_item
        S_annotation = 70                        // annotation
      };
    };

    /// (Internal) symbol kind.
    typedef symbol_kind::symbol_kind_type symbol_kind_type;

    /// The number of tokens.
    static const symbol_kind_type YYNTOKENS = symbol_kind::YYNTOKENS;

    /// A complete symbol.
    ///
    /// Expects its Base type to provide access to the symbol kind
    /// via kind ().
    ///
    /// Provide access to semantic value and location.
    template <typename Base>
    struct basic_symbol : Base
    {
      /// Alias to Base.
      typedef Base super_type;

      /// Default constructor.
      basic_symbol () YY_NOEXCEPT
        : value ()
        , location ()
      {}

#if 201103L <= YY_CPLUSPLUS
      /// Move constructor.
      basic_symbol (basic_symbol&& that)
        : Base (std::move (that))
        , value ()
        , location (std::move (that.location))
      {
        switch (this->kind ())
    {
      case symbol_kind::S_connections: // connections
      case symbol_kind::S_polymer_groups: // polymer_groups
      case symbol_kind::S_monomer_sequence: // monomer_sequence
      case symbol_kind::S_branch_monomer_group: // branch_monomer_group
        value.move< size_t > (std::move (that.value));
        break;

      case symbol_kind::S_BLOB_ID: // BLOB_ID
      case symbol_kind::S_CHEM_ID: // CHEM_ID
      case symbol_kind::S_PEPTIDE_ID: // PEPTIDE_ID
      case symbol_kind::S_RNA_ID: // RNA_ID
      case symbol_kind::S_POLYMER_ID: // POLYMER_ID
      case symbol_kind::S_SINGLE_CHARACTER_MONOMER: // SINGLE_CHARACTER_MONOMER
      case symbol_kind::S_MULTI_CHARACTER_MONOMER: // MULTI_CHARACTER_MONOMER
      case symbol_kind::S_UNKNOWN_MONOMER: // UNKNOWN_MONOMER
      case symbol_kind::S_MISSING_MONOMER: // MISSING_MONOMER
      case symbol_kind::S_MONOMER_WILDCARD: // MONOMER_WILDCARD
      case symbol_kind::S_MONOMER_RATIO: // MONOMER_RATIO
      case symbol_kind::S_INLINE_SMILES_MONOMER: // INLINE_SMILES_MONOMER
      case symbol_kind::S_ANNOTATION: // ANNOTATION
      case symbol_kind::S_REPETITIONS: // REPETITIONS
      case symbol_kind::S_UNKNOWN_SEQUENCE: // UNKNOWN_SEQUENCE
      case symbol_kind::S_VERSION_TOKEN: // VERSION_TOKEN
      case symbol_kind::S_HYDROGEN_PAIRING: // HYDROGEN_PAIRING
      case symbol_kind::S_RGROUP: // RGROUP
      case symbol_kind::S_UNDEFINED_RESIDUE_NUMBER_OR_RGROUP: // UNDEFINED_RESIDUE_NUMBER_OR_RGROUP
      case symbol_kind::S_CONNECTION_RESIDUE: // CONNECTION_RESIDUE
      case symbol_kind::S_POLYMER_GROUP_RATIO: // POLYMER_GROUP_RATIO
      case symbol_kind::S_POLYMER_GROUP_ID: // POLYMER_GROUP_ID
      case symbol_kind::S_EXTENDED_ANNOTATIONS_TOKEN: // EXTENDED_ANNOTATIONS_TOKEN
      case symbol_kind::S_extended_annotations: // extended_annotations
      case symbol_kind::S_version: // version
      case symbol_kind::S_polymer: // polymer
      case symbol_kind::S_repetitions: // repetitions
      case symbol_kind::S_monomer_list: // monomer_list
      case symbol_kind::S_monomer_and_list: // monomer_and_list
      case symbol_kind::S_monomer_or_list: // monomer_or_list
      case symbol_kind::S_monomer_list_item: // monomer_list_item
      case symbol_kind::S_smiles_monomer: // smiles_monomer
      case symbol_kind::S_monomer_id: // monomer_id
      case symbol_kind::S_blob: // blob
      case symbol_kind::S_connection_polymer: // connection_polymer
      case symbol_kind::S_attachment_point: // attachment_point
      case symbol_kind::S_connection_monomer: // connection_monomer
      case symbol_kind::S_residue_and_list: // residue_and_list
      case symbol_kind::S_residue_or_list: // residue_or_list
      case symbol_kind::S_polymer_or_list: // polymer_or_list
      case symbol_kind::S_polymer_and_list: // polymer_and_list
      case symbol_kind::S_polymer_group_item: // polymer_group_item
      case symbol_kind::S_annotation: // annotation
        value.move< std::string_view > (std::move (that.value));
        break;

      default:
        break;
    }

      }
#endif

      /// Copy constructor.
      basic_symbol (const basic_symbol& that);

      /// Constructors for typed symbols.
#if 201103L <= YY_CPLUSPLUS
      basic_symbol (typename Base::kind_type t, location_type&& l)
        : Base (t)
        , location (std::move (l))
      {}
#else
      basic_symbol (typename Base::kind_type t, const location_type& l)
        : Base (t)
        , location (l)
      {}
#endif

#if 201103L <= YY_CPLUSPLUS
      basic_symbol (typename Base::kind_type t, size_t&& v, location_type&& l)
        : Base (t)
        , value (std::move (v))
        , location (std::move (l))
      {}
#else
      basic_symbol (typename Base::kind_type t, const size_t& v, const location_type& l)
        : Base (t)
        , value (v)
        , location (l)
      {}
#endif

#if 201103L <= YY_CPLUSPLUS
      basic_symbol (typename Base::kind_type t, std::string_view&& v, location_type&& l)
        : Base (t)
        , value (std::move (v))
        , location (std::move (l))
      {}
#else
      basic_symbol (typename Base::kind_type t, const std::string_view& v, const location_type& l)
        : Base (t)
        , value (v)
        , location (l)
      {}
#endif

      /// Destroy the symbol.
      ~basic_symbol ()
      {
        clear ();
      }



      /// Destroy contents, and record that is empty.
      void clear () YY_NOEXCEPT
      {
        // User destructor.
        symbol_kind_type yykind = this->kind ();
        basic_symbol<Base>& yysym = *this;
        (void) yysym;
        switch (yykind)
        {
       default:
          break;
        }

        // Value type destructor.
switch (yykind)
    {
      case symbol_kind::S_connections: // connections
      case symbol_kind::S_polymer_groups: // polymer_groups
      case symbol_kind::S_monomer_sequence: // monomer_sequence
      case symbol_kind::S_branch_monomer_group: // branch_monomer_group
        value.template destroy< size_t > ();
        break;

      case symbol_kind::S_BLOB_ID: // BLOB_ID
      case symbol_kind::S_CHEM_ID: // CHEM_ID
      case symbol_kind::S_PEPTIDE_ID: // PEPTIDE_ID
      case symbol_kind::S_RNA_ID: // RNA_ID
      case symbol_kind::S_POLYMER_ID: // POLYMER_ID
      case symbol_kind::S_SINGLE_CHARACTER_MONOMER: // SINGLE_CHARACTER_MONOMER
      case symbol_kind::S_MULTI_CHARACTER_MONOMER: // MULTI_CHARACTER_MONOMER
      case symbol_kind::S_UNKNOWN_MONOMER: // UNKNOWN_MONOMER
      case symbol_kind::S_MISSING_MONOMER: // MISSING_MONOMER
      case symbol_kind::S_MONOMER_WILDCARD: // MONOMER_WILDCARD
      case symbol_kind::S_MONOMER_RATIO: // MONOMER_RATIO
      case symbol_kind::S_INLINE_SMILES_MONOMER: // INLINE_SMILES_MONOMER
      case symbol_kind::S_ANNOTATION: // ANNOTATION
      case symbol_kind::S_REPETITIONS: // REPETITIONS
      case symbol_kind::S_UNKNOWN_SEQUENCE: // UNKNOWN_SEQUENCE
      case symbol_kind::S_VERSION_TOKEN: // VERSION_TOKEN
      case symbol_kind::S_HYDROGEN_PAIRING: // HYDROGEN_PAIRING
      case symbol_kind::S_RGROUP: // RGROUP
      case symbol_kind::S_UNDEFINED_RESIDUE_NUMBER_OR_RGROUP: // UNDEFINED_RESIDUE_NUMBER_OR_RGROUP
      case symbol_kind::S_CONNECTION_RESIDUE: // CONNECTION_RESIDUE
      case symbol_kind::S_POLYMER_GROUP_RATIO: // POLYMER_GROUP_RATIO
      case symbol_kind::S_POLYMER_GROUP_ID: // POLYMER_GROUP_ID
      case symbol_kind::S_EXTENDED_ANNOTATIONS_TOKEN: // EXTENDED_ANNOTATIONS_TOKEN
      case symbol_kind::S_extended_annotations: // extended_annotations
      case symbol_kind::S_version: // version
      case symbol_kind::S_polymer: // polymer
      case symbol_kind::S_repetitions: // repetitions
      case symbol_kind::S_monomer_list: // monomer_list
      case symbol_kind::S_monomer_and_list: // monomer_and_list
      case symbol_kind::S_monomer_or_list: // monomer_or_list
      case symbol_kind::S_monomer_list_item: // monomer_list_item
      case symbol_kind::S_smiles_monomer: // smiles_monomer
      case symbol_kind::S_monomer_id: // monomer_id
      case symbol_kind::S_blob: // blob
      case symbol_kind::S_connection_polymer: // connection_polymer
      case symbol_kind::S_attachment_point: // attachment_point
      case symbol_kind::S_connection_monomer: // connection_monomer
      case symbol_kind::S_residue_and_list: // residue_and_list
      case symbol_kind::S_residue_or_list: // residue_or_list
      case symbol_kind::S_polymer_or_list: // polymer_or_list
      case symbol_kind::S_polymer_and_list: // polymer_and_list
      case symbol_kind::S_polymer_group_item: // polymer_group_item
      case symbol_kind::S_annotation: // annotation
        value.template destroy< std::string_view > ();
        break;

      default:
        break;
    }

        Base::clear ();
      }

#if YYDEBUG || 0
      /// The user-facing name of this symbol.
      const char *name () const YY_NOEXCEPT
      {
        return TokenParser::symbol_name (this->kind ());
      }
#endif // #if YYDEBUG || 0


      /// Backward compatibility (Bison 3.6).
      symbol_kind_type type_get () const YY_NOEXCEPT;

      /// Whether empty.
      bool empty () const YY_NOEXCEPT;

      /// Destructive move, \a s is emptied into this.
      void move (basic_symbol& s);

      /// The semantic value.
      value_type value;

      /// The location.
      location_type location;

    private:
#if YY_CPLUSPLUS < 201103L
      /// Assignment operator.
      basic_symbol& operator= (const basic_symbol& that);
#endif
    };

    /// Type access provider for token (enum) based symbols.
    struct by_kind
    {
      /// The symbol kind as needed by the constructor.
      typedef token_kind_type kind_type;

      /// Default constructor.
      by_kind () YY_NOEXCEPT;

#if 201103L <= YY_CPLUSPLUS
      /// Move constructor.
      by_kind (by_kind&& that) YY_NOEXCEPT;
#endif

      /// Copy constructor.
      by_kind (const by_kind& that) YY_NOEXCEPT;

      /// Constructor from (external) token numbers.
      by_kind (kind_type t) YY_NOEXCEPT;



      /// Record that this symbol is empty.
      void clear () YY_NOEXCEPT;

      /// Steal the symbol kind from \a that.
      void move (by_kind& that);

      /// The (internal) type number (corresponding to \a type).
      /// \a empty when empty.
      symbol_kind_type kind () const YY_NOEXCEPT;

      /// Backward compatibility (Bison 3.6).
      symbol_kind_type type_get () const YY_NOEXCEPT;

      /// The symbol kind.
      /// \a S_YYEMPTY when empty.
      symbol_kind_type kind_;
    };

    /// Backward compatibility for a private implementation detail (Bison 3.6).
    typedef by_kind by_type;

    /// "External" symbols: returned by the scanner.
    struct symbol_type : basic_symbol<by_kind>
    {
      /// Superclass.
      typedef basic_symbol<by_kind> super_type;

      /// Empty symbol.
      symbol_type () YY_NOEXCEPT {}

      /// Constructor for valueless symbols, and symbols from each type.
#if 201103L <= YY_CPLUSPLUS
      symbol_type (int tok, location_type l)
        : super_type (token_kind_type (tok), std::move (l))
#else
      symbol_type (int tok, const location_type& l)
        : super_type (token_kind_type (tok), l)
#endif
      {}
#if 201103L <= YY_CPLUSPLUS
      symbol_type (int tok, std::string_view v, location_type l)
        : super_type (token_kind_type (tok), std::move (v), std::move (l))
#else
      symbol_type (int tok, const std::string_view& v, const location_type& l)
        : super_type (token_kind_type (tok), v, l)
#endif
      {}
    };

    /// Build a parser object.
    TokenParser (TokenScanner  &scanner_yyarg, HelmParser  &helm_parser_yyarg);
    virtual ~TokenParser ();

#if 201103L <= YY_CPLUSPLUS
    /// Non copyable.
    TokenParser (const TokenParser&) = delete;
    /// Non copyable.
    TokenParser& operator= (const TokenParser&) = delete;
#endif

    /// Parse.  An alias for parse ().
    /// \returns  0 iff parsing succeeded.
    int operator() ();

    /// Parse.
    /// \returns  0 iff parsing succeeded.
    virtual int parse ();

#if YYDEBUG
    /// The current debugging stream.
    std::ostream& debug_stream () const YY_ATTRIBUTE_PURE;
    /// Set the current debugging stream.
    void set_debug_stream (std::ostream &);

    /// Type for debugging levels.
    typedef int debug_level_type;
    /// The current debugging level.
    debug_level_type debug_level () const YY_ATTRIBUTE_PURE;
    /// Set the current debugging level.
    void set_debug_level (debug_level_type l);
#endif

    /// Report a syntax error.
    /// \param loc    where the syntax error is found.
    /// \param msg    a description of the syntax error.
    virtual void error (const location_type& loc, const std::string& msg);

    /// Report a syntax error.
    void error (const syntax_error& err);

#if YYDEBUG || 0
    /// The user-facing name of the symbol whose (internal) number is
    /// YYSYMBOL.  No bounds checking.
    static const char *symbol_name (symbol_kind_type yysymbol);
#endif // #if YYDEBUG || 0


    // Implementation of make_symbol for each token kind.
#if 201103L <= YY_CPLUSPLUS
      static
      symbol_type
      make_END (location_type l)
      {
        return symbol_type (token::END, std::move (l));
      }
#else
      static
      symbol_type
      make_END (const location_type& l)
      {
        return symbol_type (token::END, l);
      }
#endif
#if 201103L <= YY_CPLUSPLUS
      static
      symbol_type
      make_YYerror (location_type l)
      {
        return symbol_type (token::YYerror, std::move (l));
      }
#else
      static
      symbol_type
      make_YYerror (const location_type& l)
      {
        return symbol_type (token::YYerror, l);
      }
#endif
#if 201103L <= YY_CPLUSPLUS
      static
      symbol_type
      make_YYUNDEF (location_type l)
      {
        return symbol_type (token::YYUNDEF, std::move (l));
      }
#else
      static
      symbol_type
      make_YYUNDEF (const location_type& l)
      {
        return symbol_type (token::YYUNDEF, l);
      }
#endif
#if 201103L <= YY_CPLUSPLUS
      static
      symbol_type
      make_BLOB_ID (std::string_view v, location_type l)
      {
        return symbol_type (token::BLOB_ID, std::move (v), std::move (l));
      }
#else
      static
      symbol_type
      make_BLOB_ID (const std::string_view& v, const location_type& l)
      {
        return symbol_type (token::BLOB_ID, v, l);
      }
#endif
#if 201103L <= YY_CPLUSPLUS
      static
      symbol_type
      make_CHEM_ID (std::string_view v, location_type l)
      {
        return symbol_type (token::CHEM_ID, std::move (v), std::move (l));
      }
#else
      static
      symbol_type
      make_CHEM_ID (const std::string_view& v, const location_type& l)
      {
        return symbol_type (token::CHEM_ID, v, l);
      }
#endif
#if 201103L <= YY_CPLUSPLUS
      static
      symbol_type
      make_PEPTIDE_ID (std::string_view v, location_type l)
      {
        return symbol_type (token::PEPTIDE_ID, std::move (v), std::move (l));
      }
#else
      static
      symbol_type
      make_PEPTIDE_ID (const std::string_view& v, const location_type& l)
      {
        return symbol_type (token::PEPTIDE_ID, v, l);
      }
#endif
#if 201103L <= YY_CPLUSPLUS
      static
      symbol_type
      make_RNA_ID (std::string_view v, location_type l)
      {
        return symbol_type (token::RNA_ID, std::move (v), std::move (l));
      }
#else
      static
      symbol_type
      make_RNA_ID (const std::string_view& v, const location_type& l)
      {
        return symbol_type (token::RNA_ID, v, l);
      }
#endif
#if 201103L <= YY_CPLUSPLUS
      static
      symbol_type
      make_POLYMER_ID (std::string_view v, location_type l)
      {
        return symbol_type (token::POLYMER_ID, std::move (v), std::move (l));
      }
#else
      static
      symbol_type
      make_POLYMER_ID (const std::string_view& v, const location_type& l)
      {
        return symbol_type (token::POLYMER_ID, v, l);
      }
#endif
#if 201103L <= YY_CPLUSPLUS
      static
      symbol_type
      make_SINGLE_CHARACTER_MONOMER (std::string_view v, location_type l)
      {
        return symbol_type (token::SINGLE_CHARACTER_MONOMER, std::move (v), std::move (l));
      }
#else
      static
      symbol_type
      make_SINGLE_CHARACTER_MONOMER (const std::string_view& v, const location_type& l)
      {
        return symbol_type (token::SINGLE_CHARACTER_MONOMER, v, l);
      }
#endif
#if 201103L <= YY_CPLUSPLUS
      static
      symbol_type
      make_MULTI_CHARACTER_MONOMER (std::string_view v, location_type l)
      {
        return symbol_type (token::MULTI_CHARACTER_MONOMER, std::move (v), std::move (l));
      }
#else
      static
      symbol_type
      make_MULTI_CHARACTER_MONOMER (const std::string_view& v, const location_type& l)
      {
        return symbol_type (token::MULTI_CHARACTER_MONOMER, v, l);
      }
#endif
#if 201103L <= YY_CPLUSPLUS
      static
      symbol_type
      make_UNKNOWN_MONOMER (std::string_view v, location_type l)
      {
        return symbol_type (token::UNKNOWN_MONOMER, std::move (v), std::move (l));
      }
#else
      static
      symbol_type
      make_UNKNOWN_MONOMER (const std::string_view& v, const location_type& l)
      {
        return symbol_type (token::UNKNOWN_MONOMER, v, l);
      }
#endif
#if 201103L <= YY_CPLUSPLUS
      static
      symbol_type
      make_MISSING_MONOMER (std::string_view v, location_type l)
      {
        return symbol_type (token::MISSING_MONOMER, std::move (v), std::move (l));
      }
#else
      static
      symbol_type
      make_MISSING_MONOMER (const std::string_view& v, const location_type& l)
      {
        return symbol_type (token::MISSING_MONOMER, v, l);
      }
#endif
#if 201103L <= YY_CPLUSPLUS
      static
      symbol_type
      make_MONOMER_WILDCARD (std::string_view v, location_type l)
      {
        return symbol_type (token::MONOMER_WILDCARD, std::move (v), std::move (l));
      }
#else
      static
      symbol_type
      make_MONOMER_WILDCARD (const std::string_view& v, const location_type& l)
      {
        return symbol_type (token::MONOMER_WILDCARD, v, l);
      }
#endif
#if 201103L <= YY_CPLUSPLUS
      static
      symbol_type
      make_MONOMER_RATIO (std::string_view v, location_type l)
      {
        return symbol_type (token::MONOMER_RATIO, std::move (v), std::move (l));
      }
#else
      static
      symbol_type
      make_MONOMER_RATIO (const std::string_view& v, const location_type& l)
      {
        return symbol_type (token::MONOMER_RATIO, v, l);
      }
#endif
#if 201103L <= YY_CPLUSPLUS
      static
      symbol_type
      make_INLINE_SMILES_MONOMER (std::string_view v, location_type l)
      {
        return symbol_type (token::INLINE_SMILES_MONOMER, std::move (v), std::move (l));
      }
#else
      static
      symbol_type
      make_INLINE_SMILES_MONOMER (const std::string_view& v, const location_type& l)
      {
        return symbol_type (token::INLINE_SMILES_MONOMER, v, l);
      }
#endif
#if 201103L <= YY_CPLUSPLUS
      static
      symbol_type
      make_ANNOTATION (std::string_view v, location_type l)
      {
        return symbol_type (token::ANNOTATION, std::move (v), std::move (l));
      }
#else
      static
      symbol_type
      make_ANNOTATION (const std::string_view& v, const location_type& l)
      {
        return symbol_type (token::ANNOTATION, v, l);
      }
#endif
#if 201103L <= YY_CPLUSPLUS
      static
      symbol_type
      make_REPETITIONS (std::string_view v, location_type l)
      {
        return symbol_type (token::REPETITIONS, std::move (v), std::move (l));
      }
#else
      static
      symbol_type
      make_REPETITIONS (const std::string_view& v, const location_type& l)
      {
        return symbol_type (token::REPETITIONS, v, l);
      }
#endif
#if 201103L <= YY_CPLUSPLUS
      static
      symbol_type
      make_UNKNOWN_SEQUENCE (std::string_view v, location_type l)
      {
        return symbol_type (token::UNKNOWN_SEQUENCE, std::move (v), std::move (l));
      }
#else
      static
      symbol_type
      make_UNKNOWN_SEQUENCE (const std::string_view& v, const location_type& l)
      {
        return symbol_type (token::UNKNOWN_SEQUENCE, v, l);
      }
#endif
#if 201103L <= YY_CPLUSPLUS
      static
      symbol_type
      make_VERSION_TOKEN (std::string_view v, location_type l)
      {
        return symbol_type (token::VERSION_TOKEN, std::move (v), std::move (l));
      }
#else
      static
      symbol_type
      make_VERSION_TOKEN (const std::string_view& v, const location_type& l)
      {
        return symbol_type (token::VERSION_TOKEN, v, l);
      }
#endif
#if 201103L <= YY_CPLUSPLUS
      static
      symbol_type
      make_HYDROGEN_PAIRING (std::string_view v, location_type l)
      {
        return symbol_type (token::HYDROGEN_PAIRING, std::move (v), std::move (l));
      }
#else
      static
      symbol_type
      make_HYDROGEN_PAIRING (const std::string_view& v, const location_type& l)
      {
        return symbol_type (token::HYDROGEN_PAIRING, v, l);
      }
#endif
#if 201103L <= YY_CPLUSPLUS
      static
      symbol_type
      make_RGROUP (std::string_view v, location_type l)
      {
        return symbol_type (token::RGROUP, std::move (v), std::move (l));
      }
#else
      static
      symbol_type
      make_RGROUP (const std::string_view& v, const location_type& l)
      {
        return symbol_type (token::RGROUP, v, l);
      }
#endif
#if 201103L <= YY_CPLUSPLUS
      static
      symbol_type
      make_UNDEFINED_RESIDUE_NUMBER_OR_RGROUP (std::string_view v, location_type l)
      {
        return symbol_type (token::UNDEFINED_RESIDUE_NUMBER_OR_RGROUP, std::move (v), std::move (l));
      }
#else
      static
      symbol_type
      make_UNDEFINED_RESIDUE_NUMBER_OR_RGROUP (const std::string_view& v, const location_type& l)
      {
        return symbol_type (token::UNDEFINED_RESIDUE_NUMBER_OR_RGROUP, v, l);
      }
#endif
#if 201103L <= YY_CPLUSPLUS
      static
      symbol_type
      make_CONNECTION_RESIDUE (std::string_view v, location_type l)
      {
        return symbol_type (token::CONNECTION_RESIDUE, std::move (v), std::move (l));
      }
#else
      static
      symbol_type
      make_CONNECTION_RESIDUE (const std::string_view& v, const location_type& l)
      {
        return symbol_type (token::CONNECTION_RESIDUE, v, l);
      }
#endif
#if 201103L <= YY_CPLUSPLUS
      static
      symbol_type
      make_POLYMER_GROUP_RATIO (std::string_view v, location_type l)
      {
        return symbol_type (token::POLYMER_GROUP_RATIO, std::move (v), std::move (l));
      }
#else
      static
      symbol_type
      make_POLYMER_GROUP_RATIO (const std::string_view& v, const location_type& l)
      {
        return symbol_type (token::POLYMER_GROUP_RATIO, v, l);
      }
#endif
#if 201103L <= YY_CPLUSPLUS
      static
      symbol_type
      make_POLYMER_GROUP_ID (std::string_view v, location_type l)
      {
        return symbol_type (token::POLYMER_GROUP_ID, std::move (v), std::move (l));
      }
#else
      static
      symbol_type
      make_POLYMER_GROUP_ID (const std::string_view& v, const location_type& l)
      {
        return symbol_type (token::POLYMER_GROUP_ID, v, l);
      }
#endif
#if 201103L <= YY_CPLUSPLUS
      static
      symbol_type
      make_EXTENDED_ANNOTATIONS_TOKEN (std::string_view v, location_type l)
      {
        return symbol_type (token::EXTENDED_ANNOTATIONS_TOKEN, std::move (v), std::move (l));
      }
#else
      static
      symbol_type
      make_EXTENDED_ANNOTATIONS_TOKEN (const std::string_view& v, const location_type& l)
      {
        return symbol_type (token::EXTENDED_ANNOTATIONS_TOKEN, v, l);
      }
#endif


  private:
#if YY_CPLUSPLUS < 201103L
    /// Non copyable.
    TokenParser (const TokenParser&);
    /// Non copyable.
    TokenParser& operator= (const TokenParser&);
#endif


    /// Stored state numbers (used for stacks).
    typedef unsigned char state_type;

    /// Compute post-reduction state.
    /// \param yystate   the current state
    /// \param yysym     the nonterminal to push on the stack
    static state_type yy_lr_goto_state_ (state_type yystate, int yysym);

    /// Whether the given \c yypact_ value indicates a defaulted state.
    /// \param yyvalue   the value to check
    static bool yy_pact_value_is_default_ (int yyvalue) YY_NOEXCEPT;

    /// Whether the given \c yytable_ value indicates a syntax error.
    /// \param yyvalue   the value to check
    static bool yy_table_value_is_error_ (int yyvalue) YY_NOEXCEPT;

    static const signed char yypact_ninf_;
    static const signed char yytable_ninf_;

    /// Convert a scanner token kind \a t to a symbol kind.
    /// In theory \a t should be a token_kind_type, but character literals
    /// are valid, yet not members of the token_kind_type enum.
    static symbol_kind_type yytranslate_ (int t) YY_NOEXCEPT;

#if YYDEBUG || 0
    /// For a symbol, its name in clear.
    static const char* const yytname_[];
#endif // #if YYDEBUG || 0


    // Tables.
    // YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
    // STATE-NUM.
    static const signed char yypact_[];

    // YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
    // Performed when YYTABLE does not specify something else to do.  Zero
    // means the default is an error.
    static const signed char yydefact_[];

    // YYPGOTO[NTERM-NUM].
    static const signed char yypgoto_[];

    // YYDEFGOTO[NTERM-NUM].
    static const unsigned char yydefgoto_[];

    // YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
    // positive, shift that token.  If negative, reduce the rule whose
    // number is the opposite.  If YYTABLE_NINF, syntax error.
    static const short yytable_[];

    static const short yycheck_[];

    // YYSTOS[STATE-NUM] -- The symbol kind of the accessing symbol of
    // state STATE-NUM.
    static const signed char yystos_[];

    // YYR1[RULE-NUM] -- Symbol kind of the left-hand side of rule RULE-NUM.
    static const signed char yyr1_[];

    // YYR2[RULE-NUM] -- Number of symbols on the right-hand side of rule RULE-NUM.
    static const signed char yyr2_[];


#if YYDEBUG
    // YYRLINE[YYN] -- Source line where rule number YYN was defined.
    static const short yyrline_[];
    /// Report on the debug stream that the rule \a r is going to be reduced.
    virtual void yy_reduce_print_ (int r) const;
    /// Print the state stack on the debug stream.
    virtual void yy_stack_print_ () const;

    /// Debugging level.
    int yydebug_;
    /// Debug stream.
    std::ostream* yycdebug_;

    /// \brief Display a symbol kind, value and location.
    /// \param yyo    The output stream.
    /// \param yysym  The symbol.
    template <typename Base>
    void yy_print_ (std::ostream& yyo, const basic_symbol<Base>& yysym) const;
#endif

    /// \brief Reclaim the memory associated to a symbol.
    /// \param yymsg     Why this token is reclaimed.
    ///                  If null, print nothing.
    /// \param yysym     The symbol.
    template <typename Base>
    void yy_destroy_ (const char* yymsg, basic_symbol<Base>& yysym) const;

  private:
    /// Type access provider for state based symbols.
    struct by_state
    {
      /// Default constructor.
      by_state () YY_NOEXCEPT;

      /// The symbol kind as needed by the constructor.
      typedef state_type kind_type;

      /// Constructor.
      by_state (kind_type s) YY_NOEXCEPT;

      /// Copy constructor.
      by_state (const by_state& that) YY_NOEXCEPT;

      /// Record that this symbol is empty.
      void clear () YY_NOEXCEPT;

      /// Steal the symbol kind from \a that.
      void move (by_state& that);

      /// The symbol kind (corresponding to \a state).
      /// \a symbol_kind::S_YYEMPTY when empty.
      symbol_kind_type kind () const YY_NOEXCEPT;

      /// The state number used to denote an empty symbol.
      /// We use the initial state, as it does not have a value.
      enum { empty_state = 0 };

      /// The state.
      /// \a empty when empty.
      state_type state;
    };

    /// "Internal" symbol: element of the stack.
    struct stack_symbol_type : basic_symbol<by_state>
    {
      /// Superclass.
      typedef basic_symbol<by_state> super_type;
      /// Construct an empty symbol.
      stack_symbol_type ();
      /// Move or copy construction.
      stack_symbol_type (YY_RVREF (stack_symbol_type) that);
      /// Steal the contents from \a sym to build this.
      stack_symbol_type (state_type s, YY_MOVE_REF (symbol_type) sym);
#if YY_CPLUSPLUS < 201103L
      /// Assignment, needed by push_back by some old implementations.
      /// Moves the contents of that.
      stack_symbol_type& operator= (stack_symbol_type& that);

      /// Assignment, needed by push_back by other implementations.
      /// Needed by some other old implementations.
      stack_symbol_type& operator= (const stack_symbol_type& that);
#endif
    };

    /// A stack with random access from its top.
    template <typename T, typename S = std::vector<T> >
    class stack
    {
    public:
      // Hide our reversed order.
      typedef typename S::iterator iterator;
      typedef typename S::const_iterator const_iterator;
      typedef typename S::size_type size_type;
      typedef typename std::ptrdiff_t index_type;

      stack (size_type n = 200) YY_NOEXCEPT
        : seq_ (n)
      {}

#if 201103L <= YY_CPLUSPLUS
      /// Non copyable.
      stack (const stack&) = delete;
      /// Non copyable.
      stack& operator= (const stack&) = delete;
#endif

      /// Random access.
      ///
      /// Index 0 returns the topmost element.
      const T&
      operator[] (index_type i) const
      {
        return seq_[size_type (size () - 1 - i)];
      }

      /// Random access.
      ///
      /// Index 0 returns the topmost element.
      T&
      operator[] (index_type i)
      {
        return seq_[size_type (size () - 1 - i)];
      }

      /// Steal the contents of \a t.
      ///
      /// Close to move-semantics.
      void
      push (YY_MOVE_REF (T) t)
      {
        seq_.push_back (T ());
        operator[] (0).move (t);
      }

      /// Pop elements from the stack.
      void
      pop (std::ptrdiff_t n = 1) YY_NOEXCEPT
      {
        for (; 0 < n; --n)
          seq_.pop_back ();
      }

      /// Pop all elements from the stack.
      void
      clear () YY_NOEXCEPT
      {
        seq_.clear ();
      }

      /// Number of elements on the stack.
      index_type
      size () const YY_NOEXCEPT
      {
        return index_type (seq_.size ());
      }

      /// Iterator on top of the stack (going downwards).
      const_iterator
      begin () const YY_NOEXCEPT
      {
        return seq_.begin ();
      }

      /// Bottom of the stack.
      const_iterator
      end () const YY_NOEXCEPT
      {
        return seq_.end ();
      }

      /// Present a slice of the top of a stack.
      class slice
      {
      public:
        slice (const stack& stack, index_type range) YY_NOEXCEPT
          : stack_ (stack)
          , range_ (range)
        {}

        const T&
        operator[] (index_type i) const
        {
          return stack_[range_ - i];
        }

      private:
        const stack& stack_;
        index_type range_;
      };

    private:
#if YY_CPLUSPLUS < 201103L
      /// Non copyable.
      stack (const stack&);
      /// Non copyable.
      stack& operator= (const stack&);
#endif
      /// The wrapped container.
      S seq_;
    };


    /// Stack type.
    typedef stack<stack_symbol_type> stack_type;

    /// The stack.
    stack_type yystack_;

    /// Push a new state on the stack.
    /// \param m    a debug message to display
    ///             if null, no trace is output.
    /// \param sym  the symbol
    /// \warning the contents of \a s.value is stolen.
    void yypush_ (const char* m, YY_MOVE_REF (stack_symbol_type) sym);

    /// Push a new look ahead token on the state on the stack.
    /// \param m    a debug message to display
    ///             if null, no trace is output.
    /// \param s    the state
    /// \param sym  the symbol (for its value and location).
    /// \warning the contents of \a sym.value is stolen.
    void yypush_ (const char* m, state_type s, YY_MOVE_REF (symbol_type) sym);

    /// Pop \a n symbols from the stack.
    void yypop_ (int n = 1) YY_NOEXCEPT;

    /// Constants.
    enum
    {
      yylast_ = 151,     ///< Last index in yytable_.
      yynnts_ = 34,  ///< Number of nonterminal symbols.
      yyfinal_ = 13 ///< Termination state number.
    };


    // User arguments.
    TokenScanner  &scanner;
    HelmParser  &helm_parser;

  };


#line 15 "../helm_parser.yy"
} // helm
#line 1686 "helm_parser.tab.hh"




#endif // !YY_YY_HELM_PARSER_TAB_HH_INCLUDED
