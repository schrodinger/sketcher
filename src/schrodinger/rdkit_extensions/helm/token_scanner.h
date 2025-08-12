#pragma once

#include <string_view>
#if !defined(yyFlexLexerOnce)
#include "schrodinger/rdkit_extensions/helm/thirdparty/FlexLexer.h"
#endif

#include "schrodinger/rdkit_extensions/helm/generated/helm_parser.tab.hh"
#include "schrodinger/rdkit_extensions/helm/generated/location.hh"

namespace helm
{

class TokenScanner : public yyFlexLexer
{
  public:
    TokenScanner(std::istream* in, const std::string_view input) :
        yyFlexLexer(in),
        ref_string_view(input){};
    virtual ~TokenScanner(){};

    int lex(helm::TokenParser::semantic_type* const lval,
            helm::TokenParser::location_type* location);

    std::string_view ref_string_view;

  private:
    helm::TokenParser::semantic_type* yylval = nullptr;
};

} /* end namespace helm */
