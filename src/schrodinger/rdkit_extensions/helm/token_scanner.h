#pragma once

#if !defined(yyFlexLexerOnce)
#include <FlexLexer.h>
#endif

#include "schrodinger/rdkit_extensions/helm/generated/helm_parser.tab.hh"
#include "schrodinger/rdkit_extensions/helm/generated/location.hh"

namespace helm
{

class TokenScanner : public yyFlexLexer
{
  public:
    TokenScanner(std::istream* in) : yyFlexLexer(in){};
    virtual ~TokenScanner(){};

    int lex(helm::TokenParser::semantic_type* const lval,
            helm::TokenParser::location_type* location);

  private:
    helm::TokenParser::semantic_type* yylval = nullptr;
};

} /* end namespace helm */
