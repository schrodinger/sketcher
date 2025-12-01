#include "schrodinger/rdkit_extensions/helm/helm_parser.h"

#include <array>
#include <fmt/format.h>
#include <fmt/ranges.h>

#include <rdkit/GraphMol/RWMol.h>
#include "schrodinger/rdkit_extensions/helm/validation.h"

namespace helm
{

namespace detail
{

/**
 * @brief Checks if a character is common whitespace
 * @note std::isspace is not constexpr
 */
constexpr bool isspace(unsigned char c) noexcept
{
    std::string_view ws{" \t\r\f\n\v"};
    return ws.find(c) != std::string_view::npos;
}

/**
 * @brief Checks if a character is a decimal digit (0-9)
 * @note std::isdigit is not constexpr
 */
constexpr bool isdigit(unsigned char c) noexcept
{
    return '0' <= c && c <= '9';
}

/**
 * @brief Standard error handler that throws a std::invalid_argument exception.
 *        Suitable for contexts where the parser must abort on the first error.
 */
struct [[nodiscard]] throw_error_handler_t {
    void operator()(std::string_view msg) const
    {
        throw std::invalid_argument(msg.data());
    }

    void operator()(std::string_view msg, size_t pos,
                    std::string_view input) const
    {
        throw std::invalid_argument(construct_error_msg(msg, pos, input));
    }

    void operator()(std::string_view msg, std::string_view bad_token,
                    std::string_view input) const
    {
        operator()(msg, std::distance(input.data(), bad_token.data()), input);
    }
};

/**
 * @brief Error handler that logs errors to stderr but allows the parser to
 *        attempt continuation (though the main parser logic may still fail).
 */
struct [[nodiscard]] log_error_handler_t {
    void operator()(std::string_view msg) const
    {
        std::cerr << msg;
    }

    void operator()(std::string_view msg, size_t pos,
                    std::string_view input) const
    {
        std::cerr << construct_error_msg(msg, pos, input);
    }

    void operator()(std::string_view msg, std::string_view bad_token,
                    std::string_view input) const
    {
        operator()(msg, std::distance(input.data(), bad_token.data()), input);
    }
};

/**
 * @brief The Lexer (Tokenizer) for the HELM format. It tokenizes the input
 *        string based on the delimiters and quoted/grouped literals.
 * @tparam error_handler_t The error handling policy to use (e.g., throw, log).
 */
template <class error_handler_t> class [[nodiscard]] lexer_t
{
  public:
    /**
     * @brief Constructs the lexer with the input HELM string
     * @param source The full HELM input
     */
    explicit lexer_t(std::string_view source) :
        m_input(source),
        m_current_pos(0)
    {
        // discard leading white space
        while (!eof() && detail::isspace(peek())) {
            ++m_current_pos;
        }
    }

    /**
     * @brief Advances the lexer position by the given number of characters
     * @param num_chars The number of characters to advance the position by.
     */
    constexpr void advance(size_t num_chars) noexcept
    {
        auto new_pos = num_chars + pos();
        m_current_pos = (new_pos > size() ? size() : new_pos);
    }

    /**
     * @brief Checks if the next character is `c` and consumes it if found.
     * @return True if consumed, false otherwise.
     */
    [[nodiscard]] bool expect(unsigned char c)
    {
        if (peek() == c) {
            ++m_current_pos;
            return true;
        }

        auto err_msg = fmt::format("Expected character '{}' here", c);
        m_error_handler(err_msg, pos(), input());
        return false;
    }

    [[nodiscard]] const std::string_view& input() const noexcept
    {
        return m_input;
    }

    [[nodiscard]] std::string_view remainder() const noexcept
    {
        return m_input.substr(m_current_pos);
    }

    [[nodiscard]] constexpr size_t size() const noexcept
    {
        return m_input.size();
    }

    [[nodiscard]] constexpr size_t pos() const noexcept
    {
        return m_current_pos;
    }

    [[nodiscard]] constexpr bool eof() const noexcept
    {
        return m_current_pos >= size();
    }

    [[nodiscard]] constexpr unsigned char peek() const noexcept
    {
        return eof() ? '\0' : m_input[m_current_pos];
    }

    /**
     * @brief Parses the next token from the input.
     * @param token Output for the parsed token
     * @return True if parsing was successful, False otherwise
     */
    [[nodiscard]] constexpr bool get_next_token(std::string_view& token)
    {
        // unlikely, but we should handle this
        if (eof()) {
            m_error_handler("unexpected end of string.", m_current_pos,
                            m_input);
            return false;
        }

        const size_t N = m_input.size();
        auto end_pos = N;
        auto seen_quote = false;
        for (size_t i = m_current_pos; i < N; ++i) {
            switch (m_input[i]) {
                case '-':
                case '.':
                case ',':
                case '(':
                case ')':
                case '[':
                case ']':
                case '|':
                case '$':
                case '{':
                case '}':
                    if (seen_quote) {
                        continue;
                    }
                    end_pos = i == m_current_pos ? i + 1 : i;
                    break;
                case '"':  // to capture inline annotations
                case '\'': // to capture monomer repetitions
                    if (i == m_current_pos) {
                        seen_quote = true;
                        continue;
                    } else if (seen_quote &&
                               m_input[i] == m_input[m_current_pos]) {
                        seen_quote = false;
                        end_pos = i + 1;
                        break;
                    } else {
                        end_pos = i == m_current_pos ? i + 1 : i;
                        break;
                    }
                default:
                    continue;
            }
            break;
        }

        // unclosed quotes should be rejected
        if (seen_quote) {
            m_error_handler("Unclosed quote", m_current_pos, m_input);
            return false;
        }

        // construct the token
        const auto length = end_pos - m_current_pos;
        token = m_input.substr(m_current_pos, length);
        m_current_pos += length;

        return true;
    }

    /**
     * @brief Parses the next grouped token from the input. A grouped
     *        token is defined as a token with unique opening and closing
     *        characters. Tokens can be nested.
     * @param token Output for the parsed token
     * @param group_start The group open character
     * @param group_end The group close character
     * @return True if parsing was successful, False otherwise
     */
    [[nodiscard]] constexpr bool get_grouped_token(std::string_view& token,
                                                   unsigned char group_start,
                                                   unsigned char group_end)
    {
        // Assume group start character is at the beginning
        auto depth = 0;
        for (size_t i = m_current_pos; i < m_input.size(); ++i) {
            if (m_input[i] == group_start) {
                ++depth;
            } else if (m_input[i] == group_end) {
                --depth;
                if (depth == 0) {
                    // construct the token
                    const auto length = i - m_current_pos + 1;
                    token = m_input.substr(m_current_pos, length);
                    m_current_pos += length;

                    return true;
                }
            }
        }

        m_error_handler("Unclosed group token", m_current_pos, m_input);
        return false;
    }

  private:
    const std::string_view m_input;
    size_t m_current_pos;
    error_handler_t m_error_handler;
};

/**
 * @brief The Parser for the HELM format. It consumes the tokens from the lexer
 *        and builds the Abstract Syntax Tree (AST) based on the HELM grammar
 *        rules.
 * @tparam error_handler_t The error handling policy to use (e.g., throw, log).
 */
template <class error_handler_t> class [[nodiscard]] parser_t
{
  public:
    parser_t() = default;

    /**
     * @brief Main entry point for parsing a HELM string.
     *        The HELM format is structured into five main sections, separated
     *        by the `$` delimiter:
     *        1. Simple Polymers
     *        2. Connections
     *        3. Polymer Groups
     *        4. Extended Annotations
     *        5. HELM Version
     * @param input The input HELM value
     * @return An optional AST structure if parsing is successful.
     */
    [[nodiscard]] std::optional<helm_info> parse(std::string_view input)
    {
        helm_info result;
        lexer_t<error_handler_t> lexer(input);

        // Enforce the strict sequential grammar for the HELM sections
        if (!(parse_simple_polymers(lexer, result.polymers) &&      //
              lexer.expect('$') &&                                  //
              parse_connections(lexer, result.connections) &&       //
              lexer.expect('$') &&                                  //
              parse_polymer_groups(lexer, result.polymer_groups) && //
              lexer.expect('$') &&                                  //
              parse_extended_annotations(lexer,
                                         result.extended_annotations) && //
              lexer.expect('$') &&                                       //
              parse_helm_version(lexer, result.helm_version))) {
            return std::optional<helm_info>{};
        }

        // Final check: ensure the entire input has been consumed
        if (!lexer.eof()) {
            m_error_handler("Unexpected tokens after HELM version string",
                            lexer.pos(), lexer.input());
            return std::optional<helm_info>{};
        }

        // at this point, the input has an acceptable syntax, but we want to
        // check for unsupported  and /or incompatible features
        std::vector<std::string> validation_errors;
        if (!validate_parsed_info(result, validation_errors, lexer.input())) {
            auto msg =
                fmt::format("Parsing failed because of the following: \n{}\n",
                            fmt::join(validation_errors, "\n"));
            m_error_handler(msg);
            return std::optional<helm_info>{};
        }

        return result;
    }

  private:
    error_handler_t m_error_handler;

    /**
     * @brief Parses a polymer ID, enforcing the prefix (e.g., PEPTIDE) and
     *        the required non-zero numeric suffix (e.g., 1).
     * @param lexer The lexer instance.
     * @param polymer_id Output for the full ID (e.g., PEPTIDE1)
     * @return True if a valid polymer ID was found and consumed.
     */
    [[nodiscard]] bool parse_polymer_id(lexer_t<error_handler_t>& lexer,
                                        std::string_view& polymer_id)
    {
        // Define all supported polymer types
        constexpr std::array<std::string_view, 4> supported_polymers{
            "BLOB", "CHEM", "PEPTIDE", "RNA"};

        auto remainder = lexer.remainder();
        auto matched_polymer =
            std::ranges::find_if(supported_polymers, [&](auto& entry) {
                return remainder.starts_with(entry);
            });

        if (matched_polymer == supported_polymers.end()) {
            m_error_handler("Expected polymer id here. Only PEPTIDE, RNA, "
                            "CHEM, and BLOB polymers are currently supported",
                            lexer.pos(), lexer.input());
            return false;
        }

        // Enforce numeric suffix
        auto suffix_start = matched_polymer->size();
        if (suffix_start == remainder.size() ||
            !detail::isdigit(remainder[suffix_start])) {
            m_error_handler(
                "Polymer ids are expected to have numeric suffixes. e.g., RNA1",
                lexer.pos(), lexer.input());
            return false;
        }

        // Reject suffixes with leading zeros.
        if (remainder[suffix_start] == '0') {
            m_error_handler("Polymer ids shouldn't have leading zeros.",
                            lexer.pos(), lexer.input());
            return false;
        }

        // construct the polymer ID
        auto suffix_end =
            remainder.find_first_not_of("0123456789", suffix_start);
        polymer_id = remainder.substr(0, suffix_end);

        lexer.advance(suffix_end); // advance the lexer position
        return true;
    }

    /**
     * @brief Attempts to parse and set an annotation string.
     * @param lexer The lexer instance
     * @param annotation Output for the annotation string
     * @return True if a annotation block was found and consumed
     */
    [[nodiscard]] bool parse_annotation(lexer_t<error_handler_t>& lexer,
                                        std::string_view& annotation)
    {
        std::string_view value;
        if (!parse_quoted_literal(lexer, value, '"')) {
            return false;
        }

        // empty values shouldn't be allowed
        if (value.size() == 2) {
            m_error_handler("Empty annotations are not supported", value,
                            lexer.input());
            return false;
        }

        annotation = value.substr(1, value.size() - 2); // strip quotes
        return true;
    }

    /**
     * @brief Parses a quoted literal, stripping the quotes. Used for parsing
     *        both repetitions('...') and general annotations("...").
     * @param lexer The lexer instance
     * @param value Output for the quoted value
     * @param quote The quoting character
     * @return True if a quoted literal was found and consumed
     */
    [[nodiscard]] bool parse_quoted_literal(lexer_t<error_handler_t>& lexer,
                                            std::string_view& value,
                                            unsigned char quote)
    {
        auto remainder = lexer.remainder();
        if (!lexer.expect(quote)) {
            return false;
        }

        // A quoted token must have a closing quote.
        auto end_pos = remainder.find(quote, 1);
        if (end_pos == std::string_view::npos) {
            m_error_handler("Unclosed quote found here", remainder,
                            lexer.input());
            return false;
        }

        value = remainder.substr(0, end_pos + 1); // construct value
        lexer.advance(end_pos); // not end_pos + 1 because of the ::expect call

        return true;
    }

    /** @brief Parses the Simple Polymers section (Section 1) */
    [[nodiscard]] bool parse_simple_polymers(lexer_t<error_handler_t>& lexer,
                                             std::vector<polymer>& polymers)
    {
        while (!lexer.eof()) {
            if (!parse_simple_polymer(lexer, polymers)) {
                return false;
            }

            if (lexer.eof() || lexer.peek() == '$') {
                break;
            }

            if (!lexer.expect('|')) {
                return false;
            }
        }

        return !polymers.empty();
    }

    /**
     * @brief Parse a Simple Polymer token.
     * @param lexer The lexer instance.
     * @param polymers Output to add parsed Simple Polymer AST to.
     * @return True if successfully parsed the Simple Polymer.
     */
    [[nodiscard]] bool parse_simple_polymer(lexer_t<error_handler_t>& lexer,
                                            std::vector<polymer>& polymers)
    {
        polymer result;

        // Enforce the strict sequential grammar for the Simple Polymer section
        if (!(parse_polymer_id(lexer, result.id) &&                //
              lexer.expect('{') &&                                 //
              parse_monomers(lexer, result.id, result.monomers) && //
              lexer.expect('}') &&
              // polymer annotations are optional
              (lexer.peek() != '"' ||
               parse_annotation(lexer, result.annotation)))) {
            return false;
        }

        polymers.push_back(std::move(result));
        return true;
    }

    /**
     * @brief Parse list of monomers in a Simple Polymer token.
     * @param lexer The lexer instance.
     * @param polymer_id The ID of the current polymer.
     * @param monomers Output to add parsed monomers to.
     * @return True if successfully parsed the list of monomers.
     */
    [[nodiscard]] bool parse_monomers(lexer_t<error_handler_t>& lexer,
                                      const std::string_view& polymer_id,
                                      std::vector<monomer>& monomers)
    {
        if (polymer_id.starts_with("BLOB")) { // needs special handling
            return parse_blob_monomer(lexer, monomers);
        }

        while (!lexer.eof() && lexer.peek() != '}' /* end of monomers */) {
            if (!monomers.empty() && !lexer.expect('.') &&
                !monomers.back().is_branch) {
                m_error_handler("Expected character '.' before this token.",
                                lexer.pos(), lexer.input());
                return false;
            }

            auto monomer = parse_monomer_unit(lexer);
            if (!monomer) { // we failed to get monomer
                return false;
            }

            monomers.push_back(std::move(*monomer));
        }

        // unlikely but we should handle it anyway
        if (monomers.size() == 0) {
            m_error_handler("Failed to parse any monomers", lexer.pos(),
                            lexer.input());
            return false;
        }

        // CHEM polymers shouldn't be multimeric
        if (polymer_id.starts_with("CHEM") && monomers.size() > 1) {
            m_error_handler("CHEM polymers can only have one monomer",
                            polymer_id, lexer.input());
            return false;
        }

        return true;
    }

    /**
     * @brief Parses a single monomer unit from the lexer
     * @param lexer The lexer instance.
     * @return an optional monomer instance if parsing is successful
     */
    std::optional<monomer> parse_monomer_unit(lexer_t<error_handler_t>& lexer)
    {
        if (lexer.eof()) {
            m_error_handler("Expected monomer unit here.", lexer.pos(),
                            lexer.input());
            return std::optional<monomer>{};
        }

        // FIXME: Add monomer list support
        if (lexer.peek() == '(') {
            m_error_handler("Monomer lists are currently unsupported.",
                            lexer.pos(), lexer.input());
            return std::optional<monomer>{};
        }

        monomer result;
        if (lexer.peek() == '[') {
            std::ignore = lexer.get_grouped_token(result.id, '[', ']');
        } else {
            std::ignore = lexer.get_next_token(result.id);
        }

        // NOTE: assume syntax error was handled by lexer
        if (result.id.empty() ||
            !is_valid_monomer_id(result.id, lexer.input())) {
            return std::optional<monomer>{};
        }

        // mark smiles monomers
        if (result.id.front() == '[' && helm::is_smiles_monomer(result.id)) {
            result.is_smiles = true;
        }

        // strip group token for multi-character ids or monomer lists
        if (result.id.front() == '[' || result.id.front() == '(') {
            result.id.remove_prefix(1);
            result.id.remove_suffix(1);
        }

        if (lexer.peek() == '"' &&
            !parse_annotation(lexer, result.annotation)) {
            return std::optional<monomer>{};
        }

        return std::make_optional<monomer>(std::move(result));
    }

    /**
     * @brief Checks if a monomer id or inline SMILES monomer follows the
     *        HELMV2.0 specification.
     * @param monomer_id the monomer id to validate
     * @param input the entire HELMV2.0 input
     * @return True if the monomer id is valid False otherwise
     */
    bool is_valid_monomer_id(std::string_view monomer_id,
                             std::string_view input)
    {
        // unlikely, but we should check this
        if (monomer_id.empty()) {
            m_error_handler("Empty monomer ids are not supported.", monomer_id,
                            input);
            return false;
        }

        // we want to handle this special case
        if (monomer_id.find('\n') != std::string_view::npos) {
            m_error_handler(
                "New-line characters are not allowed in monomer ids.",
                monomer_id, input);
            return false;
        }

        // bracket-enclosed tokens should have multiple characters.
        if (monomer_id.front() == '[') {
            if (monomer_id.size() <= 3) {
                m_error_handler("Bracket-enclosed monomer ids are only "
                                "supported for multi-character monomers.",
                                monomer_id, input);
                return false;
            }
        }
        // these should have only one character
        else if (monomer_id.size() != 1) {
            m_error_handler(
                "Multi-character monomer ids should be surrounded by brackets.",
                monomer_id, input);
            return false;
        }

        return true;
    }

    /**
     * @brief Parse unknown BLOB sequence.
     * @param lexer The lexer instance.
     * @param monomers Output to add parsed monomer to.
     * @return True if successfully parsed the unknown BLOB sequence.
     */
    [[nodiscard]] bool parse_blob_monomer(lexer_t<error_handler_t>& lexer,
                                          std::vector<monomer>& monomers)
    {
        auto remainder = lexer.remainder();
        auto end_pos = remainder.find_first_of("{}$");
        if (end_pos == std::string_view::npos) {
            m_error_handler("Missing '}' character to end BLOB polymer",
                            lexer.size(), lexer.input());
            return false;
        }

        if (remainder[end_pos] != '}') {
            m_error_handler("Expected '}' character to end BLOB polymer here",
                            lexer.pos() + end_pos, lexer.input());
            return false;
        }

        if ((end_pos == 0 && remainder[end_pos] == '}') ||
            std::ranges::all_of(remainder.substr(0, end_pos),
                                detail::isspace)) {
            m_error_handler("Empty BLOB polymers are not supported: expected "
                            "an unknown polymer sequence here. e.g. Bead.",
                            lexer.pos(), lexer.input());
            return false;
        }

        auto unknown_sequence = remainder.substr(0, end_pos);
        // prevent people from setting annotations on monomers
        if (std::ranges::count(unknown_sequence, '"') == 2 &&
            unknown_sequence.back() == '"') {
            m_error_handler("Inline annotations on BLOB monomers aren't "
                            "supported. Consider setting the annotation on the "
                            "polymer. e.g., BLOB1{Bead}\"my annotation\"",
                            unknown_sequence, lexer.input());
            return false;
        }

        monomers.push_back({.id = unknown_sequence});
        lexer.advance(end_pos);
        return true;
    }

    /** @brief Parses the Connections section (Section 2) */
    [[nodiscard]] bool parse_connections(lexer_t<error_handler_t>& lexer,
                                         std::vector<connection>& connections)
    {
        if (lexer.peek() == '$') { // connections are optional
            return true;
        }

        m_error_handler("Unsupported feature", lexer.pos(), lexer.input());
        return false;
    }

    /** @brief Parses the Polymer Groups section (Section 3) */
    [[nodiscard]] bool
    parse_polymer_groups(lexer_t<error_handler_t>& lexer,
                         std::vector<polymer_group>& polymer_groups)
    {
        if (lexer.peek() == '$') { // polymer groups are optional
            return true;
        }

        m_error_handler("Unsupported feature", lexer.pos(), lexer.input());
        return false;
    }

    /** @brief Parses the Extended Annotations section (Section 4) */
    [[nodiscard]] bool
    parse_extended_annotations(lexer_t<error_handler_t>& lexer,
                               std::string_view& extended_annotations)
    {
        if (lexer.peek() == '$') { // extended annotations are optional
            return true;
        }

        if (lexer.eof()) {
            m_error_handler("Expected extended annotations section here",
                            lexer.pos(), lexer.input());
            return false;
        }

        auto remainder = lexer.remainder();
        if (remainder.find('$') == std::string_view::npos) {
            m_error_handler("Expected a '$' character after here", lexer.pos(),
                            lexer.input());
            return false;
        }

        // NOTE: will be validated subsequently
        extended_annotations = remainder.substr(0, remainder.rfind('$'));

        lexer.advance(extended_annotations.size()); // advance lexer position
        return true;
    }

    /**
     * @brief Parses the final HELM Version string (Section 5).
     * @param lexer The lexer instance
     * @param helm_version Output for the version
     * @return True if successfully parsed the version
     */
    [[nodiscard]] bool parse_helm_version(lexer_t<error_handler_t>& lexer,
                                          std::string_view& helm_version)
    {
        // assume we're done parsing everything up to this point
        auto token = lexer.remainder();

        // remove trailing whitespace
        while (!token.empty() && detail::isspace(token.back())) {
            token.remove_suffix(1);
        }

        if (token.empty() || token == "V2.0") {
            lexer.advance(lexer.remainder().size()); // we can discard input
            helm_version = token;
            return true;
        }

        std::string_view err_msg;
        if (token.front() == 'v' || token.front() == 'V') {
            err_msg = "Only HELM and HELMV2.0 versions are currently supported";
        } else {
            err_msg = "Expected the version string or nothing here";
        }

        m_error_handler(err_msg, token, lexer.input());
        return false;
    }
};
} // namespace detail

namespace v2
{

std::optional<helm_info> parse_helm(std::string_view input, bool do_throw)
{
    using namespace helm::detail;

    return do_throw ? (parser_t<throw_error_handler_t>{}).parse(input)
                    : (parser_t<log_error_handler_t>{}).parse(input);
}

} // namespace v2
} // namespace helm
