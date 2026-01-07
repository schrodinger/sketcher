#include "schrodinger/rdkit_extensions/helm/helm_parser.h"

#include <algorithm>
#include <array>
#include <charconv>
#include <fmt/format.h>
#include <fmt/ranges.h>
#include <regex>
#include <string>

#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>

#include "schrodinger/rdkit_extensions/capture_rdkit_log.h"
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
 * @brief Checks if a value can be used as a ratio.
 * @note a ratio is defined as a non-zero unsigned number
 */
bool is_valid_ratio(const std::string_view& sv)
{
    // ratio should be non-zero value
    if (std::ranges::all_of(
            sv, [](unsigned char c) { return c == '0' || c == '.'; })) {
        return false;
    }

    static std::regex ratio_regex{"^(([0-9]+\\.?[0-9]*)|([0-9]*\\.[0-9]+))$"};
    static std::smatch match;

    std::string val{sv};
    return std::regex_match(val, match, ratio_regex);
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
 * @brief Error handler that silences errors but allows the parser to
 *        attempt continuation (though the main parser logic may still fail).
 */
struct no_op_error_handler_t {
    void operator()(std::string_view msg, size_t pos,
                    std::string_view input) const
    {
    }

    void operator()(std::string_view msg, std::string_view bad_token,
                    std::string_view input) const
    {
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

        auto err_msg = fmt::format("Expected character '{:c}' here", c);
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
            m_error_handler("Unexpected end of string.", m_current_pos,
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
                case ':':
                case '?':
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
    // allow cross-instantiation access to private members
    template <typename> friend class parser_t;

  public:
    static constexpr std::string_view BRANCH_MONOMER_GROUP_ERR_MESSAGE{
        "Only one branch monomer is allowed. If you intended to create a "
        "branched group, create a second polymer with the monomer sequence and "
        "define a custom connection with the desired linkage. E.g.: "
        "PEPTIDE1{A(C.G.K)P} could be either "
        "PEPTIDE1{A(C)P}|PEPTIDE2{G.K}$PEPTIDE1,PEPTIDE2,2:R2-1:R1$$$$V2.0 or "
        "PEPTIDE1{A.P}|PEPTIDE2{C.G.K}$PEPTIDE1,PEPTIDE2,1:R3-1:R1$$$$V2.0"};

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
        if (!(parse_simple_polymers(lexer, result.polymers) && //
              lexer.expect('$') &&                             //
              parse_connections(lexer, result.connections) &&  //
              lexer.expect('$') &&                             //
              parse_third_helm_section(lexer, result.connections,
                                       result.polymer_groups) && //
              lexer.expect('$') &&                               //
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
        std::string_view result;
        if (!lexer.get_next_token(result)) {
            return false;
        }

        if (!is_valid_polymer_id(result, lexer.input())) {
            return false;
        }

        // construct the polymer ID
        polymer_id = result;

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
     * @brief Attempts to parse and set a monomer repetition value
     * @param lexer The lexer instance
     * @param repetition_count Output for the monomer repetition value
     * @return True if a monomer repetition block was found and consumed
     */
    [[nodiscard]] bool
    parse_repetition_count(lexer_t<error_handler_t>& lexer,
                           std::string_view& repetition_count)
    {
        std::string_view value;
        if (!parse_quoted_literal(lexer, value, '\'')) {
            return false;
        }

        // empty values shouldn't be allowed
        if (value.size() == 2) {
            m_error_handler("Empty repetition counts are not supported.", value,
                            lexer.input());
            return false;
        }

        // remove surrounding quotes
        value.remove_prefix(1);
        value.remove_suffix(1);

        // make sure number of repetitions are valid
        if (auto pos = value.find('-'); pos != std::string_view::npos) {
            // leading dash
            if (pos == 0) {
                m_error_handler(
                    "Expected an unsigned non-zero integer before this.", value,
                    lexer.input());
                return false;
            }

            // trailing dash
            if (pos == value.size() - 1) {
                m_error_handler(
                    "Expected an unsigned non-zero integer after this.",
                    value.substr(pos), lexer.input());
                return false;
            }

            // multiple dashes
            auto pos2 = value.find('-', pos + 1);
            if (pos2 != std::string_view::npos) {
                m_error_handler("Repetition counts can't have multiple dashes.",
                                value.substr(pos2), lexer.input());
                return false;
            }

            // both values should be valid
            auto start = value.substr(0, pos);
            auto end = value.substr(pos + 1);
            if (!(is_valid_num_repetitions(start, lexer.input()) &&
                  is_valid_num_repetitions(end, lexer.input()))) {
                return false;
            }

            auto value1 = 0;
            auto value2 = 0;
            // assume values we can always convert values
            std::from_chars(start.data(), start.data() + start.size(), value1);
            std::from_chars(end.data(), end.data() + end.size(), value2);

            // make sure repetition makes sense
            if (value1 >= value2) {
                m_error_handler(
                    "Repetition counts should be strictly increasing. i.e. for "
                    "'N-M', N should be lesser than M.",
                    start, lexer.input());
                return false;
            }
        }
        // this isn't a range
        else if (!is_valid_num_repetitions(value, lexer.input())) {
            return false;
        }

        repetition_count = value;
        return true;
    }

    /**
     * @brief Checks if a repetition count is a positive unsigned integer
     * @param repetition_count the monomer repetition value
     * @param input the entire HELMV2.0 input
     * @return True if the repetition count is valid False otherwise
     */
    [[nodiscard]] bool
    is_valid_num_repetitions(std::string_view repetition_count,
                             std::string_view input)
    {
        if (!std::ranges::all_of(repetition_count, detail::isdigit)) {
            m_error_handler("Repetition count is expected to be an unsigned "
                            "non-zero integer.",
                            repetition_count, input);
            return false;
        }

        if (std::ranges::all_of(repetition_count,
                                [](unsigned char c) { return c == '0'; })) {
            m_error_handler("Repetition count is expected to be an unsigned "
                            "non-zero integer.",
                            repetition_count, input);
            return false;
        }

        if (repetition_count.front() == '0') {
            m_error_handler("Repetition count can't have leading zeros.",
                            repetition_count, input);
            return false;
        }

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
        while (!lexer.eof() && lexer.peek() != '$') {
            // '|' should be used to separate multiple simple polymers
            if (!polymers.empty() && !lexer.expect('|')) {
                return false;
            }

            if (!parse_simple_polymer(lexer, polymers)) {
                return false;
            }
        }

        // This is unlikely, but we should handle this
        if (polymers.empty()) {
            m_error_handler("Simple polymers are required", lexer.pos(),
                            lexer.input());
            return false;
        }

        return true;
    }

    /**
     * @brief Parse a Simple Polymer token. A structure of a polymer is as
     *        follows:
     *        * PEPTIDE and RNA polymers can consist of multiple monomer groups
     *          separated by a dot.
     *        * CHEM polymers can only have a single monomer. This monomer must
     *          conform to the monomer specifications for the PEPTIDE and RNA
     *          monomers.
     *        * BLOB polymers can only have a single monomer. This monomer can
     *          be anything. Examples are "BEAD" and "Gold particle".
     *        * Polymers can have an inline annotation
     *
     *
     *        For PEPTIDE and RNA polymers, we define a monomer group as a
     *        substructure  that can for a backbone linkage to another
     *        substructure within the same polymer. These are:
     *        * A single monomer, e.g., X, [dX], (X.X.X.X), *, _
     *        * A branch monomer group, e.g., X(X), X(X)X
     *        * A repeated monomer sequence, e.g., X'N', (X.X.X)'N'
     *
     * @param lexer The lexer instance.
     * @param polymers Output to add parsed Simple Polymer AST to.
     * @return True if successfully parsed the Simple Polymer.
     */
    [[nodiscard]] bool parse_simple_polymer(lexer_t<error_handler_t>& lexer,
                                            std::vector<polymer>& polymers)
    {
        polymer result;

        // Enforce the strict sequential grammar for the Simple Polymer section
        if (!(parse_polymer_id(lexer, result.id) && //
              lexer.expect('{') &&                  //
              parse_monomers(lexer, result.id, result.monomers,
                             result.repetitions) && //
              lexer.expect('}') &&
              // polymer annotations are optional
              (lexer.peek() != '"' ||
               parse_annotation(lexer, result.annotation)))) {
            return false;
        }

        // add extra information for validating connections. e.g. unknown
        // and wildcard monomer info
        add_ambiguous_monomer_information_to_polymer(result);

        polymers.push_back(std::move(result));
        return true;
    }

    /** @brief Records which monomers are wildcards or unknown, and records the
     *         list of seen residue names. This information will be used to
     *         validate ambiguous connections.
     */
    void add_ambiguous_monomer_information_to_polymer(polymer& result)
    {
        auto monomer_idx = 0;
        for (auto& monomer : result.monomers) {
            ++monomer_idx; // 1-indexed

            if (monomer.is_smiles) {
                continue;
            }

            auto separator = '\n';
            if (monomer.is_list) {
                bool is_monomer_union = std::ranges::count(monomer.id, '+');
                separator = (is_monomer_union ? '+' : ',');
            }

            size_t start = 0;
            while (start < monomer.id.size()) {
                size_t end = monomer.id.find(separator, start);
                if (end == std::string_view::npos) {
                    end = monomer.id.size();
                }

                auto monomer_id = monomer.id.substr(start, end - start);
                monomer_id = monomer_id.substr(0, monomer_id.find(':'));

                // strip brackets
                if (monomer_id.front() == '[') {
                    monomer_id.remove_prefix(1);
                    monomer_id.remove_suffix(1);
                }

                if (monomer_id == "*") { // this is a wildcard monomer
                    // this will be used to validate ambiguous connections
                    result.wildcard_and_unknown_residues.insert(monomer_idx);
                } else if (result.id.starts_with("PEPTIDE") &&
                           monomer_id == "X") { // this is an unknown amino acid
                    // this will be used to validate ambiguous connections
                    result.wildcard_and_unknown_residues.insert(monomer_idx);
                }
                // this is an unknown nucleotide
                else if (result.id.starts_with("RNA") && monomer_id == "N") {
                    // this will be used to validate ambiguous connections
                    result.wildcard_and_unknown_residues.insert(monomer_idx);
                }

                // this will be used to validate ambiguous connections
                result.residue_names.insert(monomer_id);

                start = end + 1; // advance location of separator
            }
        }
    }

    /**
     * @brief Parse list of monomers in a Simple Polymer token.
     * @param lexer The lexer instance.
     * @param polymer_id The ID of the current polymer.
     * @param monomers Output to add parsed monomers to.
     * @param repetitions Output to add parsed monomer repetitions to.
     * @return True if successfully parsed the list of monomers.
     */
    [[nodiscard]] bool parse_monomers(lexer_t<error_handler_t>& lexer,
                                      const std::string_view& polymer_id,
                                      std::vector<monomer>& monomers,
                                      std::vector<repetition>& repetitions)
    {
        if (polymer_id.starts_with("BLOB")) { // needs special handling
            return parse_blob_monomer(lexer, monomers);
        }

        while (!lexer.eof() && lexer.peek() != '}' /* end of monomers */) {
            if (!monomers.empty()) { // we might need a '.' to separate monomers
                // for examples like R(C)P, we don't require a '.' token.
                if (monomers.back().is_branch && lexer.peek() != '.') {
                    // nothing to do here
                }
                // Anything else should require a '.' token
                else if (!lexer.expect('.')) {
                    return false;
                }
            }

            if (next_token_is_repeated_monomer_sequence(lexer)) {
                if (!parse_repeated_monomer_sequence(lexer, polymer_id,
                                                     monomers, repetitions)) {
                    return false;
                }
            } else {
                std::string_view repetition_count;
                auto monomer = parse_monomer_unit(lexer, repetition_count);
                if (!monomer) { // we failed to get monomer
                    return false;
                }

                if (!repetition_count.empty()) {
                    repetitions.push_back({.start = monomers.size(),
                                           .size = 1,
                                           .num_repetitions = repetition_count,
                                           .annotation = {}});
                    // re-assign the annotation to the repetition since it looks
                    // like A'11'"value"
                    std::swap(repetitions.back().annotation,
                              monomer->annotation);
                }

                monomers.push_back(std::move(*monomer));
            }

            if (lexer.peek() == '(') { // assume branch monomer
                auto branch_monomer = parse_branch_monomer(lexer);
                if (!branch_monomer) {
                    return false;
                }

                // This will handle examples like PEPTIDE1{(A(C))'3'(C)}$$$$V2.0
                if (monomers.back().is_branch) {
                    m_error_handler(
                        "Branch monomers must come after backbone monomers",
                        branch_monomer->id, lexer.input());
                    return false;
                }

                monomers.push_back(std::move(*branch_monomer));
            }
        }

        // unlikely but we should handle it anyway
        if (monomers.size() == 0) {
            m_error_handler("Failed to parse any monomers", lexer.pos(),
                            lexer.input());
            return false;
        }

        if (polymer_id.starts_with("CHEM")) {
            // CHEM polymers shouldn't be multimeric
            if (monomers.size() > 1) {
                m_error_handler("CHEM polymers can only have one monomer",
                                polymer_id, lexer.input());
                return false;
            }

            // CHEM monomers can't be repeated
            if (repetitions.size() >= 1) {
                m_error_handler("CHEM polymers can only have one monomer, so "
                                "CHEM monomers can't be repeated.",
                                polymer_id, lexer.input());
                return false;
            }
        }

        return true;
    }

    /**
     * @brief Parse repeated monomer list in a Simple Polymer token. e.g.,
     *        (A.C.L)'4'
     * @param lexer The lexer instance.
     * @param polymer_id The ID of the current polymer.
     * @param monomers Output to add parsed monomers to.
     * @param repetitions Output to add parsed monomer repetitions to.
     * @return True if successfully parsed the repeated monomer sequence.
     */
    bool parse_repeated_monomer_sequence(lexer_t<error_handler_t>& lexer,
                                         const std::string_view& polymer_id,
                                         std::vector<monomer>& monomers,
                                         std::vector<repetition>& repetitions)
    {
        if (lexer.peek() != '(') {
            m_error_handler("Expected character '(' here to begin repeated "
                            "monomer sequence.",
                            lexer.pos(), lexer.input());
            return false;
        }

        std::string_view monomer_sequence;
        // this is unlikely to fail, since we make sure it's a valid monomer
        // sequence before parsing
        if (!lexer.get_grouped_token(monomer_sequence, '(', ')')) {
            return false;
        }

        // strip leading and trailing parenthesis
        monomer_sequence.remove_suffix(1);
        monomer_sequence.remove_prefix(1);

        // use a temporary list of entities so subsequent checks are easier
        lexer_t<error_handler_t> repeated_monomers_lexer(monomer_sequence);
        std::vector<monomer> repeated_monomers;
        std::vector<repetition> nested_repetitions;
        if (!parse_monomers(repeated_monomers_lexer, polymer_id,
                            repeated_monomers, nested_repetitions)) {
            return false;
        }

        if (!nested_repetitions.empty()) {
            m_error_handler("Repeated monomer lists cannot be nested",
                            nested_repetitions[0].num_repetitions,
                            lexer.input());
            return false;
        }

        if (lexer.peek() != '\'') {
            m_error_handler(
                "Expected a repetition token here. e.g. '4', '1-4'.",
                lexer.pos(), lexer.input());
            return false;
        }

        std::string_view repetition_count;
        if (!parse_repetition_count(lexer, repetition_count)) {
            return false;
        }

        std::string_view annotation;
        if (lexer.peek() == '"' && !parse_annotation(lexer, annotation)) {
            return false;
        }

        // add repetitions to output
        repetitions.push_back({.start = monomers.size(),
                               .size = repeated_monomers.size(),
                               .num_repetitions = repetition_count,
                               .annotation = annotation});

        // extend output monomer list
        monomers.insert(monomers.end(), repeated_monomers.begin(),
                        repeated_monomers.end());

        return true;
    }

    /**
     * @brief Parses a single monomer unit from the lexer
     * @param lexer The lexer instance.
     * @param repetition_count Output for the monomer repetition value
     * @return an optional monomer instance if parsing is successful
     */
    std::optional<monomer>
    parse_monomer_unit(lexer_t<error_handler_t>& lexer,
                       std::string_view& repetition_count)
    {
        if (lexer.eof()) {
            m_error_handler("Expected monomer unit here.", lexer.pos(),
                            lexer.input());
            return std::optional<monomer>{};
        }

        monomer result;
        if (lexer.peek() == '[') { // multi-character ids or inline SMILES
            std::ignore = lexer.get_grouped_token(result.id, '[', ']');
        } else if (lexer.peek() == '(') { // monomer lists
            std::ignore = lexer.get_grouped_token(result.id, '(', ')');
        } else {
            std::ignore = lexer.get_next_token(result.id);
        }

        // NOTE: assume syntax error was handled by lexer
        if (result.id.empty()) {
            return std::optional<monomer>{};
        }

        if (result.id.front() != '(' &&
            !is_valid_monomer_id(result.id, lexer.input())) {
            return std::optional<monomer>{};
        }

        if (result.id.front() == '(' &&
            !is_valid_monomer_list(result.id, lexer.input())) {
            return std::optional<monomer>{};
        }

        // mark smiles monomers
        if (result.id.front() == '[' && helm::is_smiles_monomer(result.id)) {
            result.is_smiles = true;
        }

        // mark monomer lists
        if (result.id.front() == '(') {
            result.is_list = true;
        }

        // strip group token for multi-character ids or monomer lists
        if (result.id.front() == '[' || result.id.front() == '(') {
            result.id.remove_prefix(1);
            result.id.remove_suffix(1);
        }

        if (lexer.peek() == '\'' &&
            !parse_repetition_count(lexer, repetition_count)) {
            return std::optional<monomer>{};
        }

        if (lexer.peek() == '"' &&
            !parse_annotation(lexer, result.annotation)) {
            return std::optional<monomer>{};
        }

        return std::make_optional<monomer>(std::move(result));
    }

    /**
     * @brief Parses a single branch monomer unit from the lexer
     * @param lexer The lexer instance.
     * @return an optional monomer instance if parsing is successful
     */
    std::optional<monomer> parse_branch_monomer(lexer_t<error_handler_t>& lexer)
    {
        if (!lexer.expect('(')) {
            return std::optional<monomer>{};
        }

        std::string_view repetition_count;
        if (auto result = parse_monomer_unit(lexer, repetition_count); result) {
            if (lexer.peek() != ')') {
                // NOTE: We're intentionally verbosely handling this since it's
                // a likely use case. People will find it easier to create
                // custom connections this way.
                if (lexer.peek() == '.' || lexer.peek() == '(') {
                    m_error_handler(BRANCH_MONOMER_GROUP_ERR_MESSAGE,
                                    lexer.pos(), lexer.input());
                }
                // The generic failure case
                else {
                    m_error_handler("Expected character ')' here to close "
                                    "branch monomer.",
                                    lexer.pos(), lexer.input());
                }
                return std::optional<monomer>{};
            }

            // you can't repeat branch monomers
            if (!repetition_count.empty()) {
                m_error_handler("Branch monomers can't be repeated",
                                repetition_count, lexer.input());
                return std::optional<monomer>{};
            }

            lexer.advance(1); // consume ')'

            result->is_branch = true;
            return result;
        }

        return std::optional<monomer>{};
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
        // we want to handle this special case
        if (monomer_id.find('\n') != std::string_view::npos) {
            m_error_handler(
                "New-line characters are not allowed in monomer ids.",
                monomer_id, input);
            return false;
        }

        // bracket-enclosed tokens should have multiple characters.
        if (monomer_id.front() == '[' && monomer_id.back() == ']') {
            if (monomer_id.size() <= 3) {
                m_error_handler("Bracket-enclosed monomer ids are only "
                                "supported for multi-character monomers.",
                                monomer_id, input);
                return false;
            }

            // strip surrounding parenthesis
            monomer_id.remove_prefix(1);
            monomer_id.remove_suffix(1);

            if (std::ranges::all_of(monomer_id, detail::isdigit)) {
                m_error_handler(
                    "Bracket-enclosed monomer ids shouldn't be numeric.",
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
        } else if (monomer_id != "_" && monomer_id != "*" &&
                   monomer_id != "?" && !std::isalpha(monomer_id.front())) {
            m_error_handler("The allowed single character monomer ids are "
                            "alphabets or any of '_', '?', and '*'.",
                            monomer_id, input);
            return false;
        }

        return true;
    }

    /**
     * @brief Checks if there is consistent usage of entity list separators.
     * @param entity_list the entity list to validate
     * @param input the entire HELMV2.0 input
     * @return True if the entity list has valid separators False otherwise
     */
    bool has_valid_entity_list_separators(std::string_view entity_list,
                                          std::string_view input)
    {
        auto separator_pos = entity_list.find_first_of(",+");
        if (separator_pos == std::string_view::npos) {
            // Special case for wrongly  positioned sequences
            if (is_valid_monomer_sequence(entity_list)) {
                m_error_handler("Unsupported placement of monomer sequence.",
                                entity_list, input);
                return false;
            }

            m_error_handler("Entity lists should have more than one item.",
                            entity_list, input);
            return false;
        }

        auto separator = entity_list[separator_pos];
        // leading separator
        if (entity_list.front() == separator) {
            m_error_handler("Expected list item before this", entity_list,
                            input);
            return false;
        }

        // trailing separator
        if (entity_list.back() == separator) {
            m_error_handler("Expected list item after this",
                            entity_list.substr(entity_list.size() - 1), input);
            return false;
        }

        if (separator == '+' &&
            entity_list.find(',') != std::string_view::npos) {
            m_error_handler(
                "Inconsistent use of list separators. Expected '+' here.",
                entity_list.substr(entity_list.find(',')), input);
            return false;
        } else if (separator == ',' &&
                   entity_list.find('+') != std::string_view::npos) {
            m_error_handler(
                "Inconsistent use of list separators. Expected ',' here.",
                entity_list.substr(entity_list.find('+')), input);
            return false;
        }

        size_t start = 0;
        while (start < entity_list.size()) {
            size_t end = entity_list.find(separator, start);
            if (end == start) { // assume extra delimiter e.g (A,,B)
                m_error_handler("Remove extra monomer list separator here.",
                                entity_list.substr(end), input);
                return false;
            }

            // we're done
            if (end == std::string_view::npos) {
                break;
            }

            start = end + 1; // advance location of separator
        }

        return true;
    }

    /**
     * @brief Checks if a monomer list follows the HELMV2.0 specification. This
     *        mostly checks that there is consistent usage of the monomer list
     *        separators and that monomer list items are fully formed. E.g.
     *          - A+G+C
     *          - A,[dG],L
     *          - A:1+G:2+C:3
     *          - A+G:?
     * @param monomer_list the monomer list to validate
     * @param input the entire HELMV2.0 input
     * @return True if the monomer list is valid False otherwise
     */
    bool is_valid_monomer_list(std::string_view monomer_list,
                               std::string_view input)
    {
        // it should have enough values. this will catch things like ()
        if (monomer_list.size() <= 3) {
            m_error_handler("Monomer lists should have more than one item.",
                            monomer_list, input);
            return false;
        }

        // strip surrounding parenthesis
        monomer_list.remove_prefix(1);
        monomer_list.remove_suffix(1);

        // make sure there's consistent use of + or ,
        if (!has_valid_entity_list_separators(monomer_list, input)) {
            return false;
        }

        auto separator = monomer_list[monomer_list.find_first_of(",+")];

        // Validate individual monomers
        size_t start = 0;
        while (start < monomer_list.size()) {
            size_t end = monomer_list.find(separator, start);
            if (end == std::string_view::npos) {
                end = monomer_list.size();
            }

            auto item = monomer_list.substr(start, end - start);
            start = end + 1; // advance location of separator

            // monomer list items should be monomer ids followed by an optional
            // ratio. the monomer id and the ratio should be separated by a
            // colon
            auto colon_pos = item.find(':');
            if (!is_valid_monomer_id(item.substr(0, colon_pos), input)) {
                return false;
            }

            // nothing to validate
            if (colon_pos == std::string_view::npos) {
                continue;
            }

            // it should have values after it
            if (colon_pos == item.size() - 1) {
                m_error_handler("Expected a valid ratio after this.",
                                item.substr(colon_pos), input);
                return false;
            }

            // ensure this is a numeric value
            auto ratio = item.substr(colon_pos + 1);
            if (ratio != "?" && !detail::is_valid_ratio(ratio)) {
                m_error_handler("Expected a non-zero unsigned number here.",
                                ratio, input);
                return false;
            }
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

        while (!lexer.eof() && lexer.peek() != '$') {
            // multiple connections should be separated by '|'
            if (!connections.empty() && !lexer.expect('|')) {
                return false;
            }

            // parse connection
            if (!parse_connection(lexer, connections)) {
                return false;
            }
        }

        return true;
    }

    /** @brief Parses a single connection token from the lexer. A connection
     *         token should have the following structure:
     *
     *         CONNECTION_POLYMER,CONNECTION_POLYMER,CONNECTION_MONOMER:ATTACHMENT_POINT-CONNECTION_MONOMER:ATTACHMENT_POINT
     *
     *         Some notes about HELM connections:
     *         * A custom connection can be formed between two monomers
     *         * R-groups that partake in a custom connection must be specified
     *           except in the case of BLOBs or ambiguous connections.
     *         * If multiple polymers form the same type of custom bond, you can
     *           alias the polymer ids through a polymer group and use the
     *           polymer group id as the connection polymer.
     *         * You can define multiple connection monomers if they form the
     *           same type of bond, e.g.,
     *           POLYMER,POLYMER,(MONOMER,MONOMER):R3-MONOMER:R1
     *         * Connections can have inline annotations.
     *
     */
    [[nodiscard]] bool parse_connection(lexer_t<error_handler_t>& lexer,
                                        std::vector<connection>& connections)
    {
        // parse the connection token
        connection result;
        if (!(lexer.get_next_token(result.from_id) && lexer.expect(',') &&
              lexer.get_next_token(result.to_id) && lexer.expect(',') &&
              parse_connection_monomer(lexer, result.from_res) &&
              lexer.expect(':') &&
              parse_connection_attachment_point(lexer, result.from_rgroup) &&
              lexer.expect('-') &&
              parse_connection_monomer(lexer, result.to_res) &&
              lexer.expect(':') &&
              parse_connection_attachment_point(lexer, result.to_rgroup))) {
            return false;
        }

        // parse optional inline annotation
        if (lexer.peek() == '"' &&
            !parse_annotation(lexer, result.annotation)) {
            return false;
        }

        connections.push_back(std::move(result));
        return true;
    }

    /**
     * @brief Parses a single connection monomer. A connection monomer is
     *        defined as:
     *          * connection residue
     *          * undefined residue i.e. ?
     *          * ( residue_and_list )
     *          * ( residue_or_list )
     */
    [[nodiscard]] bool
    parse_connection_monomer(lexer_t<error_handler_t>& lexer,
                             std::string_view& connection_monomer)
    {
        if (lexer.peek() != '(') {
            // this should be a connection residue or ?
            return lexer.get_next_token(connection_monomer) &&
                   (connection_monomer == "?" ||
                    is_valid_connection_residue(connection_monomer,
                                                lexer.input()));
        }

        // this should be a residue list
        std::string_view residue_list;
        if (!lexer.get_grouped_token(residue_list, '(', ')')) {
            return false;
        }

        // strip surrounding parenthesis
        residue_list.remove_prefix(1);
        residue_list.remove_suffix(1);

        // make sure there's consistent use of + or ,
        if (!has_valid_entity_list_separators(residue_list, lexer.input())) {
            return false;
        }

        auto separator = residue_list[residue_list.find_first_of(",+")];

        // Validate individual entities
        size_t start = 0;
        while (start < residue_list.size()) {
            size_t end = residue_list.find(separator, start);
            if (end == std::string_view::npos) {
                end = residue_list.size();
            }

            auto item = residue_list.substr(start, end - start);
            start = end + 1; // advance location of separator

            if (!is_valid_connection_residue(item, lexer.input())) {
                return false;
            }
        }

        connection_monomer = residue_list;
        return true;
    }

    /** @brief Checks if a token is a valid connection residue/monomer */
    [[nodiscard]] bool
    is_valid_connection_residue(std::string_view connection_residue,
                                std::string_view input)
    {
        // assume this is a residue number
        if (std::ranges::all_of(connection_residue, detail::isdigit)) {
            if (connection_residue.front() == '0') {
                m_error_handler("Connection residue numbers shouldn't have "
                                "leading zeros.",
                                connection_residue, input);
                return false;
            }

            return true;
        }

        // NOTE: only support alpha-numeric monomer ids for now.
        if (!std::ranges::all_of(connection_residue, [](unsigned char c) {
                return std::isalnum(c);
            })) {
            m_error_handler("Connection residues are expected to be residue "
                            "numbers or residue names. Only alphanumeric "
                            "residue names are currently supported.",
                            connection_residue, input);
            return false;
        }

        return true;
    }

    /** @brief Parses a single connection attachment point token from the lexer
     */
    [[nodiscard]] bool
    parse_connection_attachment_point(lexer_t<error_handler_t>& lexer,
                                      std::string_view& attachment_point)
    {
        return lexer.get_next_token(attachment_point) &&
               is_valid_connection_attachment_point(attachment_point,
                                                    lexer.input());
    }

    /** @brief Checks if a token is a valid connection attachment point */
    [[nodiscard]] bool
    is_valid_connection_attachment_point(std::string_view attachment_point,
                                         std::string_view input)
    {
        // these are special cases
        if (attachment_point == "?" || attachment_point == "pair") {
            return true;
        }

        // this has to be an r-group
        if (attachment_point.size() == 1 || attachment_point.front() != 'R') {
            m_error_handler("The supported attachment points are '?', 'pair', "
                            "and r-groups. R-groups are defined as R#, where # "
                            "is an unsigned positive number",
                            attachment_point, input);
            return false;
        }

        attachment_point.remove_prefix(1); // consume R

        if (attachment_point.front() == '0') {
            m_error_handler("R-group numbers shouldn't have leading zeros.",
                            attachment_point, input);
            return false;
        }

        if (!std::ranges::all_of(attachment_point, detail::isdigit)) {
            m_error_handler("R-groups are defined as R#, where # is an "
                            "unsigned positive number",
                            attachment_point, input);
            return false;
        }

        return true;
    }

    /** @brief Parses the Polymer Groups section (Section 3) */
    [[nodiscard]] bool
    parse_polymer_groups(lexer_t<error_handler_t>& lexer,
                         std::vector<polymer_group>& polymer_groups)
    {
        if (lexer.peek() == '$') { // polymer groups are optional
            return true;
        }

        while (!lexer.eof() && lexer.peek() != '$') {
            if (!polymer_groups.empty() && !lexer.expect('|')) {
                return false;
            }

            // parse polymer group id
            std::string_view polymer_group_id;
            if (!(lexer.get_next_token(polymer_group_id) &&
                  is_valid_polymer_group_id(polymer_group_id, lexer.input()))) {
                return false;
            }

            // parse polymer group items
            if (lexer.peek() != '(') {
                m_error_handler(
                    "Expected a polymer group list after the polymer group id.",
                    lexer.pos(), lexer.input());
                return false;
            }

            std::string_view polymer_group_list;
            if (!(lexer.get_grouped_token(polymer_group_list, '(', ')') &&
                  is_valid_polymer_group_list(polymer_group_list,
                                              lexer.input()))) {
                return false;
            }

            // strip surrounding parenthesis
            polymer_group_list.remove_prefix(1);
            polymer_group_list.remove_suffix(1);

            bool is_polymer_union = std::ranges::count(polymer_group_list, '+');
            polymer_groups.push_back({.id = polymer_group_id,
                                      .items = polymer_group_list,
                                      .is_polymer_union = is_polymer_union});
        }

        return true;
    }

    /**
     * @brief Checks if a polymer group list follows the HELMV2.0 specification.
     * This mostly checks that there is consistent usage of the list separators
     * and that list items are fully formed. E.g.
     *          - CHEM1+CHEM2
     *          - CHEM1,CHEM2
     *          - G1:0.5,CHEM2
     * @param polymer_group_list the polymer group list to validate
     * @param input the entire HELMV2.0 input
     * @return True if the polymer group list is valid False otherwise
     */
    bool is_valid_polymer_group_list(std::string_view polymer_group_list,
                                     std::string_view input)
    {
        // it should have enough values. this will catch things like ()
        if (polymer_group_list.size() <= 3) {
            m_error_handler("Monomer lists should have more than one item.",
                            polymer_group_list, input);
            return false;
        }

        // strip surrounding parenthesis
        polymer_group_list.remove_prefix(1);
        polymer_group_list.remove_suffix(1);

        // make sure there's consistent use of + or ,
        if (!has_valid_entity_list_separators(polymer_group_list, input)) {
            return false;
        }

        auto separator =
            polymer_group_list[polymer_group_list.find_first_of(",+")];

        // Validate individual monomers
        size_t start = 0;
        while (start < polymer_group_list.size()) {
            size_t end = polymer_group_list.find(separator, start);
            if (end == std::string_view::npos) {
                end = polymer_group_list.size();
            }

            auto item = polymer_group_list.substr(start, end - start);
            start = end + 1; // advance location of separator

            // polymer group list items should be monomer ids followed by an
            // optional ratio. the polymer group item and the ratio should be
            // separated by a colon
            auto colon_pos = item.find(':');
            auto id = item.substr(0, colon_pos);
            if (item.front() == 'G' && !is_valid_polymer_group_id(id, input)) {
                return false;
            }

            if (item.front() != 'G' && !is_valid_polymer_id(id, input)) {
                return false;
            }

            // nothing to validate
            if (colon_pos == std::string_view::npos) {
                continue;
            }

            // it should have values after it
            if (colon_pos == item.size() - 1) {
                m_error_handler("Expected a valid ratio after this.",
                                item.substr(colon_pos), input);
                return false;
            }

            auto ratio = item.substr(colon_pos + 1);

            // special case
            if (ratio == "?") {
                continue;
            }

            auto dash_pos = ratio.find('-');
            if (dash_pos == 0 || dash_pos == ratio.size() - 1) {
                m_error_handler("Polymer group ratios can't have leading or "
                                "trailing dashes.",
                                ratio, input);
                return false;
            }

            // ensure this is a numeric value
            if (!detail::is_valid_ratio(ratio.substr(0, dash_pos))) {
                m_error_handler("Expected a valid ratio here.", ratio, input);
                return false;
            }

            // ensure this is a numeric value
            if (dash_pos != std::string_view::npos &&
                !detail::is_valid_ratio(ratio.substr(dash_pos + 1))) {
                m_error_handler("Expected two ratios separated by a dash here. "
                                "e.g., 0.2-0.5.",
                                ratio, input);
                return false;
            }
        }

        return true;
    }

    /**
     * @brief Checks if a value is a valid polymer group id.
     * @param polymer_group_id the value to validate
     * @param input the entire HELMV2.0 input
     * @return True if the value is a valid polymer group id False otherwise
     */
    bool is_valid_polymer_group_id(std::string_view polymer_group_id,
                                   std::string_view input)
    {
        // polymer group ids should start with G
        if (polymer_group_id.front() != 'G') {
            m_error_handler("Polymer group ids should start with 'G'. e.g. G1.",
                            polymer_group_id, input);
            return false;
        }

        if (polymer_group_id.size() == 1) {
            m_error_handler("Polymer group ids are expected to have numeric "
                            "suffixes. e.g. G1.",
                            polymer_group_id, input);
            return false;
        }

        // strip leading G
        polymer_group_id.remove_prefix(1);

        if (!std::ranges::all_of(polymer_group_id, detail::isdigit)) {
            m_error_handler("Polymer group ids are expected to have numeric "
                            "suffixes. e.g. G1.",
                            polymer_group_id, input);
            return false;
        }

        // Polymer group ids with leading zeros.
        if (polymer_group_id.front() == '0') {
            m_error_handler("Polymer group ids shouldn't have leading zeros.",
                            polymer_group_id, input);
            return false;
        }

        return true;
    }

    /**
     * @brief Checks if a value is a valid polymer id.
     * @param polymer_id the value to validate
     * @param input the entire HELMV2.0 input
     * @return True if the value is a valid polymer id False otherwise
     */
    bool is_valid_polymer_id(std::string_view polymer_id,
                             std::string_view input)
    {
        // Define all supported polymer types
        constexpr std::array<std::string_view, 4> supported_polymers{
            "BLOB", "CHEM", "PEPTIDE", "RNA"};

        auto matched_polymer =
            std::ranges::find_if(supported_polymers, [&](auto& entry) {
                return polymer_id.starts_with(entry);
            });

        if (matched_polymer == supported_polymers.end()) {
            m_error_handler("Expected polymer id here. Only PEPTIDE, RNA, "
                            "CHEM, and BLOB polymers are currently supported",
                            polymer_id, input);
            return false;
        }

        // Enforce numeric suffix
        auto suffix_start = matched_polymer->size();
        if (suffix_start == polymer_id.size() ||
            !detail::isdigit(polymer_id[suffix_start])) {
            m_error_handler(
                "Polymer ids are expected to have numeric suffixes. e.g., RNA1",
                polymer_id, input);
            return false;
        }

        // Reject suffixes with leading zeros.
        if (polymer_id[suffix_start] == '0') {
            m_error_handler("Polymer id suffixes shouldn't have leading zeros.",
                            polymer_id, input);
            return false;
        }

        return true;
    }

    /** @brief Parses the third section of a HELM or HELMV2.0 input. The content
     *         of the third section is as follows:
     *              * HELM: hydrogen pairings
     *              * HELMV2.0: polymer groups
     */
    [[nodiscard]] bool
    parse_third_helm_section(lexer_t<error_handler_t>& lexer,
                             std::vector<connection>& connections,
                             std::vector<polymer_group>& polymer_groups)
    {
        if (lexer.peek() == '$') { // this section is optional
            return true;
        }

        auto remainder = lexer.remainder();

        // remove trailing whitespace
        while (!remainder.empty() && detail::isspace(remainder.back())) {
            remainder.remove_suffix(1);
        }

        // Handle the HELMV2.0 case
        if (remainder.ends_with("V2.0")) {
            return parse_polymer_groups(lexer, polymer_groups);
        }

        // This is the HELM case. The third section should only have hydrogen
        // pairings
        std::vector<connection> hydrogen_pairings;
        if (!parse_connections(lexer, hydrogen_pairings)) {
            return false;
        }

        // make sure they're only hydrogen pairings
        auto bad_connection =
            std::ranges::find_if(hydrogen_pairings, [](auto& entry) {
                return entry.from_rgroup != "pair" || entry.to_rgroup != "pair";
            });
        if (bad_connection != hydrogen_pairings.end()) {
            auto bad_token = bad_connection->from_rgroup == "pair"
                                 ? bad_connection->from_rgroup
                                 : bad_connection->to_rgroup;
            m_error_handler(
                "HELM only supports hydrogen pairings in the third section",
                bad_token, lexer.input());
            return false;
        }

        // add to result
        connections.insert(connections.end(), hydrogen_pairings.begin(),
                           hydrogen_pairings.end());
        return true;
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

    /**
     * @brief Checks if the next token in the lexer is a repeated monomer
     *        sequence. The monomer sequence should consist of more than one
     *        monomer.
     * @param lexer The lexer instance
     * @return True if the next token is a repeated monomer sequence and False
     *         otherwise.
     */
    static bool next_token_is_repeated_monomer_sequence(
        const lexer_t<error_handler_t>& lexer)
    {
        // we want to silence all errors
        lexer_t<no_op_error_handler_t> silent_lexer(lexer.remainder());

        std::string_view token;
        if (lexer.peek() == '(' &&
            silent_lexer.get_grouped_token(token, '(', ')')) {

            // strip surrounding parenthesis
            return is_valid_monomer_sequence(token.substr(1, token.size() - 2));
        }

        return false;
    }

    /**
     * @brief Checks it the input is a monomer sequence. The monomer sequence
     *        should consist of more than one monomer.
     * @param monomer_sequence The input to validate
     * @return True if the input is a monomer sequence and False otherwise.
     */
    static bool is_valid_monomer_sequence(std::string_view monomer_sequence)
    {
        // we want to silence all errors
        parser_t<no_op_error_handler_t> silent_parser;
        lexer_t<no_op_error_handler_t> silent_lexer(monomer_sequence);

        std::vector<monomer> monomers;
        std::vector<repetition> repetitions;

        // the polymer id only allows us to parse the monomers
        constexpr std::string_view dummy_polymer_id{"PEPTIDE1"};
        return silent_parser.parse_monomers(silent_lexer, dummy_polymer_id,
                                            monomers, repetitions) &&
               silent_lexer.eof() && monomers.size() > 1;
    }
};
} // namespace detail

std::optional<helm_info> parse_helm(std::string_view input, bool do_throw)
{
    using namespace helm::detail;

    return do_throw ? (parser_t<throw_error_handler_t>{}).parse(input)
                    : (parser_t<log_error_handler_t>{}).parse(input);
}

[[nodiscard]] bool is_smiles_monomer(const std::string_view& token)
{
    // Assume it's SMILES if it has whitespace
    if (token.find_first_of(" \t\n\f\r\v") != std::string_view::npos) {
        return true;
    }

    auto monomer = token.substr(1, token.size() - 2);
    // Assume it's SMILES if it has these rgroup chars
    // NOTE: token is still surrounded by [] tokens
    if (monomer.find_first_of("*:[]=+") != std::string_view::npos) {
        return true;
    }

    // Parse the SMILES literally
    static const RDKit::v2::SmilesParse::SmilesParserParams smi_opts{
        .sanitize = false, .removeHs = false, .replacements = {}};

    [[maybe_unused]] schrodinger::rdkit_extensions::CaptureRDErrorLog rdkit_log;
    const std::string smiles{monomer};
    return RDKit::v2::SmilesParse::MolFromSmiles(smiles, smi_opts) != nullptr;
}

[[nodiscard]] std::string
construct_error_msg(const std::string_view err_msg,
                    const std::string_view& failed_token,
                    const std::string_view& input)
{
    auto pos = std::distance(input.data(), failed_token.data());
    return construct_error_msg(err_msg, pos, input);
}

[[nodiscard]] std::string construct_error_msg(const std::string_view err_msg,
                                              const unsigned int pos,
                                              const std::string_view& input)
{
    // Defines the maximum length of the snippet and the prefix size before the
    // error
    static constexpr unsigned int error_size{101};
    static constexpr unsigned int prefix_size{error_size / 2};

    size_t start_idx = 0;
    size_t pointer_offset = pos;

    // Calculate the start position to center the error pointer
    if (pos >= prefix_size) {
        start_idx = pos - prefix_size;
        pointer_offset = prefix_size;
    }

    // Get truncated context snippet
    auto truncated_input = input.substr(start_idx, error_size);

    // Create pointer to line (e.g., "---^")
    std::string pointer_line(pointer_offset, '-');

    return fmt::format(
        "Malformed HELM string: check for mistakes around position {}:\n"
        "{}\n"
        "{}^\n"
        "{}\n",
        pos + 1, // 1-based index for user output
        truncated_input, pointer_line, err_msg);
}

} // namespace helm
