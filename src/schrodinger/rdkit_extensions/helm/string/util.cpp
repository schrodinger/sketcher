#include "schrodinger/rdkit_extensions/helm/string/util.h"

#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>

HelmVersion get_helm_version(const std::string& helm_string)
{
    static constexpr std::string_view V2_TOKEN{"V2.0"};

    const auto& last_dollar_sign_pos = helm_string.rfind('$');
    if (last_dollar_sign_pos == std::string::npos) {
        return HelmVersion::Unsupported;
    }

    const auto& version_string = helm_string.substr(last_dollar_sign_pos + 1);
    if (version_string.empty()) {
        return HelmVersion::V1;
    }

    return (version_string == V2_TOKEN) ? HelmVersion::V2
                                        : HelmVersion::Unsupported;
}

std::string get_helm2_from_helm(const std::string& helm_string)
{
    const auto& helm_version = get_helm_version(helm_string);
    if (helm_version == HelmVersion::V2) {
        return helm_string;
    }

    if (helm_version == HelmVersion::Unsupported) {
        throw std::invalid_argument("Unsupported HELM version.");
    }

    size_t num_seen_sections = 0, current_section_start_idx = 0;
    std::stringstream converted_helm2{};
    for (size_t idx = 0; idx < helm_string.size(); ++idx) {
        if (helm_string[idx] != '$') {
            continue;
        }
        // the simple polymers section should terminate with a '}' token.
        // there are no annotations to consider here.
        const char prev_char = (idx > 0) ? helm_string[idx - 1] : '\0';
        if (num_seen_sections == 0 && prev_char != '}') {
            continue;
        }

        if (num_seen_sections == 2) {
            const auto& has_hpairs = (idx - current_section_start_idx > 1);
            const auto& has_connections =
                (helm_string.substr(current_section_start_idx - 2, 2) != "$$");
            converted_helm2 << ((has_hpairs && has_connections) ? "|" : "");
        }
        // add section terminating token, $. We don't do this for the
        // connections section since the second and third sections are merged in
        // HELMV2.0
        converted_helm2 << helm_string.substr(current_section_start_idx,
                                              idx - current_section_start_idx)
                        << ((num_seen_sections == 1) ? "" : "$");
        // add empty third section for new polymer groups section in HELMV2.0
        converted_helm2 << ((num_seen_sections == 2) ? "$" : "");
        current_section_start_idx = idx + 1;
        ++num_seen_sections;
    }

    if (num_seen_sections != 4) {
        throw std::invalid_argument("Input HELM has missing sections.");
    }

    converted_helm2 << "V2.0";
    return converted_helm2.str();
}
