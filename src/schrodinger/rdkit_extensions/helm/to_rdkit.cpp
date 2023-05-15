#include "schrodinger/rdkit_extensions/helm/to_rdkit.h"

#include <memory>    // shared_ptr
#include <stdexcept> // invalid_argument
#include <string>
#include <string_view>

#include <GraphMol/ROMol.h>

#include "schrodinger/rdkit_extensions/helm/constants.h"
#include "schrodinger/rdkit_extensions/helm/helm_parser.h"

enum class [[nodiscard]] HelmVersion{V2, Unsupported};

[[nodiscard]] static HelmVersion
get_helm_version(const std::string& helm_string)
{
    static constexpr std::string_view V2_TOKEN{"V2.0"};

    const auto& last_dollar_sign_pos = helm_string.rfind('$');
    if (last_dollar_sign_pos == std::string::npos) {
        return HelmVersion::Unsupported;
    }
    const auto version_string = helm_string.substr(last_dollar_sign_pos + 1);
    return (version_string == V2_TOKEN) ? HelmVersion::V2
                                        : HelmVersion::Unsupported;
}

static void add_helm_specific_properties_to_mol(::RDKit::RWMol& mol)
{
    mol.setProp<bool>(HELM_MODEL_PROP_NAME, true);
}

/*
 * The structure of a HELM2 string is as follows:
 *
 *      POLYMERS$CONNECTIONS$POLYMERGROUPS$EXTENDEDANNOTATIONS$V2.0
 *
 * Other than the polymers section, all other sections are not required. This
 * means the optional sections depend on information parsed from the polymers
 * section to do validation.
 */
std::shared_ptr<::RDKit::ROMol> helm_to_rdkit(const std::string& input_helm)
{
    if (get_helm_version(input_helm) != HelmVersion::V2) {
        throw std::invalid_argument("Only HELMV2.0 is supported.");
    }

    auto mol = helm::HelmParser().parse(input_helm);
    add_helm_specific_properties_to_mol(*mol);
    return mol;
}
