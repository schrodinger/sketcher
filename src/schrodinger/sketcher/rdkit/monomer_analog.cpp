#include "schrodinger/sketcher/rdkit/monomer_analog.h"

#include <algorithm>

#include <rdkit/GraphMol/Atom.h>

namespace schrodinger
{
namespace sketcher
{

// dR has NATURAL_ANALOG=R, so it's not a key in the analog map. Pull R's
// analogs (and R itself) into dR's popup for symmetry with R.
static const std::unordered_map<std::string, std::string>
    NA_ANALOG_LOOKUP_OVERRIDES = {{"dR", "R"}};

static rdkit_extensions::opt_string_t
get_peptide_natural_analog(std::string_view sym)
{
    return rdkit_extensions::MonomerDatabase::instance().getNaturalAnalog(
        std::string(sym), rdkit_extensions::ChainType::PEPTIDE);
}

// Default DB tags sugars as RNA, but custom DBs may register them under
// DNA — try RNA first, then fall back to DNA so either works.
static rdkit_extensions::opt_string_t
get_na_natural_analog(std::string_view sym)
{
    auto& db = rdkit_extensions::MonomerDatabase::instance();
    if (auto rna = db.getNaturalAnalog(std::string(sym),
                                       rdkit_extensions::ChainType::RNA);
        rna && !rna->empty()) {
        return rna;
    }
    return db.getNaturalAnalog(std::string(sym),
                               rdkit_extensions::ChainType::DNA);
}

bool all_monomers_have_type(const std::unordered_set<const RDKit::Atom*>& atoms,
                            MonomerType type)
{
    if (atoms.empty()) {
        return false;
    }
    return std::ranges::all_of(
        atoms, [type](const auto* a) { return get_monomer_type(a) == type; });
}

bool is_d_form(std::string_view sym)
{
    if (sym.size() < 2 || sym.front() != 'd') {
        return false;
    }
    auto sym_analog = get_peptide_natural_analog(sym);
    auto stripped_analog = get_peptide_natural_analog(sym.substr(1));
    return sym_analog && stripped_analog && !sym_analog->empty() &&
           *sym_analog == *stripped_analog;
}

bool is_deoxy_sugar(std::string_view sym)
{
    if (sym.size() < 2 || sym.front() != 'd') {
        return false;
    }
    auto sym_analog = get_na_natural_analog(sym);
    auto stripped_analog = get_na_natural_analog(sym.substr(1));
    return sym_analog && stripped_analog && !sym_analog->empty() &&
           *sym_analog == *stripped_analog;
}

std::optional<std::string> toggle_d_form_symbol(std::string_view sym)
{
    if (is_d_form(sym)) {
        // Stripping the prefix is safe: is_d_form only succeeds when the
        // L-form counterpart exists with a matching analog.
        return std::string(sym.substr(1));
    }
    auto candidate = std::string("d") + std::string(sym);
    auto candidate_analog = get_peptide_natural_analog(candidate);
    auto sym_analog = get_peptide_natural_analog(sym);
    if (candidate_analog && sym_analog && !candidate_analog->empty() &&
        *candidate_analog == *sym_analog) {
        return candidate;
    }
    return std::nullopt;
}

std::optional<std::string> toggle_sugar_symbol(std::string_view sym)
{
    if (is_deoxy_sugar(sym)) {
        return std::string(sym.substr(1));
    }
    auto candidate = std::string("d") + std::string(sym);
    auto candidate_analog = get_na_natural_analog(candidate);
    auto sym_analog = get_na_natural_analog(sym);
    if (candidate_analog && sym_analog && !candidate_analog->empty() &&
        *candidate_analog == *sym_analog) {
        return candidate;
    }
    return std::nullopt;
}

bool all_residues_are_d_form(
    const std::unordered_set<const RDKit::Atom*>& atoms)
{
    return std::ranges::all_of(atoms, [](const auto* a) {
        return is_d_form(get_monomer_res_name(a));
    });
}

std::vector<MonomerMutation>
build_d_form_batch(const std::unordered_set<const RDKit::Atom*>& atoms,
                   bool target_d_form)
{
    std::unordered_map<std::string, std::unordered_set<const RDKit::Atom*>>
        by_target;
    for (const auto* a : atoms) {
        const auto current = get_monomer_res_name(a);
        if (is_d_form(current) == target_d_form) {
            continue;
        }
        if (auto next = toggle_d_form_symbol(current)) {
            by_target[*next].insert(a);
        }
    }
    std::vector<MonomerMutation> batch;
    batch.reserve(by_target.size());
    for (auto& [sym, group] : by_target) {
        batch.push_back({std::move(group), sym});
    }
    return batch;
}

std::vector<MonomerMutation>
build_sugar_toggle_batch(const std::unordered_set<const RDKit::Atom*>& atoms,
                         const std::string& target_sym)
{
    std::unordered_set<const RDKit::Atom*> group;
    for (const auto* a : atoms) {
        if (get_monomer_res_name(a) == target_sym) {
            continue;
        }
        if (auto next = toggle_sugar_symbol(get_monomer_res_name(a));
            next && *next == target_sym) {
            group.insert(a);
        }
    }
    if (group.empty()) {
        return {};
    }
    return {{std::move(group), target_sym}};
}

std::unordered_set<const RDKit::Atom*>
atoms_needing_mutation(const std::unordered_set<const RDKit::Atom*>& atoms,
                       std::string_view target)
{
    std::unordered_set<const RDKit::Atom*> out;
    for (const auto* a : atoms) {
        if (get_monomer_res_name(a) != target) {
            out.insert(a);
        }
    }
    return out;
}

AnalogMap get_merged_na_analogs()
{
    auto& db = rdkit_extensions::MonomerDatabase::instance();
    auto base = db.getMonomersByNaturalAnalog(rdkit_extensions::ChainType::RNA);
    auto other =
        db.getMonomersByNaturalAnalog(rdkit_extensions::ChainType::DNA);
    for (const auto& [key, other_list] : other) {
        auto& merged = base[key];
        std::unordered_set<std::string> seen;
        for (const auto& m : merged) {
            seen.insert(m.symbol.value_or(""));
        }
        for (const auto& m : other_list) {
            if (seen.insert(m.symbol.value_or("")).second) {
                merged.push_back(m);
            }
        }
    }
    return base;
}

std::vector<rdkit_extensions::MonomerInfo> get_analogs_for_na_button(
    const std::string& symbol, const AnalogMap& analogs_by_na,
    const std::unordered_map<std::string, std::string>& standard_names)
{
    std::vector<rdkit_extensions::MonomerInfo> analogs;
    auto override_it = NA_ANALOG_LOOKUP_OVERRIDES.find(symbol);
    const auto& lookup_key = override_it != NA_ANALOG_LOOKUP_OVERRIDES.end()
                                 ? override_it->second
                                 : symbol;
    if (override_it != NA_ANALOG_LOOKUP_OVERRIDES.end()) {
        auto parent_name_it = standard_names.find(lookup_key);
        rdkit_extensions::MonomerInfo parent_entry;
        parent_entry.symbol = lookup_key;
        parent_entry.name = parent_name_it != standard_names.end()
                                ? parent_name_it->second
                                : lookup_key;
        analogs.push_back(std::move(parent_entry));
    }
    auto analog_it = analogs_by_na.find(lookup_key);
    if (analog_it != analogs_by_na.end()) {
        for (const auto& m : analog_it->second) {
            if (m.symbol.value_or("") != symbol) {
                analogs.push_back(m);
            }
        }
    }
    return analogs;
}

} // namespace sketcher
} // namespace schrodinger
