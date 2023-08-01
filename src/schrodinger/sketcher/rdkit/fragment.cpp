#include "schrodinger/sketcher/rdkit/fragment.h"

#include <algorithm>
#include <stdexcept>

#include "schrodinger/sketcher/rdkit/rgroup.h"

namespace schrodinger
{
namespace sketcher
{

RDKit::Atom* const get_attachment_point_heavy_atom(const RDKit::ROMol& fragment)
{
    std::vector<RDKit::Atom*> all_aps;
    auto all_atoms = fragment.atoms();
    std::copy_if(all_atoms.begin(), all_atoms.end(),
                 std::back_inserter(all_aps), is_attachment_point);
    if (all_aps.size() != 1) {
        throw std::runtime_error(
            "Fragments must have exactly one attachment point.");
    }
    auto* ap_atom = all_aps.front();
    return *fragment.atomNeighbors(ap_atom).begin();
}

RDKit::Conformer translate_fragment(const RDKit::ROMol& fragment,
                                    const RDGeom::Point3D& point)
{
    auto* ap = get_attachment_point_heavy_atom(fragment);
    // note that this makes a copy of the conformer
    auto conf = fragment.getConformer();
    auto offset = point - conf.getAtomPos(ap->getIdx());
    for (auto& coords : conf.getPositions()) {
        coords += offset;
    }
    return conf;
}

RDKit::Conformer align_fragment_with_atom(const RDKit::ROMol& fragment,
                                          const RDKit::Atom* const atom)
{
    // TODO: real implementation
    return fragment.getConformer();
}

RDKit::Conformer align_fragment_with_bond(const RDKit::ROMol& fragment,
                                          const RDKit::Bond* const bond)
{
    // TODO: real implementation
    return fragment.getConformer();
}

} // namespace sketcher
} // namespace schrodinger
