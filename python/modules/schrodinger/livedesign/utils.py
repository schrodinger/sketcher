"""
Utility functions used in LiveDesign.

Copyright Schrodinger, LLC. All rights reserved.
"""
import json

from schrodinger import rdkit_extensions
from schrodinger.protein.sequence import AMINO_ACIDS_3_TO_1
from schrodinger.structutils.analyze import read_seqres_from_ct


def convert_to_rdkit_with_seqres(mol_input, st):
    """
    Convert a Structure to RDKit molecule with SEQRES information attached.

    This function extracts SEQRES (sequence records) from a PDB structure and
    attaches them as a property on the resulting RDKit molecule. This enables
    downstream conversion to monomeric representation (via `toMonomeric()`) to
    produce complete biological sequences even when residues are missing or
    incomplete in the structure coordinates.

    :param mol_input: String representation of PDB
    :type mol: string
    :param st: Input structure, typically from a PDB file
    :type st: `schrodinger.structure.Structure`

    :return: RDKit molecule with SEQRES attached as a "SEQRES" property
        containing a JSON-encoded dictionary mapping chain IDs to lists of
        1-letter amino acid codes (or "X" for unknown residues)
    :rtype: `rdkit.Chem.rdchem.Mol`

    .. note::
        The SEQRES property is stored as JSON with the format chain: list[seq]:
        ``{"A": ["M", "A", "G", "S", ...], "B": ["V", "L", ...]}``

    .. warning::
        Non-standard amino acids not in `AMINO_ACIDS_3_TO_1` are mapped to "X".
    """
    chains, sequences = read_seqres_from_ct(st)
    chain_to_seq = {}
    for chain, seqres in zip(chains, sequences):
        seqres_to_list = seqres.split()
        seqres_as_list = [
            AMINO_ACIDS_3_TO_1.get(res, "X") for res in seqres_to_list
        ]
        chain_to_seq[chain] = seqres_as_list
    mol = rdkit_extensions.to_rdkit(mol_input, rdkit_extensions.Format.PDB)
    json_dump = json.dumps(chain_to_seq)
    mol.SetProp("SEQRES", json_dump)
    return mol
