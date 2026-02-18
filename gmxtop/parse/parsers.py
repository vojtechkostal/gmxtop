from typing import List, Dict, Tuple, Any, Callable, Optional, get_type_hints, get_args

from ..sections import (
    Defaults,
    AtomType, BondType, PairType, AngleType, DihedralType, ConstraintType, NonBondParam,
    MoleculeType,
    Atom, Bond, Pair, PairNB, Angle, Dihedral,
    Exclusion, Settle, System,
    VirtualSite1, VirtualSite2, VirtualSite3, VirtualSite4, VirtualSiteN,
    PositionRestraint, DistanceRestraint, DihedralRestraint,
    OrientationRestraint, AngleRestraint, AngleRestraintZ,
    Topology
)

from .helpers import lookup_paramtype
from ..specs.interactions import (
    BONDS, PAIRS, PAIRS_NB, ANGLES, DIHEDRALS,
    CONSTRAINTS, SETTLES, NONBOND_PARAMS,
    VIRTUAL_SITES1, VIRTUAL_SITES2, VIRTUAL_SITES3, VIRTUAL_SITES4, VIRTUAL_SITESN,
    POSITION_RESTRAINTS, DISTANCE_RESTRAINTS, DIHEDRAL_RESTRAINTS,
    ORIENTATION_RESTRAINTS, ANGLE_RESTRAINTS, ANGLE_RESTRAINTS_Z
)


def resolve_input(cls, parts: List[str], top: Topology) -> None:
    """Resolve and convert input parts to appropriate
    types based on class type hints."""

    type_hints = {
        k: v
        for k, v in get_type_hints(cls).items()
        if k not in ['params', "MODIFIABLE"]
    }

    inputs = {}
    for part, (field_name, field_type) in zip(parts, type_hints.items()):
        try:
            inputs[field_name] = field_type(part)
        except Exception:
            if field_name not in ("ifdef_state", "type", "residue"):
                define_directives = [d.directive for d in top.defines]
                if part.strip("-") not in define_directives:
                    raise ValueError(
                        f"Invalid value '{part}' for field '{field_name}' "
                        f"in section {cls.__name__} in topology {top.source}."
                    )
            inputs[field_name] = part

    return inputs


def parse_defaults(
    parts: List[str],
    top: Topology,
    ifdef_state: Optional[str] = None
) -> Defaults:

    # check for defaults duplicates
    if top.defaults is not None:
        raise ValueError("Multiple [ defaults ] sections found.")

    # validate inputs
    parts += [ifdef_state]  # include ifdef_state in defaults
    inputs = resolve_input(Defaults, parts, top)

    if inputs["nbfunc"] not in (1, 2):
        raise ValueError(
            f"Invalid nbfunc in [ defaults ] in {top.source}: "
            f"expected 1 or 2, got {inputs['nbfunc']}"
        )

    if inputs["comb_rule"] not in (1, 2, 3):
        raise ValueError(
            f"Invalid comb_rule in [ defaults ] in {top.source}: "
            f"expected 1, 2, or 3, got {inputs['comb_rule']}"
        )

    if inputs["gen_pairs"] not in ('yes', 'no'):
        raise ValueError(
            f"Invalid gen_pairs in [ defaults ] in {top.source}: "
            f"expected 'yes' or 'no', got '{inputs['gen_pairs']}'"
        )

    return Defaults(**inputs)


def parse_atomtype(
    parts: List[str],
    top: Topology,
    ifdef_state: Optional[str] = None,
) -> AtomType:

    parts += [ifdef_state]  # include ifdef_state in atomtype
    inputs = resolve_input(AtomType, parts, top)

    if inputs['ptype'] not in ("A", "S", "V", "D"):
        raise ValueError(
            f"Invalid ptype of atomtype {inputs['name']} "
            f"in [ atomtypes ] in {top.source}: "
            f"expected 'A', 'S', 'V', or 'D', got '{inputs['ptype']}'"
        )

    atomtype = AtomType(**inputs)

    # Check for duplicates
    for at in top.atomtypes:
        if at.name == atomtype.name and at.ifdef_state == atomtype.ifdef_state:
            raise ValueError(
                f"Duplicate atomtype '{atomtype.name}' "
                f"found in topology {top.source}."
            )

    return atomtype


def parse_moleculetype(
    parts: List[str],
    top: Topology,
    ifdef_state: Optional[str] = None
) -> MoleculeType:

    inputs = {
        "name": parts[0],
        "nrexcl": int(parts[1]),
        'ifdef_state': ifdef_state
    }

    moleculetype = MoleculeType(**inputs)

    # check for duplicates
    for mt in top.moleculetypes:
        if mt.name == moleculetype.name and mt.ifdef_state == moleculetype.ifdef_state:
            raise ValueError(
                f"Duplicate moleculetype '{moleculetype.name}' "
                f"found in topology {top.source}."
            )

    return moleculetype


def parse_atom(
    parts: List[str],
    top: Topology,
    mol: MoleculeType,
    ifdef_state: Optional[str] = None
) -> Atom:

    parts += [ifdef_state]  # include ifdef_state in atom
    inputs = resolve_input(Atom, parts, top)

    if mol.get_atom_by_idx(inputs['nr']) is not None and ifdef_state == "free":
        raise ValueError(
            f"Duplicate atom index '{inputs['nr']}' "
            f"in molecule '{mol.name}' in topology {top.source}"
        )

    # Resolve atomtype
    atomtype_name = inputs['type']
    atomtype = next((at for at in top.atomtypes if at.name == atomtype_name), None)
    if atomtype is None:
        raise ValueError(
            f"Atom type '{atomtype_name}' "
            f"not found in topology {top.source}."
        )

    inputs["residue"] = mol
    inputs["type"] = atomtype

    return Atom(**inputs)


def parse_exclusion(
    parts: List[str],
    mol: MoleculeType,
    ifdef_state: Optional[str] = None
) -> Exclusion:

    # Resolve atoms
    excluded = []
    for a in parts:
        if a.isdigit():
            idx = int(a)
            atom = mol.get_atom_by_idx(idx)
        else:
            raise ValueError(f"Invalid atom token '{a}' in exclusions.")
        excluded.append(atom)

    return Exclusion(excluded, ifdef_state=ifdef_state)


def parse_system(
    parts: List[str],
    top: Topology,
) -> System:

    description = " ".join(parts) if parts else "system"

    if top.system is not None:
        raise ValueError("Multiple [ system ] sections found.")

    system = System(description)
    return system


def parse_molecule(
    parts: List[str],
    top: Topology,
) -> Tuple[MoleculeType, int]:

    moleculetype = next((mt for mt in top.moleculetypes if mt.name == parts[0]), None)
    if moleculetype is None:
        raise ValueError(
            f"Molecule type '{parts[0]}' "
            f"not found in topology {top.source}."
        )
    return moleculetype, int(parts[1])


def parse_interaction_type(
    parts: List[str],
    builder: object,
    section_cls: Callable,
    top: Topology,
    section: str = "",
    ifdef_state: Optional[str] = None
) -> object:

    n_atom = builder.n_atoms

    # Resolve atomtypes
    atomtypes = []
    for a in parts[:n_atom]:
        if isinstance(a, str):
            if a == "X":
                at = AtomType(
                    name="X",
                    atnum=0,
                    mass=0.0,
                    charge=0.0,
                    ptype="A",
                    sigma=0.0,
                    epsilon=0.0
                )
                atomtypes.append(at)
                continue
            for at in top.atomtypes:
                if at.name == a:
                    atomtypes.append(at)
                    break
        else:
            raise ValueError(f"Atomtype {a} in {section} not found in {top.source}.")

    # Resolve function
    func = int(parts[n_atom])

    # Resolve parameters
    param_tokens = parts[n_atom + 1:]
    params = builder.parse(
        func,
        param_tokens, ctx=f"section {section} in topology {top.source}")
    if not param_tokens:
        raise ValueError(
            f"Parameters must be provided for {section} in {top.source} "
            f"with atomtypes ({', '.join(a.name for a in atomtypes)})."
        )

    return section_cls(*atomtypes, func, params, ifdef_state=ifdef_state)


def parse_interaction(
    parts: List[str],
    builder: object,
    section_cls: object,
    top: Topology,
    mol: MoleculeType,
    section: str = "",
    ifdef_state: Optional[str] = None
) -> list[Bond, Pair, Angle, Dihedral, Exclusion]:

    # Extract atoms, function, and parameters from parts using the provided builder.
    n_atoms = builder.n_atoms

    # Resolve atoms
    atoms = []
    for a in parts[:n_atoms]:
        if a.isdigit():
            idx = int(a)
            atom = mol.get_atom_by_idx(idx)
        else:
            raise ValueError(f"Invalid atom token '{a}' in {section} of {top.source}.")
        atoms.append(atom)

    func = int(parts[n_atoms])

    # Resolve parameters
    param_tokens = parts[n_atoms + 1:]
    ctx = f"section {section} in {mol.name}"
    params = builder.parse(func, param_tokens, ctx=ctx)

    if not params:
        is_pair = section == "pairs" or section == "pairs_nb"
        if top.defaults.gen_pairs == "yes" and is_pair:
            return [section_cls(*atoms, func, params)]
        paramtypes = lookup_paramtype(top, *atoms, func=func, pair=is_pair)
        return [section_cls(*atoms, func, pt.params) for pt in paramtypes]

    return [section_cls(*atoms, func, params, ifdef_state=ifdef_state)]


INTERACTION_TYPES_MAPPINGS: Dict[str, Any] = {
    "bondtypes": (BONDS, BondType),
    "pairtypes": (PAIRS, PairType),
    "pairtypes_nb": (PAIRS_NB, PairNB),
    "angletypes": (ANGLES, AngleType),
    "dihedraltypes": (DIHEDRALS, DihedralType),
    "constrainttypes": (CONSTRAINTS, ConstraintType),
    "nonbond_params": (NONBOND_PARAMS, NonBondParam),
}

MOLECULE_INTERACTION_MAPPINGS: Dict[str, Any] = {
    "bonds": (BONDS, Bond),
    "pairs": (PAIRS, Pair),
    "pairs_nb": (PAIRS_NB, PairNB),
    "angles": (ANGLES, Angle),
    "dihedrals": (DIHEDRALS, Dihedral),
    "settles": (SETTLES, Settle),
    "virtual_sites1": (VIRTUAL_SITES1, VirtualSite1),
    "virtual_sites2": (VIRTUAL_SITES2, VirtualSite2),
    "virtual_sites3": (VIRTUAL_SITES3, VirtualSite3),
    "virtual_sites4": (VIRTUAL_SITES4, VirtualSite4),
    "virtual_sitesn": (VIRTUAL_SITESN, VirtualSiteN),
    "dummies1": (VIRTUAL_SITES1, VirtualSite1),
    "dummies2": (VIRTUAL_SITES2, VirtualSite2),
    "dummies3": (VIRTUAL_SITES3, VirtualSite3),
    "dummies4": (VIRTUAL_SITES4, VirtualSite4),
    "dummiesn": (VIRTUAL_SITESN, VirtualSiteN),
    "position_restraints": (POSITION_RESTRAINTS, PositionRestraint),
    "distance_restraints": (DISTANCE_RESTRAINTS, DistanceRestraint),
    "dihedral_restraints": (DIHEDRAL_RESTRAINTS, DihedralRestraint),
    "orientation_restraints": (ORIENTATION_RESTRAINTS, OrientationRestraint),
    "angle_restraints": (ANGLE_RESTRAINTS, AngleRestraint),
    "angle_restraints_z": (ANGLE_RESTRAINTS_Z, AngleRestraintZ),
}
