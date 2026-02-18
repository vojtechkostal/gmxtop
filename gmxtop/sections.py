from __future__ import annotations

from pathlib import Path
from dataclasses import dataclass, field
from typing import (
    Any,
    Dict,
    List,
    Optional,
    ClassVar,
    FrozenSet,
    Mapping,
    get_type_hints,
    Tuple,
    Set
)


class MixIn:
    """
    Mix-in class to provide in-place parameter updates.
    """

    MODIFIABLE: ClassVar[FrozenSet[str]] = frozenset()

    def update(self, **kwargs: Any) -> None:
        """
        Update parameters in-place.

        Policy:
        - If `self.params` exists and is a dict:
          only overwrite existing keys in that dict.
        - Otherwise: only overwrite attributes listed in `MODIFIABLE`.
        """
        params = getattr(self, "params", None)

        if isinstance(params, dict):
            unknown = kwargs.keys() - params.keys()
            if unknown:
                raise KeyError(
                    f"Cannot add new parameter(s) {sorted(unknown)}; "
                    f"allowed keys are {sorted(params.keys())}"
                )
            params.update(kwargs)
            return

        # Attribute-based update
        allowed = type(self).MODIFIABLE
        unknown = kwargs.keys() - allowed
        if unknown:
            raise KeyError(
                f"Parameter(s) {sorted(unknown)} are not modifiable "
                f"for {type(self).__name__}; allowed are {sorted(allowed)}"
            )
        for k, v in kwargs.items():
            setattr(self, k, v)


class Description:
    """Base class for all sections in a topology file.
    Provides methods for rendering section headers, legends, and values."""

    SECTION = ClassVar[str]
    W = 10
    W_FIRST = 8
    W_FLOAT = 15
    D = 4
    D_CHARGE = 10
    D_LJ = 6

    @property
    def header(self) -> str:
        """Section header string."""
        return f"[ {self.SECTION} ]"

    @property
    def args(self) -> Dict[str, Any]:
        """Extracts the fields and their values/types for rendering."""
        hints = get_type_hints(type(self))

        out: Dict[str, Tuple[Any, type]] = {}
        for k, t in hints.items():
            if k in ["MODIFIABLE", "SECTION"]:
                continue

            val = getattr(self, k)

            if t in [float, int, str]:
                out[k] = (val, t)
            elif t in [AtomType, MoleculeType]:
                out[k] = (val.name, str)
            elif t == Atom:
                out[k] = (val.nr, int)
            elif k == "excluded":
                names = " ".join(str(a.nr) for a in val)
                out[k] = (names, str)
            elif k == "params":
                for pk, pv in val.items():
                    out[pk] = (pv, type(pv))
            else:
                continue

        return out

    def _fmt(self, name, type_: type) -> int:
        """Determine format string for a given parameter based on its name and type."""
        if type_ == float:
            width = self.W_FLOAT
            if name == "charge":
                decimals = self.D_CHARGE
            elif name in ("sigma", "epsilon"):
                decimals = self.D_LJ
            else:
                decimals = self.D

            fmt = f'>{width}.{decimals}f'
        else:
            width = self.W
            fmt = f'>{width}'

        return fmt

    def _render(self, *, values: bool) -> str:
        """Render either the legend (parameter names) or values."""

        chunks: list[str] = []
        for i, (key, (val, typ)) in enumerate(self.args.items()):
            if i == 0:
                if values:
                    chunks.append(f"{val:>{self.W_FIRST}}")
                else:
                    chunks.append(f"; {key:>{self.W_FIRST - 2}}")
                continue

            fmt = self._fmt(key, typ)
            if values:
                try:
                    chunks.append(f"{val:{fmt}}")
                except ValueError:
                    # Fallback for non-numeric values
                    chunks.append(f"{str(val):>{self.W}}")
            else:
                fmt = fmt.split('.')[0]  # remove decimal part for legend
                chunks.append(f"{key:{fmt}}")

        return "".join(chunks)

    @property
    def ifdef(self) -> Optional[str]:
        """Return the ifdef directive string for this section."""
        if self.ifdef_state:
            if "ifdef" in self.ifdef_state:
                return f"#{self.ifdef_state}"
            elif "else" in self.ifdef_state:
                return "#else"
        return self.ifdef_state

    @property
    def legend(self) -> str:
        """Return the legend string (parameter names)."""
        return self._render(values=False)

    def __str__(self):
        """Return the values string (parameter values)."""
        return self._render(values=True)


@dataclass(slots=True)
class Defaults(Description):
    nbfunc: int
    comb_rule: int
    gen_pairs: str
    fudgeLJ: float
    fudgeQQ: float
    ifdef_state: Optional[str] = "free"

    SECTION = 'defaults'


@dataclass(slots=True)
class Define(MixIn):
    directive: str
    argument: Optional[str | int | float]

    MODIFIABLE: ClassVar[FrozenSet[str]] = frozenset({"argument"})

    def __hash__(self):
        return hash(self.directive)

    def __str__(self) -> str:
        if self.argument is not None:
            return f"#define {self.directive} {self.argument}"


@dataclass(slots=True)
class AtomType(MixIn, Description):
    name: str
    atnum: int
    mass: float
    charge: float
    ptype: str
    sigma: float
    epsilon: float
    ifdef_state: Optional[str] = "free"

    MODIFIABLE: ClassVar[FrozenSet[str]] = frozenset({"sigma", "epsilon"})
    SECTION = 'atomtypes'

    def __hash__(self):
        return hash(self.name)


@dataclass(slots=True)
class BondType(MixIn, Description):
    ai: AtomType
    aj: AtomType
    func: int
    params: Mapping[str, float | int | str] = field(default_factory=dict)
    ifdef_state: Optional[str] = "free"

    SECTION = 'bondtypes'


@dataclass(slots=True)
class PairType(MixIn, Description):
    ai: AtomType
    aj: AtomType
    func: int
    params: Mapping[str, float | int | str] = field(default_factory=dict)
    ifdef_state: Optional[str] = "free"

    SECTION = 'pairtypes'


@dataclass(slots=True)
class AngleType(MixIn, Description):
    ai: AtomType
    aj: AtomType
    ak: AtomType
    func: int
    params: Mapping[str, float | int | str] = field(default_factory=dict)
    ifdef_state: Optional[str] = "free"

    SECTION = 'angletypes'


@dataclass(slots=True)
class DihedralType(MixIn, Description):
    ai: AtomType
    aj: AtomType
    ak: AtomType
    al: AtomType
    func: int
    params: Mapping[str, float | int | str] = field(default_factory=dict)
    ifdef_state: Optional[str] = "free"

    SECTION = 'dihedraltypes'


@dataclass(slots=True)
class ConstraintType(MixIn, Description):
    ai: AtomType
    aj: AtomType
    func: int
    params: Mapping[str, float | int | str] = field(default_factory=dict)
    ifdef_state: Optional[str] = "free"

    SECTION = 'constrainttypes'


@dataclass(slots=True)
class NonBondParam(MixIn, Description):
    ai: AtomType
    aj: AtomType
    func: int
    params: Mapping[str, float | int | str] = field(default_factory=dict)
    ifdef_state: Optional[str] = "free"

    SECTION = 'nonbond_params'

    def __hash__(self):
        return hash((self.ai.name, self.aj.name))


@dataclass(slots=True)
class Atom(MixIn, Description):
    nr: int
    type: AtomType
    resnr: int
    residue: MoleculeType
    name: str
    cgnr: int
    charge: float
    mass: float
    ifdef_state: Optional[str] = "free"

    SECTION = 'atoms'
    MODIFIABLE: ClassVar[FrozenSet[str]] = frozenset({"charge"})


@dataclass(slots=True)
class Bond(MixIn, Description):
    ai: Atom
    aj: Atom
    func: int
    params: Mapping[str, float | int | str] = field(default_factory=dict)
    ifdef_state: Optional[str] = "free"

    SECTION = 'bonds'


@dataclass(slots=True)
class Pair(MixIn, Description):
    ai: Atom
    aj: Atom
    func: int
    params: Mapping[str, float | int | str] = field(default_factory=dict)
    ifdef_state: Optional[str] = "free"

    SECTION = 'pairs'


@dataclass(slots=True)
class PairNB(MixIn, Description):
    ai: Atom
    aj: Atom
    func: int
    params: Mapping[str, float | int | str] = field(default_factory=dict)
    ifdef_state: Optional[str] = "free"

    SECTION = 'pairs_nb'


@dataclass(slots=True)
class Angle(MixIn, Description):
    ai: Atom
    aj: Atom
    ak: Atom
    func: int
    params: Mapping[str, float | int | str] = field(default_factory=dict)
    ifdef_state: Optional[str] = "free"

    SECTION = 'angles'


@dataclass(slots=True)
class Dihedral(MixIn, Description):
    ai: Atom
    aj: Atom
    ak: Atom
    al: Atom
    func: int
    params: Mapping[str, float | int | str] = field(default_factory=dict)
    ifdef_state: Optional[str] = "free"

    SECTION = 'dihedrals'

    def __repr__(self):
        return (
            "Dihedral between atoms "
            f"{self.ai.name}, {self.aj.name}, {self.ak.name}, {self.al.name}"
        )


@dataclass(slots=True)
class Exclusion(Description):
    excluded: List[Atom]  # variable length per line
    ifdef_state: Optional[str] = "free"

    SECTION = 'exclusions'


@dataclass(slots=True)
class Settle(MixIn, Description):
    ai: Atom
    func: int
    params: Dict[str, float] = field(default_factory=dict)
    ifdef_state: Optional[str] = "free"

    SECTION = 'settles'


@dataclass(slots=True)
class VirtualSite1(MixIn, Description):
    ai: Atom
    aj: Atom
    func: int
    params: Mapping[str, float | int | str] = field(default_factory=dict)
    ifdef_state: Optional[str] = "free"

    SECTION = 'virtual_sites1'


@dataclass(slots=True)
class VirtualSite2(MixIn, Description):
    ai: Atom
    aj: Atom
    ak: Atom
    func: int
    params: Mapping[str, float | int | str] = field(default_factory=dict)
    ifdef_state: Optional[str] = "free"

    SECTION = 'virtual_sites2'


@dataclass(slots=True)
class VirtualSite3(MixIn, Description):
    ai: Atom
    aj: Atom
    ak: Atom
    al: Atom
    func: int
    params: Mapping[str, float | int | str] = field(default_factory=dict)
    ifdef_state: Optional[str] = "free"

    SECTION = 'virtual_sites3'


@dataclass(slots=True)
class VirtualSite4(MixIn, Description):
    ai: Atom
    aj: Atom
    ak: Atom
    al: Atom
    am: Atom
    func: int
    params: Mapping[str, float | int | str] = field(default_factory=dict)
    ifdef_state: Optional[str] = "free"

    SECTION = 'virtual_sites4'


@dataclass(slots=True)
class VirtualSiteN(MixIn, Description):
    ai: Atom
    func: int
    params: Mapping[str, float | int | str] = field(default_factory=dict)
    ifdef_state: Optional[str] = "free"

    SECTION = 'virtual_sitesn'


@dataclass(slots=True)
class PositionRestraint(MixIn, Description):
    ai: Atom
    func: int
    params: Mapping[str, float | int | str] = field(default_factory=dict)
    ifdef_state: Optional[str] = "free"

    SECTION = 'position_restraints'


@dataclass(slots=True)
class DistanceRestraint(MixIn, Description):
    ai: Atom
    aj: Atom
    func: int
    params: Mapping[str, float | int | str] = field(default_factory=dict)
    ifdef_state: Optional[str] = "free"

    SECTION = 'distance_restraints'


@dataclass(slots=True)
class DihedralRestraint(MixIn, Description):
    ai: Atom
    aj: Atom
    ak: Atom
    al: Atom
    func: int
    params: Mapping[str, float | int | str] = field(default_factory=dict)
    ifdef_state: Optional[str] = "free"

    SECTION = 'dihedral_restraints'


@dataclass(slots=True)
class OrientationRestraint(MixIn, Description):
    ai: Atom
    aj: Atom
    func: int
    params: Mapping[str, float | int | str] = field(default_factory=dict)
    ifdef_state: Optional[str] = "free"

    SECTION = 'orientation_restraints'


@dataclass(slots=True)
class AngleRestraint(MixIn, Description):
    ai: Atom
    aj: Atom
    ak: Atom
    al: Atom
    func: int
    params: Mapping[str, float | int | str] = field(default_factory=dict)
    ifdef_state: Optional[str] = "free"

    SECTION = 'angle_restraints'


@dataclass(slots=True)
class AngleRestraintZ(MixIn, Description):
    ai: Atom
    aj: Atom
    func: int
    params: Mapping[str, float | int | str] = field(default_factory=dict)
    ifdef_state: Optional[str] = "free"

    SECTION = 'angle_restraints_z'


@dataclass(slots=True)
class System(Description):
    description: str = "system"

    SECTION = 'system'


@dataclass(slots=True)
class MoleculeType(Description):
    name: str
    nrexcl: int
    atoms: List[Atom] = field(default_factory=list)
    bonds: List[Bond] = field(default_factory=list)
    pairs: List[Pair] = field(default_factory=list)
    pairs_nb: List[PairNB] = field(default_factory=list)
    angles: List[Angle] = field(default_factory=list)
    dihedrals: List[Dihedral] = field(default_factory=list)
    exclusions: List[Exclusion] = field(default_factory=list)
    settles: List[Settle] = field(default_factory=list)
    virtual_sites1: List[VirtualSite1] = field(default_factory=list)
    virtual_sites2: List[VirtualSite2] = field(default_factory=list)
    virtual_sites3: List[VirtualSite3] = field(default_factory=list)
    virtual_sites4: List[VirtualSite4] = field(default_factory=list)
    virtual_sitesn: List[VirtualSiteN] = field(default_factory=list)
    position_restraints: List[PositionRestraint] = field(default_factory=list)
    distance_restraints: List[DistanceRestraint] = field(default_factory=list)
    dihedral_restraints: List[DihedralRestraint] = field(default_factory=list)
    orientation_restraints: List[OrientationRestraint] = field(default_factory=list)
    angle_restraints: List[AngleRestraint] = field(default_factory=list)
    angle_restraints_z: List[AngleRestraintZ] = field(default_factory=list)
    ifdef_state: Optional[str] = "free"

    SECTION = 'moleculetype'

    def get_atom_by_idx(self, idx: int) -> Optional[Atom]:
        for atom in self.atoms:
            if atom.nr == idx:
                return atom
        return None

    @property
    def _vsite_sections(self):
        return {
            attr: getattr(self, attr)
            for attr in self.__dataclass_fields__
            if attr.startswith("virtual_sites")
        }

    @property
    def _connention_sections(self):
        return {
            "bonds", "pairs", "angles", "dihedrals",
            "position_restraints", "distance_restraints",
            "dihedral_restraints", "orientation_restraints",
            "angle_restraints", "angle_restraints_z"
        }

    def remove_vsites(self) -> None:
        from .parse.helpers import remove_vsites
        remove_vsites(self)

    def __repr__(self) -> str:
        return f"<MoleculeType name={self.name} nrexcl={self.nrexcl}>"


@dataclass(slots=True)
class Topology:
    """Representation of a molecular topology file."""
    source: Path
    defaults: Optional[Defaults] = None
    atomtypes: List[AtomType] = field(default_factory=list)
    nonbond_params: List[NonBondParam] = field(default_factory=list)
    bondtypes: List[BondType] = field(default_factory=list)
    pairtypes: List[PairType] = field(default_factory=list)
    angletypes: List[AngleType] = field(default_factory=list)
    dihedraltypes: List[DihedralType] = field(default_factory=list)
    moleculetypes: List[MoleculeType] = field(default_factory=list)

    system: Optional[System] = None
    molecules: Dict[MoleculeType, int] = field(default_factory=dict)

    defines: Set[Define] = field(default_factory=set)

    @property
    def residues(self) -> List[MoleculeType]:
        """List of all residues in the topology, expanded by their counts."""
        return [mol for mol, count in self.molecules.values() for _ in range(count)]

    @property
    def atoms(self) -> List[Atom]:
        """List of all atoms in the topology."""
        return [atom for mol in self.residues for atom in mol.atoms]

    def __post_init__(self) -> None:
        self.source = Path(self.source).resolve()
        from .read_write import read_topology
        read_topology(self.source, self)

    def write(self, fn: str | Path, **kwargs) -> None:
        from .read_write import write_topology
        write_topology(self, Path(fn).resolve(), **kwargs)

    def __repr__(self) -> str:
        return (
            f"<Topology source={self.source} "
            f"including {len(self.molecules)} molecules>"
        )
