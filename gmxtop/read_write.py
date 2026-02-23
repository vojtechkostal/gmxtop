from pathlib import Path
from typing import Optional, Dict, List

from .sections import MoleculeType, Topology, Define
from .parse.parsers import (
    parse_defaults,
    parse_atomtype, parse_moleculetype, parse_interaction_type,
    parse_atom, parse_interaction, parse_exclusion,
    parse_system, parse_molecule,
    INTERACTION_TYPES_MAPPINGS, MOLECULE_INTERACTION_MAPPINGS
)


def _group_by_func(items: List) -> Dict[int, List]:
    """Group items by their 'func' attribute."""
    groups: Dict[int, List] = {}
    for item in items:
        groups.setdefault(item.func, []).append(item)
    return groups


def _emit(items: list) -> list:
    """Emit items with proper handling of preprocessor directives."""
    lines: list[str] = []
    active_state: str | None = None  # e.g. "ifdef charge_mod", "else charge_mod"

    for item in items:
        state = item.ifdef_state

        if state == "free":
            if active_state is not None:
                lines.append("#endif")
                active_state = None
            lines.append(str(item))
            continue

        # entering a (new) conditional state
        if state != active_state:
            lines.append(item.ifdef)   # "#ifdef X" or "#else"
            active_state = state

        lines.append(str(item))

    if active_state is not None:
        lines.append("#endif")

    lines.append("")
    return lines


def read_topology(fn: Path, top: Topology) -> None:
    """
    Reads a topology file and populates the provided Topology object.

    Parameters
    ----------
    fn : Path
        Path to the topology file to be read.
    top : Topology
        An instance of the Topology class to populate with data from the file.

    Notes
    -----
    - The function processes the topology file line by line.
    - It handles conditional blocks (e.g., `#ifdef`, `#else`, `#endif`)
      and updates the `ifdef_state` accordingly.
    - Lines starting with `;` are treated as comments and ignored.
    - The function assumes the file is well-formed
      and does not handle all potential edge cases.

    Raises
    ------
    TypeError
        If the input parameters are of incorrect types.
    FileNotFoundError
        If the file specified by `fn` does not exist.
    ValueError
        If the file contains invalid or unexpected content.
    """

    # Validate input types
    if not isinstance(fn, Path):
        raise TypeError("The 'fn' parameter must be of type 'Path'.")
    if not fn.exists():
        raise FileNotFoundError(f"The file '{fn}' does not exist.")
    if not isinstance(top, Topology):
        raise TypeError(
            "The 'top' parameter must be an instance of the 'Topology' class."
        )

    top.source = fn.resolve()
    current_section: Optional[str] = None
    active_mol: Optional[MoleculeType] = None
    ifdef_block: Optional[str] = None
    ifdef_state = "free"

    with fn.open() as f:
        for raw_line in f:
            raw_line = raw_line.rstrip("\n")
            line = raw_line.split(";", 1)[0].strip()

            if not line:
                continue

            # Preprocessor directives
            if line.startswith("#"):
                directive, *rest = line.split()

                if line.startswith("#include"):
                    inc = rest[0].strip('"')
                    read_topology((fn.parent / inc).resolve(), top)

                elif line.startswith("#define"):
                    if len(rest) != 2:
                        raise ValueError(
                            f"Invalid #define directive in {top.source}: '{line}'"
                            "#define must be followed by name and value."
                        )
                    define = Define(directive=rest[0], argument=rest[1])
                    top.defines.add(define)

                elif line.startswith("#ifdef"):
                    if len(rest) != 1:
                        raise ValueError(
                            f"Invalid #ifdef directive in {top.source}: '{line}'"
                        )
                    ifdef_block = rest[0]
                    ifdef_state = f"ifdef {ifdef_block}"

                elif line.startswith("#else"):
                    if ifdef_block is None:
                        raise ValueError(
                            f"#else without matching #ifdef in {top.source}."
                        )
                    ifdef_state = f"else {ifdef_block}"

                elif line.startswith("#endif"):
                    ifdef_block = None
                    ifdef_state = "free"

                continue

            # Section header
            if line.startswith("["):
                current_section = line.strip("[]").strip()
                continue

            if not current_section:
                continue

            if current_section in MOLECULE_INTERACTION_MAPPINGS and active_mol is None:
                raise ValueError(
                    f"[ {current_section} ] section found before molecule type "
                    f"in {top.source}."
                )

            if current_section in INTERACTION_TYPES_MAPPINGS and active_mol is not None:
                raise ValueError(
                    f"[ {current_section} ] section found inside molecule type "
                    f"{active_mol.name} in {top.source}."
                    f" Parameter types must be defined globally."
                )

            # Data line
            parts = line.split()

            # ---- global definitions ----
            # defaults
            if current_section == "defaults":
                top.defaults = parse_defaults(parts, top)
                continue

            # atom types
            if current_section == "atomtypes":
                paramtype = parse_atomtype(parts, top, ifdef_state)
                top.atomtypes.append(paramtype)
                continue

            # interaction types
            if current_section in INTERACTION_TYPES_MAPPINGS:
                builder, section_cls = INTERACTION_TYPES_MAPPINGS[current_section]
                paramtype = parse_interaction_type(
                    parts,
                    builder,
                    section_cls,
                    top,
                    current_section,
                    ifdef_state
                )
                getattr(top, current_section).append(paramtype)
                continue

            # ---- molecule definitions ----
            if current_section == "moleculetype":
                active_mol = parse_moleculetype(parts, top)
                top.moleculetypes.append(active_mol)
                continue

            # atoms
            if current_section == "atoms":
                param = parse_atom(parts, top, active_mol, ifdef_state)
                active_mol.atoms.append(param)
                continue

            # interactions
            if current_section in MOLECULE_INTERACTION_MAPPINGS:

                # 'dummies3' are treated as 'virtual_sites3' in this parser
                if current_section == 'dummies3':
                    current_section = 'virtual_sites3'

                builder, section_cls = MOLECULE_INTERACTION_MAPPINGS[current_section]
                param = parse_interaction(
                    parts,
                    builder,
                    section_cls,
                    top,
                    active_mol,
                    current_section,
                    ifdef_state
                )
                getattr(active_mol, current_section).extend(param)
                continue

            if current_section == "exclusions":
                excluded = parse_exclusion(parts, active_mol, ifdef_state)
                active_mol.exclusions.append(excluded)
                continue

            # ---- system definition ----
            if current_section == "system":
                # often one-line free text
                top.system = parse_system(parts, top)
                continue

            if current_section == "molecules":
                mol, count = parse_molecule(parts, top)
                top.molecules[mol.name] = (mol, count)
                continue

            # Unsupported section in this parser: ignore.
            # raise ValueError(f"Unsupported section '{current_section}' in topology.")
            continue

        if top.defaults is None:
            raise ValueError("No [ defaults ] section found in topology.")


def write_topology(top: Topology, fn_out: Path, overwrite: bool = False) -> None:
    """Writes the Topology object to a file."""

    fn_out = Path(fn_out)
    if fn_out.exists() and not overwrite:
        raise FileExistsError(
            f"File '{fn_out}' already exists. "
            f"Use overwrite=True to overwrite."
        )

    # Collect lines to write
    lines: list[str] = []

    # --- defines ---
    if top.defines:
        for define in top.defines:
            lines.append(str(define))
        lines.append("")

    # --- defaults ---
    lines += [top.defaults.header, top.defaults.legend, str(top.defaults), ""]

    # --- atomtypes ---
    lines += [top.atomtypes[0].header, top.atomtypes[0].legend]
    used_atomtypes = list(dict.fromkeys(
        atom.type
        for mol, _ in top.molecules.values()
        for atom in mol.atoms
    ))

    lines += _emit(used_atomtypes)

    # --- nonbonded ---
    nonbond_params = getattr(top, "nonbond_params", [])
    if nonbond_params:
        lines += [nonbond_params[0].header, nonbond_params[0].legend]
        used_nonboned_params = list({
            nb
            for nb in nonbond_params
            if nb.ai in used_atomtypes and nb.aj in used_atomtypes
        })
        # groups_ifdef = _group_by_ifdef(used_nonboned_params)
        lines += _emit(used_nonboned_params)

    # --- moleculetypes ---
    for mol, _ in top.molecules.values():

        # moleculetype
        lines += [mol.header, mol.legend, str(mol), ""]

        # atoms
        lines += [mol.atoms[0].header, mol.atoms[0].legend]
        lines += _emit(mol.atoms)

        # interactions
        for section in MOLECULE_INTERACTION_MAPPINGS:
            params = getattr(mol, section, None)
            if not params:
                continue

            # group by different funcions
            groups_func = _group_by_func(params)
            for func_params in groups_func.values():
                lines += [func_params[0].header, func_params[0].legend]
                # groups_ifdef = _group_by_ifdef(func_params)
                lines += _emit(func_params)

        # exclusions are treated separately
        if mol.exclusions:
            lines += [mol.exclusions[0].header, mol.exclusions[0].legend]
            for exclusion in mol.exclusions:
                lines.append(str(exclusion))
            lines.append("")

    # --- system ---
    lines += [top.system.header, top.system.legend, str(top.system), ""]

    # --- molecules ---
    lines += ["[ molecules ]", f"; {'name':>1} {'count':>14}"]
    for mol_name, (mol, count) in top.molecules.items():
        lines.append(f"{mol_name:<10} {count:>10}")

    fn_out.write_text("\n".join(lines) + "\n")
