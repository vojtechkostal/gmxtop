# This module defines interaction specifications for molecular simulations.

from .core import Params, InteractionSpec

BONDS = InteractionSpec(
    n_atoms=2,
    funcs={
        1: Params(
            names=("b0", "kb"),
            parsers=(float, float),
            desc="harmonic bond"
        ),
        2: Params(
            names=("b0", "kb"),
            parsers=(float, float),
            desc="G96 bond"
        ),
        3: Params(
            names=("b0", "D", "beta"),
            parsers=(float, float, float),
            desc="Morse bond"
        ),
        4: Params(
            names=("b0", "C2", "C3"),
            parsers=(float, float, float),
            desc="cubic bond"
        ),
        5: Params(
            names=(),
            parsers=(),
            desc="connection (no parameters)"
        ),
        6: Params(
            names=("b0", "kb"),
            parsers=(float, float),
            desc="harmonic potential"
        ),
        7: Params(
            names=("bm", "kb"),
            parsers=(float, float),
            desc="FENE bond"
        ),
        8: Params(
            names=("table", "k"),
            parsers=(int, float),
            desc="tabulated bond"
        ),
        9: Params(
            names=("table", "k"),
            parsers=(int, float),
            desc="tabulated bond (variant)"
        ),
        10: Params(
            names=("low", "up1", "up2", "kdr"),
            parsers=(float, float, float, float),
            desc="restraint",
        ),
    },
)


PAIRS = InteractionSpec(
    n_atoms=2,
    funcs={
        1: Params(
            names=("sigma", "epsilon"),
            parsers=(float, float),
            desc="Lennard-Jones pair"
        ),
        2: Params(
            names=("fudgeQQ", "q1", "q2", "sigma", "epsilon"),
            parsers=(float, float, float, float, float),
            desc="Coulomb + LJ pair",
        ),
    },
)

PAIRS_NB = InteractionSpec(
    n_atoms=2,
    funcs={
        1: Params(
            names=("qi", "qj", "V", "W"),
            parsers=(float, float, float, float),
            desc="non-bonded pair",
        ),
    },
)

ANGLES = InteractionSpec(
    n_atoms=3,
    funcs={
        1: Params(
            names=("th0", "kth"),
            parsers=(float, float),
            desc="angle"
        ),
        2: Params(
            names=("th0", "kth"),
            parsers=(float, float),
            desc="G96 angle"
        ),
        3: Params(
            names=("r1e", "r2e", "krr"),
            parsers=(float, float, float),
            desc="cross bond-bond",
        ),
        4: Params(
            names=("r1e", "r2e", "r3e", "krth"),
            parsers=(float, float, float, float),
            desc="cross bond-angle",
        ),
        5: Params(
            names=("th0", "kth", "r13", "kub"),
            parsers=(float, float, float, float),
            desc=" Urey-Bradley",
        ),
        6: Params(
            names=("th0", "C0", "C1", "C2", "C3", "C4"),
            parsers=(float, float, float, float, float, float),
            desc="quartic angle",
        ),
        8: Params(
            names=("table", "k"),
            parsers=(int, float),
            desc="tabulated angle"
        ),
    },
)

DIHEDRALS = InteractionSpec(
    n_atoms=4,
    funcs={
        1: Params(
            names=("phi_s", "kphi", "mult"),
            parsers=(float, float, int),
            desc="proper dihedral",
        ),
        2: Params(
            names=("xi0", "kxi"),
            parsers=(float, float),
            desc="improper dihedral"
        ),
        3: Params(
            names=("C0", "C1", "C2", "C3", "C4", "C5"),
            parsers=(float, float, float, float, float, float),
            desc="RB dihedral",
        ),
        4: Params(
            names=("phi_s", "kphi", "mult"),
            parsers=(float, float, int),
            desc="periodic improper",
        ),
        5: Params(
            names=("C1", "C2", "C3", "C4", "C5"),
            parsers=(float, float, float, float, float),
            desc="Fourier dihedral",
        ),
        8: Params(
            names=("table", "k"),
            parsers=(int, float),
            desc="tabulated dihedral"
        ),
        9: Params(
            names=("phi_s", "kphi", "mult"),
            parsers=(float, float, int),
            desc="proper dihedral (multiple)",
        ),
        10: Params(
            names=("phi0", "kphi"),
            parsers=(float, float),
            desc="restricted dihedral"
        ),
        11: Params(
            names=("kphi", "a0", "a1", "a2", "a3", "a4"),
            parsers=(float, float, float, float, float, float),
            desc="combined bending-torsion",
        ),
    },
)


CONSTRAINTS = InteractionSpec(
    n_atoms=2,
    funcs={
        1: Params(
            names=("b0",),
            parsers=(float,),
            desc="constraint"
        ),
        2: Params(
            names=("b0",),
            parsers=(float,),
            desc="constraint"
        ),
    },
)


NONBOND_PARAMS = InteractionSpec(
    n_atoms=2,
    funcs={
        1: Params(
            names=("sigma", "epsilon"),
            parsers=(float, float),
            desc="Lennard-Jones parameters",
        ),
        2: Params(
            names=("a", "b", "c6"),
            parsers=(float, float, float),
            desc="Buckingham parameters",
        ),
    },
)


SETTLES = InteractionSpec(
    n_atoms=1,
    funcs={
        1: Params(
            names=("doh", "dhh"),
            parsers=(float, float),
            desc="SETTLES parameters"
        ),
    },
)

VIRTUAL_SITES1 = InteractionSpec(
    n_atoms=2,
    funcs={
        1: Params(
            names=("a",),
            parsers=(float,),
            desc="1-body virtual site"
        ),
    },
)

VIRTUAL_SITES2 = InteractionSpec(
    n_atoms=3,
    funcs={
        1: Params(
            names=("a",),
            parsers=(float,),
            desc="2-body virtual site"
        ),
        2: Params(
            names=("d",),
            parsers=(float,),
            desc="2-body virtual site (fd)"
        ),
    },
)

VIRTUAL_SITES3 = InteractionSpec(
    n_atoms=4,
    funcs={
        1: Params(
            names=("a", "b"),
            parsers=(float, float),
            desc="3-body virtual site"
        ),
        2: Params(
            names=("a", "d"),
            parsers=(float, float),
            desc="3-body virtual site (fd)"
        ),
        3: Params(
            names=("th", "d"),
            parsers=(float, float),
            desc="3-body virtual site (fad)"
        ),
        4: Params(
            names=("a", "b", "c"),
            parsers=(float, float, float),
            desc="3-body virtual site (out)",
        ),
    },
)

VIRTUAL_SITES4 = InteractionSpec(
    n_atoms=5,
    funcs={
        2: Params(
            names=("a", "b", "c"),
            parsers=(float, float, float),
            desc="4-body virtual site (fdn)",
        ),
    },
)

VIRTUAL_SITESN = InteractionSpec(
    n_atoms=1,
    funcs={
        1: Params(
            names=("from",),
            parsers=(str,),
            desc="N-body virtual site (COG)"
            ),
        2: Params(
            names=("from",),
            parsers=(str,),
            desc="N-body virtual site (COM)"
            ),
        3: Params(
            names=("from",),
            parsers=(str,),
            desc="N-body virtual site (COW)"
            ),
    },
)

POSITION_RESTRAINTS = InteractionSpec(
    n_atoms=1,
    funcs={
        1: Params(
            names=("kx", "ky", "kz"),
            parsers=(float, float, float),
            desc="position restraints",
        ),
        2: Params(
            names=("g", "r", "k"),
            parsers=(float, float, float),
            desc="position restraints variant",
        ),
    },
)


DISTANCE_RESTRAINTS = InteractionSpec(
    n_atoms=2,
    funcs={
        1: Params(
            names=("type", "label", "low", "up1", "up2", "weight"),
            parsers=(int, int, float, float, float, float),
            desc="distance restraints",
        ),
    },
)


DIHEDRAL_RESTRAINTS = InteractionSpec(
    n_atoms=4,
    funcs={
        1: Params(
            names=("phi0", "dphi", "kdihr"),
            parsers=(float, float, int),
            desc="dihedral restraints",
        ),
    },
)


ORIENTATION_RESTRAINTS = InteractionSpec(
    n_atoms=2,
    funcs={
        1: Params(
            names=("exp", "label", "alpha", "c", "obs", "weight"),
            parsers=(int, int, float, float, float, float),
            desc="orientation restraints",
        ),
    },
)


ANGLE_RESTRAINTS = InteractionSpec(
    n_atoms=4,
    funcs={
        1: Params(
            names=("theta0", "kc", "mult"),
            parsers=(float, float, int),
            desc="angle restraints",
        ),
    },
)


ANGLE_RESTRAINTS_Z = InteractionSpec(
    n_atoms=2,
    funcs={
        1: Params(
            names=("theta0", "kc", "mult"),
            parsers=(float, float, int),
            desc="angle restraints z",
        ),
    },
)
