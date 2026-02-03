# This module defines parameter specifications
# for different non-interaction sections
# of a molecular force field file.

from .core import Params


DEFAULTS = Params(
    names=("nbfunc", "comb_rule", "gen_pairs", "fudgeLJ", "fudgeQQ"),
    parsers=(int, int, str, float, float),
)

ATOMTYPES = Params(
    names=("name", "atnum", "mass", "charge", "ptype", "sigma", "epsilon"),
    parsers=(str, int, float, float, str, float, float),
)

MOLECULETYPE = Params(
    names=("name", "nrexcl"),
    parsers=(str, int),
)

SYSTEM = Params(
    names=("name",),
    parsers=(str,),
)

ATOMS = Params(
    names=("nr", "type", "resnr", "residue", "name", "cgnr", "charge", "mass"),
    parsers=(int, str, int, str, str, int, float, float),
)
