from dataclasses import dataclass
from typing import Any, Dict, Mapping, Sequence, Tuple


@dataclass(frozen=True)
class Params:
    """Specification of interaction parameters."""
    names: Tuple[str, ...]
    parsers: Tuple[type, ...]
    desc: str = ""

    def __post_init__(self):
        if len(self.names) != len(self.parsers):
            raise ValueError(
                f"Params names and parsers must have the same length: "
                f"got {len(self.names)} names and {len(self.parsers)} parsers."
            )


@dataclass(frozen=True)
class InteractionSpec:
    """Specification of an interaction type."""
    n_atoms: int
    funcs: Mapping[int, Params]

    def parse(
        self,
        func: int,
        tokens: Sequence[str],
        ctx: str = ""
    ) -> dict[str, Any]:
        """Parse interaction parameters."""

        param_spec = self.funcs.get(func)
        if param_spec is None:
            allowed = ", ".join(map(str, self.funcs))
            raise NotImplementedError(
                f"Only function(s) {allowed} are supported in {ctx}: got {func}."
            )

        params: Dict[str, Any] = {}
        for name, ptype, token in zip(param_spec.names, param_spec.parsers, tokens):
            try:
                params[name] = ptype(token)
            except ValueError as e:
                params[name] = token
                # raise ValueError(
                #     f"Invalid parameter value for '{name}' "
                #     f"(expected type {ptype.__name__}) "
                #     f"in function {func} of {ctx}: got '{token}'."
                # ) from e

        return params
