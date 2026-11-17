"""Defined mass building blocks with explicit attachment references."""
from knowledge_source import components, delta_table, molecule_table

from mass_ledger import parse_formula


# Each entry is (net formula, defined transformation). A group name alone makes no site claim.
DELTAS = delta_table()

MOLECULES = molecule_table()

SEQUENCES = components()['sequence_examples']


def group_formula(group):
    try:
        formula = DELTAS[group][0]
    except KeyError:
        raise ValueError(f"Unknown catalog group: {group!r}") from None
    return parse_formula(formula)
