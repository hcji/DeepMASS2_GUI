"""
Author: Sean Li
Date: May 1, 2023
This file contains the Formula class, as well as all classes/methods that are used to generate possible molecular formulae
from an input alphabet, mass and uncertainty. 

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as 
published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program. 
If not, see <https://www.gnu.org/licenses/>.
"""

import functools
import itertools
import re

import numpy as np

from core.annotating.subformula import constants


class Formula:
    def __init__(self, compomer):
        self.alphabet = constants.ALPHABET
        self.compomer = np.array(compomer)
        self.formula_dict = {
            element: coefficient
            for element, coefficient in zip(self.alphabet, compomer)
        }
        self.nom_mass = sum(
            map(lambda x, y: constants.INT_MASSES[x] * y, self.alphabet, self.compomer)
        )
        self.exact_mass = sum(
            map(
                lambda x, y: constants.EXACT_MASSES[x] * y, self.alphabet, self.compomer
            )
        )

    def __repr__(self):
        symbols = []
        for element in self.alphabet:
            if self.formula_dict.get(element) == 0:
                symbols.append("")
            elif self.formula_dict.get(element) == 1:
                symbols.append(element)
            else:
                symbols.append(element + str(self.formula_dict.get(element)))
        return "".join(symbols)

    def __add__(self, f2):
        new_compomer = self.compomer + f2.compomer
        return Formula(new_compomer)

    def __sub__(self, f2):
        new_compomer = self.compomer - f2.compomer
        return Formula(new_compomer)

    def __eq__(self, f2):
        return all(
            [
                self.formula_dict.get(elem) == f2.formula_dict.get(elem)
                for elem in self.alphabet
            ]
        )

    def latex_formula(self):
        symbols = []
        for element in constants.ALPHABET:
            if self.formula_dict.get(element) == 0:
                symbols.append("")
            elif self.formula_dict.get(element) == 1:
                symbols.append(r"\mathrm{" + element + "}")
            else:
                symbols.append(
                    r"\mathrm{"
                    + element
                    + "}"
                    + "_{"
                    + str(self.formula_dict.get(element))
                    + "}"
                )
        string = "".join(symbols)
        return r"${}$".format(string)

    def from_string(string):
        coefficients = [0] * len(constants.ALPHABET)
        matches = re.findall(r"([A-Z][a-z]*)(\d*)", string)
        for match in matches:
            symbol = match[0]
            coefficient = match[1] if match[1] else "1"
            coefficients[constants.ALPHABET.index(symbol)] += int(coefficient)
        return Formula(coefficients)

    def is_subformula(self, superformula):
        return min(superformula.compomer - self.compomer) >= 0

    def mass_deviation(
        self, exp_mass
    ):  # return mass deviation from experimental mass in ppm
        return 10**6 * abs(self.exact_mass - exp_mass) / exp_mass

    def dbe(self):  # ring + double bond formula: DBE = n(C) + 1 - n(H)/2 + n(N)/2
        res = 1
        for elem in self.alphabet:
            res += constants.DBES[elem] * self.compomer[self.alphabet.index(elem)]
        return res

    def get_constituent_elements(
        self,
    ):  # get an alphabet list (e.g ["C", "H", "O", "Cl"]) consisting of all elements of a given formula
        return [
            element
            for element in self.formula_dict
            if self.formula_dict.get(element) != 0
        ]


class Tree:
    def __init__(self, data=None, daughters=None):
        self.daughters = daughters
        self.data = data

    def traverse(self, index):
        return self.daughters[index]

    def sum_tree(self):
        if self.data == None:
            return 0
        elif self.daughters == None:
            return self.data
        else:
            return self.data + self.traverse(0).sum_tree()

    def expand_tree_into_formulae(
        self, alphabet, experimental_mass, tolerance, DBE_restriction
    ):
        """Expands a formula (compomer tree) into formula objects satisfying
        experimental mass and tolerance restrictions"""

        formula_list = []
        atomic_masses = np.array([constants.INT_MASSES.get(elem) for elem in alphabet])
        scaling_factor = np.array(
            [
                constants.EXACT_MASSES.get(elem) / constants.INT_MASSES.get(elem)
                for elem in alphabet
            ]
        )[
            ::-1
        ]  # e.g 1.0078 / 1

        def inner(t, r, depth):
            if depth > 0:
                r[depth - 1] = t.data

            if depth == len(alphabet):
                calc_mass = r.dot(scaling_factor)

                if (
                    10**6 * abs(calc_mass - experimental_mass) / experimental_mass
                    < tolerance
                ):
                    compomer = r.copy()[::-1] // atomic_masses
                    formula_dict = {
                        element: coefficient
                        for element, coefficient in zip(alphabet, compomer)
                    }  # converts compomer into dict
                    expanded_compomer = [
                        formula_dict.get(element, 0) for element in constants.ALPHABET
                    ]  # expands out the compomer into the full alphabet
                    formula_list.append(Formula(expanded_compomer))

            else:
                for daughter in t.daughters:
                    inner(daughter, r, depth + 1)

        inner(self, np.zeros((len(atomic_masses),), dtype=int), 0)

        if DBE_restriction:
            formula_list = filter(DBE_restriction, formula_list)

        formula_list = sorted(
            list(formula_list), key=lambda f: f.mass_deviation(experimental_mass)
        )
        return tuple(formula_list)


class FormulaGenerator:
    def __init__(self, alphabet, nominal_mass):
        self.alphabet = alphabet
        self.nominal_mass = nominal_mass

    def get_default_restrictions(self):
        """Computes the default restriction table, in the form {C: (a, b), H: (c, d)...}
        which stipulate the maximum/minimum element count in the formulae to be generated
        """

        restrictions = {}
        for element in self.alphabet:
            if element == "C":
                restrictions[element] = (
                    self.nominal_mass // (4 * constants.INT_MASSES["C"]),
                    self.nominal_mass // constants.INT_MASSES["C"],
                )
            elif element == "H":
                max_C = self.nominal_mass // constants.INT_MASSES["C"]
                restrictions[element] = (0, 2 * max_C + 2)
            else:
                restrictions[element] = (
                    0,
                    self.nominal_mass // constants.INT_MASSES[element],
                )

        return restrictions

    def get_parental_restrictions(self, formula):
        """Computes the restriction table, in the form {C: (a, b), H: (c, d)...}
        which stipulate the maximum/minimum element count in a subformula given a formula
        """

        restrictions = {}
        for element in self.alphabet:
            restrictions[element] = (0, formula.formula_dict.get(element))

        return restrictions

    def get_generating_functions(self, element_restriction=False, formula=False):
        """From an element restriction table, compute the generating functions corresponding
        to each individual element in the formula"""

        generating_functions = []

        if element_restriction:
            restriction = element_restriction
        elif formula:
            restriction = self.get_parental_restrictions(formula)
        else:
            restriction = self.get_default_restrictions()

        for element in restriction:
            interval = restriction.get(element)
            atom_mass = constants.INT_MASSES[element]
            sequence = range(
                atom_mass * interval[0], atom_mass * interval[1] + 1, atom_mass
            )
            generating_functions.append(list(sequence))

        return generating_functions

    def compute_all_formula_trees(self, generating_functions):
        """Computes the list of formula trees"""

        def build_trees_from_treelist(generating_fn, list_of_tree_lists):
            first_tree = lambda p: p[1][0]
            return [
                Tree(prod[0], prod[1])
                for prod in itertools.product(generating_fn, list_of_tree_lists)
                if first_tree(prod).sum_tree() + prod[0] < self.nominal_mass + 1
            ]

        def bin_trees(tree_list):
            res = [[] for i in range(self.nominal_mass + 1)]
            for tree in tree_list:
                if tree.sum_tree() < self.nominal_mass + 1:
                    res[tree.sum_tree()].append(tree)
            return [tuple(trees) for trees in res if len(trees) != 0]

        base_case = [[Tree()]]  # singleton list of lists
        list_of_tree_lists = functools.reduce(
            lambda t, gf: bin_trees(build_trees_from_treelist(gf, t)),
            generating_functions,
            base_case,
        )
        return list_of_tree_lists

    def compute_formula_tree(
        self, generating_functions
    ):  # returns the last formula tree, corresponding to all formulae with some given nominal mass
        all_trees = self.compute_all_formula_trees(generating_functions)
        final_subtree = all_trees[-1]
        return Tree(0, final_subtree)

    def enumerate_compomers(self, generating_functions):
        """Computes all possible formulae of a given nominal mass
        using the generating function method."""

        formula_tree = self.compute_formula_tree(generating_functions)
        return formula_tree.expand_tree_into_formulae(
            self.alphabet, self.nominal_mass, 9999999, DBE_restriction=None
        )

    def get_formula_list(
        self,
        experimental_mass,
        tolerance,
        DBE_restriction=None,
        custom_element_restriction=False,
        parent_formula=False,
    ):
        """Returns a list of formulae within a certain tolerance (ppm) and
        some DBE restriction (e.g lambda x: x.dbe() >= 0 and x.dbe().is_integer())"""

        res = []
        generating_functions = self.get_generating_functions(
            element_restriction=custom_element_restriction, formula=parent_formula
        )
        tree = self.compute_formula_tree(generating_functions)
        return tree.expand_tree_into_formulae(
            self.alphabet, experimental_mass, tolerance, DBE_restriction
        )

    def compute_all_fragment_formulae(self, exp_mass_list, tolerance, parent_formula):
        """Given a list of experimental masses and uncertainties, compute all
        possible fragment formula annotations efficiently"""

        annotations = []
        generating_functions = self.get_generating_functions(formula=parent_formula)
        tree_list = [
            Tree(0, subtree)
            for subtree in self.compute_all_formula_trees(generating_functions)
        ]
        tree_dict = {
            tree.sum_tree(): tree for tree in tree_list
        }  # create list of trees

        for mass in exp_mass_list:
            frag_nominal_mass = int(round(mass))
            to_expand = tree_dict.get(frag_nominal_mass)

            if to_expand:
                fragment_candidates = to_expand.expand_tree_into_formulae(
                    self.alphabet,
                    mass,
                    tolerance,
                    DBE_restriction=lambda f: f.dbe() >= 0,
                )
                annotations.append(fragment_candidates)

        return annotations
