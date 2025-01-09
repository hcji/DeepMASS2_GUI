'''
Author: Sean Li
Date: May 1, 2023
This file analyses parsed mass spectral data using the parent subformula graph method and 
generates the ranked lists of of whole-spectrum annotations.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as 
published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program. 
If not, see <https://www.gnu.org/licenses/>.
'''


import itertools

from core.annotating.subformula.formula import Formula, FormulaGenerator


def binary_search_closest(input_float, float_list): 
    left = 0
    right = len(float_list) - 1
    closest_index = 0
    closest_distance = abs(input_float - float_list[0])
    
    while left <= right:
        mid = (left + right) // 2
        distance = abs(input_float - float_list[mid])
        
        if distance < closest_distance:
            closest_index = mid
            closest_distance = distance
        
        if input_float < float_list[mid]:
            right = mid - 1
        
        elif input_float > float_list[mid]:
            left = mid + 1
        
        else:
            return float_list[mid]
    
    return float_list[closest_index]


class ExpMass:
    
    def __init__(self, mass, error): #error in ppm!
        self.mass = mass
        self.rounded_mass = int(round(mass,0))
        self.ppm_error = error
        self.abs_error = error*mass/10**6


    def __repr__(self):
        return "{} ± {}".format(self.mass, self.abs_error)
    
    
    def __add__(self, exp_mass2): #propagate uncertainty when adding exp masses
        new_mass = self.mass + exp_mass2.mass
        new_abs_error = (self.abs_error**2 + exp_mass2.abs_error**2)**0.5 #sum of random variables
        new_ppm_error = new_abs_error * 10**6 / new_mass
        return ExpMass(new_mass, new_ppm_error)


    def __sub__(self, exp_mass2): #propagate uncertainty when subtracting exp masses
        new_mass = self.mass - exp_mass2.mass
        new_abs_error = (self.abs_error**2 + exp_mass2.abs_error**2)**0.5 #sum of random variables
        new_ppm_error = new_abs_error * 10**6 / new_mass
        return ExpMass(new_mass, new_ppm_error)

    
    def possible_formulae(self, alphabet, element_restriction=False, parent=False, DBE_restriction = None):
        e = FormulaGenerator(alphabet, self.rounded_mass)
        return e.get_formula_list(self.mass, self.ppm_error, parent_formula=parent, custom_element_restriction=element_restriction, DBE_restriction=DBE_restriction)
    

    def possible_subformulae(formula, exp_mass_list):
        
        if len(exp_mass_list) == 0:
            return []

        generator = FormulaGenerator(formula.get_constituent_elements(), formula.nom_mass)
        tolerance = exp_mass_list[0].ppm_error
        return generator.compute_all_fragment_formulae([e.mass for e in exp_mass_list], tolerance, formula)



class MassSpectrum:

    def __init__(self, spectrum_dict, molecular_ion, error):

        '''Mass Spectrum object, which contains a dictionary of the form {mass: intensity}, along with an 
            associated ppm error. From a mass spectrum we can generate sets of formula annotations '''

        self.spectrum_dict = spectrum_dict
        self.masses = [mass for mass in spectrum_dict]
        self.intensities = [spectrum_dict[mass] for mass in spectrum_dict]
        self.error = error
        self.molecular_ion = molecular_ion
    

    def __repr__(self):
        entries = ["Mass: {}   Intensity: {}\n".format(pair[0], pair[1]) for pair in zip(self.masses, self.intensities)]
        return "".join(entries)
    
    
    def intensity_cutoff(self, factor): #removes peaks below a certain intensity
        new_dict = {mass: self.spectrum_dict[mass] for mass in self.masses if self.spectrum_dict[mass] > max(self.intensities)*factor or mass == self.molecular_ion}
        return MassSpectrum(new_dict, self.error)
    

    def get_p_filtered_mass_spectrum(self, parent_mass, intensity_cutoff = 0.1):

        '''Returns a mass spectrum with the largest mass set as the parent mass 
        and all (other) ions with an intensity below intensity_cutoff (in %!) is removed'''
        
        trimmed_spectrum = {mass: self.spectrum_dict.get(mass) for mass in self.spectrum_dict if mass <= parent_mass}
        max_intensity = max([self.spectrum_dict.get(mass) for mass in trimmed_spectrum])
        normed_spectrum = {mass: self.spectrum_dict.get(mass)*100/max_intensity for mass in trimmed_spectrum}
        filtered_spectrum = {mass: self.spectrum_dict.get(mass) for mass in normed_spectrum if self.spectrum_dict.get(mass)*100/max_intensity > intensity_cutoff or mass == parent_mass}
        return MassSpectrum(filtered_spectrum, self.error).conventional_norm()

    
    def conventional_norm(self): #normalise a mass spectrum according to convention (most intense peak = 100
        base_peak = max(self.intensities)
        new_dict = {mass: 100 * self.spectrum_dict[mass] / base_peak for mass in self.spectrum_dict}
        return MassSpectrum(new_dict, self.error)

    
    def generate_fragment_formulae(self, parent_formula, delta_frag): #generates a set of annotations for the fragments given an explicit parent formula
        exp_masses = [ExpMass(mass, delta_frag * self.error) for mass in self.masses if mass != self.molecular_ion]
        formula_annotations = [sublist[0] for sublist in ExpMass.possible_subformulae(parent_formula, exp_masses) if sublist != ()]
        formula_annotations.append(parent_formula)
        return formula_annotations
    
        
    def generate_fragment_formulae_NL(self, parent_formula, delta_frag): #generates a set of annotations for the fragments given a string molecular formula (neutral loss version)
        molecular_mass = ExpMass(self.molecular_ion, delta_frag * self.error)
        mass_differences = [molecular_mass - ExpMass(mass, delta_frag * self.error) for mass in self.masses if mass != self.molecular_ion]
        fragment_formulae = [mass_diff.possible_formulae(parent_formula.get_constituent_elements(), parent=parent_formula) for mass_diff in mass_differences]
        NL_annotation = [parent_formula - formula_list[0] for formula_list in fragment_formulae if formula_list]
        NL_annotation.append(parent_formula)
        return NL_annotation
    

    def vertex_scoring_function(self, annotation): #computes the vertex score for an annotation
        return len(annotation) 
    

    def make_edges_for_annotation(self, annotation): #returns the PSG explicitly as a list of edges
        return [(f1, f2) for f1, f2 in itertools.product(annotation, annotation) if f1 != f2 and f1.is_subformula(f2)]


    def edge_scoring_function(self, annotation): #computes the number of edges in the PSG
        return len(self.make_edges_for_annotation(annotation))
    

    def product_scoring_function(self, annotation): #computes the product score for the PSG
        vertex_score = self.vertex_scoring_function(annotation)
        
        if vertex_score == 1: #prevent division by zero 
            return 0

        return 2 * self.edge_scoring_function(annotation) / (len(self.spectrum_dict) * (vertex_score -1)) #the product is 2 |PSG|N_V / N_{peak}N_V(N_V -1), N_V cancels

    
    def get_spectral_annotations(self, alphabet, scoring_function, delta_frag = 1):

        '''Makes all spectral annotations which correspond to a
            parent-candidate subformula graph. Note that this function 
            assumes the largest mass is the molecular ion! '''

        candidate_formulae = ExpMass(self.molecular_ion, self.error * 3).possible_formulae(alphabet, DBE_restriction = lambda x: x.dbe() >= 0 and x.dbe().is_integer())
        fragment_peaks = {mass: self.spectrum_dict[mass] for mass in self.spectrum_dict if mass != self.molecular_ion}
        list_of_annotations = []
        for candidate in candidate_formulae:
            annotation = self.generate_fragment_formulae(candidate, delta_frag)
            
            if annotation: 
                list_of_annotations.append(annotation)
        
        return sorted(list_of_annotations, key=scoring_function, reverse=True) 


    def compute_most_likely_molecular_ion(self, alphabet, scoring_function, lower_bound, upper_bound, delta_frag=3, cutoff_val=0.1):

        '''Returns a list of spectral annotations ranked by a scoring function. This function will scan over
            all possible masses between lower_bound and upper_bound and generate annotations for the mass'''
        
        all_candidates = []
        candidate_masses = [mass for mass in self.masses if mass > lower_bound and mass < upper_bound]
        for mass in candidate_masses:
            new_ms = self.get_p_filtered_mass_spectrum(mass, intensity_cutoff=cutoff_val)
            annotations = new_ms.get_spectral_annotations(alphabet, scoring_function, delta_frag=delta_frag)
            annotations.sort(key=lambda a: a[-1].mass_deviation(mass)) #sorts by order of mass deviation of the parent ion, allows for tie breaking
            
            if annotations:
                all_candidates.append(annotations)
        
        all_candidates = list(itertools.chain.from_iterable(all_candidates)) #group all candidates together
        return sorted(all_candidates, key=scoring_function, reverse=True)
    

    def get_formula_dict(self, annotation):

        '''Constructs a dictionary of the form {formula: intensity} corresponding to an
            annotated mass spectrum given an input list of formulae'''
        
        formula_dict = {}
        sorted_masses = sorted([mass for mass in self.masses])
        for formula in annotation:
            target_mass = binary_search_closest(formula.exact_mass, sorted_masses)
            formula_dict[str(formula)] = self.spectrum_dict.get(target_mass)
        
        return formula_dict
            

    def compute_with_candidate_formula(self, all_candidates, scoring_function, delta_frag=3, cutoff_val=0.1):
        
        list_of_annotations = []
        all_candidates = [Formula.from_string(s) for s in all_candidates]
        for candidate in all_candidates:
            annotation = self.generate_fragment_formulae(candidate, delta_frag)
            print(annotation)
            
            if annotation: 
                list_of_annotations.append(annotation)
        
        return sorted(list_of_annotations, key=scoring_function, reverse=True) 

