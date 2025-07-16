import math
import re
# This is just a basic code designed to handle non-polyatomic molecules from Insulin (C256H381N65O77S6) to H2O
# Polyatomic molecules are written like this for example H2O(H2O) *this is just an example

# Disclaimer: the code itself is basic in nature and is meant to show coders how to implement my SASA formula into
# a model

# covalent_radii is for the bond_length_list_calculations
covalent_radii = {'H': 32, 'He': 31, 'Li': 128, 'Be': 96, 'B': 84, 'C': 77, 'N': 75, 'O': 73, 'F': 63, 'Ne': 38,
                  'Na': 166, 'Mg': 141, 'Al': 121, 'Si': 111, 'P': 110, 'S': 104, 'Cl': 99, 'Ar': 71, 'K': 203, 'Ca': 176,
                  'Sc': 162, 'Ti': 146, 'V': 134, 'Cr': 139, 'Mn': 139, 'Fe': 140, 'Co': 135, 'Ni': 124, 'Cu': 128, 'Zn': 139,
                  'Ga': 122, 'Ge': 122, 'As': 119, 'Se': 116, 'Br': 114, 'Kr': 88, 'Rb': 303, 'Sr': 249, 'Y': 211, 'Zr': 160,
                  'Nb': 146, 'Mo': 139, 'Tc': 138, 'Ru': 139, 'Rh': 137, 'Pd': 139, 'Ag': 145, 'Cd': 144, 'In': 144, 'Sn': 141,
                  'Sb': 139, 'I': 133, 'Te': 138, 'Xe': 108, 'Cs': 303, 'Ba': 253, 'La': 195, 'Ce': 198, 'Pr': 200, 'Nd': 200,
                  'Pm': 202, 'Sm': 203, 'Eu': 206, 'Gd': 207, 'Tb': 208, 'Dy': 209, 'Ho': 210, 'Er': 211, 'Tm': 212, 'Yb': 214,
                  'Lu': 215, 'Hf': 159, 'Ta': 146, 'W': 139, 'Re': 137, 'Os': 136, 'Ir': 135, 'Pt': 139, 'Au': 144, 'Hg': 150,
                  'Tl': 146, 'Pb': 154, 'Bi': 157, 'Po': 150, 'At': 150, 'Rn': 140, 'Fr': 330, 'Ra': 283, 'Ac': 221, 'Th': 206,
                  'Pa': 200, 'U': 196, 'Np': 190, 'Pu': 186, 'Am': 184, 'Cm': 183, 'Bk': 182, 'Cf': 181, 'Es': 180, 'Fm': 179,
                  'Md': 178, 'No': 177, 'Lr': 176, 'Rf': 175, 'Db': 174, 'Sg': 173, 'Bh': 172, 'Hs': 171, 'Mt': 170, 'Ds': 169,
                  'Rg': 168, 'Cn': 167, 'Nh': 166, 'Fl': 165, 'Mc': 164, 'Lv': 163, 'Ts': 162, 'Og': 161}

# atomic_radii is used for the SASA value calculations
atomic_radii = {'H': 1.10, 'He': 1.40, 'Li': 1.82, 'Be': 1.53, 'B': 1.92, 'C': 1.70, 'N': 1.55, 'O': 1.52, 'F': 0.64,
                'Ne': 1.54, 'Na': 1.86, 'Mg': 1.60, 'Al': 1.84, 'Si': 2.10, 'P': 1.80, 'S': 1.80, 'Cl': 1.75, 'Ar': 1.88,
                'K': 2.27, 'Ca': 1.97, 'Sc': 2.00, 'Ti': 2.10, 'V': 2.11, 'Cr': 2.11, 'Mn': 2.14, 'Fe': 2.16, 'Ni': 1.63,
                'Co': 1.52, 'Cu': 1.40, 'Zn': 1.39, 'Ga': 1.87, 'Ge': 2.11, 'As': 1.85, 'Se': 1.90, 'Br': 1.85, 'Kr': 2.02,
                'Rb': 2.40, 'Sr': 2.15, 'Y': 2.20, 'Zr': 2.30, 'Nb': 2.30, 'Mo': 2.30, 'Tc': 2.35, 'Ru': 2.35, 'Rh': 2.31,
                'Pd': 1.75, 'Ag': 1.60, 'Cd': 1.58, 'In': 1.93, 'Sn': 2.17, 'Sb': 2.06, 'I': 1.98, 'Xe': 2.16, 'Cs': 2.62,
                'Ba': 2.30, 'La': 2.20, 'Ce': 2.30, 'Pr': 2.30, 'Nd': 2.30, 'Pm': 2.35, 'Sm': 2.45, 'Eu': 2.50, 'Gd': 2.55,
                'Tb': 2.60, 'Dy': 2.70, 'Ho': 2.70, 'Er': 2.75, 'Tm': 2.80, 'Yb': 2.85, 'Lu': 2.90, 'Hf': 2.40, 'Ta': 2.40,
                'W': 2.40, 'Re': 2.50, 'Os': 2.50, 'Ir': 2.50, 'Pt': 1.75, 'Au': 1.44, 'Hg': 1.55, 'Tl': 1.96, 'Pb': 2.02,
                'Bi': 2.07, 'Po': 2.10, 'At': 2.10, 'Rn': 2.20, 'Fr': 2.60, 'Ra': 2.40, 'Ac': 2.40, 'Th': 2.60, 'Pa': 2.60,
                'U': 2.60, 'Np': 2.70, 'Pu': 2.70, 'Am': 2.80, 'Cm': 2.80, 'Bk': 2.90, 'Cf': 2.90, 'Es': 3.00, 'Fm': 3.00,
                'Md': 3.00, 'No': 3.10, 'Lr': 3.10, 'Rf': 2.80, 'Db': 2.80, 'Sg': 2.80, 'Bh': 2.80, 'Hs': 2.80, 'Mt': 2.80,
                'Ds': 2.80, 'Rg': 2.80, 'Cn': 2.80, 'Nh': 2.80, 'Fl': 2.80, 'Mc': 2.80, 'Lv': 2.80, 'Ts': 2.80, 'Og': 2.80}

print()
user_input = input("Please enter your molecule here: ")
probe_radius = float(input("Enter probe radius (e.g., 1.4 for water): "))
number_of_bonds = input("Please enter the Number of bonds (ionic or covalent) here: ")
# Just enter the number of electron interactions in total for example: 2 for H2O it doesn't matter if they are ionic or
# covalent
print('---------------------------------------------------------------------------------------')

# This will store the final parsed list of element–subscript dictionaries
# Example: [{'C': 6}, {'H': 12}, {'O': 6}]
elements = []

# Regular expression to find element symbols and optional numbers
# ([A-Z][a-z]?) matches element symbols (e.g., H, He, O, Cl, etc.)
# (\d*) matches the number after the element symbol (subscript); optional
pattern = r'([A-Z][a-z]?)(\d*)'

# Find all matches in the formula, e.g., C6H12O6 -> [('C', '6'), ('H', '12'), ('O', '6')]
matches = re.findall(pattern, user_input)

# Loop through each match and build the elements list
for symbol, count in matches:
    # If there's no number after the element (e.g., H in H2O), count defaults to 1
    count = int(count) if count else 1

    # Append a dictionary with the element and its subscript
    elements.append({symbol: count})

bond_length_list = []
element_counts = {}
for ele_dict in elements:
    for ele, count in ele_dict.items():
        element_counts[ele] = element_counts.get(ele, 0) + count
unique_elements = list(element_counts.keys())
for i in range(len(unique_elements)):
    for j in range(i, len(unique_elements)):
        e1 = unique_elements[i]
        e2 = unique_elements[j]
        if e1 == e2 and element_counts[e1] < 2:
            continue
        if e1 in covalent_radii and e2 in covalent_radii:
            bond_length = covalent_radii[e1] + covalent_radii[e2]
            bond_name = f"{e1}-{e2}"
            bond_length_list.append({bond_name: bond_length})
# All bond values are in "pm"

# 1. Uncorrected SASA (sum of 4πr² for each atom)
total_sasa_uncorrected = sum(
    count * (4 * math.pi * atomic_radii[atom] ** 2)
    for item_ele in elements
        for atom, count in item_ele.items())

total_exclusion = 0.0
average_bond_length = 0.0
average_reff = 0.0
for bond in bond_length_list:
    for key_bond, value_bond in bond.items():
        key_bond_split = key_bond.split("-")
        ele_1 = key_bond_split[0]
        ele_2 = key_bond_split[1]
        ra = atomic_radii[ele_1]
        rn = atomic_radii[ele_2]
        d = value_bond
        reff = (ra + rn) / 2 + probe_radius
        average_bond_length += value_bond
        average_reff += reff
average_bond_length = average_bond_length / len(bond_length_list)
converted_bond_length = average_bond_length / 100
average_reff = average_reff / len(bond_length_list)
R_exclusion = math.sqrt((2 * average_reff) ** 2 - converted_bond_length ** 2) / 2
total_exclusion += R_exclusion
total_exclusion *= float(number_of_bonds)

# 3. Final corrected SASA
total_sasa_corrected = total_sasa_uncorrected - total_exclusion

print(bond_length_list, 'Estimated_bond_lengths')
print(total_sasa_uncorrected, 'total_sasa_uncorrected')
print(total_sasa_corrected, 'total_sasa_corrected')
print('---------------------------------------------------------------------------------------')


