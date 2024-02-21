''' This file contains the allowed elements along with some of their properties.
If one wishes to add more (allowed) elements, simply change the ALPHABET list. Be sure to
update the INT_MASSES, EXACT_MASSES and DBES for said element as well!'''

#these are all for the most abundant isotope!
INT_MASSES = {"H": 1,
                  "C": 12,
                  "N": 14,
                  "O": 16,
                  "F": 19,
                  "Si": 28, 
                  "P": 31,
                  "S": 32,
                  "Cl": 35,
                  "Br": 79,
                  "I": 127}


EXACT_MASSES =   {"H": 1.00782503207,
                  "C": 12,
                  "N": 14.0030740048,
                  "O": 15.99491461956,
                  "F": 18.99840322,
                  "Si": 27.9769265,
                  "P": 30.97376163,
                  "S": 31.97207100,
                  "Cl": 34.96885268,
                  "Br": 78.9183361,
                  "I": 126.904477}


#This is the default list of elements, it can be changed if other lists/dictionaries are modified accordingly
ALPHABET = ["C",
"H",
"N",
"O",
"F",
"Si",
"P",
"S", 
"Cl",
"Br",
"I"
]

#Contribution of each atom to the RDBE (ring plus double bond equivalent) value
DBES ={"C": 1,
    "O": 0,
    "S": 0,
    "Si": 1,
    "H": -0.5,
    "N": 0.5,
    "Cl": -0.5,
    "Br": -0.5,
    "F": -0.5,
    "P": 0.5,
    "I": -0.5}



#unicode for subscripts!
SUBSCRIPTS = {"1": u"\u2081",
              "2": u"\u2082",
              "3": u"\u2083",
              "4": u"\u2084",
              "5": u"\u2085",
              "6": u"\u2086",
              "7": u"\u2087",
              "8": u"\u2088",
              "9": u"\u2089",
              "0": u"\u2080"
}