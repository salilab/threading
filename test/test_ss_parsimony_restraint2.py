from __future__ import print_function
import IMP
import numpy
import IMP.pmi.tools
import IMP.atom
import IMP.core
import IMP.threading
import os

def sort_ses(ses):
    # Given a list of structural elements, sort them by increasing first residue
    res = sorted([(s.get_first_residue_number(), s) for s in ses], key=lambda x: x[0])
    #print(res)
    return [x[1] for x in res]

def setup_structural_element(root_hier, element, max_translation=1):
    ## SETUP STRUCTURAL ELEMENTS
    # Get Particle Indexes for CA atoms 
    # from residue range element[0]:element[1]
    # From chain X
    #print(range(element[0],element[0]+element[1]))
    pis = IMP.atom.Selection(root_hier, chain_id="A", 
                residue_indexes=range(element[0],element[0]+element[1]),
                atom_type=IMP.atom.AT_CA).get_selected_particles()

    pi = IMP.Particle(root_hier.get_model())

    h = IMP.atom.Hierarchy.setup_particle(pi)
    # Get XYZs
    xyz = []
    i = 0
    for p in pis:
        np = IMP.Particle(root_hier.get_model())
        hp = IMP.atom.Hierarchy.setup_particle(np)
        xyz = IMP.core.XYZ.setup_particle(np)
        #m = IMP.atom.Mass.setup_particle(root_hier.get_model(), np)
        xyz.set_coordinates(IMP.core.XYZ(p).get_coordinates())
        h.add_child(hp)


#    IMP.threading.StructureElement.setup_particle(root_hier.get_model(), pi.get_index(), element[0]-4, 1, element[1], 0, 'H')
    IMP.threading.StructureElement.setup_particle(root_hier.get_model(), pi.get_index(), element[0], 1, element[1], 0, 'A')
    #se = se.setup_particle(p, element[0], 1, element[1], 0)

    # Set up this element as a helix
    if element[-1] == 'H':
        IMP.atom.SecondaryStructureResidue.setup_particle(pi, 1, 0, 0)
    elif element[-1] == 'S':
        IMP.atom.SecondaryStructureResidue.setup_particle(pi, 0, 1, 0)
    elif element[-1] == 'C':
        IMP.atom.SecondaryStructureResidue.setup_particle(pi, 0, 0, 1)
    else:
        IMP.atom.SecondaryStructureResidue.setup_particle(pi, 0.33, 0.33, 0.33)
        print("secondary structure type is not defined for the corresponding element")

    se = IMP.threading.StructureElement(pi)
    se.set_keys_are_optimized(True)

    return se


# The "true" sequence
#  MET at 11 and 37
seq = "SQPAKKTYTWNTKEEAKQAFKEALKEKRVPSNASWEQAMKMIINDPRYSALAKLSEKKQAFNAYKVQTEK"
# A toy system consisting of the two helices and 50 residues.

# Restraint Weights
length_slope = 0.1
xl_slope = 0.1
semet_slope = 0.1
psipred_slope = 0.1


m = IMP.Model()
######################################
DATADIR = os.path.abspath(os.path.join(os.path.dirname(__file__), 'input'))
pdbfile = os.path.join(DATADIR, "pdb2lv9_A.ent")
root_hier = IMP.atom.read_pdb(pdbfile, m)
#IMP.atom.show_molecular_hierarchy(root_hier)

# Define Structural Elements
# These are the residues in the PDB that correspond to structural elements.
# (start_res, length, SSID)
#-------------------------
#elements=[(3,8,'H'),(19,19,'H')]#, (56,8,'H')]
elements=[(15,11,'H'),(13,16,'C'), (48,12,'H')]

se = []

for e in elements:
    se.append(setup_structural_element(root_hier, e))

se = sort_ses(se)

#######################
# Set up Sequence
#######################
seq_chain = IMP.atom.Chain.setup_particle(IMP.Particle(m), "A")
seq_chain.set_name(seq_chain.get_id())
root_hier.add_child(seq_chain)
for i in range(len(seq)):
    res = IMP.atom.Residue.setup_particle(IMP.Particle(m),
                                                    IMP.pmi.alphabets.amino_acid.get_residue_type_from_one_letter_code(seq[i]),
                                                    i+1)
    IMP.atom.Mass.setup_particle(res.get_particle(), IMP.atom.get_mass(res.get_residue_type()))
    IMP.core.XYZR.setup_particle(res.get_particle())
    IMP.atom.SecondaryStructureResidue.setup_particle(res.get_particle(), 0.8, 0, 0.2)
    seq_chain.add_child(res)

#######################
# Set up Scoring Function
#######################

# SS Parsimony Restraint
ssp_weight = 1.0

# Restraint evaluated at the SE level? Pass SE and seq_chain.  
# Have restraint look up residue from the seq_chain and get SS

a = IMP.atom.Hierarchy(seq_chain.get_particle()).get_children()
se_id = [s.get_particle_index() for s in se]
print(len(se_id))
print("get first sequence residue's all probabilities", IMP.atom.SecondaryStructureResidue(m, a[0].get_particle_index()).get_all_probabilities())
print("get first SE residue's all probabilities", IMP.atom.SecondaryStructureResidue(m, se_id[0]).get_all_probabilities())

par_res = IMP.threading.SecondaryStructureParsimonyRestraint(m, se_id, seq_chain.get_particle_index(), 1.0)
a_value  = par_res.unprotected_evaluate(None)
print(a_value)



