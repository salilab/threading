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
        #print(IMP.core.XYZ(p).get_coordinates())

    #print(xyz)

    IMP.threading.StructureElement.setup_particle(root_hier.get_model(), pi.get_index(), element[0]-4, 1, element[1], 0)
    #se = se.setup_particle(p, element[0], 1, element[1], 0)
    
    # Set up this element as a helix
    IMP.atom.SecondaryStructureResidue.setup_particle(pi, 1, 0, 0)

    se = IMP.threading.StructureElement(pi)
    se.set_keys_are_optimized(True)

    print(se.get_max_res())

    #print("XXXxx", IMP.core.XYZ(root_hier.get_model(), h.get_children()[0].get_particle_index()), h.get_children())
    return se

def setup_conditional_pair_restraint(p1, p2, length, constant):
    #dps = IMP.core.DistancePairScore(IMP.core.HarmonicUpperBound(length, xl_slope))
    r = IMP.threading.ConditionalPairRestraint(m, IMP.core.HarmonicUpperBound(length, xl_slope), 
        p1, p2, constant)
    return r

def setup_pair_restraint(p1, p2, length):
    #dps = IMP.core.DistancePairScore(IMP.core.HarmonicUpperBound(length, xl_slope))
    r = IMP.core.DistanceRestraint(m, IMP.core.HarmonicUpperBound(length, xl_slope), 
        p1, p2)
    return r

def setup_length_restraint(s):
    # Setup a restraint that biases the structural element towards
    # the length of the number of coordinates.

    # score = -log(#)
    uf = IMP.core.Linear(s.get_number_of_coordinates(), -1*length_slope)
    sf = IMP.core.AttributeSingletonScore(uf, IMP.FloatKey("length"))
    r = IMP.core.SingletonRestraint(m, sf, s.get_particle())
    print("SSR", r)
    return r

def add_SECR(p1, p2, slope=1, dpr=3.4):
    r = IMP.threading.StructureElementConnectivityRestraint(m, IMP.core.HarmonicUpperBound(0, slope), p1, p2, dpr, "")
    return r

def add_all_SECR(se_list, slope=1, dpr=3.4):
    SECR_restraints = []
    for i in range(len(se_list-1)):
        p1 = se_list[i].get_particle_index()
        p2 = se_list[i+1].get_particle_index()
        rs.append(IMP.threading.StructureElementConnectivityRestraint(m, IMP.core.HarmonicUpperBound(0, slope), p1, p2, dpr, ""))

    return SECR_restraints

def modify_all_SECR(se_list, rst_list):

    for i in range(len(rst_list)):
        p1 = se_list[i].get_particle_index()
        p2 = se_list[i+1].get_particle_index()
        rst_list.assign_particles(p1, p2)

    return SECR_restraints

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
pdbfile = os.path.join(DATADIR, "pdb2lv9.ent")
root_hier = IMP.atom.read_pdb(pdbfile, m)
#IMP.atom.show_molecular_hierarchy(root_hier)

# Define Structural Elements
# These are the residues in the PDB that correspond to structural elements.
# (start_res, length, SSID)
#-------------------------
#elements=[(3,8,'H'),(19,19,'H')]#, (56,8,'H')]
elements=[(35,11,'H'),(13,16,'H'), (48,12,'H')]

se = []

for e in elements:
    se.append(setup_structural_element(root_hier, e))

se = sort_ses(se)
#######################
# Set up Sequence
#######################
seq_chain = IMP.atom.Chain.setup_particle(IMP.Particle(m), "S")
root_hier.add_child(seq_chain)
for i in range(len(seq)):
    res = IMP.atom.Residue.setup_particle(IMP.Particle(m),
                                                    IMP.pmi.tools.get_residue_type_from_one_letter_code(seq[i]),
                                                    i+1)

    IMP.atom.Mass.setup_particle(res.get_particle(), IMP.atom.get_mass(res.get_residue_type()))
    IMP.core.XYZR.setup_particle(res.get_particle())
    IMP.atom.SecondaryStructureResidue.setup_particle(res.get_particle(), 0.8, 0, 0.2)
    seq_chain.add_child(res)

#######################
# Set up Scoring Function
#######################

# XL Restraint
# ---------
# (res1, res2, length)
xl_constant = 1 # slope of line
xls = [(14, 59, 7), (24, 49, 8), (25, 39, 7)]
rests = []
for xl in xls:
    p1 = IMP.atom.Selection(root_hier, chain_id='S', residue_index=xl[0]).get_selected_particle_indexes()[0]
    p2 = IMP.atom.Selection(root_hier, chain_id='S', residue_index=xl[1]).get_selected_particle_indexes()[0]

    r = setup_conditional_pair_restraint(p1, p2, xl[2], xl_constant)
    print("XL Setup between residues: ", xl[0], xl[1], " at a length of ", xl[2])
    rests.append(r)

# Completeness Restraint
# ---------
for s in se:
    r = setup_length_restraint(s)
    rests.append(r)


# SeMet Position Restraint
# ---------
#Se_pos = [(25.018, -40.381, 32.070), (43.517, -26.296, 33.279)]
Se_pos = [(2.382,  -3.664, -12.726), (6.511, -12.492, -16.446)]

Se_dist = 3.4  # real value is ~3.28 from xtal structure
met_particles = IMP.atom.Selection(root_hier, chain_id='S', residue_type=IMP.atom.MET).get_selected_particles()

for s in Se_pos:
    p = IMP.Particle(m)
    xyzr = IMP.core.XYZR.setup_particle(p)
    xyzr.set_coordinates(s)
    xyzr.set_radius(1)
    rs = []
    for met in met_particles:
        r = setup_pair_restraint(met, p, Se_dist)
        rs.append(r)
    mr = IMP.core.MinimumRestraint(1, rs, name="SeMetRestraint %1%") 
    rests.append(mr)

# MODELLER CA-CA Restraint
# ---------

# Loop Length distances
# ---------

# Structural domain : sequence of SEs.  Use the sequence to 
# add the E2E restraint.  model distance = last coord of SE1 - first coord of SE2.
# evaluated distance = (first resid of SE1 - last resid of SE2) * res_dist
# Scoring function is persistance length of the random coil?
se_rests = []
se_pairs = []
for i in range(len(se)-1):
    print("RESIS: ", se[i].get_last_residue_number(), se[i+1].get_first_residue_number())
    se_pairs.append((se[i].get_particle_index(), se[i+1].get_particle_index()))

for s in se_pairs:
    r = add_SECR(s[0], s[1])
    print(s, "NRES:", r.get_number_of_residues(), r.get_model_distance(), r.get_max_distance(), r.unprotected_evaluate(None))
    rests.append(r)
    se_rests.append(r)

print([r.unprotected_evaluate(None) for r in se_rests])


for i in range(22):
    se[0].set_start_res_key(se[0].get_start_res()+1)
    r = se_rests[0]
    print("SECR:", se[0].get_start_res(), se[1].get_start_res(), r.get_number_of_residues(), r.get_model_distance(), r.get_max_distance(), r.unprotected_evaluate(None))



exit()



