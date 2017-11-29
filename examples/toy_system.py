import IMP
import numpy
import IMP.pmi.tools
import IMP.atom
import IMP.core
import IMP.threading



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
        #print IMP.core.XYZ(p).get_coordinates() 

    #print xyz

    IMP.threading.StructureElement.setup_particle(root_hier.get_model(), pi.get_index(), element[0]+10, 1, element[1], 0)
    #se = se.setup_particle(p, element[0], 1, element[1], 0)

    se = IMP.threading.StructureElement(pi)
    se.set_keys_are_optimized(True)

    print se.get_max_res()

    #print "XXXxx", IMP.core.XYZ(root_hier.get_model(), h.get_children()[0].get_particle_index()), h.get_children()
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
    print "SSR", r
    return r

# The "true" sequence
#  MET at 11 and 37
seq = "ALPLGRMALYGAAVVATYEPISDHERMPNAATSDFEDELVVAGEFYIKKE"
# A toy system consisting of the two helices and 50 residues.

# Restraint Weights
length_slope = 2
xl_slope = 1
semet_slope = 1
psipred_slope = 1


m = IMP.Model()
######################################
pdbfile = "./data/toy_system.pdb"
root_hier = IMP.atom.read_pdb(pdbfile, m)
#IMP.atom.show_molecular_hierarchy(root_hier)

# Define Structural Elements
#-------------------------
#elements=[(3,8,'H'),(19,19,'H')]#, (56,8,'H')]
elements=[(3,8,'H'),(19,19,'H')]#, (56,8,'H')]




se = []

for e in elements:
    se.append(setup_structural_element(root_hier, e))


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
    seq_chain.add_child(res)

#######################
# Set up Scoring Function
#######################

# XL Restraint
# ---------
# (res1, res2, length)
xl_constant = 10
xls = [(4, 32, 8), (10, 28, 10), (8, 29, 9)]
rests = []
for xl in xls:
    p1 = IMP.atom.Selection(root_hier, chain_id='S', residue_index=xl[0]).get_selected_particle_indexes()[0]
    p2 = IMP.atom.Selection(root_hier, chain_id='S', residue_index=xl[1]).get_selected_particle_indexes()[0]

    r = setup_conditional_pair_restraint(p1, p2, xl[2], xl_constant)
    print "XL Setup between residues: ", xl[0], xl[1], " at a length of ", xl[2]
    rests.append(r)

# Completeness Restraint
# ---------
for s in se:
    r = setup_length_restraint(s)
    rests.append(r)


# SeMet Position Restraint
# ---------
Se_pos = [(25.018, -40.381, 32.070), (43.517, -26.296, 33.279)]
Se_dist = 3.4  # real value is ~3.28 from xtal structure
met_particles = IMP.atom.Selection(root_hier, chain_id='S', residue_type=IMP.atom.MET).get_selected_particles()
print met_particles

for s in Se_pos:
    p = IMP.Particle(m)
    xyzr = IMP.core.XYZR.setup_particle(p)
    xyzr.set_coordinates(s)
    xyzr.set_radius(1)
    rs = []
    for met in met_particles:
        r = setup_pair_restraint(met, p, Se_dist)
        rs.append(r)
    mr = IMP.core.MinimumRestraint(1, rs) 
    rests.append(mr)

# MODELLER CA-CA Restraint
# ---------

# End to End distances
# ---------
# Structural domain : sequence of SEs.  Use the sequence to 
# add the E2E restraint.  model distance = last coord of SE1 - first coord of SE2.
# evaluated distance = (first resid of SE1 - last resid of SE2) * res_dist
# Scoring function is persistance length of the random coil?


# PSIPRED Restraint
# ---------



sf = IMP.core.RestraintsScoringFunction(rests)

print "EVALUATE - NO COORDS", sf.evaluate(False)

#########################
# Set up Structure Elements Samplers
#########################
# Add structural hierarchy
struct_hier = IMP.atom.Chain.setup_particle(IMP.Particle(m), "X")
IMP.core.XYZR.setup_particle(struct_hier)
IMP.atom.Mass.setup_particle(struct_hier, 1)
root_hier.add_child(struct_hier)

p = IMP.Particle(m)
f = IMP.atom.Fragment.setup_particle(p, range(6,13))

# The "correct" alignment results from Key values:
#  [(3,8,1,0), (19,19,1,0)]


#########################
# Move structure to Chain S
#########################

mc = IMP.core.MonteCarlo(m)
mc.set_scoring_function(sf)

for s in se:
    #ssee = IMP.threading.StructureElement(s)
    #resis = IMP.atom.Selection(seq_chain, residue_indexes=s.get_resindex_list()).get_selected_particles()
    coordinates = s.get_coordinates()
    sem = IMP.threading.StructureElementMover(m, s.get_particle_index(), seq_chain.get_particle())
    #for i in range(len(resis)):
    #    r = resis[i]
    #    IMP.core.XYZ(r).set_coordinates(coordinates[i])

    #sem.zero_coordinates()

    sem.transform_coordinates()

    #for i in seq_chain.get_children():
    #    print i, IMP.core.XYZ(i.get_particle())
    #sem.propose()
    #for i in seq_chain.get_children():
    #    print i, IMP.core.XYZ(i.get_particle())
    #sem.reject()
    #print s.get_all_key_values()
    #for i in seq_chain.get_children():
    #    print i, IMP.core.XYZ(i.get_particle())
    #print IMP.atom.Selection(seq_chain, residue_indexes=s.get_resindex_list()).get_selected_particles()

    #print sf.evaluate(False)
    mc.add_mover(sem)

for r in rests:
    print r, r.evaluate(False)

for s in se:
    print "KEY VALUES BEFORE:", s.get_all_key_values()

print "INIT EVALUATE", sf.evaluate(False)
print mc.get_kt()
mc.set_kt(100)

nframes = 100
for f in range(nframes):
    mc.optimize(10)
    print se[0].get_all_key_values(), se[1].get_all_key_values(), [(r, r.evaluate(False)) for r in rests]

for s in se:
    print "KEY VALUES AFTER:", s.get_all_key_values()
print "FINAL EVAL", sf.evaluate(False), mc.get_last_accepted_energy()
s

print mc.get_number_of_accepted_steps(), mc.get_number_of_proposed_steps()



for r in rests:
    print r, r.evaluate(False),

IMP.atom.destroy(root_hier)
print "destroyed"




