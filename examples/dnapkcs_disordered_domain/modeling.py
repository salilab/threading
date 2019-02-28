import IMP
import numpy
import IMP.pmi.tools
import IMP.atom
import IMP.core
import IMP.threading
import IMP.pmi.samplers


# NMR structure of the FF domain L24A mutant's folding transition state

def randomize_ses_start_residues(ses, seq_length):
    # given a number of SEs, randomize the start residues.
    import random

    avail_seq = range(1,seq_length+2)

    random.shuffle(ses)

    subseqs = [avail_seq]

    for s in ses:
        s_len = s.get_length()
        available_spots = []
        for ss in subseqs:
            if len(ss) > s_len:
                available_spots += ss[:-(s_len+2)]

        #print available_spots

        # pick an available spot
        spot = random.choice(available_spots)

        s.set_start_res_key(float(spot))

        for i in range(len(subseqs)):
            if spot in subseqs[i]:
                ss1 = range(subseqs[i][0],spot)
                ss2 = range(spot+s_len,subseqs[i][-1])

                prev_si = subseqs[i]
                del subseqs[i]
                if len(ss1) > 0: 
                    subseqs.append(ss1)
                if len(ss2) > 0:
                    subseqs.append(ss2)
                break


        #print spot, subseqs        



def read_psipred_ss2(ss2_file, ps):

    f = open(ss2_file, 'r')

    ss_probs = []
    for line in f.readlines():
        fields = line.strip().split()
        try: 
            int(fields[0])
            ss_probs.append((fields[4], fields[5], fields[3]))
        except:
            continue

    if len(ss_probs) != len(ps):
        raise Exception("SS file does not have same # of residues as protein")

    for i in range(len(ss_probs)):
        #print ps[i], ss_probs[i][0], ss_probs[i][1], ss_probs[i][2], ss_probs[i][0]
        IMP.atom.SecondaryStructureResidue.setup_particle(ps[i], float(ss_probs[i][0]), float(ss_probs[i][1]), float(ss_probs[i][2]))

    return ps 

def sort_ses(ses):
    # Given a list of structural elements, sort them by increasing first residue
    res = sorted([(s.get_first_residue_number(), s) for s in ses], key=lambda x: x[0])
    #print res
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
        #print IMP.core.XYZ(p).get_coordinates() 

    #print xyz

    IMP.threading.StructureElement.setup_particle(root_hier.get_model(), pi.get_index(), element[0], 1, element[1], 0)
    #se = se.setup_particle(p, element[0], 1, element[1], 0)

    # Set up this element as a helix
    IMP.atom.SecondaryStructureResidue.setup_particle(pi, 1, 0, 0)

    se = IMP.threading.StructureElement(pi)
    se.set_keys_are_optimized(True)
    se.set_max_res(200)

    #print "XXXxx", IMP.core.XYZ(root_hier.get_model(), h.get_children()[0].get_particle_index()), h.get_children()
    return se

def setup_conditional_pair_restraint(p1, p2, length, xl_slope, constant):
    #dps = IMP.core.DistancePairScore(IMP.core.HarmonicUpperBound(length, xl_slope))
    r = IMP.threading.ConditionalPairRestraint(m, IMP.core.HarmonicUpperBound(length, xl_slope), 
        p1, p2, constant)
    return r

def setup_pair_restraint(p1, p2, length):
    #dps = IMP.core.DistancePairScore(IMP.core.HarmonicUpperBound(length, xl_slope))
    r = IMP.core.DistanceRestraint(m, IMP.core.HarmonicUpperBound(length, xl_slope), 
        p1, p2)
    return r

def setup_length_restraint(s, length_slope):
    # Setup a restraint that biases the structural element towards
    # the length of the number of coordinates.

    # score = -log(#)
    uf = IMP.core.Linear(s.get_number_of_coordinates(), -1*length_slope)
    sf = IMP.core.AttributeSingletonScore(uf, IMP.FloatKey("length"))
    r = IMP.core.SingletonRestraint(m, sf, s.get_particle())
    print "SSR", r
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

#  Sequence of the disordered region
#  2601 to 2800
seq = "VLTPMFVETQASQGTLQTRTQEGSLSARWPVAGQIRATQQQHDFTLTQTADGRSSFDWLTGSSTDPLVDHTSPSSDSLLFAHKRSERLQRAPLKSVGPDFGKKRLGLPGDEVDNKVKGAAGRTDLLRLRRRFMRDQEKLSLMYARKGVAEQKREKEIKSELKMKQDAQVVLYRSYRHGDLPDIQIKHSSLITPLQAVAQR"

# Restraint Weights
length_slope = 0.1
xl_slope = 0.1
semet_slope = 0.1
psipred_slope = 0.1


m = IMP.Model()
######################################
pdbfile = "./5luq_unknown_A.pdb"
root_hier = IMP.atom.read_pdb(pdbfile, m)
#IMP.atom.show_molecular_hierarchy(root_hier)

# Define Structural Elements
# These are the residues in the PDB that correspond to structural elements.
# (start_res, length, SSID)
#-------------------------
#elements=[(3,10,'H'),(19,19,'H')], (56,8,'H')]
elements=[(2,14,'H'),(23,25,'H'), (55,10,'H')]
#elements=[(25,11,'H'),(3,16,'H'), (44,12,'H')]
se = []

for e in elements:
    se.append(setup_structural_element(root_hier, e))

randomize_ses_start_residues(se, len(seq)+1)
se = sort_ses(se)
#######################
# Set up Sequence
#######################
seq_chain = IMP.atom.Chain.setup_particle(IMP.Particle(m), "S")
root_hier.add_child(seq_chain)
res_particles = []
for i in range(len(seq)):
    res = IMP.atom.Residue.setup_particle(IMP.Particle(m),
                                                    IMP.pmi.tools.get_residue_type_from_one_letter_code(seq[i]),
                                                    i+1)
    res_particles.append(res.get_particle())
    IMP.atom.Mass.setup_particle(res.get_particle(), IMP.atom.get_mass(res.get_residue_type()))
    IMP.core.XYZR.setup_particle(res.get_particle())
    #IMP.atom.SecondaryStructureResidue.setup_particle(res.get_particle(), 0.8, 0, 0.2)
    seq_chain.add_child(res)
rests = []
#######################
# Set up Scoring Function
#######################
'''
# XL Restraint
# ---------
# (res1, res2, length)
xl_constant = 1 # Penalty for unmodeled crosslink
xl_slope = 0.1
xls = [(13, 58, 7), (20, 51, 8), (24, 38, 7)]

for xl in xls:
    p1 = IMP.atom.Selection(root_hier, chain_id='S', residue_index=xl[0]).get_selected_particle_indexes()[0]
    p2 = IMP.atom.Selection(root_hier, chain_id='S', residue_index=xl[1]).get_selected_particle_indexes()[0]

    r = setup_conditional_pair_restraint(p1, p2, xl[2], xl_slope, xl_constant)
    print "XL Setup between residues: ", xl[0], xl[1], " at a length of ", xl[2]
    rests.append(r)
'''
# Completeness Restraint
# ---------
for s in se:
    r = setup_length_restraint(s, length_slope=length_slope)
    rests.append(r)

'''
# SeMet Position Restraint
# ---------
#Se_pos = [(25.018, -40.381, 32.070), (43.517, -26.296, 33.279)]
#Se_pos = [(2.382,  -3.664, -12.726), (6.511, -12.492, -16.446)]
Se_pos = []
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
    mr = IMP.core.MinimumRestraint(1, rs, "SeRestraint%1%") 
    rests.append(mr)
'''
# MODELLER CA-CA Restraint
# ---------

# Loop Length distances
# ---------

# Structural domain : sequence of SEs.  Use the sequence to 
# add the E2E restraint.  model distance = last coord of SE1 - first coord of SE2.
# evaluated distance = (first resid of SE1 - last resid of SE2) * res_dist
# Scoring function is persistance length of the random coil?
se_pairs = []
for i in range(len(se)-1):
    print "RESIS: ", se[i].get_last_residue_number(), se[i+1].get_first_residue_number()
    se_pairs.append((se[i].get_particle_index(), se[i+1].get_particle_index()))


secrs = []
for s in se_pairs:
    r = add_SECR(s[0], s[1])
    secrs.append(r)
    print s, "NRES:", r.get_number_of_residues(), r.get_model_distance(), r.get_max_distance(), r.unprotected_evaluate(None)
    rests.append(r)


# PSIPRED Restraint
# ---------

# Need a SS restraint for each residue in the SEs. 
# If the residue is 
read_psipred_ss2("./data/psipred_dnapkcs_disordered.txt", res_particles)

#for i in res_particles:
#    print IMP.atom.SecondaryStructureResidue(i)

se_part_indexes = [s.get_particle_index() for s in se]

r = IMP.threading.SecondaryStructureParsimonyRestraint(m, se_part_indexes, seq_chain.get_particle_index(), 1.0)

print r.get_baseline_probabilities(), r.unprotected_evaluate(None)

rests.append(r)

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
mc.set_return_best(False)

sems = []

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
    #for i in seq_chain.get_children():
    #    print i, IMP.core.XYZ(i.get_particle())
    #sem.reject()
    #print s.get_all_key_values()
    #for i in seq_chain.get_children():
    #    print i, IMP.core.XYZ(i.get_particle())
    #print IMP.atom.Selection(seq_chain, residue_indexes=s.get_resindex_list()).get_selected_particles()
    sems.append(sem)

    mc.add_mover(sem)

def coolate_output(se, sf):
    outlist = [[s.get_all_key_values() for s in se], sf.evaluate(False)]
    return outlist

def run_one_sim(mc, n_equil_frames, n_prod_frames):
    mc.set_return_best(False)
    for f in range(n_equil_frames):
        mc.optimize(1)
    mc.set_return_best(True)
    for f in range(n_prod_frames):
        mc.optimize(10)


for r in rests:
    print r, r.evaluate(False)

for s in se:
    print "KEY VALUES BEFORE:", s.get_all_key_values()


#exit()

print "INIT EVALUATE", sf.evaluate(False), [(r.get_name(), r.evaluate(False)) for r in rests]
#init_score = sf.evaluate(False)

'''
for r in secrs:
    print r.get_name(), r.evaluate(False), r.get_model_distance(), r.get_number_of_residues()

for s in se:
    print s.get_coordinates()[0], s.get_coordinates()[-1], s.get_first_residue_number(), s.get_last_residue_number(), s.get_resindex_list()
    print seq[s.get_first_residue_number():s.get_last_residue_number()+1]

for mp in met_particles:
    print "MP", IMP.core.XYZ(mp)

for p in seq_chain.get_children():
    print IMP.atom.Residue(p.get_particle()), IMP.core.XYZ(p.get_particle())
'''


start_positions = []
for s in se:
    start_positions.append((s.get_start_res(), s.get_length(), 'H'))
#for sem in sems:
#    sem.transform_coordinates()

all_output=[]
print "NEW EVAL", sf.evaluate(False), [(r.get_name(), r.evaluate(False)) for r in rests]
#print mc.get_kt()
#mc.set_kt(100)


#rex = IMP.pmi.samplers.ReplicaExchange(m, 0.5, 2, [mc])
#mpivs = IMP.pmi.samplers.MPI_values(rex)


output_file = open("output2_json."+ str(rex.get_my_index()) + ".dat", "w")

print [start_positions, rex.get_my_index(), rex.get_my_temp()]
#output_file.write('{} {} {}\n'.format(start_positions, rex.get_my_index(), rex.get_my_temp()))



def enumerate_start_resis_3(seq_len, se_lengths, force=False):
    # Given three SE's, enumerate the potential start residues
    import itertools

    if se_lengths > 5 and force is True:
        print "Are you sure you want to enumerate this many SEs?"
        print "If so, set force=True"
        exit()

    all_se_perms = list(itertools.permutations(se_lengths, len(se_lengths)))
    start_resis = []

    for i in range(1, seq_len+1-se_lengths[0]-se_lengths[1]-se_lengths[2]):
        for j in range(se_lengths[0]+i, seq_len+1-se_lengths[1]-se_lengths[2]):
            for k in range(se_lengths[1]+j, seq_len+1-se_lengths[2]):
                start_resis.append((i,j,k))

    return start_resis


se_lengths=(14,25,10)

def enumerate_start_resis_3_eval(start_resis, ses, sems, sf):

    out_dict={}
    out_dict['models'] = []
    print "Let's perform ", len(start_resis), " calculations!"
    for i in range(len(start_resis)):
        if i%1000==0:
            print i,"th step"
        for n in range(len(start_resis[i])):
            sr = start_resis[i][n]
            ste = ses[n]
            sem = sems[n]

            ste.set_start_res_key(sr)
            sem.transform_coordinates()

        output = coolate_output(ses, sf)
        new_mod = {}
        new_mod['frame'] = i
        new_mod['model'] = output[0]
        new_mod['score'] = output[-1]
        new_mod['restraints'] = {}
        for r in rests:
            new_mod['restraints'][r.get_name()] = r.evaluate(False)
        out_dict['models'].append(new_mod)   
        
    return out_dict 

srs = enumerate_start_resis_3(len(seq), se_lengths)

out_dict = enumerate_start_resis_3_eval(srs, se, sems, sf)

json.dump(out_dict, "enumerate_disordered_1_2_3.dat")

exit()

out_dict={}
out_dict['models'] = []
all_output=[]
import json
rex.rem.set_was_used(True)
for i in range(10000):

    for so in rex.samplerobjects:
        so.optimize(10)
    #run_one_sim(mc, 5, 10)
    output = coolate_output(se, sf)
    new_mod = {}
    new_mod['frame'] = i
    new_mod['model'] = output[0]
    new_mod['score'] = output[-1]
    new_mod['restraints'] = {}
    for r in rests:
        new_mod['restraints'][r.get_name()] = r.evaluate(False)
    #all_output.append(output)
    rex.swap_temp(i, output[-1])

    out_dict['models'].append(new_mod)
    #mpivs.set_value("score",score)
    #output_file.write('{} {} {}\n'.format(i, output, [(r.get_name(), r.evaluate(False)) for r in rests]))

    #output_file.writelines([i, output, [(r.get_name(), r.evaluate(False)) for r in rests], mc.get_number_of_upward_steps(), mc.get_number_of_downward_steps()])
    if i%10==0:
        print i, output
#print all_output

json.dump(out_dict, output_file)
print start_positions

exit()


for i in range(1000):
    run_one_sim(mc, 5, 5)
    op = coolate_output(se, sf)
    all_output.append(op)

    print i, op, [(r.get_name(), r.evaluate(False)) for r in rests], mc.get_number_of_upward_steps(), mc.get_number_of_downward_steps()




print all_output

print start_positions


'''
for i in range(135):
    for s in se:
        s.set_start_res_key(s.get_start_res()+1)
    for sem in sems:
        sem.transform_coordinates()
    all_output.append(coolate_output(se, sf))
    print i, coolate_output(se, sf), [(r.get_name(), r.evaluate(False)) for r in rests], mc.get_number_of_upward_steps(), mc.get_number_of_downward_steps()


print all_output
'''

#print init_score


for s in se:
    print "KEY VALUES AFTER:", s.get_all_key_values()
print "FINAL EVAL", sf.evaluate(False), mc.get_best_accepted_energy(), i


print mc.get_number_of_accepted_steps(), mc.get_number_of_proposed_steps()



for r in rests:
    print r, r.evaluate(False),

IMP.atom.destroy(root_hier)
print "destroyed"




