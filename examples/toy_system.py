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
        print IMP.core.XYZ(p).get_coordinates() 

    #print xyz

    IMP.threading.StructureElement.setup_particle(root_hier.get_model(), pi.get_index(), element[0], 1, element[1], 0)
    #se = se.setup_particle(p, element[0], 1, element[1], 0)

    se = IMP.threading.StructureElement(pi)
    #se.set_coordinates(xyz)

    #print "LIKE", IMP.atom.Hierarchy(pi).get_children()

    #print "XXXxx", IMP.core.XYZ(root_hier.get_model(), h.get_children()[0].get_particle_index()), h.get_children()
    print "YYY", pi
    print se.get_coordinates(), "XXX"
    return se



# The "true" sequence
#  MET at 11 and 37
seq = "MLPLGRDALYMAAVVATYEPISDHERRPNAATSDFEDELVVAGEFYIKKE"
# A toy system consisting of the two helices and 50 residues.
m = IMP.Model()
######################################
pdbfile = "./data/toy_system.pdb"
root_hier = IMP.atom.read_pdb(pdbfile, m)
#IMP.atom.show_molecular_hierarchy(root_hier)

elements=[(5,10,'H'),(19,15,'H')]

#print root_hier.get_children(), len(root_hier.get_children())
se = []

for e in elements:
    se.append(setup_structural_element(root_hier, e))

# The "correct" alignment is [(6,19), (29,45)]

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

#########################
# Set up Structure Elements
#########################
# Add structural hierarchy
struct_hier = IMP.atom.Chain.setup_particle(IMP.Particle(m), "X")
IMP.core.XYZR.setup_particle(struct_hier)
IMP.atom.Mass.setup_particle(struct_hier, 1)
root_hier.add_child(struct_hier)

p = IMP.Particle(m)
f = IMP.atom.Fragment.setup_particle(p, range(6,13))


#########################
# Move structure to Chain S
#########################

for s in se:
    #ssee = IMP.threading.StructureElement(s)
    #resis = IMP.atom.Selection(seq_chain, residue_indexes=s.get_resindex_list()).get_selected_particles()
    coordinates = s.get_coordinates()
    sem = IMP.threading.StructureElementMover(m, s.get_particle_index(), seq_chain.get_particle())
    #for i in range(len(resis)):
    #    r = resis[i]
    #    IMP.core.XYZ(r).set_coordinates(coordinates[i])

    print "DHJHFJK", s.get_all_key_values(), coordinates
    #sem.zero_coordinates()

    sem.transform_coordinates()

    for i in seq_chain.get_children():
        print i, IMP.core.XYZ(i.get_particle())
    print s.get_all_key_values()
    sem.propose()
    print s.get_all_key_values()

    for i in seq_chain.get_children():
        print i, IMP.core.XYZ(i.get_particle())

    print IMP.atom.Selection(seq_chain, residue_indexes=s.get_resindex_list()).get_selected_particles()





#se.set_transformation_key((5,1,10,0))








