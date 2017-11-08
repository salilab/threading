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
    print(range(element[0],element[0]+element[1]))
    pis = IMP.atom.Selection(root_hier, chain_id="A", 
                residue_indexes=range(element[0],element[0]+element[1]),
                atom_type=IMP.atom.AT_CA).get_selected_particles()

    # Get XYZs
    xyz = []
    for p in pis:
        xyz.append(IMP.core.XYZ(p).get_coordinates())

    print xyz

    IMP.threading.StructureElement.setup_particle(p, element[0], 1, element[1], 0)
    #se = se.setup_particle(p, element[0], 1, element[1], 0)

    se = IMP.threading.StructureElement(p)
    se.set_coordinates(xyz)

    print "transform+1", se.get_resindex_list()
    se.flip_polarity_key()
    print "transform-1", se.get_resindex_list()
    return p



# The "true" sequence
#  MET at 11 and 37
seq = "MLPLGRDALYMAAVVATYEPISDHERRPNAATSDFEDELVVAGEFYIKKE"
# A toy system consisting of the two helices and 50 residues.
m = IMP.Model()
######################################
pdbfile = "./data/toy_system.pdb"
root_hier = IMP.atom.read_pdb(pdbfile, m)
IMP.atom.show_molecular_hierarchy(root_hier)

elements=[(1,13,'H'),(19,15,'H')]

print root_hier.get_children(), len(root_hier.get_children())
se = []

for e in elements:
    se.append(setup_structural_element(root_hier, e))

# The "correct" alignment is [(6,19), (29,45)]

#######################
# Set up Sequence
#######################
ch = IMP.atom.Chain.setup_particle(IMP.Particle(m), "S")
root_hier.add_child(ch)
for i in range(len(seq)):
    res = IMP.atom.Residue.setup_particle(IMP.Particle(m),
                                                    IMP.pmi.tools.get_residue_type_from_one_letter_code(seq[i]),
                                                    i+1)
    ch.add_child(res)

#########################
# Set up Structure Elements
#########################
# Add structural hierarchy
struct_hier = IMP.atom.Chain.setup_particle(IMP.Particle(m), "X")
root_hier.add_child(struct_hier)

p = IMP.Particle(m)
f = IMP.atom.Fragment.setup_particle(p, range(6,13))



#se.set_transformation_key((5,1,10,0))








