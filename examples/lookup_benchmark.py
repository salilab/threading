import IMP
import numpy
import IMP.pmi.tools
import IMP.atom
import IMP.core
import cProfile


seq = "MLPLGRDALYMAAVVATYEPISDHERRPNAATSDFEDELVVAGEFYIKKE"
for i in range(10):
    seq += seq
# A toy system consisting of the two helices and 50 residues.
m = IMP.Model()
root_hier = IMP.atom.Hierarchy.setup_particle(IMP.Particle(m))

elements=[(1,13,'H'),(19,15,'H')]

######################################
pdbfile = "./data/toy_system.pdb"


hier = IMP.atom.read_pdb(pdbfile, root_hier.get_model())
chain = IMP.atom.Chain.setup_particle(m, hier.get_particle_index(), "X")
root_hier.add_child(hier)

ch = IMP.atom.Chain.setup_particle(IMP.Particle(m), "S")
root_hier.add_child(ch)
for i in range(len(seq)):
    res = IMP.atom.Residue.setup_particle(IMP.Particle(m),
                                                    IMP.pmi.tools.get_residue_type_from_one_letter_code(seq[i]),
                                                    i+1)
    xyzr = IMP.core.XYZR.setup_particle(res.get_particle())
    mass = IMP.atom.Mass.setup_particle(res.get_particle(), 1)
    xyzr.set_coordinates((-10,-10,-10))
    xyzr.set_radius(0.1)
    ch.add_child(res)
    #print(res, res.get_index())



for i in range(1000):
    res = IMP.atom.Selection(root_hier, chain_id="X", residue_index=12, atom_type=IMP.atom.AT_CA).get_selected_particles()
    #print(res)
    xyz = IMP.core.XYZ(res[0]).get_coordinates()

    sres = IMP.atom.Selection(root_hier, chain_id="S", residue_index=12).get_selected_particles()
    #print(sres)
    IMP.core.XYZ(sres[0]).set_coordinates(xyz)

