from __future__ import print_function
import IMP
import IMP.test
import IMP.pmi.tools
import IMP.atom
import IMP.core
import IMP.threading
import os

class Tests(IMP.test.TestCase):

    def sort_ses(self, ses):
        # Given a list of structural elements, sort them by increasing first residue
        res = sorted([(s.get_first_residue_number(), s) for s in ses], key=lambda x: x[0])
        #print(res)
        return [x[1] for x in res]

    def setup_structural_element(self, root_hier, element, max_translation=1):
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
            coord = IMP.core.XYZ(p).get_coordinates()
            xyz.set_coordinates(coord)
            h.add_child(hp)


        IMP.threading.StructureElement.setup_particle(root_hier.get_model(), pi.get_index(), element[0], 1, element[1], 0, 'A')

        # Set up this element as a helix or strand or a coil

        probmap = {'H': (1, 0, 0), 'S': (0, 1, 0), 'C': (0, 0, 1)}
        prob = probmap.get(element[-1], (0.33, 0.33, 0.33))
        IMP.atom.SecondaryStructureResidue.setup_particle(pi, *prob)

        se = IMP.threading.StructureElement(pi)

        se.set_keys_are_optimized(True)

        return se

    def add_SECR(self, m, p1, p2, slope=1, dpr=3.4):

        r = IMP.threading.StructureElementConnectivityRestraint(m, IMP.core.HarmonicUpperBound(0, slope), p1, p2, dpr, "")
        return r

    def setup_test_system(self, pdb_name, seq, elements):

        m = IMP.Model()
        ######################################
        DATADIR = os.path.abspath(os.path.join(os.path.dirname(__file__), 'input'))
        pdbfile = os.path.join(DATADIR, pdb_name)
        root_hier = IMP.atom.read_pdb(pdbfile, m)
        
        se = []

        for e in elements:
            se_ = self.setup_structural_element(root_hier, e)
            se.append(se_)
            
        se = self.sort_ses(se)

        #######################
        # Set up Sequence
        #######################
        
        seq_chain = IMP.atom.Chain.setup_particle(IMP.Particle(m), 'A')
        root_hier.add_child(seq_chain)

        for i in range(len(seq)):
            res = IMP.atom.Residue.setup_particle(IMP.Particle(m),
                                                            IMP.pmi.alphabets.amino_acid.get_residue_type_from_one_letter_code(seq[i]),
                                                            i+1)

            IMP.atom.Mass.setup_particle(res.get_particle(), IMP.atom.get_mass(res.get_residue_type()))
            IMP.core.XYZR.setup_particle(res.get_particle())
            IMP.atom.SecondaryStructureResidue.setup_particle(res.get_particle(), 0.8, 0, 0.2)
            seq_chain.add_child(res)

        return m, se, seq_chain

    def test_SECRS(self):
        pdb_name = 'pdb2lv9_A.ent'
        seq = 'SQPAKKTYTWNTKEEAKQAFKEALKEKRVPSNASWEQAMKMIINDPRYSALAKLSEKKQAFNAYKVQTEK'

        # Define Structural Elements
        # These are the residues in the PDB that correspond to structural elements.
        # (start_res, length, SSID)
        #-------------------------
        
        elements=[(35,5,'S'),(23,5,'S')]
    
        m, se, seq_chain = self.setup_test_system(pdb_name, seq, elements)
        
        #######################
        # Set up Scoring Function
        #######################

        rests = []
        se_rests = []
        se_pairs = []
        
        num_res = []
        for i in range(len(se)-1):
            num_res.append(int(se[i+1].get_first_residue_number()) - int(se[i].get_last_residue_number()))
            se_pairs.append((se[i].get_particle_index(), se[i+1].get_particle_index()))

        for s in se_pairs:
            r = self.add_SECR(m, s[0], s[1])

            self.assertEqual(num_res[se_pairs.index(s)], r.get_number_of_residues())
            rests.append(r)
            se_rests.append(r)


        for i in range(10):
            se[0].set_start_res_key(se[0].get_start_res()+1)
            r = se_rests[0]
            max_dist = r.get_max_distance()
            print("SECR:", se[0].get_start_res(), se[1].get_start_res(), r.get_number_of_residues(), r.get_model_distance(), r.get_max_distance(), r.unprotected_evaluate(None))
            # maximum distance between two SSEs should not be zero, which means no residues between these two secondary structure elements, this should be penalized by having a high score
            if round(max_dist) == 0.0:
                self.assertGreater(r.unprotected_evaluate(None), 1000.0)




if __name__ == '__main__':
    IMP.test.main()
