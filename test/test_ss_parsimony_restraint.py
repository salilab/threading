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

    def setup_structural_element(self, root_hier, element, sse_chain_id):
        ## SETUP STRUCTURAL ELEMENTS
        # Get Particle Indexes for CA atoms 
        # from residue range element[0]:element[1]
        # From chain X
        #print(range(element[0],element[0]+element[1]))
        pis = IMP.atom.Selection(root_hier, chain_id='A', 
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
            xyz.set_coordinates(IMP.core.XYZ(p).get_coordinates())
            h.add_child(hp)

        # Setting up structure element, element[0] is start_res, 1 for polarity, element[1] for length of the element, 0 for offset and finally chain id
        IMP.threading.StructureElement.setup_particle(root_hier.get_model(), pi.get_index(), element[0], 1, element[1], 0, sse_chain_id)

        # Set up this element as a helix or strand or a coil

        probmap = {'H': (1, 0, 0), 'S': (0, 1, 0), 'C': (0, 0, 1)}
        prob = probmap.get(element[-1], (0.33, 0.33, 0.33))
        IMP.atom.SecondaryStructureResidue.setup_particle(pi, *prob)


        se = IMP.threading.StructureElement(pi)
        se.set_keys_are_optimized(True)

        return se
    def set_up_system_cal_score(self, pdb_name, seq, elements, seq_chain_id, sse_chain_id):
#       seq = 'SQPAKKTYTWNTKEEAKQAFKEALKEKRVPSNASWEQAMKMIINDPRYSALAKLSEKKQAFNAYKVQTEK'

        m = IMP.Model()
        ######################################
        DATADIR = os.path.abspath(os.path.join(os.path.dirname(__file__), 'input'))
        pdbfile = os.path.join(DATADIR, pdb_name)
        root_hier = IMP.atom.read_pdb(pdbfile, m)



        se = []

        for e in elements:
            se.append(self.setup_structural_element(root_hier, e, sse_chain_id))

        se = self.sort_ses(se)

        #######################
        # Set up Sequence
        #######################

        seq_chain = IMP.atom.Chain.setup_particle(IMP.Particle(m), seq_chain_id)
        seq_chain.set_name(seq_chain.get_id())

        root_hier.add_child(seq_chain)

        for i in range(len(seq)):
            res = IMP.atom.Residue.setup_particle(IMP.Particle(m),
                                                            IMP.pmi.alphabets.amino_acid.get_residue_type_from_one_letter_code(seq[i]),
                                                            i+1)
            IMP.atom.Mass.setup_particle(res.get_particle(), IMP.atom.get_mass(res.get_residue_type()))
            IMP.core.XYZR.setup_particle(res.get_particle())
            IMP.atom.SecondaryStructureResidue.setup_particle(res.get_particle(), 1.0, 0.0, 0.0)
            seq_chain.add_child(res)

        #######################
        # Set up Scoring Function
        #######################

        # SS Parsimony Restraint
        psipred_slope = 1.0

        # Restraint evaluated at the SE level, Pass SE and seq_chain.  
        # Have restraint look up residue from the seq_chain and get SS 
        chain_par = IMP.atom.Hierarchy(seq_chain.get_particle()).get_children()

        # get the secondary structure probabilities of the first chain residue
#        print(IMP.atom.SecondaryStructureResidue(m, chain_par[0].get_particle_index()).get_all_probabilities())

        par_res = IMP.threading.SecondaryStructureParsimonyRestraint(m, [s.get_particle_index() for s in se], seq_chain.get_particle(), psipred_slope)
        score  = par_res.unprotected_evaluate(None)
        return score

    def test_parsimony_res(self):
        pdb_name = 'pdb2lv9_A.ent'
        seq = 'SQPAKKTYTWNTKEEAKQAF'
        seq_chain_id = 'B'
        sse_chain_id = 'B'
        # Define Structural Elements
        # These are the residues in the PDB that correspond to structural elements.
        # (start_res, length, SSID)
        #-------------------------

        # all helix

        elements=[(3, 10,'H')]

        score = self.set_up_system_cal_score(pdb_name, seq, elements, seq_chain_id, sse_chain_id)

        # since 10 residues should be coiled ans our sequence length is 20 then (-0.1*10)/20
        self.assertEqual(round(score, 2), -0.05)

        # check what happens when the chain ids are not matching
        sse_chain_id = 'A'
        score = self.set_up_system_cal_score(pdb_name, seq, elements, seq_chain_id, sse_chain_id)
        # since chain ids are not matching, then (-0.1*20)/20 
        self.assertEqual(round(score, 2), -0.1)

        # check when sse type is not helix, while the sequence based sse assignment suggests it is a helix 
        
        elements=[(3, 10,'S')]
        score = self.set_up_system_cal_score(pdb_name, seq, elements, seq_chain_id, sse_chain_id)
        # since chain ids are not matching, then (-0.1*20)/20
        self.assertEqual(round(score, 2), -0.1)



if __name__ == '__main__':
    IMP.test.main()


