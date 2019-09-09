from __future__ import print_function
import IMP
import IMP.test
import IMP.threading
import IMP.core
import IMP.atom
import os
import time
import io



class Tests(IMP.test.TestCase):

    def setup_structure_element(self, p, start_res, polarity, length, offset, nres=10):
        model = p.get_model()
        IMP.threading.StructureElement.setup_particle(p, start_res, polarity, length, offset)
        h = IMP.atom.Hierarchy.setup_particle(p)
        se = IMP.threading.StructureElement(p) 

        for c in self.make_coordinates(nres):
            np = IMP.Particle(model)
            hp = IMP.atom.Hierarchy.setup_particle(np)
            xyz = IMP.core.XYZR.setup_particle(np)
            xyz.set_coordinates(c)
            IMP.atom.Mass.setup_particle(np, 1.0)
            h.add_child(hp)
        
        return se   

    def make_coordinates(self, ncoord=3):
        coords = []
        for i in range(ncoord):
            coords.append(IMP.algebra.Vector3D(1.1,0.1,i*2))
        return coords

    def test_setup_particle(self):
        m = IMP.Model()
        p = IMP.Particle(m)

        se = self.setup_structure_element(p, 2, 1, 5, 0)

        sem = IMP.threading.StructureElementMolecule.setup_particle(p)
        
        self.assertTrue(IMP.threading.StructureElementMolecule.get_is_setup(p))


    @IMP.test.expectedFailure
    def test_fail_if_SE_not_setup(self):
        m = IMP.Model()
        p = IMP.Particle(m)
        sem = IMP.threading.StructureElementMolecule.setup_particle(p)

    def test_semol_lists(self):
        m = IMP.Model()
        p = IMP.Particle(m)

        se = self.setup_structure_element(p, 2, 1, 5, 0)

        semol = IMP.threading.StructureElementMolecule.setup_particle(p)

        semol.set_molname_list(["Mol1", "Mol2"])

        self.assertEqual(semol.get_molnames()[0], "Mol1")
        self.assertEqual(semol.get_maxcopy_list()[0], -1)
        self.assertEqual(len(semol.get_maxcopy_list()), 2)

        semol.set_maxcopy_list([1,1])
        self.assertEqual(semol.get_maxcopy_list()[0], 1)

if __name__ == '__main__':
    IMP.test.main()
