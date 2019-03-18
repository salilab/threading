from __future__ import print_function
import IMP
import IMP.test
import IMP.threading
import IMP.core
import IMP.atom
import IMP.algebra
import os
import time
import io



class Tests(IMP.test.TestCase):

    def setup_structure_element(self, p, start_res, polarity, length, offset):
        IMP.threading.StructureElement.setup_particle(p, start_res, polarity, length, offset)
        se = IMP.threading.StructureElement(p)  
        return se  

    def make_coordinates(self, ncoord=3):
        coords = []
        for i in range(ncoord):
            coords.append(IMP.algebra.Vector3D(1.1,0.1,i*2))
        return coords

    def test_setup_particle(self):
        m = IMP.Model()
        p = IMP.Particle(m)
        se = IMP.threading.StructureElement(p, coords)
        print(se.get_polarity(), se.length(), se.get_coordinates())

    def test_flip_polarity_key(self):
        m = IMP.Model()
        p = IMP.Particle(m)
        se = self.setup_structure_element(p, 2, 1, 5, 0)

        self.assertEqual(se.get_polarity(), 1)  
        se.flip_polarity_key()
        self.assertEqual(se.get_polarity(), -1) 

    def test_set_length_key(self):
        m = IMP.Model()
        p = IMP.Particle(m)
        se = self.setup_structure_element(p, 2, 1, 5, 0)

        self.assertEqual(se.get_length(), 5)  
        se.set_length_key(7)
        self.assertEqual(se.get_length(), 7)       

    def test_set_offset_key(self):
        m = IMP.Model()
        p = IMP.Particle(m)
        se = self.setup_structure_element(p, 2, 1, 5, 0)

        self.assertEqual(se.get_offset(), 0)  
        se.set_offset_key(1)
        self.assertEqual(se.get_offset(), 1)  

    def test_set_start_res_key(self):
        m = IMP.Model()
        p = IMP.Particle(m)
        se = self.setup_structure_element(p, 2, 1, 5, 0)

        self.assertEqual(se.get_start_res(), 2)  
        se.set_start_res_key(11)
        self.assertEqual(se.get_start_res(), 11) 

    def test_resindex_list(self):
        m = IMP.Model()
        p = IMP.Particle(m)
        se = self.setup_structure_element(p, 2, 1, 5, 0)

        self.assertEqual(se.get_resindex_list(), [2,3,4,5,6])
        se.flip_polarity_key()
        self.assertEqual(se.get_resindex_list(), [6,5,4,3,2])

if __name__ == '__main__':
    IMP.test.main()
