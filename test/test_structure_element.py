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
        IMP.threading.StructureElement.setup_particle(p, start_res, polarity, length, offset, 'A')
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

        self.assertEqual(se.get_start_res(), 2)
        self.assertEqual(se.get_offset(), 0)
        self.assertEqual(se.get_polarity(), 1)
        self.assertEqual(se.get_length(), 5)

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

    def test_get_coordinates(self):
        m = IMP.Model()
        p = IMP.Particle(m)
        se = self.setup_structure_element(p, 2, 1, 5, 0) 
        coordinates = self.make_coordinates(ncoord=10)
        self.assertEqual(len(se.get_coordinates()), len(coordinates[0:5]))
        self.assertEqual(se.get_coordinates()[1][2], coordinates[0:5][1][2])
        se.set_offset_key(2)   
        self.assertEqual(len(se.get_coordinates()), len(coordinates[2:7]))
        self.assertEqual(se.get_coordinates()[1][2], coordinates[2:7][1][2])

        se.set_length_key(3)
        self.assertEqual(len(se.get_coordinates()), len(coordinates[2:5]))
        self.assertEqual(se.get_coordinates()[1][2], coordinates[2:5][1][2])

    def test_coordinate_resindex_compatability(self):
        # Make sure the resindex list and coordinate list have the same number
        # of elements.
        m = IMP.Model()
        p = IMP.Particle(m)
        se = self.setup_structure_element(p, 2, 1, 5, 0) 

        self.assertEqual(len(se.get_resindex_list()), len(se.get_coordinates()))

      

if __name__ == '__main__':
    IMP.test.main()
