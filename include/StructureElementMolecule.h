/**
 *  \file IMP/threading/StructureElementMolecule.h
 *  \brief A special type of IMP.atom.Fragment that contains
 *   keys for specifying sequence assignment
 *
 *  Copyright 2007-2017 IMP Inventors. All rights reserved.
 *
 */

#ifndef IMPTHREADING_STRUCTURE_ELEMENT_MOLECULE_H
#define IMPTHREADING_STRUCTURE_ELEMENT_MOLECULE_H

#include <IMP/threading/threading_config.h>

#include <IMP/base_types.h>
//#include <IMP/exception.h>
#include <IMP/decorator_macros.h>
#include <IMP/Decorator.h>
#include <IMP/algebra/Vector3D.h>
#include <IMP/atom/Hierarchy.h>
#include <IMP/core/XYZ.h>
IMPTHREADING_BEGIN_NAMESPACE

//!
/* This Decorator holds a list of molecule names and list of integers of number of copies
   and two FloatKeys that act as integer indices to those lists. It is used in conjunction with 
   threading::StructureElement to define the sequence assignment of the coordinates in the 
   StructureElement.

   StructureElementMolecule REQUIRES that the particle already be set up as a StructureElement

    "molname"   : string :: The name of the molecule to which this SE is assigned
    "copy"      : int    :: the copy number of the molecule to which it is assigned (-1 means no copy number) 
 
 
    Lists of molecule names and, optionally, copy numbers must be added.

    E.G.
    molnames = ["Mol0", "Mol1", "Mol2", "Mol3"]
    maxcopy_list = [1, 1, 3, -1]

    means two copies of "Mol0" and "Mol1", four copies of "Mol2"
    and one copy of "Mol3". 

    -1 is used to indicate that copy number should NOT be used during selection (see StructureElementMover.cpp::transform_coordinates())

    Copy number begins at zero, so 3, means that there exist "Mol3.0", "Mol3.1", "Mol3.2", "Mol3.3"

    PYTHON USAGE:

    # Using the lists defined above, we want to initialize this SEMOL to "Mol2" copy 1
    se = IMP.threading.StructureElement(model, pi, start_res, polarity, length, offset)
    semol = IMP.threading.StructureElementMolecule(model, pi)
    semol.set_molname_list(molnames, maxcopy_list_)
    semol.set_molname_index_key[2]
    semol.set_copy_index_key[1]
 
 */

class IMPTHREADINGEXPORT StructureElementMolecule : public Decorator {
  ParticleIndex pi_;

  static void do_setup_particle(Model *m, ParticleIndex pi) {

    if(StructureElement::get_is_setup(m, pi)==false){ 
      IMP_FAILURE("Particle must be set up as a StructureElement");        
    };

    ParticleIndex pi_ = pi;
    m->add_attribute(get_molname_index_key(), pi_, 0);
    m->add_attribute(get_copy_key(), pi_, 0);
  
  }

 public:

  Strings molnames_;
  Ints maxcopy_list_;

  static bool get_is_setup(Model *m, ParticleIndex pi) {
    return m->get_has_attribute(get_molname_index_key(), pi);
  }

  IMP_DECORATOR_METHODS(StructureElementMolecule, Decorator);
  // Decorator setup works as:  ( NAME, v1type, v1name, v2type, v2name ...)
  IMP_DECORATOR_SETUP_0(StructureElementMolecule);

  static FloatKey get_molname_index_key();
  static FloatKey get_copy_key();

  // Return all key values in a nice little vector
  //Floats get_all_key_values(){
  //  Floats vals;
  //  vals.push_back(get_molname());
  //  vals.push_back(get_copy());
  //  return vals;
  //}

  std::string get_molname() {
    std::string i = molnames_[int(get_model()->get_attribute(get_molname_index_key(), get_particle_index()))];
    return i;
  }

  // A list of strings corresponding to molecule names
  void set_molname_list(Strings molnames, Ints copies={}){
    molnames_ = molnames;
    
    // If there are no copies given, set all copies to zero
    // (We are assuming that all these molecules are singular)
    if(copies.size()==0){
       for(int i=0;i<molnames.size();i++){copies.push_back(-1);};
    }
    set_maxcopy_list(copies);
  }

  // maxcopy_list is the number of copies of each molecule.  
  void set_maxcopy_list(Ints copies){
    IMP_USAGE_CHECK(molnames_.size() == copies.size(), "Length of copy list must equal the length of the molecule name list");
    maxcopy_list_ = copies;
  }

  Ints get_maxcopy_list(){
    return maxcopy_list_;
  }

  Int get_maxcopy(){
    return maxcopy_list_[get_molname_index()];
  }

  Strings get_molnames(){
    return molnames_;
  }

  Int get_max_mol_index(){
    return molnames_.size();
  }
  
  int get_molname_index() {
    int i = int(get_model()->get_attribute(get_molname_index_key(), get_particle_index()));
    return i;
  }

  int get_copy() {
    int i = int(get_model()->get_attribute(get_copy_key(), get_particle_index()));
    return i;
  }

  // Return all of the coordinates
  algebra::Vector3Ds get_all_coordinates(){
    algebra::Vector3Ds coords;
    atom::Hierarchies hs=atom::Hierarchy(get_particle()).get_children();

    for (unsigned int i=0; i<hs.size(); i++){
      algebra::Vector3D c=core::XYZ(get_model(), hs[i].get_particle_index()).get_coordinates();
      coords.push_back(c);
    }
    return coords;
  };

  //! Set which keys are optimized
  void set_keys_are_optimized(bool tf) const {
    get_particle()->set_is_optimized(get_molname_index_key(), tf);
    get_particle()->set_is_optimized(get_copy_key(), tf);
  }

  void set_molname_is_optimized(bool tf) const {
    get_particle()->set_is_optimized(get_molname_index_key(), tf);
  }

  void set_copy_is_optimized(bool tf) const {
    get_particle()->set_is_optimized(get_copy_key(), tf);
  }

  bool get_molname_is_optimized() {
    get_particle()->get_is_optimized(get_molname_index_key());
  }
  bool get_copy_is_optimized() {
    get_particle()->get_is_optimized(get_copy_key());
  }

  // **********************************************
  // Functions for modifying keys

  void set_all_keys(Floats vals){
    IMP_USAGE_CHECK( (2 == vals.size()), "Must pass exactly two values");
    set_molname_index_key(vals[0]);
    set_copy_key(vals[1]);
  }

  // molname_index_key must be less than the size of the molnames_ list
  void set_molname_index_key(float i){
    IMP_USAGE_CHECK( (i < molnames_.size()), "Molname index too large");
    get_particle()->set_value(get_molname_index_key(), i);
  };
  
  // copy number must be less than the max copy
  void set_copy_key(float i){ 
    IMP_USAGE_CHECK( (i <= get_maxcopy()), "Copy number is too large");
    get_particle()->set_value(get_copy_key(), i);
  };

  // *********************************************************
  // Functions for determining available move sets

  int get_number_of_molecule_indices(){return molnames_.size();}
  int get_max_copy_number_indices(){return get_maxcopy();}

  //algebra::Vector3D get_last_coordinate(){}
};

IMP_DECORATORS(StructureElementMolecule, StructureElementMolecules, Decorator);

IMPTHREADING_END_NAMESPACE

#endif /* IMPTHREADING_STRUCTURE_ELEMENT_MOLECULE_H */
