/**
 *  \file IMP/threading/StructureElement.h
 *  \brief A special type of IMP.atom.Fragment that contains
 *   keys for specifying sequence assignment
 *
 *  Copyright 2007-2017 IMP Inventors. All rights reserved.
 *
 */

#ifndef IMPTHREADING_STRUCTURE_ELEMENT_H
#define IMPTHREADING_STRUCTURE_ELEMENT_H

#include <IMP/threading/threading_config.h>

//#include <IMP/base_types.h>
//#include <IMP/exception.h>
#include <IMP/decorator_macros.h>
#include <IMP/Decorator.h>
#include <IMP/algebra/Vector3D.h>
#include <IMP/atom/Hierarchy.h>
#include <IMP/core/XYZ.h>
IMPTHREADING_BEGIN_NAMESPACE

//!
/* This Decorator defines a set of contiguous residues with a single coordinate per residue
  (IMP::core::XYZ). Internally, these residues have residue indexes from 1-N. 

  IntKey decorators then assign these residues and their associated structure to a larger
  sequence. 
    "start_res" : >=1    :: the starting residue of the sequence it is assigned to
    "polarity"  : (1;-1) :: the direction at which the sequence is assigned (+1 for forward, -1 for reverse)
    "length"    : (0,N)  :: the number of residues with which to assign to the sequence
    "offset"    : (0,N-1):: the offset from the first residue of the Fragment 
       (i.e. offset = 2; polarity = 1; Fragment residue 3 is assigned to "start_res")
 */

class IMPTHREADINGEXPORT StructureElement : public Decorator {
  ParticleIndex pi_;

  static void do_setup_particle(Model *m, ParticleIndex pi,
                                double start_res = 0, double polarity = 1,
                                double length = 0, double offset = 0) {

    ParticleIndex pi_ = pi;

    m->add_attribute(get_start_res_key(), pi_, start_res);
    IMP_USAGE_CHECK( (polarity == -1 || polarity == 1), "Polarity must be 1 or -1" );
    m->add_attribute(get_polarity_key(), pi_, polarity);
    m->add_attribute(get_length_key(), pi_, length);
    m->add_attribute(get_offset_key(), pi_, offset);
  
  }
  int n_coordinates_;
 public:

  static bool get_is_setup(Model *m, ParticleIndex pi) {
    return m->get_has_attribute(get_length_key(), pi);
  }

  IMP_DECORATOR_METHODS(StructureElement, Decorator);
  // Decorator setup works as:  ( NAME, v1type, v1name, v2type, v2name ...)
  IMP_DECORATOR_SETUP_4(StructureElement, 
                        int, start_res, 
                        int, polarity, 
                        int, length, 
                        int, offset);

  static FloatKey get_start_res_key();
  static FloatKey get_polarity_key();
  static FloatKey get_length_key();
  static FloatKey get_offset_key();

  // Return all key values in a nice little vector
  Floats get_all_key_values(){
    Floats vals;
    vals.push_back(get_start_res());
    vals.push_back(get_polarity());
    vals.push_back(get_length());
    vals.push_back(get_offset());
    return vals;
  }

  int get_start_res() {
    int i = int(get_model()->get_attribute(get_start_res_key(), get_particle_index()));
    return i;
  }

  int get_polarity() {
    float p = get_model()->get_attribute(get_polarity_key(), get_particle_index());
    int i = int(p);
    return i;
  }

  int get_offset() {
    int i = int(get_model()->get_attribute(get_offset_key(), get_particle_index()));
    return i;
  }

  int get_length() {
    int i = int(get_model()->get_attribute(get_length_key(), get_particle_index()));
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

  // Return the [offset, offset+length] coordinates
  algebra::Vector3Ds get_coordinates(){
    int offset = get_offset();
    int length = get_length();

    algebra::Vector3Ds coords;
    atom::Hierarchies hs=atom::Hierarchy(get_particle()).get_children();
    //std::cout << "PIx " << pi_ << " " << hs << " " << offset << " " << length << " " << hs <<std::endl;
    for (int i=0; i<length; i++){
      int ix = i + offset;
      //std::cout << "GET_COORDS " << i << " " << ix << " " << hs[ix] << " " << hs[ix].get_particle_index() << " " << core::XYZ(get_model(), hs[ix].get_particle_index()).get_coordinates() <<std::endl;
      algebra::Vector3D c=core::XYZ(get_model(), hs[ix].get_particle_index()).get_coordinates();
      coords.push_back(c);
    }
    return coords;
  }


  //! Set which keys are optimized
  void set_keys_are_optimized(bool tf) const {
    get_particle()->set_is_optimized(get_offset_key(), tf);
    get_particle()->set_is_optimized(get_start_res_key(), tf);
    get_particle()->set_is_optimized(get_length_key(), tf);
    get_particle()->set_is_optimized(get_polarity_key(), tf);
  }

  void set_polarity_is_optimized(bool tf) const {
    get_particle()->set_is_optimized(get_polarity_key(), tf);
  }

  void set_start_res_is_optimized(bool tf) const {
    get_particle()->set_is_optimized(get_start_res_key(), tf);
  }

  void set_length_is_optimized(bool tf) const {
    get_particle()->set_is_optimized(get_length_key(), tf);
  }

  void set_offset_is_optimized(bool tf) const {
    get_particle()->set_is_optimized(get_offset_key(), tf);
  }
  
  bool get_start_res_is_optimized() {
    return get_particle()->get_is_optimized(get_start_res_key());
  }
  bool get_polarity_is_optimized() {
    return get_particle()->get_is_optimized(get_polarity_key());
  }
  bool get_length_is_optimized() {
    return get_particle()->get_is_optimized(get_length_key());
  }
  bool get_offset_is_optimized() {
    return get_particle()->get_is_optimized(get_offset_key());
  }

  // Get the list of residue indexes that these
  // coordinates map to.
  Ints get_resindex_list() {
    int length = get_length();
    Ints resindex_transform;
    // The polarity key determines whether we build the
    // residue list forwards or backwards
    if (get_polarity() == 1){
      for (int i=0; i<length; i++){
        resindex_transform.push_back(i + get_start_res());
      }
    } else {
      for (int i=length; i>0; i--){
        resindex_transform.push_back(i + get_start_res()-1);
      }      
    }
    return resindex_transform;
  };

  // **********************************************
  // Functions for modifying keys

  void set_all_keys(Floats vals){
    IMP_USAGE_CHECK( (4 == vals.size()), "Must pass exactly four values");
    set_start_res_key(vals[0]);
    set_polarity_key(vals[1]);
    set_length_key(vals[2]);
    set_offset_key(vals[3]);
  }

  // start_res_key can be anything in the structure
  void set_start_res_key(float i){ 
    get_particle()->set_value(get_start_res_key(), i);
  };

  // polarity_key must be 1 or -1.
  void set_polarity_key(float i){ 
    IMP_USAGE_CHECK( (i == -1 || i == 1), "Polarity must be 1 or -1" );
    get_particle()->set_value(get_polarity_key(), i);
  };

  void flip_polarity_key(){ 
    int o = get_particle()->get_value(get_polarity_key());
    get_particle()->set_value(get_polarity_key(), o * -1.);
  };

  // Length key cannot be greater than the number of elements in the 
  void set_length_key(float i){ 
    atom::Hierarchies hs=atom::Hierarchy(get_particle()).get_children();
    n_coordinates_ = hs.size();
    //std::cout << "set_length_key: " << i << " " << n_coordinates_ << std::endl;
    IMP_USAGE_CHECK( (i <= n_coordinates_), "Length key cannot be greater than number of coordinates in StructureElement");
    get_particle()->set_value(get_length_key(), i);
  };

  // Length plus offset cannot be greater than the number of coordinates
  void set_offset_key(float i) { 
    atom::Hierarchies hs=atom::Hierarchy(get_particle()).get_children();
    n_coordinates_ = hs.size();
    //std::cout << "L+soC+i " << i << " + " << get_length() << "=" << i + get_length() << " " << n_coordinates_ << std::endl;
    IMP_USAGE_CHECK( (i + get_length() <= n_coordinates_), "Length plus offset cannot be greater than number of coordinates in StructureElement");
    get_particle()->set_value(get_offset_key(), i);
  };

  // *********************************************************
  // Functions for determining available move sets

  int get_max_offset() {
    // The offset can range from 0 to sizeof(coordinates_) - length
    atom::Hierarchies hs=atom::Hierarchy(get_particle()).get_children();
    n_coordinates_ = hs.size();
    //std::cout << n_coordinates_ << " ||| " << get_length() << std::endl;
    return n_coordinates_ - get_length();
  };

  int get_number_of_coordinates() {
    atom::Hierarchies hs=atom::Hierarchy(get_particle()).get_children();
    return hs.size();    
  };

  int get_last_residue_number(){ return get_start_res() + get_offset() + get_length() - 1; };
  int get_first_residue_number(){ return get_start_res() + get_offset() ; }

  algebra::Vector3D get_first_coordinate(){
    atom::Hierarchies hs=atom::Hierarchy(get_particle()).get_children();
    algebra::Vector3D c;
    int i = get_offset();
    if (get_polarity() == 1){
    c=core::XYZ(get_model(), hs[i].get_particle_index()).get_coordinates();
    } else {
    int j = get_length();
    c=core::XYZ(get_model(), hs[i+j].get_particle_index()).get_coordinates();    
    }
    return c;
  }

  //algebra::Vector3D get_last_coordinate(){}
};

IMP_DECORATORS(StructureElement, StructureElements, Decorator);

IMPTHREADING_END_NAMESPACE

#endif /* IMPTHREADING_STRUCTURE_ELEMENT_H */
