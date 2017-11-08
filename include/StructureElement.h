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


IMPTHREADING_BEGIN_NAMESPACE

//!
/* This Decorator defines a set of contiguous residues with a single coordinate per residue
  (IMP::core::XYZ). Internally, these residues have residue indexes from 1-N. 

  IntKey decorators then assign these residues and their associated structure to a larger
  sequence. 
    "start_res" : >=1    :: the starting residue of the sequence it is assigned to
    "polarity"  : (1;-1) :: the direction at which the sequence is assigned (+1 for forward, -1 for reverse)
    "length"    : (0,N)  :: the number of residues with which to assign to the sequence
    "offset"    : (0,N-1):: the foffset from the first residue of the Fragment 
        (i.e. offset = 2; polarity = 1; Fragment residue 3 is assigned to "start_res")
 */

//class StructureElement;

class IMPTHREADINGEXPORT StructureElement : public Decorator {

  std::vector<int> resindex_transform_;


  static void do_setup_particle(Model *m, ParticleIndex pi,
                                double start_res = 0, double polarity = 1,
                                double length = 0, double offset = 0) {
    //std::cout << "DSP " << m->get_particle(pi)->has_attribute(get_transformation_key(0)) << " " << m->get_particle(pi)->has_attribute(get_transformation_key(1)) << std::endl;
    m->add_attribute(get_start_res_key(), pi, start_res);
    IMP_USAGE_CHECK( (polarity == -1 || polarity == 1), "Polarity must be 1 or -1" );
    m->add_attribute(get_polarity_key(), pi, polarity);
    m->add_attribute(get_length_key(), pi, length);
    m->add_attribute(get_offset_key(), pi, offset);
  }

 public:

  algebra::Vector3Ds coordinates;

  static bool get_is_setup(Model *m, ParticleIndex pi) {
    /*
    IMP_USAGE_CHECK(
        ( m->get_has_attribute(get_start_res_key(), pi) &&
          m->get_has_attribute(get_polarity_key(), pi) &&
          m->get_has_attribute(get_offset_key(), pi) &&
          m->get_has_attribute(get_length_key(), pi) 
        ),
        "Particle must have all four transformation attributes");
    */
    return m->get_has_attribute(get_length_key(), pi);

  }


  //bool get_contains_sequence_residue(int rindex) const;
  //IMP_DECORATOR_SETUP_0(StructureElement);
  // Decorator setup works as:  ( NAME, v1type, v1name, v2type, v2name ...)
  IMP_DECORATOR_METHODS(StructureElement, Decorator);
  //IMP_DECORATOR_SETUP_0(StructureElement);
  IMP_DECORATOR_SETUP_4(StructureElement, 
                        int, start_res, 
                        int, polarity, 
                        int, length, 
                        int, offset);

  static FloatKey get_start_res_key();
  static FloatKey get_polarity_key();
  static FloatKey get_length_key();
  static FloatKey get_offset_key();

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


  void set_coordinates(algebra::Vector3Ds c){
    coordinates = c;
  };

  // Return all of the coordinates
  algebra::Vector3Ds get_all_coordinates(){
    return coordinates;
  };

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


  // Return the [offset, offset+length] coordinates
  //algebra::Vector3Ds get_coordinates();

  // Get the list of residue indexes that these
  // coordinates map to.
  Ints get_resindex_list() {
    int length = get_length();
    Ints resindex_transform;
    // The polarity key determines whether we build the
    // residue list forwards or backwards
    if (get_polarity() == 1){
      for (unsigned int i=0; i<length; i++){
        resindex_transform.push_back(i + get_start_res());
      }
    } else {
      for (unsigned int i=length; i>0; i--){
        resindex_transform.push_back(i + get_start_res()-1);
      }      
    }
    return resindex_transform;
  };

  // Functions for modifying keys
  void set_start_res_key(float i){ get_particle()->set_value(get_start_res_key(), i);};
  void set_polarity_key(float i){ get_particle()->set_value(get_polarity_key(), i);};
  void flip_polarity_key(){ int o = get_particle()->get_value(get_polarity_key());
                            get_particle()->set_value(get_polarity_key(), o * -1.);
                          }
  void set_length_key(float i){ get_particle()->set_value(get_length_key(), i);};
  void set_offset_key(float i) { get_particle()->set_value(get_offset_key(), i);};

  //ParticleIndexes get_structure_particle_indexes(){return pis_;};
};
IMP_DECORATORS(StructureElement, StructureElements, Decorator);

IMPTHREADING_END_NAMESPACE

#endif /* IMPTHREADING_STRUCTURE_ELEMENT_H */
