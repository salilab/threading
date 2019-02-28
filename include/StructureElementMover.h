/**
 *  \file IMP/threading/StructureElementMover.h
 *  \brief A modifier which moves the sequence in a fragment.
 *
 *  Copyright 2007-2016 IMP Inventors. All rights reserved.
 *
 */

#ifndef IMPTHREADING_STRUCTURE_ELEMENT_MOVER_H
#define IMPTHREADING_STRUCTURE_ELEMENT_MOVER_H

#include <IMP/threading/threading_config.h>

#include <IMP/core/MonteCarloMover.h>
//#include <IMP/threading/StructureElement.h>
#include <IMP/base_types.h>
#include <IMP/exception.h>
#include <IMP/atom/Hierarchy.h>
//#include <IMP/atom/Residue.h>
//#include <IMP/core/XYZ.h>
IMPTHREADING_BEGIN_NAMESPACE

//! A mover for modifying the sequence assignment keys 
//! from an IMP.threading.StructureElement.

class IMPTHREADINGEXPORT StructureElementMover : public IMP::core::MonteCarloMover {
  ParticleIndex pi_;
  ParticleIndex s_hier_pi_; // sequence root hierarchy node particle index (for right now, just a single chain)
  Floats orig_key_values_;
  int pct_flip_ = 0; // percentage of time to try a polarity flip move
  int pct_offset_ = 30; // percentage of time to try an offset move
  int pct_length_ = 30; // percentage of time to try a length move

  void initialize(Model *m, ParticleIndex pi, ParticleIndex s_hier_pi);

 public:

 StructureElementMover(Model *m, ParticleIndex pi, ParticleIndex s_hier_pi);

 void transform_coordinates();
 void zero_coordinates();

 void set_flip_rate(float i){
  bool j = check_if_one_to_hundred_(i);
  IMP_USAGE_CHECK((j == true), "StructureElementMover::set_flip_rate - Value must be between 0 and 100");
  pct_flip_ = int(i);
 }

 void set_offset_rate(float i){
  bool j = check_if_one_to_hundred_(i);
  IMP_USAGE_CHECK((j == true), "StructureElementMover::set_offset_rate - Value must be between 0 and 100");
  pct_offset_ = int(i);
 }

 void set_length_rate(float i){
  bool j = check_if_one_to_hundred_(i);
  IMP_USAGE_CHECK((j == true), "StructureElementMover::set_length_rate - Value must be between 0 and 100");
  pct_length_ = int(i);
 }

 protected:
  virtual ModelObjectsTemp do_get_inputs() const IMP_OVERRIDE;

  //! Move particle attributes within a ball, as specified in constructor
  virtual IMP::core::MonteCarloMoverResult do_propose() IMP_OVERRIDE;

  //! restore original attributes from before do_propose
  virtual void do_reject() IMP_OVERRIDE;
  IMP_OBJECT_METHODS(StructureElementMover);

  bool check_if_one_to_hundred_(float i){
    if (i>=0 and i<100){
      return true;
    } else {
      return false;
    }
  }


};

IMPTHREADING_END_NAMESPACE

#endif /* IMPTHREADING_STRUCTURE_ELEMENT_MOVER_H */
