/**
 *  \file IMP/threading/SSEThreadMover.h
 *  \brief This mover updates the sequence assignment of a set of StructureElements
 *  with respect to a given sequence hierarchy. This class is meant to be a base-class
 *  for specific move types that modify specific SE keys.  
 *  
 *  THis object is not meant to be used as an actual mover.  do_propose() 
 *  and do_reject() must be formulated for the specific move type in an inherited object
 *
 *  Copyright 2007-2020 IMP Inventors. All rights reserved.
 *
 */

#ifndef IMPTHREADING_SSETHREAD_MOVER_H
#define IMPTHREADING_SSETHREAD_MOVER_H

#include <IMP/threading/threading_config.h>
#include <IMP/threading/StructureElement.h>
#include <IMP/core/MonteCarloMover.h>
#include <IMP/base_types.h>
#include <IMP/exception.h>
#include <IMP/atom/Hierarchy.h>
#include <boost/random.hpp>
//#include <IMP/statistics/statistics_config.h>
//#include <IMP/atom/Residue.h>
//#include <IMP/core/XYZ.h>
IMPTHREADING_BEGIN_NAMESPACE

//! A mover for modifying sequence assignment keys for structure elements in an SSETHread system
//! from an IMP.threading.StructureElement.
//! 
//! Input is a list of structure_element indexes and the structural hierarchy particle index
//!
//! The class allows for easy random selection of a StructureElement and the modification of the 
//! sequence hierarchy based on the keys from the StructureElement
//!
//! ####  This class cannot be used as is. ####  
//! 
//! Please inherit from it and define do_propose(), do_reject() and do_accept()
//! to create the mover behavior you want

class IMPTHREADINGEXPORT SSEThreadMover : public IMP::core::MonteCarloMover {
  ParticleIndexes se_pis_;     // Structure Element particle indexes
  ParticleIndex s_hier_pi_;  // Sequence root hierarchy node particle index

  void initialize(ParticleIndexes se_pis, ParticleIndex s_hier_pi);

 public:

 SSEThreadMover(Model *m, ParticleIndexes se_pis, ParticleIndex s_hier_pi);

 void transform_coordinates(ParticleIndex se_pi);
 void zero_coordinates(ParticleIndex se_pi);

StructureElement get_random_structure_element_particle() {
  // Get a random particle index from the list
  ::boost::uniform_int<> rand(0, se_pis_.size() - 1);
  StructureElement ret = StructureElement(get_model(), se_pis_[int(rand(random_number_generator))]);  
  return ret;
}

 protected:

  // These functions must be defined in classes that inherit this object.
  virtual ModelObjectsTemp do_get_inputs() const IMP_OVERRIDE;
  virtual IMP::core::MonteCarloMoverResult do_propose() IMP_OVERRIDE;
  virtual void do_reject() IMP_OVERRIDE;
  
  IMP_OBJECT_METHODS(SSEThreadMover);

};

IMPTHREADING_END_NAMESPACE

#endif /* IMPTHREADING_SSETHREAD_MOVER_H */
