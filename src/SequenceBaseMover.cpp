/**
 *  \file FragmentBaseMover.cpp  \brief A modifier to a Fragment object that perturbs 
 *     the start residue by an integer value.
 *
 *  Copyright 2007-2016 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/threading/SequenceBaseMover.h>


//#include <IMP/atom/Fragment.h>
//#include <IMP/atom/Hierarchy.h>
//#include <IMP/core/XYZ.h>
//#include <IMP/atom/Residue.h>
//#include <IMP/particle_index.h>
#include <IMP/random.h>
#include <IMP/algebra/vector_generators.h>
//#include <IMP/core/MonteCarloMover.h>
#include <boost/random/uniform_int_distribution.hpp>
#include <random>

IMPTHREADING_BEGIN_NAMESPACE


namespace {
std::string get_sequence_base_mover_name(Model *m, ParticleIndex fpi) {
  return "SequenceBaseMover-" + m->get_particle(fpi)->get_name();
}
}

void SequenceBaseMover::initialize(ParticleIndex fpi) {
  fpi_ = fpi; //this should eventually be a StructureElement 
  start_res_ = IntKey("start_res");
  polarity_ = IntKey("polarity");
  original_ = 0;

}

SequenceBaseMover::SequenceBaseMover(Model *m, ParticleIndex fpi)
    : IMP::core::MonteCarloMover(m, get_sequence_base_mover_name(m, fpi)) {
  initialize(fpi);
}

IMP::core::MonteCarloMoverResult SequenceBaseMover::do_propose() {
  IMP_OBJECT_LOG;
  // 1) get the original value of start_res
  original_ = get_model()->get_attribute(start_res_, fpi_);
  
  // 2) get available moves from the associated StructureElement 

  // 3) Chose a move from this set
  int new_sr = original_ - 1;

  // 4) Set the IntKey
  get_model()->set_attribute(start_res_, fpi_, new_sr);

  // 5) Run the macro within the StructureElement to move the coordinates
  //  from structure to sequence
  get_model()->get_particle(fpi_);

  //std::cout << resis;
  //Ints nrr = std::irange(nv, nv + pis_[i]->get_residue_indexes())
  ///IMP::core::Fragment(pis_[i])->set_residue_indexes()
  //std::cout << "FRAGS" << i << originals_[i] << trans << get_model()->get_attribute(key_, pis_[i]) << std::endl;
  IMP_LOG_TERSE("SequenceBaseMover " << original_  
              << get_model()->get_attribute(start_res_, fpi_) << std::endl);

  ParticleIndexes out(1);
  out[0] = fpi_;
  return IMP::core::MonteCarloMoverResult(out, 1.0);
}


void SequenceBaseMover::do_reject() {
  IMP_OBJECT_LOG;
  get_model()->set_attribute(start_res_, fpi_, original_);

}

ModelObjectsTemp SequenceBaseMover::do_get_inputs() const {
  ModelObjectsTemp ret(1);
  ret[0] = get_model()->get_particle(fpi_);
  //ModelObjectsTemp ret = {inp};
  return ret;
}

IMPTHREADING_END_NAMESPACE
