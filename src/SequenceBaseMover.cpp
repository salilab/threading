/**
 *  \file FragmentBaseMover.cpp  \brief A modifier to a Fragment object that perturbs 
 *     the start residue by an integer value.
 *
 *  Copyright 2007-2016 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/threading/SequenceBaseMover.h>


//#include <IMP/atom/Fragment.h>
#include <IMP/atom/Hierarchy.h>
//#include <IMP/core/XYZ.h>
//#include <IMP/atom/Residue.h>
#include <IMP/Model.h>
#include <IMP/random.h>
#include <IMP/algebra/vector_generators.h>
#include <IMP/core/MonteCarloMover.h>
#include <boost/random/uniform_int_distribution.hpp>
#include <random>

IMPTHREADING_BEGIN_NAMESPACE

namespace {
std::string get_sequence_base_mover_name(Model *m, ParticleIndex fpi) {
  return "SequenceBaseMover-" + m->get_particle(fpi)->get_name();
}
}

void SequenceBaseMover::initialize(ParticleIndex fpi, ParticleIndexes rpis) {
  fpi_ = fpi;
  rpis_ = rpis;
  start_res_ = IntKey("start_res");
  polarity_ = IntKey("polarity");
  original_ = 0;

  // This statement can be compiled
  f_hier_ = IMP::atom::Hierarchy();

  // This statement cannot be compiled
  f_hier_ = IMP::atom::Hierarchy(get_model(), fpi_);


  //f_resis_ = f_hier_.get_children();
  //std::cout << "FRG" << f_resis_ <<std::endl;
  //xyzs_.resize(resis.size());
  //rtypes_.resize(resis.size());

  //for (unsigned int i = 0; i < resis.size(); ++i) {
  //  xyzs_[i] = XYZ(resis[i].get_particle());
    //rtypes_[i] = IMP::atom::Residue(resis[i].get_particle()).get_residue_type();
  //}

  //std::cout << xyzs_ << rtypes_ <<std::endl;

}

SequenceBaseMover::SequenceBaseMover(Model *m, ParticleIndex fpi, ParticleIndexes rpis)
    : IMP::core::MonteCarloMover(m, get_sequence_base_mover_name(m, fpi)) {
  initialize(fpi, rpis);
}

IMP::core::MonteCarloMoverResult SequenceBaseMover::do_propose() {
  IMP_OBJECT_LOG;
  // First, get the original value of start_res
  original_ = get_model()->get_attribute(start_res_, fpi_);
  //Ints f_resis = IMP::atom::Fragment(get_model(), fpi_).get_residue_indexes();
/*
  for (unsigned int i = 0; i < xyzs_.size(); ++i){
    std::cout << "DO_PROPOSE :: " << i << xyzs_[i] << std::endl;
  }
*/
  //std::cout << "DO_PROPOSE :: " << f_resis << original_ << std::endl;
  //int trans = 0;
  //int nv = 0;

    //int polarity = 0;
    //IntKey ikp;
    //ikp = IntKey("polarity");
    //polarity = get_model()->get_attribute(ikp, pis_[i]);

    //std::cout << "TRANS_prop " << trans << " " << nv << " Polarity " << ikp << polarity <<std::endl;     
   
  IMP_USAGE_CHECK(
      get_model()->get_is_optimized(start_res_, fpi_),
      "SequenceBaseMover can't move non-optimized attribute. "
          << "particle: " << get_model()->get_particle_name(fpi_)
          << "attribute: " << start_res_);
  get_model()->set_attribute(start_res_, fpi_, 10);

  //std::cout << resis;
  //Ints nrr = std::irange(nv, nv + pis_[i]->get_residue_indexes())
  ///IMP::core::Fragment(pis_[i])->set_residue_indexes()
  //std::cout << "FRAGS" << i << originals_[i] << trans << get_model()->get_attribute(key_, pis_[i]) << std::endl;
  IMP_LOG_TERSE("FRAGS " i << original_  
              << get_model()->get_attribute(start_res_, fpi_) << std::endl);

  //ParticleIndexes out = fpi_;
  return IMP::core::MonteCarloMoverResult(rpis_, 1.0);
}


void SequenceBaseMover::do_reject() {
  IMP_OBJECT_LOG;
  get_model()->set_attribute(start_res_, fpi_, original_);

}

ModelObjectsTemp SequenceBaseMover::do_get_inputs() const {
  ModelObjectsTemp ret(rpis_.size());
  for (unsigned int i = 0; i < rpis_.size(); ++i) {
    ret[i] = get_model()->get_particle(rpis_[i]);
  }

  return ret;
}

IMPTHREADING_END_NAMESPACE
