/**
 *  \file StructureElementMover.cpp  \brief A modifier to a Fragment object that perturbs 
 *     the start residue by an integer value.
 *
 *  Copyright 2007-2016 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/threading/StructureElementMover.h>



IMPTHREADING_BEGIN_NAMESPACE


namespace {
std::string get_structure_element_mover_name(Model *m, ParticleIndex pi) {
  return "StructureElementMover-" + m->get_particle(pi)->get_name();
}
}

void StructureElementMover::initialize(Model *m, ParticleIndex pi, ParticleIndex s_hier_pi) {

  IMP_USAGE_CHECK(StructureElement().get_is_setup(get_model(), pi), "Particle is not set up as a StructureElement");
  pi_ = pi;
  s_hier_pi_ = s_hier_pi;

}

StructureElementMover::StructureElementMover(Model *m, ParticleIndex pi, ParticleIndex s_hier_pi)
    : IMP::core::MonteCarloMover(m, get_structure_element_mover_name(m, pi)) {
  initialize(m, pi, s_hier_pi);
}

IMP::core::MonteCarloMoverResult StructureElementMover::do_propose() {
  IMP_OBJECT_LOG;

  // Initialize the structure element
  StructureElement se(get_model(), pi_);
  atom::Hierarchy h(get_model(), s_hier_pi_);
  // 1) get the original key values
  orig_key_values_ = se.get_all_key_values();

  int r;

  zero_coordinates();

  if (se.get_start_res_is_optimized()) {
    r = statistics::internal::random_int(2) * 2 - 1;
    se.set_start_res_key(orig_key_values_[0] - r);
    std::cout << "start_res: " << orig_key_values_[0] << " " << r << " " << se.get_start_res() << std::endl;
  };

  if (se.get_polarity_is_optimized() && rand() % 100 < pct_flip_) {
    se.flip_polarity_key();
  };

  if (se.get_length_is_optimized() && rand() % 100 < pct_length_) {
    // TODO: find available moves from the StructureElement class
    r = statistics::internal::random_int(2) * 2 - 1;
    se.set_length_key(orig_key_values_[2] - r);
    std::cout << "length: " << orig_key_values_[2] << " " << r << " " << se.get_length() << std::endl;
  };

  if (se.get_offset_is_optimized() && rand() % 100 < pct_offset_) {
    // TODO: find available moves from the StructureElement class
    r = statistics::internal::random_int(2) * 2 - 1;
    se.set_offset_key(orig_key_values_[3] - r);
    std::cout << "offset: " << orig_key_values_[3] << " " << r << " " << se.get_offset() << std::endl;
  };

  transform_coordinates();

  return IMP::core::MonteCarloMoverResult(h.get_children(), 1.0);
}


void StructureElementMover::zero_coordinates(){
  // Zero out the coordinates
  StructureElement se(get_model(), pi_);
  atom::Hierarchy h(get_model(), s_hier_pi_);
  atom::Selection sel(h);
  sel.set_residue_indexes(se.get_resindex_list());
  ParticlesTemp oldsel = sel.get_selected_particles();
  for (unsigned int i=0; i<se.get_resindex_list().size(); i++) {
    core::XYZ coord(oldsel[i]);
    coord.set_coordinates(algebra::Vector3D(0,0,0));
  }; 

}

void StructureElementMover::transform_coordinates(){
  // First, zero out the old coordiantes using the orig_keys_
  StructureElement se(get_model(), pi_);
  atom::Hierarchy h(get_model(), s_hier_pi_);
  atom::Selection sel(h);
  sel.set_residue_indexes(se.get_resindex_list());
  ParticlesTemp oldsel = sel.get_selected_particles();
  algebra::Vector3Ds new_coords = se.get_coordinates(); // does not take into account polarity!!
  for (unsigned int i=0; i<se.get_resindex_list().size(); i++) {
    core::XYZ coord(oldsel[i]);
    coord.set_coordinates(new_coords[i]);
  }; 
}

void StructureElementMover::do_reject() {
  IMP_OBJECT_LOG;

  StructureElement se(get_model(), pi_);
  zero_coordinates();
  //Return StructureElement key values to their original position
  std::cout << orig_key_values_ << std::endl;
  se.set_all_keys(orig_key_values_);
  //Return the coordinates of the ThreadingSequence to their original value.
  transform_coordinates();

}

ModelObjectsTemp StructureElementMover::do_get_inputs() const {
  ModelObjectsTemp ret(2);
  ret[0] = get_model()->get_particle(pi_);
  ret[0] = get_model()->get_particle(s_hier_pi_);
  //ModelObjectsTemp ret = {inp};
  return ret;
}

IMPTHREADING_END_NAMESPACE
