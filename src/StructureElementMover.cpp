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
  threading::StructureElement se(get_model(), pi_);
  atom::Hierarchy h(get_model(), s_hier_pi_);
  
  // 1) get the original key values
  orig_key_values_ = se.get_all_key_values();
  int seq_len = atom::Hierarchy(get_model(), s_hier_pi_).get_children().size();
  int r;
  
  //  Zero out the coordinates in the Sequence hierarchy
  zero_coordinates();
  // Propose move to start_res
  if (se.get_start_res_is_optimized()) {
    r = statistics::internal::random_int(2) * 2 - 1;
    int new_start_res = orig_key_values_[0] - r;
      //std::cout << "new_start_res + length : " << (new_start_res + orig_key_values_[2]) << std::endl;
    if (new_start_res > 0 && (new_start_res + orig_key_values_[2]) <= seq_len-2) {
        // new_start_res = orig_key_values_[0];
        // std::cout << "new_start_res: " << new_start_res << std::endl;
        se.set_start_res_key(new_start_res);
    };
  };

  // Propose move to polarity
  if (se.get_polarity_is_optimized() && rand() % 100 < pct_flip_) {
      //std::cout << "tried polarity change and flipped it"<< std::endl;
    se.flip_polarity_key();
  };

  // Propose move to length key
  if (se.get_length_is_optimized() && rand() % 100 < pct_length_) {
    r = statistics::internal::random_int(2) * 2 - 1;
    int new_length = orig_key_values_[2] - r;
    if ( new_length >= 1 && (new_length + se.get_offset()) <= se.get_number_of_coordinates() && (new_length + se.get_start_res()) <= seq_len-2) {
        //std::cout << "proposed new length "<< new_length << " current start res " << se.get_start_res() << " where seq_len is "  << seq_len << std::endl;
        se.set_length_key(new_length);
    }
      //std::cout << "changed total length "<< (new_length + se.get_start_res()) << " where seq_len is "  << seq_len << std::endl;
  };

  // Propose move to offset
  if (se.get_offset_is_optimized() && rand() % 100 < pct_offset_) {
    r = statistics::internal::random_int(2) * 2 - 1;
    int new_offset = orig_key_values_[3] - r;
    if ( new_offset >= 0 && new_offset <= se.get_max_offset() ) {
        //std::cout << "offset changed and set to less than max_offset"<< std::endl;
        se.set_offset_key(new_offset);
    }
  };
  // With the new key values, move the SE coordiantes to the sequence hierarchy

  transform_coordinates();
  return IMP::core::MonteCarloMoverResult(h.get_children(), 1.0);
}


void StructureElementMover::zero_coordinates(){
  // Zero out the coordinates
  StructureElement se(get_model(), pi_);
  atom::Hierarchy h(get_model(), s_hier_pi_);
  atom::Selection sel(h);
  sel.set_residue_indexes(se.get_resindex_list());
  sel.set_chain(se.get_chain());
  
  ParticlesTemp oldsel = sel.get_selected_particles();
  for (unsigned int i=0; i<sel.get_selected_particles().size(); i++) {
    core::XYZ coord(oldsel[i]);
    coord.set_coordinates(algebra::Vector3D(0,0,0));
    //Set this as a flag for evaluation or not
    coord.set_coordinates_are_optimized(false);
  }; 

}

void StructureElementMover::transform_coordinates(){
  StructureElement se(get_model(), pi_);

  // If the SR has been changed to zero, then zero out these coordinates
  if (se.get_start_res()==0){
    return;
  }

  // We want the sequence root hierarchy
  atom::Hierarchy h(get_model(), s_hier_pi_);

  // Create a selection to find the residues we wish to 
  // add these coordinates to.
  atom::Selection sel(h);

  sel.set_residue_indexes(se.get_resindex_list());
  sel.set_chain(se.get_chain());

  // Select these particles
  ParticlesTemp new_resis = sel.get_selected_particles();

  algebra::Vector3Ds new_coords = se.get_coordinates();

  // Flip the order of newresis if polarity is -1
  if(se.get_polarity()==-1)
   std::reverse(new_resis.begin(), new_resis.end());

  // For polarity, assign forward, for polarity = -1, assign backwards


  // TO AVOID SEG FAULT Dibyendu replaced se.get_resindex_list() to new_resis.size()
  //for (unsigned int i=0; i<se.get_resindex_list().size(); i++) {
    for (unsigned int i=0; i<new_resis.size(); i++) {

    core::XYZ(new_resis[i]).set_coordinates(new_coords[i]);

    // Set this as a flag for evaluation or not
    core::XYZ(new_resis[i]).set_coordinates_are_optimized(true);
  }; 
}

void StructureElementMover::do_reject() {
  IMP_OBJECT_LOG;
  
  StructureElement se(get_model(), pi_);
  
  // Reset newly assigned coordinates to zero
  zero_coordinates();
  
  //Return StructureElement key values to their original position
  se.set_all_keys(orig_key_values_);

  //Return the coordinates of the sequence hierarchy  to their original value.
  transform_coordinates();
}

ModelObjectsTemp StructureElementMover::do_get_inputs() const {
  ModelObjectsTemp ret(2);
  ret[0] = get_model()->get_particle(pi_);
  ret[1] = get_model()->get_particle(s_hier_pi_);
  return ret;
}

IMPTHREADING_END_NAMESPACE
