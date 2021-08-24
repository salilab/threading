/**
 *  \file SSEThreadMover.cpp  \brief A modifier to a Fragment object that perturbs 
 *     the start residue by an integer value.
 *
 *  Copyright 2007-2020 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/threading/SSEThreadMover.h>
#include <IMP/threading/StructureElement.h>

IMPTHREADING_BEGIN_NAMESPACE

namespace {
std::string get_structure_element_mover_name() {
  return "SSEThreadMover-BASE";
}
}

void SSEThreadMover::initialize(Model *m, ParticleIndexes se_pis, ParticleIndex s_hier_pi) {
  se_pis_ = se_pis;
  s_hier_pi_ = s_hier_pi;
}

SSEThreadMover::SSEThreadMover(Model *m, ParticleIndexes se_pis, ParticleIndex s_hier_pi)
    : IMP::core::MonteCarloMover(m, get_structure_element_mover_name()) {
  initialize(m, se_pis, s_hier_pi);
}

IMP::core::MonteCarloMoverResult SSEThreadMover::do_propose() {
      return IMP::core::MonteCarloMoverResult(atom::Hierarchy(get_model(), s_hier_pi_).get_children(),1.0);
}


void SSEThreadMover::do_reject() {
  IMP_OBJECT_LOG;
}

void SSEThreadMover::transform_coordinates(ParticleIndex se_pi){
  IMP_USAGE_CHECK(StructureElement().get_is_setup(get_model(), se_pi), "Particle is not set up as a StructureElement");
  threading::StructureElement se(get_model(), se_pi);

  // If the SR has been changed to zero, then zero out these coordinates
  if (se.get_start_res()==0){
    return;
  }
  
  // We want the sequence root hierarchy
  atom::Hierarchy h(get_model(), s_hier_pi_);

  // Create a selection to find the residues we wish to 
  // add these coordinates to.
  // get the chain id of the SecondaryElement
 std::string se_chain = se.get_chain();

 atom::Selection sel(h);
 
  sel.set_residue_indexes(se.get_resindex_list());
  // assign the change to the correct chain 
  sel.set_chain_id(se_chain);
  // Select these particles
  ParticlesTemp new_resis = sel.get_selected_particles();
  
  algebra::Vector3Ds new_coords = se.get_coordinates();

  // Flip the order of newresis if polarity is -1
  if(se.get_polarity()==-1)
   std::reverse(new_resis.begin(), new_resis.end());

  // For polarity, assign forward, for polarity = -1, assign backwards
  for (unsigned int i=0; i<se.get_resindex_list().size(); i++) {
    core::XYZ(new_resis[i]).set_coordinates(new_coords[i]);
    // Set this as a flag for evaluation or not
    core::XYZ(new_resis[i]).set_coordinates_are_optimized(true);    
  }; 
}


void SSEThreadMover::zero_coordinates(ParticleIndex se_pi){
  // Zero out the coordinates from this StructureElement
  IMP_USAGE_CHECK(StructureElement().get_is_setup(get_model(), se_pi), "Particle is not set up as a StructureElement");
  threading::StructureElement se(get_model(), se_pi);
  atom::Hierarchy h(get_model(), s_hier_pi_);
  atom::Selection sel(h);
  std::string se_chain = se.get_chain();
  
  sel.set_residue_indexes(se.get_resindex_list());

  // assign the change to the correct chain
  sel.set_chain_id(se_chain);

  ParticlesTemp oldsel = sel.get_selected_particles();
  for (unsigned int i=0; i<oldsel.size();i++){
    core::XYZ coord(oldsel[i]);
    
    coord.set_coordinates(algebra::Vector3D(0,0,0));
    //Set this as a flag for evaluation or not
    coord.set_coordinates_are_optimized(false);
  };
}

ModelObjectsTemp SSEThreadMover::do_get_inputs() const {
  ModelObjectsTemp ret;
 
  ret.push_back(get_model()->get_particle(s_hier_pi_));
  for (unsigned int i=0; i<se_pis_.size(); i++){
    ret.push_back(get_model()->get_particle(se_pis_[i]));
  } 
  
  return ret;
}

IMPTHREADING_END_NAMESPACE
