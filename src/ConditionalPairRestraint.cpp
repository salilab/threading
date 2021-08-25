/**
 *  \file example/ExampleRestraint.cpp
 *  \brief Restrain a list of particle pairs.
 *
 *  Copyright 2007-2017 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/threading/ConditionalPairRestraint.h>
#include <IMP/core/XYZ.h>
#include <IMP/atom/Residue.h>
#include <IMP/atom/Hierarchy.h>
#include <IMP/atom/Selection.h>
#include <IMP/algebra/Vector3D.h>

IMPTHREADING_BEGIN_NAMESPACE

ConditionalPairRestraint::ConditionalPairRestraint(Model *m, UnaryFunction *score_func,                     
                    ParticleIndex a,
                    ParticleIndex b,
                    double dpr,
                    std::string name)
    : Restraint(m, name), score_func_(score_func), a_(a), b_(b), dpr_(dpr) {}

/* Apply the pair score to each particle pair listed in the container.
 */

Particles ConditionalPairRestraint::get_closest_built_residue_particles(ParticleIndex pi) const
{
  ParticlesTemp endpoints;
  atom::Residue r = atom::Residue(get_model(), pi);
  int ri = r.get_index();
  std::cout << "Finding endpoints for: " << ri << "| "; 
  // Get the chain of this residue
  atom::Hierarchy res(get_model(), pi);
  atom::Hierarchy seq_chain = res.get_parent();
  // Go forward in sequence space looking for a coordinate that is optimized
  for (unsigned int i=1; i<100; i++){
    atom::Selection sel(seq_chain);
    sel.set_residue_index(ri + i);
    ParticlesTemp nextres = sel.get_selected_particles();
    // If you hit the end of the chain, then just leave
    if (nextres.size() == 0){ break;}
    if (core::XYZ(nextres[0]).get_coordinates_are_optimized()){
      endpoints.push_back(nextres[0]);
      std::cout << " " << ri+i;
      break;
    } 
  }

  // Go backward in sequence space looking for a coordinate that is optimized
  for (int i=-1; i>-100; i--){
    atom::Selection sel(seq_chain);
    sel.set_residue_index(ri + i);
    ParticlesTemp nextres = sel.get_selected_particles();
    // If you hit the end of the chain, then just leave
    if (nextres.size() == 0){ std::cout << std::endl; break;}
    if (core::XYZ(nextres[0]).get_coordinates_are_optimized()){
      endpoints.push_back(nextres[0]);
      std::cout << " " << ri+i << std::endl;   
      break;
    } 
  }
  return endpoints;
}

double ConditionalPairRestraint::unprotected_evaluate(DerivativeAccumulator *accum)
    const {
  //IMP_CHECK_OBJECT(a_.get());
  //IMP_CHECK_OBJECT(b_.get());
  IMP_CHECK_OBJECT(score_func_);
  IMP_OBJECT_LOG;
  // Find if either or both residues have coordinates
  bool a_opt = core::XYZ(get_model(), a_).get_coordinates_are_optimized();
  bool b_opt = core::XYZ(get_model(), b_).get_coordinates_are_optimized();
  double score;

  // Determine XL evaluation distance in all circumstances...this could be written much better
  float distance = 1000;
  // If both residues are structured
  if (a_opt==true and b_opt==true) {
    distance = core::get_distance(core::XYZ(get_model(), a_), core::XYZ(get_model(), b_));

  // Only a is structured
  } else if (a_opt==true and b_opt==false) {
    core::XYZ a_coord(get_model(), a_);
    // Find particles of closest built residues to b_
    Particles ps = get_closest_built_residue_particles(b_);

    for (unsigned int i=0; i<ps.size(); i++) {
      // get residue index
      int r_index = atom::Residue(ps[i]).get_index();  
      float r_dist = abs(r_index - atom::Residue(get_model(), b_).get_index()) * dpr_;
      float new_dist = core::get_distance(core::XYZ(get_model(), a_), core::XYZ(ps[i])) - r_dist;
      if (new_dist < distance){ distance = new_dist; }
    }

  // Only b is structured
  } else if (a_opt==false and b_opt==true) {
    core::XYZ b_coord(get_model(), b_);
    // Find particles of closest built residues to b_
    Particles ps = get_closest_built_residue_particles(a_);

    for (unsigned int i=0; i<ps.size(); i++) {
      // get residue index
      int r_index = atom::Residue(ps[i]).get_index();  
      float r_dist = abs(r_index - atom::Residue(get_model(), a_).get_index()) * dpr_;
      float new_dist = core::get_distance(core::XYZ(get_model(), b_), core::XYZ(ps[i])) - r_dist;
      if (new_dist < distance){ distance = new_dist; }
    }    

  // both are unstructured
  } else {
    Particles ps_a = get_closest_built_residue_particles(a_);
    Particles ps_b = get_closest_built_residue_particles(b_);

    for (unsigned int i=0; i<ps_a.size(); i++) {
      int r_a_index = atom::Residue(ps_a[i]).get_index();  
      float r_a_dist = abs(r_a_index - atom::Residue(get_model(), a_).get_index()) * dpr_;

      for (unsigned int j=0; j<ps_b.size(); j++) {
        int r_b_index = atom::Residue(ps_b[j]).get_index();  
        float r_b_dist = abs(r_b_index - atom::Residue(get_model(), b_).get_index()) * dpr_;

        float new_dist = core::get_distance(core::XYZ(ps_a[i]), core::XYZ(ps_b[j])) - r_a_dist - r_b_dist;
        std::cout << "New dist: " << new_dist << " " << core::get_distance(core::XYZ(ps_a[i]), core::XYZ(ps_b[j])) << " Indices:" << r_a_index << " " << r_b_index << " | " << r_a_dist << " " << r_b_dist << std::endl; 
        if (new_dist < distance){ distance = new_dist; }
      }
    }
  }

  std::cout << "Distance: " << distance << std::endl;

  score = score_func_->evaluate(distance);
  return score;
}

/* Return all particles whose attributes are read by the restraints. To
   do this, ask the pair score what particles it uses.*/
ModelObjectsTemp ConditionalPairRestraint::do_get_inputs() const {
  return ModelObjectsTemp(1, get_model()->get_particle(a_));
}

IMPTHREADING_END_NAMESPACE
