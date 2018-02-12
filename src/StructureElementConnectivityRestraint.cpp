/**
 *  \file example/ExampleRestraint.cpp
 *  \brief Restrain a list of particle pairs.
 *
 *  Copyright 2007-2017 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/threading/StructureElementConnectivityRestraint.h>
#include <IMP/core/XYZ.h>
#include <IMP/algebra/Vector3D.h>
#include <IMP/threading/StructureElement.h>
#include <IMP/algebra/VectorD.h>
IMPTHREADING_BEGIN_NAMESPACE

StructureElementConnectivityRestraint::StructureElementConnectivityRestraint(Model *m, UnaryFunction *score_func,                     
                    ParticleIndex a,
                    ParticleIndex b,
                    double dpr,
                    std::string name)
    : Restraint(m, "SEConnectivityRestraint%1%"), score_func_(score_func), a_(a), b_(b), dpr_(dpr) {}

/* Apply the pair score to each particle pair listed in the container.
 */


double StructureElementConnectivityRestraint::unprotected_evaluate(DerivativeAccumulator *accum)
    const {
  //IMP_CHECK_OBJECT(a_.get());
  //IMP_CHECK_OBJECT(b_.get());
  IMP_CHECK_OBJECT(score_func_);
  IMP_OBJECT_LOG;

  // Get the coordinates of the first and last residues for a_ and b_
  algebra::Vector3D a_coords = threading::StructureElement(get_model(), a_).get_coordinates().back();
  algebra::Vector3D b_coords = threading::StructureElement(get_model(), b_).get_coordinates().front();
  
  // calculate distance between termini
  float model_distance = IMP::algebra::get_distance(a_coords, b_coords);

  // Get number of residues between last residue of SE a_ and SE b_
  // calculate the # of residues (should have a check to make sure not negative)
  int residues = threading::StructureElement(get_model(), b_).get_first_residue_number() - threading::StructureElement(get_model(), a_).get_last_residue_number();
  //std::cout << "RESIS:" << threading::StructureElement(get_model(), b_).get_first_residue_number() << " " << threading::StructureElement(get_model(), a_).get_last_residue_number() << std::endl;
  float theor_max_dist = residues * dpr_;

  double distance = model_distance - theor_max_dist;

  double score;
  if (residues <= 1) {
    score = 10000;
  } else if (distance <= 0) {
    score = 0;
  } else {
    score = score_func_->evaluate(distance);
  }
  return score;
}

/* Return all particles whose attributes are read by the restraints. To
   do this, ask the pair score what particles it uses.*/
ModelObjectsTemp StructureElementConnectivityRestraint::do_get_inputs() const {
  return ModelObjectsTemp(1, get_model()->get_particle(a_));
}

IMPTHREADING_END_NAMESPACE
