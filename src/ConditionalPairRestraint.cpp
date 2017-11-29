/**
 *  \file example/ExampleRestraint.cpp
 *  \brief Restrain a list of particle pairs.
 *
 *  Copyright 2007-2017 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/threading/ConditionalPairRestraint.h>
#include <IMP/core/XYZ.h>
#include <IMP/algebra/Vector3D.h>

IMPTHREADING_BEGIN_NAMESPACE

ConditionalPairRestraint::ConditionalPairRestraint(Model *m, UnaryFunction *score_func,                     
                    ParticleIndex a,
                    ParticleIndex b,
                    double c,
                    std::string name)
    : Restraint(m, "ConditionalPairRestraint%1%"), score_func_(score_func), a_(a), b_(b), c_(c) {}

/* Apply the pair score to each particle pair listed in the container.
 */


double ConditionalPairRestraint::unprotected_evaluate(DerivativeAccumulator *accum)
    const {
  //IMP_CHECK_OBJECT(a_.get());
  //IMP_CHECK_OBJECT(b_.get());
  IMP_CHECK_OBJECT(score_func_);
  IMP_OBJECT_LOG;
  bool o1 = core::XYZ(get_model(), a_).get_coordinates_are_optimized();
  bool o2 = core::XYZ(get_model(), b_).get_coordinates_are_optimized();
  double score;
  if (o1==false or o2==false) {
    score = c_;
  } else {
    score = score_func_->evaluate(core::get_distance(core::XYZ(get_model(), a_), core::XYZ(get_model(), b_)));
  }
  return score;
}

/* Return all particles whose attributes are read by the restraints. To
   do this, ask the pair score what particles it uses.*/
ModelObjectsTemp ConditionalPairRestraint::do_get_inputs() const {
  return ModelObjectsTemp(1, get_model()->get_particle(a_));
}

IMPTHREADING_END_NAMESPACE
