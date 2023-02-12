/**
 *  \file IMP/threading/ConditionalPairRestraint.h
 *  \brief A restraint on a pair of particles
 *  \that is either evaluated using a UnaryFunction
 *  \if both particles have XYZ decorators with 
 *  \get_is_optimized() = True. Otherwise a constant
 *  \value is returned
 *
 *  Copyright 2007-2017 IMP Inventors. All rights reserved.
 *
 */

#ifndef IMPTHREADING_CONDITIONAL_PAIR_RESTRAINT_H
#define IMPTHREADING_CONDITIONAL_PAIR_RESTRAINT_H

#include <IMP/threading/threading_config.h>
#include <IMP/Restraint.h>
#include <IMP/base_types.h>
#include <IMP/exception.h>
#include <IMP/UnaryFunction.h>

IMPTHREADING_BEGIN_NAMESPACE

//! Restrain a particle to be in the x,y plane
/** \note Be sure to check out the swig wrapper file and how it
    wraps this class.

    The source code is as follows:
    \include ExampleRestraint.h
    \include ExampleRestraint.cpp
*/
class IMPTHREADINGEXPORT ConditionalPairRestraint : public Restraint {
  IMP::PointerMember<UnaryFunction> score_func_;
  ParticleIndex a_;
  ParticleIndex b_;
  double dpr_;

 public:
  //! Create the restraint.
  /** Restraints should store the particles they are to act on,
      preferably in a Singleton or PairContainer as appropriate.
   */
  ConditionalPairRestraint(Model *m, UnaryFunction *score_func,                     
                    ParticleIndex a,
                    ParticleIndex b,
                    double dpr,
                    std::string name = "ConditionalPairRestraint %1%");


  Particles get_closest_built_residue_particles(ParticleIndex pi) const;
  double unprotected_evaluate(DerivativeAccumulator *accum) const override;

  ModelObjectsTemp do_get_inputs() const override;
  IMP_OBJECT_METHODS(ConditionalPairRestraint);
};

IMPTHREADING_END_NAMESPACE

#endif /* IMPTHREADING_CONDITIONAL_PAIR_RESTRAINT_H */
