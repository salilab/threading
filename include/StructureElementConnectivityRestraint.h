/**
 *  \file IMP/threading/StructureElementConnectivityRestraint.h
 *  \brief A restraint on a pair of particles
 *  \that is either evaluated using a UnaryFunction
 *  \if both particles have XYZ decorators with 
 *  \get_is_optimized() = True. Otherwise a constant
 *  \value is returned
 *
 *  Copyright 2007-2017 IMP Inventors. All rights reserved.
 *
 */

#ifndef IMPTHREADING_STRUCTURE_ELEMENT_CONECTIVITY_RESTRAINT_H
#define IMPTHREADING_STRUCTURE_ELEMENT_CONECTIVITY_RESTRAINT_H

#include <IMP/threading/threading_config.h>
#include <IMP/Restraint.h>
#include <IMP/UnaryFunction.h>

IMPTHREADING_BEGIN_NAMESPACE

//! Restraint on the sequence distance between two StructureElements
// Particles a and b refer to structural elements
// Particle a always stays constant.
// Particle b can be updated.

class IMPTHREADINGEXPORT StructureElementConnectivityRestraint : public Restraint {
  ParticleIndex a_;
  ParticleIndex b_;
  int n_sds_; // The distance per residue
  IMP::PointerMember<UnaryFunction> score_func_;
  
 public:
  //! Create the restraint.
  /** Restraints should store the particles they are to act on,
      preferably in a Singleton or PairContainer as appropriate.
   */
  StructureElementConnectivityRestraint(Model *m, UnaryFunction *score_func,                     
                    ParticleIndex a,
                    ParticleIndex b,
                    float n_sds,
                    std::string name = "EndToEndRestraint %1%");

  double unprotected_evaluate(DerivativeAccumulator *accum) const
      IMP_OVERRIDE;

  void assign_particles(ParticleIndex a, ParticleIndex b) {
    a_ = a;
    b_ = b;
  };

  //void set_distance_per_residue(float dpr) {
  //  dpr_ = dpr;
  //};

  //Currently only works for helices

  ParticleIndex get_pia() { return a_; };
  ParticleIndex get_pib() { return b_; };

  float get_mean_distance_per_residue() const;

  float get_sd_distance_per_residue() const;

  int get_number_of_residues() const;

  float get_mean_distance() const;

  float get_max_distance() const;

  float get_model_distance() const;



  ModelObjectsTemp do_get_inputs() const;
  IMP_OBJECT_METHODS(StructureElementConnectivityRestraint);
};

IMPTHREADING_END_NAMESPACE

#endif /* IMPTHREADING_STRUCTURE_ELEMENT_CONNECTIVITY_RESTRAINT_H */
