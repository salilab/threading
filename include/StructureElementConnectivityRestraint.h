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
  double dpr_; // The distance per residue
  IMP::PointerMember<UnaryFunction> score_func_;
  
 public:
  //! Create the restraint.
  /** Restraints should store the particles they are to act on,
      preferably in a Singleton or PairContainer as appropriate.
   */
  StructureElementConnectivityRestraint(Model *m, UnaryFunction *score_func,                     
                    ParticleIndex a,
                    ParticleIndex b,
                    double dpr,
                    std::string name = "EndToEndRestraint %1%");

  double unprotected_evaluate(DerivativeAccumulator *accum) const
      IMP_OVERRIDE;

  void assign_particles(ParticleIndex a, ParticleIndex b) {
    a_ = a;
    b_ = b;
  };

  void set_distance_per_residue(float dpr) {
    dpr_ = dpr;
  };

  int get_number_of_residues() {
    // Get resindex lists for each SE
    int residues = threading::StructureElement(get_model(), b_).get_first_residue_number() - threading::StructureElement(get_model(), a_).get_last_residue_number();
    return residues;
  };

  float get_max_distance() {

    int residues = get_number_of_residues();
    return residues * dpr_;
  };

  float get_model_distance() {
    // Get the coordinates of the first and last residues for a_ and b_
    algebra::Vector3D a_coords = threading::StructureElement(get_model(), a_).get_coordinates().back();
    algebra::Vector3D b_coords = threading::StructureElement(get_model(), b_).get_coordinates().front();
    
    float model_distance = IMP::algebra::get_distance(a_coords, b_coords);

    return model_distance;
  };

  ModelObjectsTemp do_get_inputs() const;
  IMP_OBJECT_METHODS(StructureElementConnectivityRestraint);
};

IMPTHREADING_END_NAMESPACE

#endif /* IMPTHREADING_STRUCTURE_ELEMENT_CONNECTIVITY_RESTRAINT_H */
