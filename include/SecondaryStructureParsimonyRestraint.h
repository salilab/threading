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

#ifndef IMPTHREADING_SECONDARY_STRUCTURE_PARSIMONY_RESTRAINT_H
#define IMPTHREADING_SECONDARY_STRUCTURE_PARSIMONY_RESTRAINT_H

#include <IMP/threading/threading_config.h>
#include <IMP/Restraint.h>
#include <IMP/base_types.h>

IMPTHREADING_BEGIN_NAMESPACE
/** 
*/
class IMPTHREADINGEXPORT SecondaryStructureParsimonyRestraint : public Restraint {
  ParticleIndexes ses_; // The structural elements
  ParticleIndex b_; // The sequence residue
  Float c_;
  Float helix_prob_;
  Float coil_prob_;
  Float sheet_prob_;
 public:
  /* ! Compare Secondary Structure Residue values for the 
   */
  SecondaryStructureParsimonyRestraint(Model *m,                    
                    ParticleIndexes ses,
                    ParticleIndex b,
                    float c,
                    std::string name = "SecondaryStructureParsimonyRestraint %1%");
  
  //ParticleIndexes ses_; // The structural elements
  //ParticleIndex b_; // The sequence residue
  //Float c_;

  double unprotected_evaluate(DerivativeAccumulator *accum) const
      IMP_OVERRIDE;
  ModelObjectsTemp do_get_inputs() const;
  IMP_OBJECT_METHODS(SecondaryStructureParsimonyRestraint);

void set_baseline_ss_probabilities(Float a, Float b, Float c) {
  float norm = a + b + c;

  helix_prob_ = a / norm;
  coil_prob_ = c / norm;
  sheet_prob_ = b / norm;
};

Floats get_baseline_probabilities() {
  Floats probs;
  probs.push_back(helix_prob_);
  probs.push_back(sheet_prob_);
  probs.push_back(coil_prob_);
  return probs;
};

void set_weight(Float c) {
  c_ = c;
};

void set_structural_elements(ParticleIndexes a){
  ses_ = a;
};

void set_sequence_residue(ParticleIndex b){
  b_ = b;
};

/*
std::vector<std::vector<float> > get_se_sequence_probabilities(){
  return ses_ss_probs;
};
*/
};
IMPTHREADING_END_NAMESPACE

#endif /* IMPTHREADING_SECONDARY_STRUCTURE_PARSIMONY_RESTRAINT_H */
