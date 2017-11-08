/**
 *  \file IMP/core/FragmentMover.h
 *  \brief A modifier which moves the sequence in a fragment.
 *
 *  Copyright 2007-2016 IMP Inventors. All rights reserved.
 *
 */

#ifndef IMPTHREADING_SEQUENCE_BASE_MOVER_H
#define IMPTHREADING_SEQUENCE_BASE_MOVER_H

#include <IMP/threading/threading_config.h>

#include <IMP/core/MonteCarloMover.h>

#include <IMP/base_types.h>
#include <IMP/exception.h>
#include <IMP/atom/Hierarchy.h>
//#include <IMP/atom/Residue.h>
//#include <IMP/core/XYZ.h>
IMPTHREADING_BEGIN_NAMESPACE

//! A mover for modifying the assignment of residues to a sequ
//! perturbing them within a ball.
/** The variables are perturbed within a ball of the
    given radius.
    \see MonteCarlo
 */
class IMPTHREADINGEXPORT SequenceBaseMover : public IMP::core::MonteCarloMover {
  ParticleIndex fpi_;
  IntKey start_res_;
  IntKey polarity_;
  Int original_;

  //IMP::core::XYZs xyzs_;
  //IMP::atom::ResidueTypes rtypes_;

  void initialize(ParticleIndex fpi);

 public:
  //! Move specified variables of particle pi within a ball of specified radius
  /** Construct a mover that in each move, perturbs the specified
      variables (attributes) of particle pi in model m, within a ball
      of specified radius, whose dimensionality is the total number of
      attributes
  */
  //FragmentBaseMover(Model *m, ParticleIndexes pis, const IntKey &var,
  //          int max_trans);
 SequenceBaseMover(Model *m, ParticleIndex fpi);


 protected:
  virtual ModelObjectsTemp do_get_inputs() const IMP_OVERRIDE;

  //! Move particle attributes within a ball, as specified in constructor
  virtual IMP::core::MonteCarloMoverResult do_propose() IMP_OVERRIDE;

  //! restore original attributes from before do_propose
  virtual void do_reject() IMP_OVERRIDE;
  IMP_OBJECT_METHODS(SequenceBaseMover);
};

IMPTHREADING_END_NAMESPACE

#endif /* IMPTHREADING_SEQUENCE_BASE_MOVER_H */
