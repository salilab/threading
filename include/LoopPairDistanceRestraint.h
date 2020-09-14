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

#ifndef IMPTHREADING_LOOP_PAIR_DISTANCE_RESTRAINT_H
#define IMPTHREADING_LOOP_PAIR_DISTANCE_RESTRAINT_H

#include <IMP/threading/threading_config.h>
#include <IMP/Restraint.h>
#include <IMP/base_types.h>
#include <IMP/exception.h>
#include <IMP/UnaryFunction.h>
#include <boost/unordered_map.hpp>

IMPTHREADING_BEGIN_NAMESPACE

//! Restrain a distance between two residues that may or may not have coordinates
/** When evaluating a crosslink to a residue with no coordinates, in disordered loop, we need to infer where it is.

    For each unstructured residue, the coordinates of the nearest structured residue(s) are used to 
    evaluate the XL distance.  A tolerance is then subtracted
    from this distance to account for the loop length.  This length calculated based on the number
    of residues. In most cases, there is more than one endpoint, so all permutations are evaluated and the
    minimum distance returned.

    ##
    The average and SD of the distance per residue
    ##

    This distance is then evaluated as "normal", given the input scoring function. (upperharmonic, sigmoidal, etc...)

    For example:

    An XL is to be evaluated between residues A and B.  B is unstructured.  It is residue 6 in a loop of 
    15 residue (5 residues between B and its N-terminal structured residue (B_NTS) and 8 residues between it and 
    the C-terminal structured residue (B_CTS)). We decide to allow a tolerance of +1SD from the loop length mean (n_sds=1)

    Say the distance between A and B_NTS is 45.0 Ang. The loop length is 6 * (2.178 + 1 * 0.499) = 16.1 Ang
    A -- B_NTS = 28.9 Angstroms  (45.0 - 16.1)

    Say the distance between A and B_CTS is 50.0 Ang. The loop length is 8 * (1.876 + 1 * 0.537) = 19.0 Ang
    A -- B_CTS = 31.0 Angstroms  (50.0 - 19.0)

    The restraint is evaluated at the smallest distance, 28.9 angstroms, using whatever scoring function you want

*/
class IMPTHREADINGEXPORT LoopPairDistanceRestraint : public Restraint {
  ParticleIndex a_;
  ParticleIndex b_;
  ParticleIndexes as_;
  ParticleIndexes bs_;

  bool amb_rest_ = false;

  double n_sds_;
  double distance_threshold_;
  // Mapping chain, residue to endpoints
  mutable boost::unordered_map<ParticleIndex, ParticlesTemp> map_endpoints_;
  
  IMP::PointerMember<UnaryFunction> score_func_;

 public:
  //! Create the restraint.
  /** Restraints should store the particles they are to act on,
      preferably in a Singleton or PairContainer as appropriate.
   */
  LoopPairDistanceRestraint(Model *m, UnaryFunction *score_func,                     
                    ParticleIndex a,
                    ParticleIndex b,
                    double n_sds,
                    std::string name = "LoopPairDistanceRestraint %1%");

  LoopPairDistanceRestraint(Model *m, UnaryFunction *score_func,                     
                    ParticleIndexes as,
                    ParticleIndexes bs,
                    double n_sds,
                    double distance_threshold,        
                    std::string name = "LoopPairDistanceRestraint %1%");

  double calc_sphere_cap_distance(float R, float d, float alpha) const;

  Particles get_closest_built_residue_particles(ParticleIndex pi) const;

  ParticleIndexes get_sequence_residue_particles() const;

  algebra::Vector3D get_sphere_cap_center(Particle* P0, Particle* P1, float R0, float R1) const;

  double get_loop_distance(int a, int b) const;

  double get_pair_distance(ParticleIndex pa, ParticleIndex pb) const;

  algebra::Vector3D get_cap_center(ParticleIndex p_sel) const;

  double unprotected_evaluate(DerivativeAccumulator *accum) const
      IMP_OVERRIDE;

  ModelObjectsTemp do_get_inputs() const;
  IMP_OBJECT_METHODS(LoopPairDistanceRestraint);
};

IMPTHREADING_END_NAMESPACE

#endif /* IMPTHREADING_LOOP_PAIR_DISTANCE_RESTRAINT_H */
