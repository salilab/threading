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
                    float n_sds,
                    std::string name)
    : Restraint(m, "SEConnectivityRestraint%1%"), score_func_(score_func), a_(a), b_(b), n_sds_(n_sds) {}

/* Apply the pair score to each particle pair listed in the container.
 */

static float nhh_means_[30] = {3.81, 3.036, 2.836, 2.511, 2.275, 2.178, 2.026, 1.876, 1.835, 1.669, 1.658, 1.666, 1.625, 1.53,
                1.445, 1.374, 1.292, 1.212, 1.164, 1.133, 1.049, 1.043, 1.074, 0.977, 0.965, 0.938, 0.868, 0.824, 0.805,
                0.788};

static float nhh_sds_[30] = {0.027, 0.284, 0.397, 0.441, 0.483, 0.499, 0.504, 0.537, 0.534, 0.538, 0.545, 0.507, 0.494, 0.468,
              0.447, 0.428, 0.439, 0.415, 0.432, 0.392, 0.382, 0.38, 0.401, 0.381, 0.38, 0.317, 0.328, 0.304, 0.318,
              0.273};

float StructureElementConnectivityRestraint::get_sd_distance_per_residue() const {
  int nres = get_number_of_residues();
  float sd;
  if (nres > 29) { sd = nhh_sds_[29];
  } else { sd = nhh_sds_[nres];};
  return sd;
};

int StructureElementConnectivityRestraint::get_number_of_residues()  const{
  int residues = threading::StructureElement(get_model(), b_).get_first_residue_number() - threading::StructureElement(get_model(), a_).get_last_residue_number();
  return residues;
};

float StructureElementConnectivityRestraint::get_mean_distance() const{
  int residues = get_number_of_residues();
  return residues * get_mean_distance_per_residue();
};

float StructureElementConnectivityRestraint::get_max_distance() const {
  float meandist = get_mean_distance();
  float sddist = get_sd_distance_per_residue();
  //std::cout << " > " << meandist << " " << sddist << std::endl;
  return meandist + sddist * n_sds_;
};

float StructureElementConnectivityRestraint::get_model_distance() const {
  // Get the coordinates of the first and last residues for a_ and b_
  algebra::Vector3D a_coords = threading::StructureElement(get_model(), a_).get_coordinates().back();
  algebra::Vector3D b_coords = threading::StructureElement(get_model(), b_).get_coordinates().front();
  float model_distance = IMP::algebra::get_distance(a_coords, b_coords);
  return model_distance;
};


float StructureElementConnectivityRestraint::get_mean_distance_per_residue() const {
  int nres = get_number_of_residues();
  float mean;
  if (nres > 29) { mean = nhh_means_[29];
  } else { mean = nhh_means_[nres];};
  return mean;
};

double StructureElementConnectivityRestraint::unprotected_evaluate(DerivativeAccumulator *accum)
    const {
  //IMP_CHECK_OBJECT(a_.get());
  //IMP_CHECK_OBJECT(b_.get());
  IMP_CHECK_OBJECT(score_func_);
  IMP_OBJECT_LOG;

  // Get the coordinates of the first and last residues for a_ and b_
  //algebra::Vector3D a_coords = threading::StructureElement(get_model(), a_).get_coordinates().back();
  //algebra::Vector3D b_coords = threading::StructureElement(get_model(), b_).get_coordinates().front();
  
  // calculate distance between termini
  //float model_distance = get_model_distance(); //IMP::algebra::get_distance(a_coords, b_coords);

  // Get number of residues between last residue of SE a_ and SE b_
  // calculate the # of residues (should have a check to make sure not negative)
  //int residues = get_number_of_residues(); //threading::StructureElement(get_model(), b_).get_first_residue_number() - threading::StructureElement(get_model(), a_).get_last_residue_number();
  //std::cout << "RESIS:" << threading::StructureElement(get_model(), b_).get_first_residue_number() << " " << threading::StructureElement(get_model(), a_).get_last_residue_number() << std::endl;
  //float theor_max_dist = get_max_distance();

  double distance = get_model_distance() - get_max_distance();

  double score;
  if (get_number_of_residues() <= 1) {
    score = 10000;
  } else if (distance <= 0) {
    score = 0;
  } else {
    score = score_func_->evaluate(distance);
  }
  //std::cout << "Connectivity: " << get_name() << score << " " << get_number_of_residues() << " | " << threading::StructureElement(get_model(), b_).get_first_residue_number() << " " << threading::StructureElement(get_model(), a_).get_last_residue_number() << std::endl;
  return score;
}


/* Return all particles whose attributes are read by the restraints. To
   do this, ask the pair score what particles it uses.*/
ModelObjectsTemp StructureElementConnectivityRestraint::do_get_inputs() const {
  return ModelObjectsTemp(1, get_model()->get_particle(a_));
}

IMPTHREADING_END_NAMESPACE
