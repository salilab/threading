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

/* 
 * Apply the pair score to each particle pair listed in the container.
 * The mean and standard deviation of the distance-per-residue for loops based on number of residues and bounding secondary structure element was determined from the CATH S100 database as described 
 * in Section 2.4.1 in https://doi.org/10.1016/j.pbiomolbio.2019.09.003. 
 * Helix-loop-helix and sheet-loop-sheet values were computed using 1,668,774 and 614,152 Smotif elements, respectively. 
 * Sheet-loop-helix and helix-loop-sheet observations were combined, with a total of 1,416,244 Smotif elements (702,557 sheet-loop-helix and 713,687 helix-loop-sheet) used to compute the values. 
 * The values are taken from  Table 1 in https://doi.org/10.1016/j.pbiomolbio.2019.09.003
 */

static float ll_means_[90] = {3.81, 3.036, 2.836, 2.511, 2.275, 2.178, 2.026, 1.876, 1.835, 1.669, 1.658, 1.666, 1.625, 1.53,
                              1.445, 1.374, 1.292, 1.212, 1.164, 1.133, 1.049, 1.043, 1.074, 0.977, 0.965, 0.938, 0.868, 0.824, 0.805,
                              0.788, 3.81, 3.19, 1.846, 1.607, 1.274, 1.14, 1.139, 1.198, 1.177, 1.115, 1.029, 1.048, 0.935, 0.91, 0.908,
                              0.85, 0.83, 0.852, 0.849, 0.761, 0.722, 0.742, 0.684, 0.677, 0.611, 0.587, 0.596, 0.565, 0.576, 0.532,
                              3.809, 3.137, 2.818, 2.482, 2.154, 1.928, 1.749, 1.67, 1.531, 1.428, 1.377, 1.282, 1.261, 1.203,
                              1.135, 1.045, 1.004, 1.02, 0.977, 0.928, 0.865, 0.834, 0.811, 0.756, 0.761, 0.749, 0.777, 0.74, 0.655,
                              0.648};


static float ll_sds_[90] = {0.027, 0.284, 0.397, 0.441, 0.483, 0.499, 0.504, 0.537, 0.534, 0.538, 0.545, 0.507, 0.494, 0.468,
                            0.447, 0.428, 0.439, 0.415, 0.432, 0.392, 0.382, 0.38, 0.401, 0.381, 0.38, 0.317, 0.328, 0.304, 0.318,
                            0.273, 0.027, 0.313, 0.293, 0.469, 0.419, 0.474, 0.49, 0.505, 0.447, 0.501, 0.475, 0.479, 0.417, 0.451, 0.416,
                            0.373, 0.395, 0.47, 0.418, 0.36, 0.349, 0.359, 0.312, 0.302, 0.281, 0.279, 0.264, 0.259, 0.346, 0.257,
                            0.067, 0.278, 0.361, 0.418, 0.45, 0.448, 0.455, 0.436, 0.452, 0.438, 0.416, 0.407, 0.402, 0.411, 0.405,
                            0.381, 0.378, 0.373, 0.36, 0.372, 0.338, 0.322, 0.308, 0.285, 0.289, 0.296, 0.298, 0.294, 0.286, 0.208};

float StructureElementConnectivityRestraint::get_sd_distance_per_residue() const {
  int nres = get_number_of_residues();
  float sd;
  int sset;
  Floats ses_ss_prob_a = atom::SecondaryStructureResidue(get_model(), a_).get_all_probabilities();
  Floats ses_ss_prob_b = atom::SecondaryStructureResidue(get_model(), b_).get_all_probabilities();
  if (ses_ss_prob_a[0]==1 && ses_ss_prob_b[0]==1) {
        sset = 0;
  } else if (ses_ss_prob_a[1]==1 && ses_ss_prob_b[1]==1){
        sset = 1;
  } else if ((ses_ss_prob_a[1]==0 && ses_ss_prob_b[1]==1) || (ses_ss_prob_a[1]==1 && ses_ss_prob_b[1]==0)){
        sset = 2;
  } else {sset = 0;}

 if (nres > 29) {
      sd = ll_sds_[29+(sset*30)];
  }else { sd = ll_sds_[nres+(sset*30)];}

  return sd;
};

int StructureElementConnectivityRestraint::get_number_of_residues()  const{
  int residues = threading::StructureElement(get_model(), b_).get_first_residue_number() - threading::StructureElement(get_model(), a_).get_last_residue_number();

  return residues;
};

float StructureElementConnectivityRestraint::get_mean_distance() const{
  int residues = get_number_of_residues();
  //std::cout<<"n_residues "<< residues<<" "<<get_mean_distance_per_residue()<<" "<<residues * get_mean_distance_per_residue()<<std::endl;
  
  return residues * get_mean_distance_per_residue();
};

float StructureElementConnectivityRestraint::get_max_distance() const {
  float meandist = get_mean_distance();
  float sddist = get_sd_distance_per_residue();
  return meandist + sddist * n_sds_;
};

float StructureElementConnectivityRestraint::get_model_distance() const {
  // Get the coordinates of the last and first residues for a_ and b_

  float model_distance;
  
  std::string chain_1 = threading::StructureElement(get_model(), a_).get_chain();
  std::string chain_2 = threading::StructureElement(get_model(), b_).get_chain();  
  //std::cout << "chain 1: " << chain_1 << " & chain 2: " << chain_2 << std::endl;
  // Compute and return the distance
  if (chain_1 == chain_2){
 
     algebra::Vector3D a_coords = threading::StructureElement(get_model(), a_).get_coordinates().back();
      //std::cout<<" a_coords "<<a_coords<<std::endl;
     algebra::Vector3D b_coords = threading::StructureElement(get_model(), b_).get_coordinates().back();
     //std::cout<<" b_coords "<<b_coords<<std::endl;
    
     model_distance = IMP::algebra::get_distance(a_coords, b_coords);
     //std::cout<<"model_distance "<<model_distance<<std::endl;
  }
  else {
    model_distance = 0;
  }
  //std::cout<<"CHAINS "<<chain_1<<" "<<chain_2<<std::endl;
  
  return model_distance;
};


float StructureElementConnectivityRestraint::get_mean_distance_per_residue() const {
  int nres = get_number_of_residues();
 
  float mean;

    int sset;
    Floats ses_ss_prob_a = atom::SecondaryStructureResidue(get_model(), a_).get_all_probabilities();
    Floats ses_ss_prob_b = atom::SecondaryStructureResidue(get_model(), b_).get_all_probabilities();
    if (ses_ss_prob_a[0]==1 && ses_ss_prob_b[0]==1) {
      sset = 0;
    } else if (ses_ss_prob_a[1]==1 && ses_ss_prob_b[1]==1){
      sset = 1;
    } else if ((ses_ss_prob_a[1]==0 && ses_ss_prob_b[1]==1) || (ses_ss_prob_a[1]==1 && ses_ss_prob_b[1]==0)){
      sset = 2;
    } else {sset = 0;}

    if (nres > 29) {
      mean = ll_means_[29+(sset*30)];
    }else {mean = ll_means_[nres+(sset*30)];};

  return mean;
};

double StructureElementConnectivityRestraint::unprotected_evaluate(DerivativeAccumulator *accum)
    const {
  //IMP_CHECK_OBJECT(a_.get());
  //IMP_CHECK_OBJECT(b_.get());
  IMP_CHECK_OBJECT(score_func_);
  IMP_OBJECT_LOG;

  // If particle a_ is not set up as StructureElements, this is a dummy restraint that evaluates to zero
  if(threading::StructureElement().get_is_setup(get_model(), a_)==false){
	  return 0.0;
  }


  // Compute distance between C-term of SE a_ and N-term of SE b_
  double mod_dist = get_model_distance();
  
  // Compute "max" distance per restraint
  double max_dist = get_max_distance();
  
  // Evaluate the difference between our distance and max distance
  double distance = mod_dist - max_dist;

  //std::cout<<"model distance "<<mod_dist<<" max distance "<<max_dist<<std::endl;
  
  double score;

  // If the number of residues is less than zero, then we have overlap between SEs.
  // This should not happen in practice, but if checks are not implemented at the system
  // level, an overlap will return a score of 10000
  if (get_number_of_residues() <= 1) {
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
