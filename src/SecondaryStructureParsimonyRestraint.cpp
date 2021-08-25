/**
 *  \file example/ExampleRestraint.cpp
 *  \brief Restrain a list of particle pairs.
 *
 *  Copyright 2007-2017 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/threading/SecondaryStructureParsimonyRestraint.h>
#include <IMP/atom/Residue.h>
#include <IMP/atom/Hierarchy.h>
#include <IMP/threading/StructureElement.h>
#include <IMP/atom/SecondaryStructureResidue.h>
#include <math.h>
#include <boost/algorithm/string.hpp>

IMPTHREADING_BEGIN_NAMESPACE

SecondaryStructureParsimonyRestraint::SecondaryStructureParsimonyRestraint(Model *m,                    
                    ParticleIndexes ses,
                    ParticleIndex b,
                    float c,
                    std::string name)
    : Restraint(m, name), ses_(ses), b_(b), c_(c){
      helix_prob_ = 0;
      coil_prob_ = 1;
      sheet_prob_ = 0;
    }

double SecondaryStructureParsimonyRestraint::unprotected_evaluate(DerivativeAccumulator *accum)
    const {

  // Ge the number of residues
  int nres = atom::Hierarchy(get_model(), b_).get_children().size();
  std::string prot_chain = atom::Hierarchy(get_model(), b_)->get_name();
  
  IMP_OBJECT_LOG;


  // Get the SE sequence probabilities
  //std::vector<std::vector<float> > ses_ss_probs = get_se_sequence_probabilities();

  // Initialize the SE secondary structure probability array to the default value
  std::vector<std::vector<double> > ses_ss_probs(nres, std::vector<double>(3,1)={helix_prob_,sheet_prob_,coil_prob_});

  // Populate the SE sec st array with covered SE residues
  for (unsigned nse = 0; nse < ses_.size(); nse++){

    // Get chain
    std::string chain = threading::StructureElement(get_model(),ses_[nse]).get_chain();
    if (boost::algorithm::contains(chain, prot_chain)) {
    
        // First, get the residue indexes for this SE.
        Ints resinds = threading::StructureElement(get_model(), ses_[nse]).get_resindex_list();
        // Second, get the probabilities for this SE
        Floats ses_ss_prob = atom::SecondaryStructureResidue(get_model(), ses_[nse]).get_all_probabilities();

        for(unsigned i = 0; i < resinds.size(); i++){
      
          for (int k = 0; k < 3; k++){
            ses_ss_probs[resinds[i]-1][k] = ses_ss_prob[k];
          }
        }
    }
  }

  atom::Residues resis = atom::Hierarchy(get_model(), b_).get_children();
  Float score = 0;

  // Loop over all the residues
  // Compare to sequence ss_probs
  for (int j = 0; j < nres; j++){

    Floats seq_prob = atom::SecondaryStructureResidue(get_model(), resis[j].get_particle_index()).get_all_probabilities();
    std::vector<double> se_prob = ses_ss_probs[j];
    float dot = 0;
    for (int k = 0; k < 3; k++){
        //std::cout << "seq_prob " << seq_prob[k] << " se_prob " << se_prob[k] << std::endl;
      dot += seq_prob[k] * se_prob[k];
    }
    score += -1 * log(dot+0.0001); // to avoid log(0) problem
    //std::cout << "RES# " << j << " (" << se_prob[0] << se_prob[1] << se_prob[2] << "), (" << seq_prob[0] << seq_prob[1] << seq_prob[2] << ") " << dot << std::endl;

  }
  return score / (nres * 1.0) * c_;
}


/* Return all particles whose attributes are read by the restraints. To
   do this, ask the pair score what particles it uses.*/
ModelObjectsTemp SecondaryStructureParsimonyRestraint::do_get_inputs() const {
  return ModelObjectsTemp(1, get_model()->get_particle(b_));
}

IMPTHREADING_END_NAMESPACE
