/**
 *  \file threading/LoopPairDistanceRestraint.cpp
 *  \brief A distance restraint that can be evaluated on residues with no structure
 *
 *  Copyright 2007-2019 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/threading/LoopPairDistanceRestraint.h>
#include <IMP/core/XYZ.h>
#include <IMP/atom/Residue.h>
#include <IMP/atom/Hierarchy.h>
#include <IMP/atom/Selection.h>
#include <IMP/algebra/Vector3D.h>

IMPTHREADING_BEGIN_NAMESPACE

LoopPairDistanceRestraint::LoopPairDistanceRestraint(Model *m, UnaryFunction *score_func,                     
                    ParticleIndex a,
                    ParticleIndex b,
                    double n_sds,
                    std::string name)
    : Restraint(m, "LoopPairDistanceRestraint%1%"), score_func_(score_func), a_(a), b_(b), n_sds_(n_sds) {}

/* Find the closest built residues to particle pi. (Should put this as a helper function)
 */

static float hh_means_[30] = {3.81, 3.036, 2.836, 2.511, 2.275, 2.178, 2.026, 1.876, 1.835, 1.669, 1.658, 1.666, 1.625, 1.53,
                1.445, 1.374, 1.292, 1.212, 1.164, 1.133, 1.049, 1.043, 1.074, 0.977, 0.965, 0.938, 0.868, 0.824, 0.805,
                0.788};

static float hh_sds_[30] = {0.027, 0.284, 0.397, 0.441, 0.483, 0.499, 0.504, 0.537, 0.534, 0.538, 0.545, 0.507, 0.494, 0.468,
              0.447, 0.428, 0.439, 0.415, 0.432, 0.392, 0.382, 0.38, 0.401, 0.381, 0.38, 0.317, 0.328, 0.304, 0.318,
              0.273};


double LoopPairDistanceRestraint::calc_sphere_cap_distance(float R, float d, float alpha) const {
  // Given the radius of the sphere, R, distance from the center of the sphere to the chord, r, 
  // and an angle from the sphere cap plane, alpha = 0-pi/2, return the distnace from the center
  // of the sphere cap circle to the points on the spherical segment at angle alpha

  // The solution to this is quadratic and for values of alpha from -0.5 - pi/2, the
  // plus form will always be the positive distance
  float sina = sin(alpha);
  
  float dist = -d * sina + sqrt(abs(d*d*(sina*sina - 1 ) + R*R));
  if (dist < 0) { dist = -1*dist;}
  
  // For debugging.  If DIST == NaN, this trips.
  if (dist!=dist){
   int i=1;  
   std::cout << "THe sphere cap distance is undefined.  There is some error in the LoopPairDistanceRestraint calculation" << std::endl;
   //std::cout << "    NANNAN-R,d,a,RR,dd: " << R << ", " << d << ", " << alpha << " " << R*R << " " << d*d*(sina*sina - 1 ) << " " << sina << " || " << dist <<std::endl;
  }
  return dist;
}

ParticleIndexes LoopPairDistanceRestraint::get_sequence_residue_particles() const {

  ParticleIndexes parts;

  parts.push_back(a_);
  parts.push_back(b_);
  
  return parts;
}

Particles LoopPairDistanceRestraint::get_closest_built_residue_particles(ParticleIndex pi) const
{
  // TODO: Ensure that input particle is a residue
  
  ParticlesTemp endpoints;
  atom::Residue r = atom::Residue(get_model(), pi);
  atom::Hierarchy nextres;
  atom::Hierarchy prevres;
  int ri = r.get_index();
  
  // Get the chain of this residue
  atom::Hierarchy res(get_model(), pi);
  atom::Hierarchy seq_chain = res.get_parent();
  
  // Go forward in sequence space looking for a coordinate that is optimized (max = 1000)
  // This is very slow, especially for sparsely built models  and should be optimized
  for (unsigned int i=1; i<1000; i++){
    // Get the next residue particle
    nextres = atom::get_next_residue(r);
    
    // If you hit the end of the chain, then just leave
    if (nextres.get_particle()==0){ break;}
    if (core::XYZ(nextres.get_particle()).get_coordinates_are_optimized()){
      endpoints.push_back(nextres.get_particle());
      break;
    }
    r = atom::Residue(nextres.get_particle()); 
  }
  
  r = atom::Residue(get_model(), pi);
  // Go backward in sequence space looking for a coordinate that is optimized
  for (int i=-1; i>-1000; i--){
    // Get the previous residue
    prevres = atom::get_previous_residue(r);
    
    // If you hit the end of the chain, then just leave
    if (prevres.get_particle()==0){ break;}
    if (core::XYZ(prevres.get_particle()).get_coordinates_are_optimized()){
      endpoints.push_back(prevres.get_particle());
      break;
    } 
    r = atom::Residue(prevres.get_particle());
  }
  return endpoints;
}

algebra::Vector3D LoopPairDistanceRestraint::get_sphere_cap_center(Particle* P0, Particle* P1, float R0, float R1) const {
  // First, define the unit vector connecting the two centers
  algebra::Vector3D Dvec = core::XYZ(get_model(), P0->get_index()).get_coordinates() - core::XYZ(get_model(), P1->get_index()).get_coordinates();
  //algebra::Vector3D uvec = Dvec / Dvec.get_magnitude();
  
  // Then compute the distance needed to travel along the unit vector
  float h = 1/2 + (R1*R1-R0*R0)/(2 * Dvec.get_magnitude()*Dvec.get_magnitude());

  // First, if the loop lengths cannot span the two endpoints, return the center of the distance between them
  algebra::Vector3D center = core::XYZ(get_model(), P1->get_index()).get_coordinates() + h * Dvec;

  return center;
}

double LoopPairDistanceRestraint::get_loop_distance(int a, int b) const {
  // Get the statistical potential difference for the maximum loop length
  // between residues a and b
  int res_dist = abs(a - b);
  float mean_a;
  float sd_a;

  if(res_dist < 30) { 
    mean_a = hh_means_[res_dist-1]; sd_a = hh_sds_[res_dist-1];
  } else { 
    // If loop is >30 residues, use the dist/residue for 30 residues
    mean_a = hh_means_[29]; sd_a = hh_sds_[29];
  }

  // Max dist is mean plus a number of standard deviations above (hard-coded). 
  float dist = res_dist * (mean_a + n_sds_ * sd_a);

  return dist;
}

double LoopPairDistanceRestraint::unprotected_evaluate(DerivativeAccumulator *accum)
    const {
  //IMP_CHECK_OBJECT(a_.get());
  //IMP_CHECK_OBJECT(b_.get());
  IMP_CHECK_OBJECT(score_func_);
  IMP_OBJECT_LOG;
  algebra::Vector3D a_center;
  algebra::Vector3D b_center;
  float offset_a;
  float offset_b;
  float center_dist;

  // Find if either or both residues have coordinates (using are_optimized as a proxy for this)
  bool a_opt = core::XYZ(get_model(), a_).get_coordinates_are_optimized();
  bool b_opt = core::XYZ(get_model(), b_).get_coordinates_are_optimized();
  double score;

  // Determine XL evaluation distance in all circumstances...this could be written much better
  float distance = 1000;  // model distance
  float loop_dist;
  float new_dist;

  // If both residues are structured, computing the model distance is easy
  
  if (a_opt==true and b_opt==true) {
    distance = core::get_distance(core::XYZ(get_model(), a_), core::XYZ(get_model(), b_));
  
    // If not this gets very ugly
  } else {

    algebra::Vector3D D_vec_a;
    algebra::Vector3D D_vec_b;
    Particles ps_a;
    Particles ps_b;
    float R0_a;
    float R0_b;
    float R1_a;
    float R1_b;
    float alpha;
    float avec;
    double pi_over_2 = IMP::PI * 0.5;

    // If A is structured, then its center is 
    if (a_opt==true) {
      a_center = core::XYZ(get_model(), a_).get_coordinates();
    } else {

      // Find particles of closest built residues to b_
      ps_a = get_closest_built_residue_particles(a_);

      if (ps_a.size()==1){
	  // If there is only one endpoint, then the center of the uncertainty sphere is this atom.
	a_center = core::XYZ(get_model(), ps_a[0]->get_index()).get_coordinates();
      } else if (ps_a.size()==0) {
	  // If there are no endpoints found, then the crosslink endpoint is on a completely disordered chain.
	  // What does this mean?  For now, return a constant value, however this is not optimal
	  // However, the likelihood of a completely disordered chain is small.
          return 10;
      } else {
      // Get sphere cap parameters
        D_vec_a = core::XYZ(get_model(), ps_a[0]->get_index()).get_coordinates() - core::XYZ(get_model(), ps_a[1]->get_index()).get_coordinates(); // Vector between structured endpoints
        
        R0_a = get_loop_distance(atom::Residue(ps_a[0]).get_index(), atom::Residue(get_model(), a_).get_index());
        R1_a = get_loop_distance(atom::Residue(ps_a[1]).get_index(), atom::Residue(get_model(), a_).get_index());
      
        // Compute Center of sphere cap
        if (R0_a > R1_a + D_vec_a.get_magnitude()){
          a_center = core::XYZ(get_model(), ps_a[0]->get_index()).get_coordinates();
        } else if (R1_a > R0_a + D_vec_a.get_magnitude()) {
          a_center = core::XYZ(get_model(), ps_a[1]->get_index()).get_coordinates();
        } else {
          a_center = get_sphere_cap_center(ps_a[0], ps_a[1], R0_a, R1_a);        
        }

      }

    }

    //----------------------------------------
    // Now we do the same as above for Particle B
    if (b_opt==true) {
      b_center = core::XYZ(get_model(), b_).get_coordinates();
    } else {

      // Find particles of closest built residues to b_
      ps_b = get_closest_built_residue_particles(b_);
      //if both endpoints are in the same loop, it is satisfied - Do we really want to do this?
      //if (ps_b[0]==ps_a[0] && ps_b[1]==ps_a[1]){return 0;}
      //if (ps_b[1]==ps_a[0] && ps_b[0]==ps_a[1]){return 0;}


      // Get sphere cap parameters
      if (ps_b.size()==1) {
	  b_center = core::XYZ(get_model(), ps_b[0]->get_index()).get_coordinates();
      } else if (ps_b.size()==0) {
	  return 10;
      } else {
        D_vec_b = core::XYZ(get_model(), ps_b[0]->get_index()).get_coordinates() - core::XYZ(get_model(), ps_b[1]->get_index()).get_coordinates(); // Vector between structured endpoints
      
        R0_b = get_loop_distance(atom::Residue(ps_b[0]).get_index(), atom::Residue(get_model(), b_).get_index());
        R1_b = get_loop_distance(atom::Residue(ps_b[1]).get_index(), atom::Residue(get_model(), b_).get_index());
      
        // Center of uncertainty volume
        if (R0_b > R1_b + D_vec_b.get_magnitude()){
          b_center = core::XYZ(get_model(), ps_b[0]->get_index()).get_coordinates();
        } else if (R1_b > R0_b + D_vec_b.get_magnitude()) {
          b_center = core::XYZ(get_model(), ps_b[1]->get_index()).get_coordinates();
        } else {
          b_center = get_sphere_cap_center(ps_b[0], ps_b[1], R0_b, R1_b);
        }
      }
    }

    // ------------------------------------------------------
    // Now that we've computed the sphere cap centers, we can compute the model distance
    // Get vector from sc_center to other endpoint
    algebra::Vector3D xl_vec = a_center - b_center;
    center_dist = algebra::get_distance(a_center, b_center);
    
    // Now, subtract the distance from the center to the intersections of the sphere caps
    // We also consider the loop lengths when a_ is not ordered
    if (a_opt==true){
      offset_a = 0;
      // If a_ only has one endpoint (is on a terminal loop) the loop distance is simple to compute
    } else if (ps_a.size()==1){
      offset_a = get_loop_distance(atom::Residue(ps_a[0]).get_index(), atom::Residue(get_model(), a_).get_index());
    } else  {
      //If a_ has two endpoints, then we must compute the angle between the XL vector and D vector to find where it
      // intersects the sphere cap  
      avec = (D_vec_b * xl_vec) / (D_vec_b.get_magnitude() * xl_vec.get_magnitude());
      
      // Sometimes floating point errors cause problems at the bounds
      if (avec > 1.0){alpha=0;}
      else if (avec < -1.0){alpha=3.14159;}
      else{alpha = std::acos(avec);}

      if (D_vec_a.get_magnitude() > R1_a + R0_a){
        offset_a = 0;
      } else if (R0_a > R1_a + D_vec_a.get_magnitude()){
        offset_a = R0_a;
      } else if (R1_a > R0_a + D_vec_a.get_magnitude()){
        offset_a = R1_a;
      } else if (alpha > pi_over_2) {
        offset_a = calc_sphere_cap_distance(R1_a, algebra::get_distance(a_center, core::XYZ(get_model(), ps_a[1]->get_index()).get_coordinates()), alpha-pi_over_2);
      } else {
        offset_a = calc_sphere_cap_distance(R0_a, algebra::get_distance(a_center, core::XYZ(get_model(), ps_a[0]->get_index()).get_coordinates()), alpha);
      }
    }

    // Now, we do the same for B
    if (b_opt==true){
      offset_b = 0;
    } else if (ps_b.size()==1){
      offset_b = get_loop_distance(atom::Residue(ps_b[0]).get_index(), atom::Residue(get_model(), b_).get_index());
    } else  {
      // Get sphere cap angle
      avec = (D_vec_b * xl_vec) / (D_vec_b.get_magnitude() * xl_vec.get_magnitude());
      if (avec > 1.0){alpha=0;}
      else if (avec < -1.0){alpha=3.14159;}
      else{alpha = std::acos(avec);}
      
      if (D_vec_b.get_magnitude() > R1_b + R0_b){
        offset_b = 0;
      } else if (R0_b > R1_b + D_vec_b.get_magnitude()){
        offset_b = R0_b;
      } else if (R1_b > R0_b + D_vec_b.get_magnitude()){
        offset_b = R1_b;
      } else if (alpha < pi_over_2) {
        offset_b = calc_sphere_cap_distance(R1_b, algebra::get_distance(b_center, core::XYZ(get_model(), ps_b[1]->get_index()).get_coordinates()), alpha);
      } else {
        offset_b = calc_sphere_cap_distance(R0_b, algebra::get_distance(b_center, core::XYZ(get_model(), ps_b[0]->get_index()).get_coordinates()), alpha-pi_over_2);
      }
    }
    
    // The model distance is finally computed as the distance between the sphere cap centers minus offsets due to the sphere cap surfaces
    distance = center_dist - offset_a - offset_b;
  }

  // Evaluate the score based on this distance
  score = score_func_->evaluate(distance);
  //std::cout << " || Score: " << score << " " << distance << " " << center_dist << " " << offset_a << " " << offset_b << " | " << a_opt << " " << b_opt << std::endl; 
  return score;
}

/* Return all particles whose attributes are read by the restraints. To
   do this, ask the pair score what particles it uses.*/
ModelObjectsTemp LoopPairDistanceRestraint::do_get_inputs() const {
  ModelObjectsTemp ret;
  ret.push_back(get_model()->get_particle(a_));
  ret.push_back(get_model()->get_particle(b_));
  return ret; //ModelObjectsTemp(1, get_model()->get_particle(a_));
}

IMPTHREADING_END_NAMESPACE
