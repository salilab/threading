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

  float dist = -d * sina + sqrt(d*d*(sina*sina - 1 ) + R*R);
  if (dist < 0) { dist = -1*dist;}
  //std::cout << "     R,d,a,RR,dd: " << R << ", " << d << ", " << alpha << " " << R*R << " " << d*d*(sina*sina - 1 ) << " || " << dist <<std::endl;
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

  // create vector for the endpoint particles
  ParticlesTemp endpoints;
  atom::Residue r = atom::Residue(get_model(), pi);
  int ri = r.get_index();
  //std::cout << "Closest built resis for: " << ri << "| "; 
  // Get the chain of this residue
  atom::Hierarchy res(get_model(), pi);
  atom::Hierarchy seq_chain = res.get_parent();
  // Go forward in sequence space looking for a coordinate that is optimized (max = 1000)
  for (unsigned int i=1; i<1000; i++){
    atom::Selection sel(seq_chain);
    sel.set_residue_index(ri + i);
    ParticlesTemp nextres = sel.get_selected_particles();
    // If you hit the end of the chain, then just leave
    if (nextres.size() == 0){ break;}
    if (core::XYZ(nextres[0]).get_coordinates_are_optimized()){
      endpoints.push_back(nextres[0]);
      //std::cout << " " << ri+i;
      break;
    } 
  }

  // Go backward in sequence space looking for a coordinate that is optimized
  for (int i=-1; i>-1000; i--){
    atom::Selection sel(seq_chain);
    sel.set_residue_index(ri + i);
    ParticlesTemp nextres = sel.get_selected_particles();
    // If you hit the end of the chain, then just leave
    if (nextres.size() == 0){ std::cout << std::endl; break;}
    if (core::XYZ(nextres[0]).get_coordinates_are_optimized()){
      endpoints.push_back(nextres[0]);
      //std::cout << " " << ri+i << std::endl;   
      break;
    } 
  }
  return endpoints;
}

algebra::Vector3D LoopPairDistanceRestraint::get_sphere_cap_center(Particle* P0, Particle* P1, float R0, float R1) const {
  // First, define the unit vector

  algebra::Vector3D Dvec = core::XYZ(get_model(), P0->get_index()).get_coordinates() - core::XYZ(get_model(), P1->get_index()).get_coordinates();
  //algebra::Vector3D uvec = Dvec / Dvec.get_magnitude();

  float h = 1/2 + (R1*R1-R0*R0)/(2 * Dvec.get_magnitude()*Dvec.get_magnitude());

  //std::cout <<"H: " << h; //" | R0, R1 " << R0 << ", " << R1 << " |Dvec " << Dvec.get_magnitude() << std::endl;

  // First, if the loop lengths cannot span the two endpoints, return the center of the distance between them
  algebra::Vector3D center = core::XYZ(get_model(), P1->get_index()).get_coordinates() + h * Dvec;

  //std::cout << " Res Endpoints: " << atom::Residue(get_model(), P0->get_index()) << " | " << atom::Residue(get_model(), P1->get_index()) << std::endl;

  //std::cout << " P0;P1;cen: " << core::XYZ(get_model(), P0->get_index()).get_coordinates() << " | " << core::XYZ(get_model(), P1->get_index()).get_coordinates() << " || " << center << std::endl;

  return center;
}

double LoopPairDistanceRestraint::get_loop_distance(int a, int b) const {
  
  int res_dist = abs(a - b);
  float mean_a;
  float sd_a;

  if(res_dist < 30) { 
    mean_a = hh_means_[res_dist-1]; sd_a = hh_sds_[res_dist-1];
  } else { 
    mean_a = hh_means_[29]; sd_a = hh_sds_[29];
  }

  float dist = res_dist * (mean_a + n_sds_ * sd_a);

  //std::cout << a << " " << b << " Res dist: " << res_dist << " Dist:" << dist << " Mean/SD: " << mean_a << "/" << sd_a << std::endl; 

  return dist;
}

double LoopPairDistanceRestraint::unprotected_evaluate(DerivativeAccumulator *accum)
    const {
  //IMP_CHECK_OBJECT(a_.get());
  //IMP_CHECK_OBJECT(b_.get());
  IMP_CHECK_OBJECT(score_func_);
  IMP_OBJECT_LOG;

  // Find if either or both residues have coordinates (using are_optimized as a proxy for this)
  bool a_opt = core::XYZ(get_model(), a_).get_coordinates_are_optimized();
  bool b_opt = core::XYZ(get_model(), b_).get_coordinates_are_optimized();
  double score;

  // Determine XL evaluation distance in all circumstances...this could be written much better
  float distance = 1000;
  float loop_dist;
  float new_dist;
  //std::cout << "||||---- resis:" << atom::Residue(get_model(), a_).get_index() << " " << atom::Residue(get_model(), b_).get_index() << std::endl;

  // If both residues are structured, this is easy
  if (a_opt==true and b_opt==true) {
    distance = core::get_distance(core::XYZ(get_model(), a_), core::XYZ(get_model(), b_));

  } else {

    algebra::Vector3D a_center;
    algebra::Vector3D b_center;
    algebra::Vector3D D_vec_a;
    algebra::Vector3D D_vec_b;
    Particles ps_a;
    Particles ps_b;
    float offset_a;
    float offset_b;
    float R0_a;
    float R0_b;
    float R1_a;
    float R1_b;
    double pi_over_2 = IMP::PI * 0.5;


    if (a_opt==true) {
      a_center = core::XYZ(get_model(), a_).get_coordinates();
    } else {

      // Find particles of closest built residues to b_
      ps_a = get_closest_built_residue_particles(a_);

      // Get sphere cap parameters
      D_vec_a = core::XYZ(get_model(), ps_a[0]->get_index()).get_coordinates() - core::XYZ(get_model(), ps_a[1]->get_index()).get_coordinates(); // Vector between structured endpoints

      R0_a = get_loop_distance(atom::Residue(ps_a[0]).get_index(), atom::Residue(get_model(), a_).get_index());
      R1_a = get_loop_distance(atom::Residue(ps_a[1]).get_index(), atom::Residue(get_model(), a_).get_index());

      // Center of sphere cap
      if (R0_a > R1_a + D_vec_a.get_magnitude()){
        a_center = core::XYZ(get_model(), ps_a[0]->get_index()).get_coordinates();
      } else if (R1_a > R0_a + D_vec_a.get_magnitude()) {
        a_center = core::XYZ(get_model(), ps_a[1]->get_index()).get_coordinates();
      } else {
        a_center = get_sphere_cap_center(ps_a[0], ps_a[1], R0_a, R1_a);        
      }


    }

    if (b_opt==true) {
      b_center = core::XYZ(get_model(), b_).get_coordinates();
    } else {

      // Find particles of closest built residues to b_
      ps_b = get_closest_built_residue_particles(b_);

      // Get sphere cap parameters
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

    // **********************************************************
    // Get vector from sc_center to other endpoint
    algebra::Vector3D xl_vec = a_center - b_center;
    float center_dist = algebra::get_distance(a_center, b_center);

    // Now, subtract the distance from the center to the intersections of the sphere caps
    if (a_opt==true){
      offset_a = 0;
    } else {
      // Get sphere cap angle
      float alpha = std::acos((D_vec_a * xl_vec) / (D_vec_a.get_magnitude() * xl_vec.get_magnitude()));
      //std::cout << "ACEN: " << a_center << " a: " << alpha << " " << D_vec_a << " " << xl_vec << " " << pi_over_2 << std::endl;
      //std::cout << " ---dveca, xlvec " << D_vec_a.get_magnitude() << " " << xl_vec.get_magnitude() << " " << alpha << std::endl;
      // If less than PI/2, use the R1 sphere cap
      if (D_vec_a.get_magnitude() > R1_a + R0_a){
        offset_a = 0;
      } else if (R0_a > R1_a + D_vec_a.get_magnitude()){
        offset_a = R0_a;
      } else if (R1_a > R0_a + D_vec_a.get_magnitude()){
        offset_a = R1_a;
      } else if (alpha > pi_over_2) {
        //std::cout << core::XYZ(get_model(), ps_a[1]->get_index()).get_coordinates() << std::endl;
        offset_a = calc_sphere_cap_distance(R1_a, algebra::get_distance(a_center, core::XYZ(get_model(), ps_a[1]->get_index()).get_coordinates()), alpha-pi_over_2);
      } else {
        offset_a = calc_sphere_cap_distance(R0_a, algebra::get_distance(a_center, core::XYZ(get_model(), ps_a[0]->get_index()).get_coordinates()), alpha);
      }
    }

    // Now, subtract the distance from the center to the intersections of the sphere caps
    if (b_opt==true){
      offset_b = 0;
    } else {
      // Get sphere cap angle
      float alpha = std::acos((D_vec_b * xl_vec) / (D_vec_b.get_magnitude() * xl_vec.get_magnitude()));
      //std::cout <<  "BCEN: " << b_center << " a: " << alpha << " " << D_vec_b << " " << xl_vec << std::endl;
      //std::cout << " ---dvecb, xlvec " << D_vec_b.get_magnitude() << " " << xl_vec.get_magnitude() << " " << alpha << " " << R0_b << " " << R1_b << std::endl;
      // If less than PI/2, use the R1 sphere cap
      if (D_vec_b.get_magnitude() > R1_b + R0_b){
        offset_b = 0;
      } else if (R0_b > R1_b + D_vec_b.get_magnitude()){
        offset_b = R0_b;
      } else if (R1_b > R0_b + D_vec_b.get_magnitude()){
        offset_b = R1_b;
      } else if (alpha < pi_over_2) {
        offset_b = calc_sphere_cap_distance(R1_b, algebra::get_distance(b_center, core::XYZ(get_model(), ps_b[1]->get_index()).get_coordinates()), alpha);
        //std::cout << "Radius: " << R1_b << " Chord dist: " << algebra::get_distance(b_center, core::XYZ(get_model(), ps_b[1]->get_index()).get_coordinates())  << std::endl;
      } else {
        offset_b = calc_sphere_cap_distance(R0_b, algebra::get_distance(b_center, core::XYZ(get_model(), ps_b[0]->get_index()).get_coordinates()), alpha-pi_over_2);
        //std::cout << "Radius: " << R0_b << " Chord dist: " << algebra::get_distance(b_center, core::XYZ(get_model(), ps_b[0]->get_index()).get_coordinates())  << std::endl;
      }
    }
    distance = center_dist - offset_a - offset_b;
    //std::cout << "---Cendist: " << center_dist << " offA: " << offset_a << " offB: " << offset_b << std::endl;   
  }

  // Evaluate the score based on the shortest distance
  score = score_func_->evaluate(distance);
  //std::cout << " || Score: " << score << std::endl; 
  return score;
}

/* Return all particles whose attributes are read by the restraints. To
   do this, ask the pair score what particles it uses.*/
ModelObjectsTemp LoopPairDistanceRestraint::do_get_inputs() const {
  return ModelObjectsTemp(1, get_model()->get_particle(a_));
}

IMPTHREADING_END_NAMESPACE
