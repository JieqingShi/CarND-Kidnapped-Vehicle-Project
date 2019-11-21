/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

#include <limits>

using std::string;
using std::vector;

using std::normal_distribution;
static std::default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  num_particles = 100;  // TODO: Set the number of particles
  double std_x = std[0];  // 0.3
  double std_y = std[1];  // 0.3
  double std_theta = std[2];  // 0.01

  // Generate normal distributions according to respective means and stds
  normal_distribution<double> dist_x(x, std_x);
  normal_distribution<double> dist_y(y, std_y);
  normal_distribution<double> dist_theta(theta, std_theta);

  
  for(int i=0; i < num_particles; i++){
    Particle pa;  // current particle
    pa.id = i;
    pa.x = dist_x(gen);
    pa.y = dist_y(gen);
    pa.theta = dist_theta(gen);
    pa.weight = 1.0;

    particles.push_back(pa);  // push back to particles vector
    weights.push_back(pa.weight);
  }
  is_initialized = true;  // set flag
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  double std_x = std_pos[0];
  double std_y = std_pos[1];
  double std_theta = std_pos[2];
  
  // Generate normal distributions for noise according to mean 0 and respective stds
  normal_distribution<double> dist_x(0, std_x);
  normal_distribution<double> dist_y(0, std_y);
  normal_distribution<double> dist_theta(0, std_theta);

  // todo: iterate through vector instead of index
  for(int i=0; i<num_particles; i++){
    double x = particles[i].x;
    double y = particles[i].y;
    double theta = particles[i].theta;

    // Prediction step
    // todo: add if else mechanism to avoid division by zero if yaw_rate very small
    if (fabs(yaw_rate) > 0.00001) {
      x += velocity/yaw_rate*(sin(theta+yaw_rate*delta_t)-sin(theta));
      y += velocity/yaw_rate*(cos(theta)-cos(theta+yaw_rate*delta_t));
      theta += yaw_rate*delta_t;
    }

    else{
      x += velocity*delta_t*cos(theta);
      y += velocity*delta_t*sin(theta);
    }

    // Add random gaussian noise
    particles[i].x = x + dist_x(gen);
    particles[i].y = y + dist_y(gen);
    particles[i].theta = theta + dist_theta(gen);
    
	}
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  for(unsigned int i=0; i<observations.size(); i++){
    double obs_x = observations[i].x;
    double obs_y = observations[i].y;
    double distance;
    double min_dist = std::numeric_limits<double>::infinity();
    int nearestId = -1;

    for(unsigned int j=0; j<predicted.size(); j++){
      // double preds_id = predicted[j].id;  // not needed
      double preds_x = predicted[j].x;
      double preds_y = predicted[j].y;
      distance = dist(obs_x, obs_y, preds_x, preds_y);
      if(distance < min_dist){
        min_dist = distance;
        nearestId = predicted[j].id;  // or is it j?
      }
    }
    observations[i].id = nearestId;
  }
}


void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {

  double std_ldm_x = std_landmark[0];
  double std_ldm_y = std_landmark[1];
  

  for(int i=0; i<num_particles; i++){
   
    double xp = particles[i].x;
    double yp = particles[i].y;
    double thetap = particles[i].theta;
    vector<LandmarkObs> landmarkWithinSensorRange;

    // For each particle, consider only landmarks in vicinity of sensor_range
    for(unsigned int j=0; j<map_landmarks.landmark_list.size();++j){
      double xm = map_landmarks.landmark_list[j].x_f;
      double ym = map_landmarks.landmark_list[j].y_f;
      int lm_id = map_landmarks.landmark_list[j].id_i; 
      double dist_particle2lm = dist(xp, yp, xm, ym);

      if(dist_particle2lm<=sensor_range){
        // push measurement to vector of immediate landmarks
        landmarkWithinSensorRange.push_back(LandmarkObs{lm_id, xm, ym});
      }
    }
    // Transform observation coordinates to particles coordinate system
    vector<LandmarkObs> transformedObs;
    for(unsigned int o=0; o<observations.size(); o++){
      int obs_id = observations[o].id;
      double xo = observations[o].x;
      double yo = observations[o].y;

      // double xt, yt;
      double xt = xp + cos(thetap) * xo - sin(thetap) * yo;
      double yt = yp + sin(thetap) * xo + cos(thetap) * yo;
      transformedObs.push_back(LandmarkObs{obs_id, xt, yt});
    }
    
    // Associate transformed observations to nearest landmark
    dataAssociation(landmarkWithinSensorRange, transformedObs);


    // compute weights using multivariate Gaussian
    //particles[i].weight = 1.0;  // todo: does it have to be re-initialized?
    for(unsigned int t=0; t<transformedObs.size(); t++){
      double x = transformedObs[t].x;
      double y = transformedObs[t].y;
      int associatedLmId = transformedObs[t].id;
      

      // find x and y of nearest landmark
      double mu_x, mu_y;
      for(unsigned int l=0; l<landmarkWithinSensorRange.size(); ++l){
        if(landmarkWithinSensorRange[l].id == associatedLmId){
          mu_x = landmarkWithinSensorRange[l].x;
          mu_y = landmarkWithinSensorRange[l].y;
          break;
        }
      }

      double weight_obs = multiv_prob(std_ldm_x, std_ldm_y, x, y, mu_x, mu_y);
      particles[i].weight *= weight_obs;
      weights[i] = particles[i].weight;  // use this vector for the resampling step
    }
  }
}

  

/*
void ParticleFilter::resample() {

  vector<Particle> sampledParticles(num_particles);
  //std::default_random_engine gen;
  
  // get all of the current weights
  vector<double> weights;
  for (int i = 0; i < num_particles; i++) {
    weights.push_back(particles[i].weight);
  }

  double beta = 0.0;
  double mw = *std::max_element(weights.begin(), weights.end());
  std::cout<<"max weight: "<<mw<<"\n";
  std::uniform_real_distribution<double> distribution(0.0, mw);

  int index = rand() % num_particles;  // todo: to num_particles or num_particles-1; can also use std::uniform_int_distribution
  std::cout<<"Generating random index: "<<index<<"\n";
  for(int i=0; i<num_particles; ++i){
    beta += distribution(gen)*2.0;
    while(weights[index] < beta){
      beta -= weights[index];
      index = (index+1) % num_particles;
    }
    sampledParticles.push_back(particles[index]);
    //std::cout<<"Sampling particle: "<<particles[index]<<"\n";
  }
  std::cout<<"Sampled "<<sampledParticles.size()<<"particles \n";
  particles = sampledParticles;
}
*/
  

void ParticleFilter::resample() {
  // TODO: Resample particles with replacement with probability proportional to their weight. 
  // NOTE: You may find std::discrete_distribution helpful here.
  //   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
  
  vector<Particle> new_particles;

  // get all of the current weights
  /*
  vector<double> weights;
  for (int i = 0; i < num_particles; i++) {
    weights.push_back(particles[i].weight);
  }
  */

  // generate random starting index for resampling wheel
  std::uniform_int_distribution<int> uniintdist(0, num_particles-1);
  auto index = uniintdist(gen);

  // get max weight
  double max_weight = *std::max_element(weights.begin(), weights.end());

  // uniform random distribution [0.0, max_weight)
  std::uniform_real_distribution<double> unirealdist(0.0, max_weight);

  double beta = 0.0;

  // spin the resample wheel!
  for (int i = 0; i < num_particles; i++) {
    beta += unirealdist(gen) * 2.0;
    while (beta > weights[index]) {
      beta -= weights[index];
      index = (index + 1) % num_particles;
    }
    new_particles.push_back(particles[index]);
  }

  particles = new_particles;
}


void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}