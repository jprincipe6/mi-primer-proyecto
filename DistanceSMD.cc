
#include "DistanceSMD.hh"
#include <numeric/xyzVector.fwd.hh>
#include <numeric/xyzVector.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/conformation/Conformation.hh>
#include <vector>
#include <string>
#include <sstream>
#include <map>
#include <utility>

#define NORM_MULTIPLIED_FACTOR 0.5f

DistanceSMD::DistanceSMD(const core::pose::PoseOP& p,  std::string ss_in) {
  ss = ss_in;
  pose_ind = p->clone();
  bool init_found = false;
  bool end_found = false;
  int init_residue = 0, end_residue = 0;
  for (int i = 1; i < ss.size(); i++) {
    if ((ss[i] != 'L') && (ss[i -1] != ss[i])) {
      std::cout << "ss at " << i << " val " << ss[i] << std::endl;
      init_residue = i + 1;
      init_found = true;
    }
    if ((ss[i] == 'L') && (ss[i -1] != ss[i])) {
      std::cout << "ss at " << i - 1 << " val " << ss[i - 1] << std::endl;
      end_residue = i - 1;
      end_found = true;
    }
    if (init_found && end_found) {
      init_found = false;
      end_found = false;
      int half_residue = std::floor( (end_residue + init_residue) / 2 );
      selected_residues_for_rmsd.push_back( half_residue );
    }
  }

  build_inter_distances_of_straight();
  build_inter_distance_for_an_individual(*p, inter_distance_target);
}

void
DistanceSMD::build_inter_distances_of_straight() {
  core::pose::PoseOP straight = pose_ind->clone();
  for (int i = 1; i < straight->total_residue(); i++) {
    straight->set_phi(i, 180.0);
    straight->set_psi(i, 180.0);
  }
  int half_res = straight->total_residue() / 2;
  core::PointPosition calpha1_pos_aux  = straight->residue(1).xyz("CA");
  core::PointPosition calpha2_pos_aux = straight->residue(half_res).xyz("CA");
  double max_dist = calpha1_pos_aux.distance(calpha2_pos_aux);

  for ( int i = 0; i < selected_residues_for_rmsd.size(); ++i)  {
    for (int j  = i + 1; j < selected_residues_for_rmsd.size(); j++) {
      core::Size res_num_1 = core::Size(selected_residues_for_rmsd[i]);
      core::Size res_num_2 = core::Size(selected_residues_for_rmsd[j]);
      core::PointPosition calpha1_pos  = straight->residue(res_num_1).xyz("CA");
      core::PointPosition calpha2_pos = straight->residue(res_num_2).xyz("CA");
      double dist = calpha1_pos.distance(calpha2_pos);
      // falta normalizar esta distancia dividiendo por al distancia maxima ( la misma dist para la proteina estirada)
      //inter_distance_individual_1.push_back( std::sqrt( std::pow(dist, 2) / 2)  );
      inter_dist_norm_max[std::make_pair<int, int>(res_num_1, res_num_2)] = max_dist ;
    }
  }
}


void
DistanceSMD::build_inter_distance_for_an_individual(core::pose::PoseOP pose_ind, std::vector<double>& inter_distance_individual) {
  inter_distance_individual.resize(0);
  for ( int i = 0; i < selected_residues_for_rmsd.size(); ++i)  {
    for (int j  = i + 1; j < selected_residues_for_rmsd.size(); j++) {
      core::Size res_num_1 = core::Size(selected_residues_for_rmsd[i]);
      core::Size res_num_2 = core::Size(selected_residues_for_rmsd[j]);
      core::PointPosition calpha1_pos  = pose_ind->residue(res_num_1).xyz("CA");
      core::PointPosition calpha2_pos = pose_ind->residue(res_num_2).xyz("CA");
      double dist = calpha1_pos.distance(calpha2_pos);
      // normalizar esta distancia dividiendo por al distancia maxima ( la misma dist para la proteina estirada)
      //inter_distance_individual.push_back( std::sqrt( std::pow(dist, 2) / 2)  );
      if (inter_dist_norm_max[std::pair<int, int>(res_num_1, res_num_2)] == 0) {
	inter_distance_individual.push_back(0);
      } else {
	double normalized_dist = dist / ( NORM_MULTIPLIED_FACTOR * inter_dist_norm_max[std::pair<int, int>(res_num_1, res_num_2)] );
	//std::cout << " dist " << dist << " norm " << normalized_dist << " max " << inter_dist_norm_max[std::pair<int, int>(res_num_1, res_num_2)] << std::endl;
	inter_distance_individual.push_back( normalized_dist );
      }
    }
  }
}

double
DistanceSMD::distance_calculation(core::pose::Pose& pose) {
  double sum = 0.0;
  std::vector<double> inter_distance_individual;
  build_inter_distance_for_an_individual(*pose, inter_distance_individual);

  for ( int i = 0; i < inter_distance_individual_1.size(); ++i)  {
    sum += std::pow( inter_distance_target[i] - inter_distance_individual[i], 2);
  }

  double result = std::sqrt(sum);
  result = result / inter_distance_target.size();

  return result;
}

double
DistanceSMD::current_distance_calculation(core::pose::Pose& pose_1, core::pose::Pose& pose_2) {
  double sum = 0.0;
  std::vector<double> inter_distance_individual_1, inter_distance_individual_2;
  build_inter_distance_for_an_individual(*pose_1, inter_distance_individual_1);
  build_inter_distance_for_an_individual(*pose_2, inter_distance_individual_2);

  for ( int i = 0; i < inter_distance_individual_1.size(); ++i)  {
    sum += std::pow( inter_distance_individual_1[i] - inter_distance_individual_2[i], 2);
  }

  double result = std::sqrt(sum);
  result = result / inter_distance_individual_1.size();

  return result;
}

