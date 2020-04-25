
/*
  Calculo de euclidean distance de las interdistancias usando la definicion de Mario et.al..

  Ejemplo de uso:

  Se necesita un pose "native" y una estructura secundaria.

  Puede ser una cadena introducida manualmente:
  std::string secondary_structure = "LHHHLHHHLHHHL"

  o la cadena que viene con la estructura nativa:
  std::string secondary_structure = native->secstruct();

  Estanciar con:
  DistanceSMDPtr calculo_smd = DistanceSMDPtr( new DistanceSMD(native, secondary_structure)  );

  Es importante iniciar la clase con la proteina native.

  para calcular la distancia de 1 individuo "pose", 2 optionces:

  1 ) resultado = calculo_smd->distance_calculation(pose);
  2 ) resultado = calculo_smd->current_distance_calculation(native, pose);

  Esta segunda opción permite hacer el calculo con 2 individuos cualesquiera, ejemplo:

  resultado = calculo_smd->current_distance_calculation(pose1, pose2);

 */

#ifndef CALCULATEDISTANCEPOPULATION_H
#define CALCULATEDISTANCEPOPULATION_H

#include <map>
#include <boost/shared_ptr.hpp>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/rms_util.hh>

class DistanceSMD
{
public:
  core::pose::PoseOP pose_ind;
  std::string ss;
  std::vector<double> inter_distance_target;
  std::vector<int> selected_residues_for_rmsd;
  std::vector<std::vector<double> > inter_distances_per_ind;
  std::map<std::pair<int, int>, double > inter_dist_norm_max;

  DistanceSMD() {}

  DistanceSMD(const core::pose::PoseOP& p , std::string ss_in) ;

  double
  distance_calculation(core::pose::Pose& pose);


  double
  current_distance_calculation(core::pose::Pose& pose_1, core::pose::Pose& pose_2);

  void
  build_inter_distances_of_straight();

  void
  build_inter_distance_for_an_individual(core::pose::Pose pose_ind_1, std::vector<double>& inter_distance_individual_1);

};

typedef boost::shared_ptr<DistanceSMD> DistanceSMDPtr;


/* Other class members that could be included
  // core::pose::PoseOP pose_ind, pose_other;
  // std::vector<core::pose::PoseOP> popul_pdb;
  // core::scoring::ScoreFunctionOP scorefxn;
  // std::string ss;
  // core::pose::PoseOP native_;
*/



#endif /* CALCULATEDISTANCEPOPULATION_H */
