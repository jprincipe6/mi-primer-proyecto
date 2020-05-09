// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file TrialMover
/// @brief performs a move and accepts it according to Monte Carlo accept/reject criterion.
/// @author Monica Berrondo

#include <protocols/moves/TrialMover.hh>

// Rosetta Headers
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <core/pose/PDBInfo.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/rms_util.hh>

#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <basic/options/option.hh>
#include <basic/options/util.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

// Utility headers

// C++ headers
#include <string>
#include <basic/datacache/DataMap.hh>
#include <utility/vector0.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>
#include <vector>
#include <sys/types.h>
#include <dirent.h>
#include <fstream>      // std::ofstream
#include <iostream>
#include <boost/numeric/conversion/cast.hpp>

std::vector<std::string> get_paths_pdbs_from_dir(const char* path){
    std::vector<std::string> results;
    DIR* dirFile = opendir(path);
    if (dirFile){
        struct dirent* hFile;
        errno = 0;
        while ((hFile = readdir(dirFile)) != NULL) {
            if( !strcmp(hFile->d_name, ".")) continue;
            if( !strcmp(hFile->d_name, "..")) continue;
            
            // IN LINUX HIDDEN FILES ALL START WITH '.'
            //if ( gIgnoreHidden && ( hFile->d_name[0] == '.' )) continue;
            // dirFile.name is the name of the file. Do whatever string comparison
            // you want here. Something like:
            if ( strstr( hFile->d_name, ".pdb" )){
                printf( " found an .pdb file: %s \n", hFile->d_name);
                results.push_back(hFile->d_name);
            }
        }
        closedir(dirFile);
    }
    return results;
}

namespace protocols {
namespace moves {

static basic::Tracer tr( "protocols.TrialMover" );

using namespace core;
const char *path_input = "./soluciones_1elwA";
MonteCarloUtil::MonteCarloUtil() : mc_(/* 0 */)
{
    
}

MonteCarloUtil::MonteCarloUtil( protocols::moves::MonteCarloOP mc) : mc_(std::move(mc))
{
    
}

MonteCarloUtil::~MonteCarloUtil() = default;

void MonteCarloUtil::apply(Pose & pose)
{
    if ( mode_ == "reset" ) {
        mc_->reset(pose);
    } else if ( mode_ == "recover_low" ) {
        mc_->recover_low(pose);
    } else {
        utility_exit_with_message("MonteCarloUtil mode must be 'reset' or 'recover_low', this should have been caught earlier, dying");
    }
    
}

void MonteCarloUtil::parse_my_tag(
                                  utility::tag::TagCOP tag,
                                  basic::datacache::DataMap & data,
                                  protocols::filters::Filters_map const & /* filters */,
                                  protocols::moves::Movers_map const & /* movers */,
                                  core::pose::Pose const & /* pose */)
{
    if ( !tag->hasOption("mode") ) {
        throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "you must specify option mode in MonteCarloUtil");
    }
    if ( !tag->hasOption("montecarlo") ) {
        throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "you must specify the option montecarlo in MonteCarloUtil");
    }
    
    mode_ = tag->getOption<std::string>("mode");
    std::string const mc_name(tag->getOption<std::string>("montecarlo"));
    
    if ( mode_ != "reset" && mode_ != "recover_low" ) {
        throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "the option mode must be set to either reset or recover_low in MonteCarloUtil");
    }
    
    mc_ = data.get_ptr<protocols::moves::MonteCarlo>( "montecarlos",mc_name);
}

protocols::moves::MoverOP MonteCarloUtil::clone() const
{
    return utility::pointer::make_shared< MonteCarloUtil >(mc_);
}

protocols::moves::MoverOP MonteCarloUtil::fresh_instance() const
{
    return utility::pointer::make_shared< MonteCarloUtil >();
}


std::string MonteCarloUtil::get_name() const
{
    return "MonteCarloUtil";
}

TrialMover::TrialMover() :
start_weight_( 0.0 ),
original_weight( 0.0 ),
ramp_weight( 0.0 ),
delta( 0.0 ),
stats_type_( all_stats )
{}

std::vector<std::string> tokenize(const std::string& s, char c) {
    auto end = s.cend();
    auto start = end;

    std::vector<std::string> v;
    for( auto it = s.cbegin(); it != end; ++it ) {
        if( *it != c ) {
            if( start == end )
                start = it;
            continue;
        }
        if( start != end ) {
            v.emplace_back(start, it);
            start = end;
        }
    }
    if( start != end )
        v.emplace_back(start, end);
    return v;
}


core::Real getUmbralLimite(){
    std::ifstream myfile;
    std::string line;
    std::string value;
    std::vector<std::string> elems;
    double umbral = 0;
    myfile.open ("./umbral/umbral.txt", std::ios::in);
    while (getline(myfile, line)) {
        std::string str2 = line.substr (8,line.size());
        value = str2;
    }
    elems = tokenize(value,';');
//    std::cout << "Valores: " << elems[3] << std::endl;
    umbral  = std::stod(elems[3]);
    return umbral;
}
// constructor with arguments
TrialMover::TrialMover( MoverOP mover_in, MonteCarloOP mc_in ) :
start_weight_( 0.0 ),
original_weight( 0.0 ),
ramp_weight( 0.0 ),
delta( 0.0 ),
stats_type_( all_stats )
{
    
    using namespace core;
    using namespace core::import_pose;
    using namespace pose;
    using namespace basic::options;
    using namespace basic::options::OptionKeys;
    mover_ = mover_in;
    mc_ = mc_in;
    
    rmsd_vs_actual_acc = 0;
    cont_total_rmsd_vs_actual_acc = 0;
    paths_soluciones_pdbs = get_paths_pdbs_from_dir(path_input);
    acomuladorDeAceptadosCustom = 0;
    acomuladorDeAceptadosNormal = 0;
    countApplys = 0;
                                                                                                
}

// Copy constructor
TrialMover::TrialMover(TrialMover const & object_to_copy) : Mover(object_to_copy),
mover_(object_to_copy.mover_),
mc_(object_to_copy.mc_),
start_weight_(object_to_copy.start_weight_),
original_weight(object_to_copy.original_weight),
ramp_weight(object_to_copy.ramp_weight),
delta(object_to_copy.delta),
stats_type_(object_to_copy.stats_type_)
{}

TrialMover::~TrialMover() = default;

// set the weights for the score_type for ramping
void TrialMover::initialize_weights(
                                    Real const start_weight,
                                    Real const end_weight,
                                    core::scoring::ScoreType score_type,
                                    int const ramp_cycles
                                    ) {
    original_weight = mc_->score_function().get_weight( score_type );
    delta = ( end_weight - start_weight ) / ramp_cycles;
    ramp_weight = start_weight;
    start_weight_ = start_weight;
}

void TrialMover::set_mc( MonteCarloOP mc_in ) {
    mc_ = mc_in;
}

/// Devuelve el tamaño de un fichero
int get_file_size(std::string filename) // path to file
{
    FILE *p_file = NULL;
    p_file = fopen(filename.c_str(),"rb");
    fseek(p_file,0,SEEK_END);
    int size = ftell(p_file);
    fclose(p_file);
    return size;
}
///@brief:
///the apply print funcion for stats
///@details:
/// this function show the statics
/// for each cycle.
/// @author: Jean 2P. Príncipe

struct Estadisticas_Acc {
    double acomuladorDeAceptadosNormal;
    double acomuladorDeAceptadosCustom;
    double cont_total_rmsd_vs_actual_acc;
    double media;
    double numApply;
} estadisticasStage4[3];

void TrialMover::setEstadisticasStage4(int index, int numApplys){
    double media = 0;
    double normal = 0;
    double custom = 0;
    double totalRmsdVsActual = 0;
    if (cont_total_rmsd_vs_actual_acc > 0) {
        totalRmsdVsActual = cont_total_rmsd_vs_actual_acc;
        media = (rmsd_vs_actual_acc / cont_total_rmsd_vs_actual_acc);
    }else{
        media = 0;
    }
    if (acomuladorDeAceptadosNormal > 0){
        normal = (acomuladorDeAceptadosNormal);
    } else {
        normal = 0;
    }
    if (acomuladorDeAceptadosCustom > 0){
        custom = (acomuladorDeAceptadosCustom);
    } else {
        custom = 0;
    }

    estadisticasStage4[index].acomuladorDeAceptadosNormal = normal;
    estadisticasStage4[index].acomuladorDeAceptadosCustom = custom;
    estadisticasStage4[index].media = media;
    estadisticasStage4[index].cont_total_rmsd_vs_actual_acc = totalRmsdVsActual;
    estadisticasStage4[index].numApply = numApplys;

}
void TrialMover::imprimirEstadisticasStage4(int numPdb){
    int tamanoStage4 = 3;
    double media = 0;
    double normal = 0;
    double custom = 0;
    double numApplys = 0;
    double totalRmsdVsActual = 0;
    std::ofstream myfile;
    myfile.open ("salida_4.txt", std::ios::out | std::ios::app );
    if (myfile.is_open())
    {
        myfile << "============================================" << std::endl;
        myfile << "================ FINAL STATS "<< numPdb<<" ===============" << std::endl;
        myfile << "================ STAGE - 4 =================" << std::endl;
        
        for (int index=0; index < tamanoStage4; index++) {
            media = estadisticasStage4[index].cont_total_rmsd_vs_actual_acc;
            totalRmsdVsActual = estadisticasStage4[index].cont_total_rmsd_vs_actual_acc;
            normal = estadisticasStage4[index].acomuladorDeAceptadosNormal;
            custom = estadisticasStage4[index].acomuladorDeAceptadosCustom;
            numApplys = estadisticasStage4[index].numApply;
            myfile << "================ ITERATE - " << (index + 1) <<" ===============" << std::endl;
            if (totalRmsdVsActual > 0) {
                myfile << "checkeos totales rmsd vs actual: " << totalRmsdVsActual << std::endl;
                myfile << "media del rmsd vs actual: " << media << std::endl;
            }else{
                myfile << "media del rmsd vs actual: 0" << std::endl;
            }
            myfile <<"Número de veces que se llama Apply: " << numApplys<< std::endl;
            if (normal > 0){
                myfile <<"Porcentaje total de aceptados (Normal): "<< (normal * 100)/numApplys<<"%" << std::endl;
                myfile <<"Numero total de aceptados (Normal): "<< normal << std::endl;
            } else {
                myfile <<"Porcentaje total de aceptados (Normal): 0%"<< std::endl;
                myfile <<"Numero total de aceptados (Normal): 0" << std::endl;
            }
            if (custom > 0){
                myfile <<"Porcentaje total de aceptados (Custom): "<< (custom*100)/numApplys <<"%" << std::endl;
                myfile <<"Numero total de aceptados (Custom): "<< custom << std::endl;
            } else {
                myfile <<"Porcentaje total de aceptados (Custom): 0%"<< std::endl;
                myfile <<"Numero total de aceptados (Custom): 0" << std::endl;
            }
        }
        myfile << "============================================" << std::endl;
    }else {
        std::cout << "Unable to open file";
    }

}


void TrialMover::imprimir_estadisticas(int numApplys, int stage, int numPdb)
{
    std::ofstream myfile;
    myfile.open ("salida_"+std::to_string(stage)+".txt", std::ios::out | std::ios::app );
    
    if (myfile.is_open())
    {
        myfile << "============================================" << std::endl;
        myfile << "================ FINAL STATS "<< numPdb << " ===============" << std::endl;
        myfile << "================ STAGE - " << stage <<" =================" << std::endl;
        
        if (cont_total_rmsd_vs_actual_acc > 0) {
            myfile << "checkeos totales rmsd vs actual: " << (cont_total_rmsd_vs_actual_acc) << std::endl;
            myfile << "media del rmsd vs actual: " << (rmsd_vs_actual_acc / cont_total_rmsd_vs_actual_acc) << std::endl;
        }else{
            myfile << "media del rmsd vs actual: 0" << std::endl;
        }
        myfile <<"Número de veces que se llama Apply: "<< numApplys << std::endl;
        std::cout << "NUM APPLYS: " << numApplys << " Aceptados-Normal "<< acomuladorDeAceptadosNormal <<" Aceptados-Custom "<< acomuladorDeAceptadosCustom << std::endl;
        if (acomuladorDeAceptadosNormal > 0){
            myfile <<"Porcentaje total de aceptados (Normal): " << (acomuladorDeAceptadosNormal * 100)/numApplys <<"%" << std::endl;
             myfile <<"Numero total de aceptados (Normal): " << acomuladorDeAceptadosNormal<< std::endl;
        } else {
            myfile <<"Porcentaje total de aceptados (Normal): 0%" << std::endl;
            myfile <<"Numero total de aceptados (Normal): 0" << std::endl;
        }
        if (acomuladorDeAceptadosCustom > 0){
            myfile <<"Porcentaje total de aceptados (Custom): " << (acomuladorDeAceptadosCustom * 100)/numApplys <<"%" << std::endl;
            myfile <<"Numero total de aceptados (Custom): " << acomuladorDeAceptadosCustom << std::endl;
        } else {
            myfile <<"Porcentaje total de aceptados (Custom): 0%" << std::endl;
            myfile <<"Numero total de aceptados (Custom): 0" << std::endl;
        }
        myfile << "============================================" << std::endl;
    }else {
        std::cout << "Unable to open file";
    }

    
}

void TrialMover::resetAcomuladores(){
    acomuladorDeAceptadosNormal = 0;
    acomuladorDeAceptadosCustom = 0;
    cont_total_rmsd_vs_actual_acc = 0;
    countApplys = 0;

}


/// @brief:
///  the apply function for a trial
/// @details:
///  the trial object is created with an mc object
///  the mover is applied before doing an mc.boltzmann
///
/// @author: Monica Berrondo
void TrialMover::apply( pose::Pose & pose )
{
    using scoring::total_score;
    
    /// get the initial scores
    if ( keep_stats_type() == all_stats ) {
        stats_.add_score( mc_->last_accepted_score() ); ///< initial_last_accepted_score
        stats_.add_score( pose.energies().total_energy() ); ///< initial_pose_score
    }
    core::pose::Pose pose_anterior = pose;
    /// make the move
    mover_->apply( pose );
    
    // if ( keep_stats_type() == all_stats ) { //// score and get residue energies
    // Stupid and wasteful.  The structure will be scored inside mc_->boltzman.  mc_->score_function()( pose );
    // Unneccessary since revision 23846 --- mc_->score_function().accumulate_residue_total_energies( pose );
    // WAIT FOR IT. stats_.add_score( pose.energies().total_energy() ); ///< score_after_move
    // }
    using namespace core;
    using namespace core::import_pose;
    using namespace pose;
    std::vector<PoseOP>::iterator it_pose;
    /// test if MC accepts or rejects it
    //  CODIGO ANTERIOR:  bool accepted_move = mc_->boltzmann( pose, mover_->type() );
    
    
    bool accepted_move = false;
    int inicio;
    bool reemplazo_rechazado = false;
    std::vector<std::string>::iterator it;
    //std::cout << umbralLimite<< std::endl;
    //Contador de Applys
    countApplys++;
    std::cout << "===== Satage " << "- " << stage << " =====" << std::endl;
    std::cout << "Umbral límite: " << umbral_apply << " soluciones_anteriores "<< soluciones_anteriores.size()<< std::endl;
    std::string fase4 = "Stage 4";
    std::cout << "==== apply fragment boltzmann ==== " << "Umbral límite: " << umbral_apply << "===== Satage " << "- " << stage<< std::endl;
    if(umbral_apply != 0){
        
//        if ((!stage.compare(fase4)) && (soluciones_anteriores.size() > 0)) {
//             std::cout << "===== DENTRO " << "- " << stage << " =====" << std::endl;
//            mc_->set_temperature(1000000000); //aumentamos la temperatura
//            accepted_move = mc_->boltzmann( pose, mover_->type() );
//        }else{
//            accepted_move = mc_->boltzmann( pose, mover_->type() );
//        }
        
//        if (accepted_move == 1) {
//            acomuladorDeAceptadosNormal += 1;
//        }
        
        if (soluciones_anteriores.size() > 0) {
            accepted_move = mc_->boltzmann( pose, mover_->type() );
            std::cout << "boltzmann 1 " << std::endl;
            if (accepted_move == 1) {
                acomuladorDeAceptadosNormal += 1;
            }
            if (boost::numeric_cast<int>(soluciones_anteriores.size()) > 10){
                inicio = boost::numeric_cast<int>(soluciones_anteriores.size()) - 10;
            }else{
                inicio = 0;
            }
            std::cout << "Umbral (>0): " << umbral_apply << " soluciones_anteriores "<< soluciones_anteriores.size()<< " Estado " << stage<<std::endl;
            
            for(int idx_pose = inicio; idx_pose < boost::numeric_cast<int>(soluciones_anteriores.size()) && !reemplazo_rechazado; idx_pose++){
                PoseOP current_pose = soluciones_anteriores[idx_pose];

                //core::Real SMD_vs_actual = calculo_smd->current_distance_calculation(*current_pose, pose);
                
                core::Real rmsd_vs_actual = core::scoring::CA_rmsd(*current_pose, pose);
                //std::cout << "VALOR SMD " << SMD_vs_actual << std::endl;
                std::cout << "FASE-"<< stage << " VALOR RMSD " << rmsd_vs_actual <<std::endl;
                //rmsd_vs_actual_acc = rmsd_vs_actual_acc + SMD_vs_actual;
                rmsd_vs_actual_acc = rmsd_vs_actual_acc + rmsd_vs_actual;
                cont_total_rmsd_vs_actual_acc += 1;

                //if (accepted_move == 1 && SMD_vs_actual < umbral_apply) {
                if (accepted_move == 1 && rmsd_vs_actual < umbral_apply) {
                    reemplazo_rechazado = true;
                    break;
                }
            }
            //Solo se acepta el reemplazo si todas las soluciones son mayor
            //que el actual
            if (accepted_move == 1 && reemplazo_rechazado == false) {
                acomuladorDeAceptadosCustom +=1;
                std::cout << "Acomulador custom " << acomuladorDeAceptadosCustom << " Umbral apply: " << umbral_apply << std::endl;
            }else {
                std::cout<< "deshace boltzmann "<< std::endl;
                pose = pose_anterior;
            }
        }else{//empieza el caso soluciones-anteriores =0; Rosetta
            accepted_move = mc_->boltzmann( pose, mover_->type() );
            std::cout << "boltzmann 2 " << std::endl;
            if (accepted_move == 1) {
                acomuladorDeAceptadosNormal += 1;
            }
        }// fin para el caso soluciones_anteriores = 0
        
    } else { // empieza el caso umbral == 0; Rosetta
        accepted_move = mc_->boltzmann( pose, mover_->type() );
        std::cout << "boltzmann 3 " << std::endl;
        if (accepted_move == 1) {
            acomuladorDeAceptadosNormal += 1;
        }
    }// fin (umbral_apply > 0)
    
    if ( keep_stats_type() == all_stats ) {
        stats_.add_score( mc_->total_score_of_last_considered_pose() );
    }
    
    if ( stats_type_ <= accept_reject ) {
        stats_.accepted( accepted_move );
    }
    if ( keep_stats_type() > no_stats ) {
        stats_.print( mc_, mover_->type() );
    }
}

std::string
TrialMover::get_name() const {
    return "TrialMover";
}

protocols::moves::MoverOP
TrialMover::clone() const
{
    return utility::pointer::make_shared< TrialMover >(*this);
}

protocols::moves::MoverOP
TrialMover::fresh_instance() const
{
    return utility::pointer::make_shared< TrialMover >();
}

Real TrialMover::acceptance_rate() const
{
    tr << "Acceptance rate: " << stats_.acceptance_rate() << std::endl;
    return stats_.acceptance_rate();
}

int TrialMover::num_accepts() const
{
    return stats_.num_accepted();
}

// sets the input pose also for the contained mover (barak)
void TrialMover::set_input_pose( PoseCOP pose )
{
    this->Mover::set_input_pose( pose );
    if ( mover_ ) {
        mover_->set_input_pose( pose );
    }
}

// sets the native pose also for the contained mover (barak)
void TrialMover::set_native_pose( PoseCOP pose )
{
    this->Mover::set_native_pose( pose );
    if ( mover_ ) {
        mover_->set_native_pose( pose );
    }
}

void TrialMover::parse_my_tag(
                              TagCOP const tag,
                              basic::datacache::DataMap & data,
                              Filters_map const &,
                              Movers_map const & movers,
                              Pose const &
                              )
{
    // 1. MonteCarlo object's name
    std::string const mc_name( tag->getOption< std::string > ( "montecarlo", "" ));
    if ( mc_name == "" ) {
        throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,  "TrialMover requires the 'montecarlo' option which was not provided" );
    }
    
    // 2. Mover
    std::string const movername( tag->getOption< std::string > ( "mover", "" ));
    if ( movername == "" ) {
        throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,  "TrialMover requires the 'mover' option which was not provided" );
    }
    auto  find_mover ( movers.find( movername ));
    if ( find_mover == movers.end() && movername != "" ) {
        throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,  "TrialMover was not able to find the mover named '" + movername + "' in the Movers_map" );
    }
    mover_ = find_mover->second;
    mc_ = data.get_ptr< protocols::moves::MonteCarlo>( "montecarlos", mc_name );
    
    // 3. stats_type.
    std::string const statstype( tag->getOption< std::string > ( "keep_stats", "no_stats" ));
    if ( statstype != "no_stats" && statstype != "accept_reject" && statstype != "all_stats" ) {
        throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,  "TrialMover keep_stats may only be given the values:\n'no_stats', 'accept_reject', and 'all_stats'.\nRead value '" + statstype + "'" );
    }
    if ( statstype == "no_stats" ) {
        stats_type_ = no_stats;
    } else if ( statstype == "accept_reject" ) {
        stats_type_ = accept_reject;
    } else {
        stats_type_ = all_stats;
    }
    
}

std::ostream &operator<< (std::ostream &os, TrialMover const &mover)
{
    //moves::operator<<(os, mover); // this line causes unexpected segmentation fault.
    os << "Mover name: " << mover.get_name() << ", Mover type: " << mover.get_type() << ", Mover current tag: " << mover.get_current_tag() << std::endl <<
    "Mover being tried:   " << mover.mover() << std::endl <<
    "Moves were accepted: " << mover.num_accepts() << " times." << std::endl <<
    "Acceptance rate:     " << mover.stats_.acceptance_rate() << std::endl;
    os << "MonteCarlo:          ";
    if ( mover.mc_ != nullptr ) { os << mover.mc_ << std::endl; }
    else { os << "none" << std::endl; }
    
    return os;
}

}  // namespace moves
}  // namespace protocols

