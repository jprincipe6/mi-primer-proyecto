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

#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>
#include <basic/datacache/DataMap.hh>

// Utility headers

// C++ headers
#include <string>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/rms_util.hh>
#include <utility/vector0.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>
#include <vector>
#include <sys/types.h>
#include <dirent.h>

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
const char *path_input = "/Users/principe/Documents/Rosetta/rosetta_bin_mac_2019.35.60890_bundle/demos/public/abinitio/soluciones_1elwA";
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
    std::vector<PoseOP>::iterator it_pose;
    std::vector<std::string>::iterator it;
    
    mover_ = mover_in;
    mc_ = mc_in;
    
    
    
    rmsd_vs_actual_acc = 0;
    cont_total_rmsd_vs_actual_acc = 0;
    paths_soluciones_pdbs = get_paths_pdbs_from_dir(path_input);
    acomuladorDeAceptadosCustom = 0;
    acomuladorDeAceptadosNormal = 0;
    
    for (it= paths_soluciones_pdbs.begin(); it < paths_soluciones_pdbs.end(); it++) {
        std::string path_file_pdb = std::string (path_input) + std::string("/") +std::string(*it);
        PoseOP ejecucion_previa = pose_from_file(path_file_pdb);
        soluciones_anteriores.push_back(ejecucion_previa);
    }
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
    double acumuladorNumApplys;
};

/* POSIBILIDAD UTILIZANDO UN ACUMULADOR EXTERNO EN CLASSIC AB INITIO
 Estadisticas_Acc TrialMover::imprimir_estadisticas(Estadisticas_Acc stats_anteriores)
 {
 acomuladorDeAceptadosNormal = acomuladorDeAceptadosNormal + stats_anteriores.acomuladorDeAceptadosNormal;
 
 std::cout << "============================================" << std::endl;
 std::cout << "================ FINAL STATS ===============" << std::endl;
 
 if (cont_total_rmsd_vs_actual_acc > 0) {
 std::cout << "media del rmsd vs actual: " << (rmsd_vs_actual_acc / cont_total_rmsd_vs_actual_acc) << std::endl;
 }
 std::cout <<"Número de veces que se llama Apply: "<< numApplys << std::endl;
 if (acomuladorDeAceptadosNormal > 0){
 std::cout <<"Porcentaje total de aceptados (Normal): " << (acomuladorDeAceptadosNormal * 100)/numApplys <<"%" << std::endl;
 }
 if (acomuladorDeAceptadosCustom > 0){
 std::cout <<"Porcentaje total de aceptados (Custom): " << (acomuladorDeAceptadosCustom * 100)/numApplys <<"%" << std::endl;
 //::cout <<"número total de no aceptados (Custom): " << (soluciones_anteriores.size() - acomuladorDeAceptadosCustom) << std::endl;
 }
 std::cout << "============================================" << std::endl;
 
 Estadisticas_Acc actuales;
 actuales.acomuladorDeAceptadosNormal = acomuladorDeAceptadosNormal;
 return actuales;
 }*/

void TrialMover::imprimir_estadisticas(int numApplys, int stage)
{
    
    std::cout << "============================================" << std::endl;
    std::cout << "================ FINAL STATS ===============" << std::endl;
    std::cout << "================ STAGE - " << stage <<" =================" << std::endl;
    
    if (cont_total_rmsd_vs_actual_acc > 0) {
        std::cout << "checkeos totales rmsd vs actual: " << (cont_total_rmsd_vs_actual_acc) << std::endl;
        std::cout << "media del rmsd vs actual: " << (rmsd_vs_actual_acc / cont_total_rmsd_vs_actual_acc) << std::endl;
    }else{
        std::cout << "media del rmsd vs actual: 0" << std::endl;
    }
    std::cout <<"Número de veces que se llama Apply: "<< numApplys << std::endl;
    if (acomuladorDeAceptadosNormal > 0){
        std::cout <<"Porcentaje total de aceptados (Normal): " << (acomuladorDeAceptadosNormal * 100)/numApplys <<"%" << std::endl;
    } else {
        std::cout <<"Porcentaje total de aceptados (Normal): 0%" << std::endl;
    }
    if (acomuladorDeAceptadosCustom > 0){
        std::cout <<"Porcentaje total de aceptados (Custom): " << (acomuladorDeAceptadosCustom * 100)/numApplys <<"%" << std::endl;
    } else {
        std::cout <<"Porcentaje total de aceptados (Custom): 0%" << std::endl;
    }
    std::cout << "============================================" << std::endl;
    
}

void TrialMover::resetAcomuladores(){
    acomuladorDeAceptadosNormal = 0;
    acomuladorDeAceptadosCustom = 0;
    cont_total_rmsd_vs_actual_acc = 0;

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
    
    
    core::Real umbral_limite = 100;
    bool accepted_move = false;
    std::vector<std::string>::iterator it;
        
    if (soluciones_anteriores.size() > 0){
        
        accepted_move = mc_->boltzmann( pose, mover_->type() );
        bool reemplazo_rechazado = false;
        core::pose::Pose pose_anterior = pose;
        
        if (accepted_move == 1) {
            acomuladorDeAceptadosNormal += 1;
        }
        
        for(it_pose = soluciones_anteriores.begin(); it_pose < soluciones_anteriores.end() && !reemplazo_rechazado; it_pose++){
            //Para cada POSE de entrada se calcula la distancia RMSD a la actual
            core::Real rmsd_vs_actual = core::scoring::CA_rmsd(**it_pose, pose);
            rmsd_vs_actual_acc = rmsd_vs_actual_acc + rmsd_vs_actual;
            cont_total_rmsd_vs_actual_acc += 1;
            
            if (accepted_move == 1 && rmsd_vs_actual < umbral_limite) {
                reemplazo_rechazado = true;
            }
        }
        //Solo se acepta el reemplazo si todas las soluciones son mayor
        //que el actual
        if (reemplazo_rechazado == false) {
            acomuladorDeAceptadosCustom +=1;
        }else {
            pose = pose_anterior;
            std::cout << "Valor del normal: " << acomuladorDeAceptadosNormal<< std::endl;
        }
    } else {
        accepted_move = mc_->boltzmann( pose, mover_->type() );
        if (accepted_move == 1) {
            acomuladorDeAceptadosNormal += 1;
            std::cout << "No hay soluciones anteriores (.pdb)"<< std::endl;
        }
        
    }
    
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

