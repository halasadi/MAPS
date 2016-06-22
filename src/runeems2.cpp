
#include "eems2.hpp"
#include "draw.hpp"

// The distance metric is a global variable, so that
// the pairwise_distance function can see it
// Choose 'euclidean' (the default) or 'greatcirc' (great circle distance)
string dist_metric;


void proposeMove(EEMS2 &eems2, Proposal &proposal, MCMC &mcmc){
    switch ( eems2.choose_move_type( ) ) {
        case Q_VORONOI_BIRTH_DEATH:
            eems2.propose_birthdeath_qVoronoi(proposal);
            break;
        case M_VORONOI_BIRTH_DEATH:
            eems2.propose_birthdeath_mVoronoi(proposal);
            break;
        case Q_VORONOI_POINT_MOVE:
            eems2.propose_move_one_qtile(proposal);
            break;
        case M_VORONOI_POINT_MOVE:
            eems2.propose_move_one_mtile(proposal);
            break;
        case Q_VORONOI_RATE_UPDATE:
            eems2.propose_rate_one_qtile(proposal);
            break;
        case M_VORONOI_RATE_UPDATE:
            eems2.propose_rate_one_mtile(proposal);
            break;
        case M_MEAN_RATE_UPDATE:
            eems2.propose_overall_mrate(proposal);
            break;
        case Q_MEAN_RATE_UPDATE:
            eems2.propose_overall_qrate(proposal);
            break;
        case DF_UPDATE:
            eems2.propose_df(proposal,mcmc);
            break;
        default:
            cerr << "[RunEEMS2] Unknown move type" << endl;
            exit(EXIT_FAILURE);
    }
    
}

int main(int argc, char** argv)
{
    try {
        
        long seed_from_command_line = 1L;
        string params_file; bool error;
        
        po::options_description options("EEMS2 options from command line");
        po::variables_map vm;
        options.add_options()
        ("help", "Produce this help message")
        ("seed", po::value<long>(&seed_from_command_line)->default_value(time(NULL)), "Set the random seed")
        ("params", po::value<string>(&params_file)->required(), "Specify input parameter file") ;
        
        po::store(po::parse_command_line(argc, argv, options), vm);
        po::notify(vm);
        
        if(vm.count("help")) {
            cerr << options << endl; return(EXIT_FAILURE);
        }
        
        Params params(params_file,seed_from_command_line);
        error = params.check_input_params( );
        if (error) {
            cerr << "[RunEEMS2] Error parametrizing EEMS2." << endl;
            return(EXIT_FAILURE);
        }
        
        // Specify the distance metric in the params.ini file
        dist_metric = params.distance;
        
        EEMS2 eems2_coldchain(params);
        MCMC mcmc_coldchain(params);
        
        EEMS2 eems2_hotchain(params);
        MCMC mcmc_hotchain(params);
        
        boost::filesystem::path dir(eems2_coldchain.prevpath().c_str());
        if (exists(dir)) {
            cerr << "Load final EEMS2 state from " << eems2_coldchain.prevpath() << endl << endl;
            eems2_coldchain.load_final_state();
            eems2_hotchain.load_final_state();
        } else {
            cerr << "Initialize EEMS2 random state" << endl << endl;
            eems2_coldchain.initialize_state();
            eems2_hotchain.initialize_state();
        }
        
        error = eems2_coldchain.start_eems(mcmc_coldchain);
        eems2_hotchain.start_eems(mcmc_hotchain);
        if (error) {
            cerr << "[RunEEMS2] Error starting EEMS2." << endl;
            return(EXIT_FAILURE);
        }
        
        Proposal proposal_coldchain;
        Proposal proposal_hotchain;
        
        
        double Temperature = 2;
        double s = 0.5;
        Draw draw; // Random number generator

        
        while (!mcmc_coldchain.finished) {
            
            if (draw.runif() < s) {
                
                // update both hot and cold chain
                
                proposeMove(eems2_coldchain, proposal_coldchain, mcmc_coldchain);
                mcmc_coldchain.add_to_total_moves(proposal_coldchain.move);
                if (eems2_coldchain.accept_proposal(proposal_coldchain, 1)) { mcmc_coldchain.add_to_okay_moves(proposal_coldchain.move); }
                
                eems2_coldchain.update_hyperparams( );
                mcmc_coldchain.end_iteration( );
                
                // Check whether to save the current parameter state,
                // as the thinned out iterations are not saved
                // only save iterations from cold chain
                int iter = mcmc_coldchain.to_save_iteration( );
                if (iter>=0) {
                    eems2_coldchain.print_iteration(mcmc_coldchain);
                    eems2_coldchain.save_iteration(mcmc_coldchain);
                    eems2_coldchain.writePopSizes();
                }
                if (error) { cerr << "[RunEEMS2] Error saving eems results to " << eems2_coldchain.mcmcpath() << endl; }
            
                proposeMove(eems2_hotchain, proposal_hotchain, mcmc_hotchain);
                mcmc_hotchain.add_to_total_moves(proposal_hotchain.move);
                if (eems2_hotchain.accept_proposal(proposal_hotchain, Temperature)) { mcmc_hotchain.add_to_okay_moves(proposal_hotchain.move); }
                
                eems2_hotchain.update_hyperparams( );
                mcmc_hotchain.end_iteration( );
                
            } else{
                
                // hot-chain parameters transfers to cold-chain
                
                double loga = eems2_hotchain.getLogPosterior() + (eems2_coldchain.getLogPosterior())/Temperature - eems2_coldchain.getLogPosterior() - (1/Temperature)*eems2_hotchain.getLogPosterior();
                cout << exp(loga) << endl;
                double u = draw.runif();
                if (log(u) < min(0.0,loga)){
                    cout << "Making a switch" << endl;
                    eems2_coldchain.nowmtiles = eems2_hotchain.nowmtiles;
                    eems2_coldchain.nowqtiles = eems2_hotchain.nowqtiles;
                    eems2_coldchain.nowmSeeds = eems2_hotchain.nowmSeeds;
                    eems2_coldchain.nowmEffcts = eems2_hotchain.nowmEffcts;
                    eems2_coldchain.nowmrateMu = eems2_hotchain.nowmrateMu;
                    eems2_coldchain.nowqSeeds = eems2_hotchain.nowqSeeds;
                    eems2_coldchain.nowqEffcts = eems2_hotchain.nowqEffcts;
                    eems2_coldchain.nowqrateS2 = eems2_hotchain.nowqrateS2;
                    eems2_coldchain.nowmrateS2 = eems2_hotchain.nowmrateS2;
                    eems2_coldchain.nowqrateMu = eems2_hotchain.nowqrateMu;
                    eems2_coldchain.nowpi = eems2_hotchain.nowpi;
                    eems2_coldchain.nowll = eems2_hotchain.nowll;
                    eems2_coldchain.nowdf = eems2_hotchain.nowdf;
                    eems2_coldchain.nowqColors = eems2_hotchain.nowqColors;
                    eems2_coldchain.nowmColors = eems2_hotchain.nowmColors;
                    
                }
                
                
            }
      
            
        }
        error = eems2_coldchain.output_results(mcmc_coldchain);
        
    } catch(exception& e) {
        cerr << e.what() << endl;
        return(EXIT_FAILURE);
    }
    
    return(0);
}
