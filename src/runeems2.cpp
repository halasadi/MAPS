
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
        
        EEMS2 eems2(params);
        
        MCMC mcmc(params);
        
        boost::filesystem::path dir(eems2.prevpath().c_str());
        if (exists(dir)) {
            cerr << "Load final EEMS2 state from " << eems2.prevpath() << endl << endl;
            eems2.load_final_state();
        } else {
            cerr << "Initialize EEMS2 random state" << endl << endl;
            eems2.initialize_state();
        }
        
        error = eems2.start_eems(mcmc);
        if (error) {
            cerr << "[RunEEMS2] Error starting EEMS2." << endl;
            return(EXIT_FAILURE);
        }
        
        
        Proposal proposal;
        
        
        double Temperature = 10;
        double s = 0.5;
        Draw draw; // Random number generator
	
        // now adjust to multiple chains 6/12/2016
        while (!mcmc.finished) {
            
            if (draw.runif() < s) {
                
                proposeMove(eems2, proposal, mcmc);
                mcmc.add_to_total_moves(proposal);
                if (eems2.accept_proposal(proposal, 1)) { mcmc.add_to_okay_moves(proposal.move); }
                
                eems2.update_hyperparams( );
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
                
                double loga = eems2_hotchain.nowll - eems2_coldchain.nowll + (eems2_coldchain.nowll - eems2_hotchain.nowll) * (1/Temperature);
                double u = draw.runif();
                if (log(u) < min(0.0,loga)){
                    // transfer parameters
                }
                
                
            }
      
            
        }
        error = eems2.output_results(mcmc_coldchain);
        
    } catch(exception& e) {
        cerr << e.what() << endl;
        return(EXIT_FAILURE);
    }
    
    return(0);
}
