
#include "eems2.hpp"
#include "draw.hpp"

// The distance metric is a global variable, so that
// the pairwise_distance function can see it
// Choose 'euclidean' (the default) or 'greatcirc' (great circle distance)
string dist_metric;


void proposeMove(EEMS2 &eems2, Proposal &proposal, MCMC &mcmc, int chain_no){
    switch ( eems2.choose_move_type(chain_no) ) {
        case TRANSFER_FROM_HOT_CHAIN:
            proposal.move = TRANSFER_FROM_HOT_CHAIN;
            break;
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
        eems2.initialize_state();
        
        if (error) {
            cerr << "[RunEEMS2] Error starting EEMS2." << endl;
            return(EXIT_FAILURE);
        }
        
        Proposal proposal;
        
        double s = 0.5;
        Draw draw; // Random number generator
        
        // start with running the hottest chain
        
        int nChains = 5;
        double Temperatures [5] = { 50, 25, 12, 5, 1};
        MCMC mcmc(params);
        error = eems2.start_eems(mcmc.num_iters_to_save( ) );
        

        for (int chain_no = 0; chain_no < nChains; chain_no++){
            
            double Temperature = Temperatures[chain_no];
            eems2.setTemperature(Temperature);
            mcmc.restart(params, Temperature);
            
            
            if (chain_no > 0){
                eems2.prev_stored_mcmc_states = eems2.now_stored_mcmc_states;
                eems2.now_stored_mcmc_states.clear();
            }
            
            cout << "chain no :" << chain_no << endl;
            
            while (!mcmc.finished) {
                
                proposeMove(eems2, proposal, mcmc, chain_no);

                mcmc.add_to_total_moves(proposal.move);
                
                if (proposal.move == TRANSFER_FROM_HOT_CHAIN){

                    if (eems2.accept_transfer(proposal, Temperatures[chain_no], Temperatures[chain_no-1])){
                        mcmc.add_to_okay_moves(proposal.move);
                    }
                } else{
                    
                    if (eems2.accept_proposal(proposal)) { mcmc.add_to_okay_moves(proposal.move); }
                }
                mcmc.end_iteration( );
                eems2.update_hyperparams( );
                
                int iter = mcmc.to_save_iteration( );
                if (iter>=0) {
                    eems2.print_iteration(mcmc);
                    eems2.store_iteration();
                }
                
                // Check whether to save the current parameter state,
                // as the thinned out iterations are not saved
                // only save iterations from cold chain
                if (Temperature == 1 && iter >= 0 ){
                    eems2.save_iteration(mcmc);
                    eems2.writePopSizes();
                    
                }
                
                
            }
            
            eems2.print_iteration(mcmc);
            
            if (Temperature == 1){
                error = eems2.output_results(mcmc);
                if (error) { cerr << "[RunEEMS2] Error saving eems results to " << eems2.mcmcpath() << endl; }
                
            }
            
        }
        
        
    } catch(exception& e) {
        cerr << e.what() << endl;
        return(EXIT_FAILURE);
    }
    
    return(0);
    
}
    
