
#include "eems2.hpp"
#include "draw.hpp"

// The distance metric is a global variable, so that
// the pairwise_distance function can see it
// Choose 'euclidean' (the default) or 'greatcirc' (great circle distance)
string dist_metric;


void proposeMove(EEMS2 &eems2, Proposal &proposal, int chainIndex){
    switch ( eems2.choose_move_type( ) ) {
        case Q_VORONOI_BIRTH_DEATH:
            eems2.propose_birthdeath_qVoronoi(proposal, chainIndex);
            break;
        case M_VORONOI_BIRTH_DEATH:
            eems2.propose_birthdeath_mVoronoi(proposal, chainIndex);
            break;
        case Q_VORONOI_POINT_MOVE:
            eems2.propose_move_one_qtile(proposal, chainIndex);
            break;
        case M_VORONOI_POINT_MOVE:
            eems2.propose_move_one_mtile(proposal, chainIndex);
            break;
        case Q_VORONOI_RATE_UPDATE:
            eems2.propose_rate_one_qtile(proposal, chainIndex);
            break;
        case M_VORONOI_RATE_UPDATE:
            eems2.propose_rate_one_mtile(proposal, chainIndex);
            break;
        case M_MEAN_RATE_UPDATE:
            eems2.propose_overall_mrate(proposal, chainIndex);
            break;
        case Q_MEAN_RATE_UPDATE:
            eems2.propose_overall_qrate(proposal, chainIndex);
            break;
        case DF_UPDATE:
            eems2.propose_df(proposal, chainIndex);
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
        
        MCMC mcmc_coldchain(params);
        
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
        double s = 0.5;
        const int nChains = eems2.getnChains();
        Draw draw; // Random number generator
        
        while (!mcmc_coldchain.finished) {
            // update the cold chains
            for (int chain = 0; chain < nChains; chain++){
                
                // re-write code.. a little confusing
                double u = draw.runif();
                if (chain == (nChains-1)){
                    u = 1.1;
                }
                
                // update parameters
                if (u < s){
                    proposeMove(eems2, proposal, chain);
                    if (eems2.accept_proposal(chain)) {
                        // if it's the cold chain, save and print
                        if (chain == 0){
                            mcmc_coldchain.add_to_okay_moves(proposal.move);
                            if (iter>=0) {
                                eems2.print_iteration(chain);
                                eems2.save_iteration(chain);
                                eems2.writePopSizes(chain);
                            }
                        }
                        
                    }
                    eems2.update_hyperparams(chain);
                
                    if (error) { cerr << "[RunEEMS2] Error saving eems results to " << eems2.mcmcpath() << endl; }
                    
                 // transfer parameters from the hot(ter) chain to the cold(er) chain
                } else{
                    double temperature = eems2.getTemperature(chain);
                    double loga = ((temperature-1)/temperature) * (eems2.getll(chain+1) - eems2.getll(chain));
                    double u = draw.runif();
                    if (log(u) < min(0.0,loga)){
                        // transfer parameters
                        eems2.transferChain(chain+1, chain);
                    }
                    
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
