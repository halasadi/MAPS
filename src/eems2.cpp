#include "eems2.hpp"

EEMS2::EEMS2(const Params &params) {
    this->params = params;
    draw.initialize(params.seed);
    habitat.generate_outer(params.datapath);
    habitat.dlmwrite_outer(params.mcmcpath);
    graph.generate_grid(params.datapath,params.gridpath,
                        habitat,params.nDemes,params.nIndiv);
    graph.dlmwrite_grid(params.mcmcpath);
    
    
    o = graph.get_num_obsrv_demes();
    d = graph.get_num_total_demes();
    n = params.nIndiv;
    initialize_sims();
}
EEMS2::~EEMS2( ) { }
string EEMS2::datapath( ) const { return params.datapath; }
string EEMS2::mcmcpath( ) const { return params.mcmcpath; }
string EEMS2::prevpath( ) const { return params.prevpath; }
string EEMS2::gridpath( ) const { return params.gridpath; }
// Draw points randomly inside the habitat: the habitat is two-dimensional, so
// a point is represented as a row in a matrix with two columns
void EEMS2::randpoint_in_habitat(MatrixXd &Seeds) {
    for (int i = 0 ; i < Seeds.rows() ; i++ ) {
        bool in = false;
        double x,y;
        while (!in) {
            x = habitat.get_xmin() + habitat.get_xspan() * draw.runif();
            y = habitat.get_ymin() + habitat.get_yspan() * draw.runif();
            in = habitat.in_point(x,y);
        }
        Seeds(i,0) = x;
        Seeds(i,1) = y;
    }
}

void EEMS2::rnorm_effects(const double lowerBnd, const double upperBnd, const double rateS2, VectorXd &Effcts) {
    for (int i = 0 ; i < Effcts.rows() ; i++ ) {
        Effcts(i) = draw.rtrnorm(0.0,rateS2,lowerBnd, upperBnd);
    }
}

void EEMS2::initialize_sims( ) {
    cerr << "[Sims::initialize]" << endl;
    MatrixXd sims = readMatrixXd(params.datapath + ".sims");
    MatrixXi Sims = sims.cast <int> ();
    if ((Sims.rows()!=n)||(Sims.cols()!=n)) {
        cerr << "  Error reading similarities matrix " << params.datapath + ".sims" << endl
        << "  Expect a " << n << "x" << n << " matrix of pairwise similarities" << endl; exit(1);
    }
    cerr << "  Loaded similarities matrix from " << params.datapath + ".sims" << endl;
    
    // the number of comparisons for each deme.
    
    
    cMatrix = MatrixXd::Zero(o,o);
    cvec = VectorXd::Zero(o);
    observedIBD = MatrixXd::Zero(o, o);
    maxCnt = Sims.maxCoeff();
    
    // counts the number of IBD segments that are 0, 1, 2, etc.
    // pre-computation
    cClasses = VectorXd::Zero(maxCnt+1);
    
    //observedIBD(i,j) is the sum of number of blocks shared between individuals in deme i and deme j that are greater than u cM.
    
    int demei;
    int demej;
    // all (n choose 2) comparisons at the individual level
    // TO DO: there is a more simple way to fill in cMatrix (or get it from cvec)
    
    
    for ( int i = 0 ; i < n ; i ++ ) {
        for (int j = (i+1); j < n; j++){
            demei = graph.get_deme_of_indiv(i);
            demej = graph.get_deme_of_indiv(j);
            cMatrix(demei, demej) += 1;
            cMatrix(demej, demei) = cMatrix(demei, demej);
            
            observedIBD(demei, demej) += Sims(i,j);
            observedIBD(demej, demei) = observedIBD(demei, demej);
            
            cClasses(Sims(i,j)) += 1;
        }
    }
    
    for ( int i = 0 ; i < n ; i ++ ) {
        cvec(graph.get_deme_of_indiv(i)) += 1;
    }
    
    JtDhatJ = MatrixXd::Zero(o,o);
    expectedIBD = MatrixXd::Zero(o,o);
    
    
}

void EEMS2::initialize_state( ) {
    cerr << "[EEMS2::initialize_state]" << endl;
    
    for (int i = 0; i < nChains; i++){
        
        Chain chain = chains[i];
        
        // first chain is the cold chain
        chain.Temperature = i * 10;
        
        chain.df = params.dfmin + draw.runif() * (params.dfmax - params.dfmin);
        // Initialize the two Voronoi tessellations
        chain.qtiles = draw.rnegbin(2*o,0.5); // o is the number of observed demes
        chain.mtiles = draw.rnegbin(2*o,0.5);

        // Draw the Voronoi centers Coord uniformly within the habitat
        chain.qSeeds = MatrixXd::Zero(chain.qtiles,2); randpoint_in_habitat(chain.qSeeds);
        chain.mSeeds = MatrixXd::Zero(chain.mtiles,2); randpoint_in_habitat(chain.mSeeds);
        chain.mrateS2 = draw.rinvgam(0.5,0.5);
        chain.qrateS2 = draw.rinvgam(0.5,0.5);
        
        // Assign migration rates to the Voronoi tiles
        chain.mrateMu = params.mrateMuLowerBound + draw.runif() * (params.mrateMuUpperBound - params.mrateMuLowerBound);
        chain.qrateMu = params.qrateMuLowerBound + draw.runif() * (params.qrateMuUpperBound - params.qrateMuLowerBound);
        
        // Assign rates to the Voronoi tiles
        chain.qEffcts = VectorXd::Zero(chain.qtiles); rnorm_effects(params.qEffctLowerBound, params.qEffctUpperBound, chain.qrateS2, chain.qEffcts);
        chain.mEffcts = VectorXd::Zero(chain.mtiles); rnorm_effects(params.mEffctLowerBound, params.mEffctUpperBound, chain.mrateS2, chain.mEffcts);
        // Initialize the mapping of demes to qVoronoi tiles
        graph.index_closest_to_deme(chain.qSeeds,chain.qColors);
        // Initialize the mapping of demes to mVoronoi tiles
        graph.index_closest_to_deme(chain.mSeeds,chain.mColors);
        
    }

    cerr << "[EEMS2::initialize_state] Done." << endl << endl;
}
void EEMS2::load_final_state( ) {
    
    // first chain is always the cold chain
    // load parameters into the cold chain
    Chain chain = chains[0]
    
    cerr << "[EEMS2::load_final_state]" << endl;
    MatrixXd tempi; bool error = false;
    tempi = readMatrixXd(params.prevpath + "/lastqtiles.txt");
    if ((tempi.rows()!=1) || (tempi.cols()!=1)) { error = true; }
    chain.qtiles = tempi(0,0);
    tempi = readMatrixXd(params.prevpath + "/lastmtiles.txt");
    if ((tempi.rows()!=1) || (tempi.cols()!=1)) { error = true; }
    chain.mtiles = tempi(0,0);
    cerr << "  EEMS starts with " << chain.qtiles << " qtiles and " << chain.mtiles << " mtiles" << endl;
    tempi = readMatrixXd(params.prevpath + "/lastthetas.txt");
    if ((tempi.rows()!=1) || (tempi.cols()!=1)) { error = true; }
    chain.df = tempi(0,0);
    tempi = readMatrixXd(params.prevpath + "/lastdfpars.txt");
    if ((tempi.rows()!=1) || (tempi.cols()!=2)) { error = true; }
    chain.dfmin = tempi(0,0);
    chain.dfmax = tempi(0,1);
    tempi = readMatrixXd(params.prevpath + "/lastmhyper.txt");
    if ((tempi.rows()!=1) || (tempi.cols()!=2)) { error = true; }
    chain.mrateMu = tempi(0,0);
    chain.mrateS2 = tempi(0,1);
    tempi = readMatrixXd(params.prevpath + "/lastqhyper.txt");
    if ((tempi.rows()!=1) || (tempi.cols()!=2)) { error = true; }
    chain.qrateMu = tempi(0,0);
    chain.qrateS2 = tempi(0,1);
    tempi = readMatrixXd(params.prevpath + "/lastqeffct.txt");
    if ((tempi.rows()!=chain.qtiles) || (tempi.cols()!=1)) { error = true; }
    chain.qEffcts = tempi.col(0);
    tempi = readMatrixXd(params.prevpath + "/lastmeffct.txt");
    if ((tempi.rows()!=chain.mtiles) || (tempi.cols()!=1)) { error = true; }
    chain.mEffcts = tempi.col(0);
    chain.qSeeds = readMatrixXd(params.prevpath + "/lastqseeds.txt");
    if ((chain.qSeeds.rows()!=chain.qtiles) || (chain.qSeeds.cols()!=2)) { error = true; }
    chain.mSeeds = readMatrixXd(params.prevpath + "/lastmseeds.txt");
    if ((chain.mSeeds.rows()!=chain.mtiles) || (chain.mSeeds.cols()!=2)) { error = true; }
    // Initialize the mapping of demes to qVoronoi tiles
    graph.index_closest_to_deme(chain.mSeeds,chain.mColors);
    // Initialize the mapping of demes to mVoronoi tiles
    graph.index_closest_to_deme(chain.qSeeds,chain.qColors);
    if (error) {
        cerr << "  Error loading MCMC state from " << params.prevpath << endl; exit(1);
    }
    cerr << "[EEMS::load_final_state] Done." << endl << endl;
    
    
}
bool EEMS2::start_eems(const MCMC &mcmc) {
    bool error = false;
    
    // The deviation of move proposals is scaled by the habitat range
    params.mSeedsProposalS2x = params.mSeedsProposalS2 * habitat.get_xspan();
    params.mSeedsProposalS2y = params.mSeedsProposalS2 * habitat.get_yspan();
    params.qSeedsProposalS2x = params.qSeedsProposalS2 * habitat.get_xspan();
    params.qSeedsProposalS2y = params.qSeedsProposalS2 * habitat.get_yspan();
    // MCMC draws are stored in memory, rather than saved to disk,
    // so it is important to thin
    
    int niters = mcmc.num_iters_to_save();
    mcmcmhyper = MatrixXd::Zero(niters,2);
    mcmcqhyper = MatrixXd::Zero(niters,2);
    mcmcpilogl = MatrixXd::Zero(niters,2);
    mcmcmtiles = VectorXd::Zero(niters);
    mcmcqtiles = VectorXd::Zero(niters);
    mcmcthetas = VectorXd::Zero(niters);
    mcmcmRates.clear();
    mcmcqRates.clear();
    mcmcxCoord.clear();
    mcmcyCoord.clear();
    mcmcwCoord.clear();
    mcmczCoord.clear();
    
    for (int i = 0; i < nChains; i++){
        Chain chain = chains[i];
        chain.pi = eval_prior(chain.mSeeds,chain.mEffcts,chain.mrateMu,chain.mrateS2,
                           chain.qSeeds, chain.qEffcts, chain.qrateMu, chain.qrateS2,
                           chain.df);
        chain.ll = eems2_likelihood(chain.mSeeds, chain.mEffcts, chain.mrateMu, chain.qSeeds, chain.qEffcts, chain.qrateMu, chain.df, true);
        
    }
    
    Chain cold_chain = chains[0];
    cerr << "Input parameters: " << endl << params << endl
    << "Initial log prior of cold chain: " << cold_chain.pi << endl
    << "Initial log llike of cold chain: " << cold_chain.ll << endl << endl;
    if ((cold_chain.pi==-Inf) || (cold_chain.pi==Inf) || (cold_chain.ll==-Inf) || (cold_chain.ll==Inf)) { error = true; }
    return(error);
}
MoveType EEMS2::choose_move_type( ) {
    double u1 = draw.runif( );
    double u2 = draw.runif( );
    // There are 4 types of proposals:
    // * birth/death (with equal probability)
    // * move a tile (chosen uniformly at random)
    // * update the rate of a tile (chosen uniformly at random)
    // * update the mean migration rate or the mean coalescent rate (with equal probability)
    MoveType move = UNKNOWN_MOVE_TYPE;
    if (u1 < 0.25) {
        // Propose birth/death to update the Voronoi tessellation of the effective diversity,
        // with probability params.qVoronoiPr (which is 0.05 by default). Otherwise,
        // propose birth/death to update the Voronoi tessellation of the effective migration.
        if (u2 < params.qVoronoiPr) {
            move = Q_VORONOI_BIRTH_DEATH;
        } else {
            move = M_VORONOI_BIRTH_DEATH;
        }
    } else if (u1 < 0.5) {
        if (u2 < params.qVoronoiPr) {
            move = Q_VORONOI_POINT_MOVE;
        } else {
            move = M_VORONOI_POINT_MOVE;
        }
    } else if (u1 < 0.75) {
        if (u2 < params.qVoronoiPr) {
            move = Q_VORONOI_RATE_UPDATE;
        } else {
            move = M_VORONOI_RATE_UPDATE;
        }
    } else {
        if (u2 < 0.333) {
            move = M_MEAN_RATE_UPDATE;
        } else if (u2 < 0.6666) {
            move = Q_MEAN_RATE_UPDATE;
        } else{
            move = DF_UPDATE;
        }
    }
    return(move);
}

void EEMS2::propose_df(Proposal &proposal, int chainIndex) {
    
    Chain chain = chains[chainIndex];
    proposal.move = DF_UPDATE;
    proposal.newpi = -Inf;
    proposal.newll = -Inf;
    double newdf = draw.rnorm(chain.df,params.dfProposalS2);
    if ( (newdf>params.dfmin) && (newdf<params.dfmax) ) {
      proposal.newdf = newdf;
      proposal.newpi = eval_prior(chain.mSeeds,chain.mEffcts,chain.mrateMu,chain.mrateS2,
				  chain.qSeeds,chain.qEffcts,chain.qrateMu,chain.qrateS2,
				  newdf);
      proposal.newll = eems2_likelihood(chain.mSeeds, chain.mEffcts, chain.mrateMu, chain.qSeeds, chain.qEffcts, chain.qrateMu, newdf, true);
    }
}

void EEMS2::propose_rate_one_qtile(Proposal &proposal, int chainIndex) {
    Chain chain = chains[chainIndex];
    // Choose a tile at random to update
    int qtile = draw.runif_int(0,chain.qtiles-1);
    // Make a random-walk proposal, i.e., add small offset to current value
    double curqEffct = chain.qEffcts(qtile);
    double newqEffct = draw.rnorm(curqEffct,params.qEffctProposalS2);
    proposal.move = Q_VORONOI_RATE_UPDATE;
    proposal.newqEffcts = chain.qEffcts;
    proposal.newqEffcts(qtile) = newqEffct;
    // The prior distribution on the tile effects is truncated normal
    // So first check whether the proposed value is in range
    // Then update the prior and evaluate the new likelihood
    if ( (newqEffct > params.qEffctLowerBound) && (newqEffct < params.qEffctUpperBound) ) {
        proposal.newpi = eval_prior(chain.mSeeds,chain.mEffcts,chain.mrateMu,chain.mrateS2,
                                    chain.qSeeds,proposal.newqEffcts,chain.qrateMu,chain.qrateS2,
                                    chain.df);
        proposal.newll = eems2_likelihood(chain.mSeeds, chain.mEffcts, chain.mrateMu, chain.qSeeds, proposal.newqEffcts, chain.qrateMu, chain.df, false);
    } else {
        proposal.newpi = -Inf;
        proposal.newll = -Inf;
    }
}
void EEMS2::propose_rate_one_mtile(Proposal &proposal, int chainIndex) {
    // Choose a tile at random to update
    int mtile = draw.runif_int(0,chain.mtiles-1);
    // Make a random-walk proposal, i.e., add small offset to current value
    double curmEffct = chain.mEffcts(mtile);
    double newmEffct = draw.rnorm(curmEffct,params.mEffctProposalS2);
    proposal.move = M_VORONOI_RATE_UPDATE;
    proposal.newmEffcts = chain.mEffcts;
    proposal.newmEffcts(mtile) = newmEffct;
    // The prior distribution on the tile effects is truncated normal
    // So first check whether the proposed value is in range
    // Then update the prior and evaluate the new likelihood
    if ( (newmEffct > params.mEffctLowerBound) && (newmEffct < params.mEffctUpperBound)  ) {
        proposal.newpi = eval_prior(chain.mSeeds,proposal.newmEffcts,chain.mrateMu,chain.mrateS2,
                                    chain.qSeeds,chain.qEffcts,chain.qrateMu,chain.qrateS2,
                                    chain.df);
        proposal.newll = eems2_likelihood(chain.mSeeds, proposal.newmEffcts, chain.mrateMu, chain.qSeeds, chain.qEffcts, chain.qrateMu, chain.df, true);
        
    } else {
        proposal.newpi = -Inf;
        proposal.newll = -Inf;
    }
}
void EEMS2::propose_overall_mrate(Proposal &proposal, int &chainIndex) {
    // Make a random-walk Metropolis-Hastings proposal
    double newmrateMu = draw.rnorm(chain.mrateMu,params.mrateMuProposalS2);
    proposal.move = M_MEAN_RATE_UPDATE;
    proposal.newmrateMu = newmrateMu;
    // If the proposed value is in range, the prior probability does not change
    // as the prior distribution on mrateMu is uniform
    // Otherwise, setting the prior and the likelihood to -Inf forces a rejection
    if (newmrateMu <= params.mrateMuUpperBound && newmrateMu >= params.mrateMuLowerBound) {
        proposal.newpi = chain.pi;
        proposal.newll = eems2_likelihood(chain.mSeeds, chain.mEffcts, proposal.newmrateMu, chain.qSeeds, chain.qEffcts, chain.qrateMu, chain.df, true);

    } else {
        proposal.newpi = -Inf;
        proposal.newll = -Inf;
    }
}


void EEMS2::propose_overall_qrate(Proposal &proposal, int chainIndex) {
    
    Chain chain = chains[chainIndex];
    // Make a random-walk Metropolis-Hastings proposal
    double newqrateMu = draw.rnorm(chain.qrateMu,params.qrateMuProposalS2);
    proposal.move = Q_MEAN_RATE_UPDATE;
    proposal.newqrateMu = newqrateMu;
    // If the proposed value is in range, the prior probability does not change
    // as the prior distribution on qrateMu is uniform
    // Otherwise, setting the prior and the likelihood to -Inf forces a rejection
    if (newqrateMu <= params.qrateMuUpperBound && newqrateMu >= params.qrateMuLowerBound) {
        proposal.newpi = chain.pi;
        proposal.newll = eems2_likelihood(chain.mSeeds, chain.mEffcts, chain.mrateMu, chain.qSeeds, chain.qEffcts, proposal.newqrateMu, chain.df, false);
    } else {
        proposal.newpi = -Inf;
        proposal.newll = -Inf;
    }
}
void EEMS2::propose_move_one_qtile(Proposal &proposal, int chainIndex) {
    // Choose a tile at random to move
    Chain chain = chains[chainIndex];
    int qtile = draw.runif_int(0,chain.qtiles-1);
    // Make a random-walk proposal, i.e., add small offset to current value
    // In this case, there are actually two values -- longitude and latitude
    double newqSeedx = draw.rnorm(chain.qSeeds(qtile,0),params.qSeedsProposalS2x);
    double newqSeedy = draw.rnorm(chain.qSeeds(qtile,1),params.qSeedsProposalS2y);
    proposal.move = Q_VORONOI_POINT_MOVE;
    proposal.newqSeeds = chain.qSeeds;
    proposal.newqSeeds(qtile,0) = newqSeedx;
    proposal.newqSeeds(qtile,1) = newqSeedy;
    if (habitat.in_point(newqSeedx,newqSeedy)) {
        proposal.newpi = chain.pi;
        proposal.newll = eems2_likelihood(chain.mSeeds, chain.mEffcts, chain.mrateMu, proposal.newqSeeds, chain.qEffcts, chain.qrateMu, chain.df, false);
    } else {
        proposal.newpi = -Inf;
        proposal.newll = -Inf;
    }
}
void EEMS2::propose_move_one_mtile(Proposal &proposal, int chainIndex) {
    Chain chain = chains[chainIndex];
    // Choose a tile at random to move
    int mtile = draw.runif_int(0,chain.mtiles-1);
    // Make a random-walk proposal, i.e., add small offset to current value
    // In this case, there are actually two values -- longitude and latitude
    double newmSeedx = draw.rnorm(chain.mSeeds(mtile,0),params.mSeedsProposalS2x);
    double newmSeedy = draw.rnorm(chain.mSeeds(mtile,1),params.mSeedsProposalS2y);
    proposal.move = M_VORONOI_POINT_MOVE;
    proposal.newmSeeds = chain.mSeeds;
    proposal.newmSeeds(mtile,0) = newmSeedx;
    proposal.newmSeeds(mtile,1) = newmSeedy;
    if (habitat.in_point(newmSeedx,newmSeedy)) {
        proposal.newpi = chain.pi;
        proposal.newll = eems2_likelihood(proposal.newmSeeds, chain.mEffcts, chain.mrateMu, chain.qSeeds, chain.qEffcts, chain.qrateMu, chain.df, true);
    } else {
        proposal.newpi = -Inf;
        proposal.newll = -Inf;
    }
}
void EEMS2::propose_birthdeath_qVoronoi(Proposal &proposal, int chainIndex) {
    Chain chain = chains[chainIndex];
    int newqtiles = chain.qtiles,r;
    double u = draw.runif();
    double pBirth = 0.5;
    double pDeath = 0.5;
    proposal.newqEffcts = chain.qEffcts;
    proposal.newqSeeds = chain.qSeeds;
    // If there is exactly one tile, rule out a death proposal
    if ((chain.qtiles==1) || (u<0.5)) { // Propose birth
        if (chain.qtiles==1) { pBirth = 1.0; }
        newqtiles++;
        MatrixXd newqSeed = MatrixXd::Zero(1,2);
        randpoint_in_habitat(newqSeed);
        pairwise_distance(chain.qSeeds,newqSeed).col(0).minCoeff(&r);
        // The new tile is assigned a rate by perturbing the current rate at the new seed
        double chain.qEffct = chain.qEffcts(r);
        double newqEffct = draw.rtrnorm(chain.qEffct,params.qEffctProposalS2,params.qEffctLowerBound, params.qEffctUpperBound);
        insertRow(proposal.newqSeeds,newqSeed.row(0));
        insertElem(proposal.newqEffcts,newqEffct);
        // Compute log(proposal ratio) and log(prior ratio)
        proposal.newratioln = log(pDeath/pBirth)
        - dtrnormln(newqEffct,chain.qEffct,params.qEffctProposalS2,params.qEffctLowerBound, params.qEffctUpperBound);
        
        proposal.newpi = eval_prior(chain.mSeeds,chain.mEffcts,chain.mrateMu,chain.mrateS2,
                                    proposal.newqSeeds,proposal.newqEffcts,chain.qrateMu,chain.qrateS2,
                                    chain.df)
        + log((chain.qtiles+params.qnegBiSize)/(newqtiles/params.qnegBiProb));
    } else {                      // Propose death
        if (chain.qtiles==2) { pBirth = 1.0; }
        newqtiles--;
        int qtileToRemove = draw.runif_int(0,newqtiles);
        MatrixXd oldqSeed = chain.qSeeds.row(qtileToRemove);
        removeRow(proposal.newqSeeds,qtileToRemove);
        removeElem(proposal.newqEffcts,qtileToRemove);
        pairwise_distance(proposal.newqSeeds,oldqSeed).col(0).minCoeff(&r);
        double chain.qEffct = proposal.newqEffcts(r);
        double oldqEffct = chain.qEffcts(qtileToRemove);
        // Compute log(prior ratio) and log(proposal ratio)
        proposal.newratioln = log(pBirth/pDeath)
        + dtrnormln(oldqEffct,chain.qEffct,params.qEffctProposalS2,params.qEffctLowerBound, params.qEffctUpperBound);
        
        proposal.newpi = eval_prior(chain.mSeeds,chain.mEffcts,chain.mrateMu,chain.mrateS2,
                                    proposal.newqSeeds,proposal.newqEffcts,chain.qrateMu,chain.qrateS2,
                                    chain.df)
        + log((chain.qtiles/params.qnegBiProb)/(newqtiles+params.qnegBiSize));
    }
    proposal.move = Q_VORONOI_BIRTH_DEATH;
    proposal.newqtiles = newqtiles;
    proposal.newll = eems2_likelihood(chain.mSeeds, chain.mEffcts, chain.mrateMu, proposal.newqSeeds, proposal.newqEffcts, chain.qrateMu, chain.df, false);
}
void EEMS2::propose_birthdeath_mVoronoi(Proposal &proposal, int chainIndex) {
    int newmtiles = chain.mtiles,r;
    double u = draw.runif();
    double pBirth = 0.5;
    double pDeath = 0.5;
    proposal.newmEffcts = chain.mEffcts;
    proposal.newmSeeds = chain.mSeeds;
    if ((chain.mtiles==1) || (u<0.5)) { // Propose birth
        if (chain.mtiles==1) { pBirth = 1.0; }
        newmtiles++;
        MatrixXd newmSeed = MatrixXd::Zero(1,2);
        randpoint_in_habitat(newmSeed);
        pairwise_distance(chain.mSeeds,newmSeed).col(0).minCoeff(&r);
        double chain.mEffct = chain.mEffcts(r);
        double newmEffct = draw.rtrnorm(chain.mEffct,params.mEffctProposalS2,params.mEffctLowerBound, params.mEffctUpperBound);
        insertRow(proposal.newmSeeds,newmSeed.row(0));
        insertElem(proposal.newmEffcts,newmEffct);
        // Compute log(prior ratio) and log(proposal ratio)
        proposal.newratioln = log(pDeath/pBirth)
        - dtrnormln(newmEffct,chain.mEffct,params.mEffctProposalS2,params.mEffctLowerBound, params.mEffctUpperBound);
        
        proposal.newpi = eval_prior(proposal.newmSeeds,proposal.newmEffcts,chain.mrateMu,chain.mrateS2,
                                    chain.qSeeds,chain.qEffcts,chain.qrateMu,chain.qrateS2,
                                    chain.df)
        + log((chain.mtiles+params.mnegBiSize)/(newmtiles/params.mnegBiProb));
    } else {                      // Propose death
        if (chain.mtiles==2) { pBirth = 1.0; }
        newmtiles--;
        int mtileToRemove = draw.runif_int(0,newmtiles);
        MatrixXd oldmSeed = chain.mSeeds.row(mtileToRemove);
        removeRow(proposal.newmSeeds,mtileToRemove);
        removeElem(proposal.newmEffcts,mtileToRemove);
        pairwise_distance(proposal.newmSeeds,oldmSeed).col(0).minCoeff(&r);
        double chain.mEffct = proposal.newmEffcts(r);
        double oldmEffct = chain.mEffcts(mtileToRemove);
        // Compute log(prior ratio) and log(proposal ratio)
        proposal.newratioln = log(pBirth/pDeath)
        + dtrnormln(oldmEffct,chain.mEffct,params.mEffctProposalS2,params.mEffctLowerBound, params.mEffctUpperBound);
        
        proposal.newpi = eval_prior(proposal.newmSeeds,proposal.newmEffcts,chain.mrateMu,chain.mrateS2,
                                    chain.qSeeds,chain.qEffcts,chain.qrateMu,chain.qrateS2,
                                    chain.df)
        + log((chain.mtiles/params.mnegBiProb)/(newmtiles+params.mnegBiSize));
    }
    proposal.move = M_VORONOI_BIRTH_DEATH;
    proposal.newmtiles = newmtiles;
    proposal.newll = eems2_likelihood(proposal.newmSeeds, proposal.newmEffcts, chain.mrateMu, chain.qSeeds, chain.qEffcts, chain.qrateMu, chain.df, true);
}
void EEMS2::update_hyperparams(int chainIndex) {
    
    Chain chain = chains[chainIndex];
    
    double SSq = chain.qEffcts.squaredNorm();
    double SSm = chain.mEffcts.squaredNorm();
    
    chain.qrateS2 = draw.rinvgam(params.qrateShape_2 + 0.5 * chain.qtiles, params.qrateScale_2 + 0.5 * SSq);
    chain.mrateS2 = draw.rinvgam(params.mrateShape_2 + 0.5 * chain.mtiles, params.mrateScale_2 + 0.5 * SSm);
    chain.pi = eval_prior(chain.mSeeds,chain.mEffcts,chain.mrateMu,chain.mrateS2,
                       chain.qSeeds,chain.qEffcts,chain.qrateMu,chain.qrateS2,
                       chain.df);
}
bool EEMS2::accept_proposal(Proposal &proposal, int chainIndex) {
    
    Chain chain = chains[chainIndex]
    double u = draw.runif( );
    // The proposal cannot be accepted because the prior is 0
    // This can happen if the proposed value falls outside the parameter's support
    if ( proposal.newpi == -Inf ) {
        proposal.newpi = chain.pi;
        proposal.newll = chain.ll;
        return false;
    }
    double ratioln = (proposal.newpi - chain.pi + proposal.newll - chain.ll)/chain.Temperature;
    // If the proposal is either birth or death, add the log(proposal ratio)
    if (proposal.move==Q_VORONOI_BIRTH_DEATH ||
        proposal.move==M_VORONOI_BIRTH_DEATH) {
        ratioln += proposal.newratioln;
    }

    if ( log(u) < min(0.0,ratioln) ) {
        switch (proposal.move) {
            case Q_VORONOI_RATE_UPDATE:
                chain.qEffcts = proposal.newqEffcts;
                break;
            case Q_VORONOI_POINT_MOVE:
                chainqSeeds = proposal.newqSeeds;
                // Update the mapping of demes to qVoronoi tiles
                graph.index_closest_to_deme(chain.qSeeds,chain.qColors);
                break;
            case Q_VORONOI_BIRTH_DEATH:
                chain.qSeeds = proposal.newqSeeds;
                chain.qEffcts = proposal.newqEffcts;
                chain.qtiles = proposal.newqtiles;
                // Update the mapping of demes to qVoronoi tiles
                graph.index_closest_to_deme(chain.qSeeds,chain.qColors);
                break;
            case M_VORONOI_RATE_UPDATE:
                chain.mEffcts = proposal.newmEffcts;
                break;
            case M_MEAN_RATE_UPDATE:
                chain.mrateMu = proposal.newmrateMu;
                break;
            case M_VORONOI_POINT_MOVE:
                chain.mSeeds = proposal.newmSeeds;
                // Update the mapping of demes to mVoronoi tiles
                graph.index_closest_to_deme(chain.mSeeds,chain.mColors);
                break;
            case M_VORONOI_BIRTH_DEATH:
                chain.mSeeds = proposal.newmSeeds;
                chain.mEffcts = proposal.newmEffcts;
                chain.mtiles = proposal.newmtiles;
                // Update the mapping of demes to mVoronoi tiles
                graph.index_closest_to_deme(chain.mSeeds,chain.mColors);
                break;
            case Q_MEAN_RATE_UPDATE:
                chain.qrateMu = proposal.newqrateMu;
                break;
            case DF_UPDATE:
                chain.df = proposal.newdf;
                break;
            default:
                cerr << "[RJMCMC] Unknown move type" << endl;
                exit(1);
        }
        chain.pi = proposal.newpi;
        chain.ll = proposal.newll;
        return true;
    } else {
        proposal.newpi = chain.pi;
        proposal.newll = chain.ll;
        return false;
    }
}
///////////////////////////////////////////
void EEMS2::print_iteration(const MCMC &mcmc, int chainIndex) const {
    Chain chain = chains[chainIndex]
    cerr << " Ending iteration " << mcmc.currIter
    << " with acceptance proportions:" << endl << mcmc
    << " and over-dispersion parameter = " << pow(10, chain.df) << setprecision(4) << endl
    << "         number of qVoronoi tiles = " << chain.qtiles << endl
    << "         number of mVoronoi tiles = " << chain.mtiles << endl
    << "          Log prior = " << chain.pi << setprecision(4) << endl
    << "          Log llike = " << chain.ll << setprecision(4) << endl;
}
void EEMS2::save_iteration(const MCMC &mcmc, int chainIndex) {
    Chain chain = chains[chainIndex];
    int iter = mcmc.to_save_iteration( );
    mcmcqhyper(iter,0) = chain.qrateMu;
    mcmcqhyper(iter,1) = chain.qrateS2;
    mcmcmhyper(iter,0) = chain.mrateMu;
    mcmcmhyper(iter,1) = chain.mrateS2;
    mcmcpilogl(iter,0) = chain.pi;
    mcmcpilogl(iter,1) = chain.ll;
    mcmcqtiles(iter) = chain.qtiles;
    mcmcmtiles(iter) = chain.mtiles;
    mcmcthetas(iter) = chain.df;
    for ( int t = 0 ; t < chain.qtiles ; t++ ) {
        mcmcqRates.push_back(pow(10.0,chain.qEffcts(t) + chain.qrateMu));
    }
    for ( int t = 0 ; t < chain.qtiles ; t++ ) {
        mcmcwCoord.push_back(chain.qSeeds(t,0));
    }
    for ( int t = 0 ; t < chain.qtiles ; t++ ) {
        mcmczCoord.push_back(chain.qSeeds(t,1));
    }
    for ( int t = 0 ; t < chain.mtiles ; t++ ) {
        mcmcmRates.push_back(pow(10.0,chain.mEffcts(t) + chain.mrateMu));
    }
    for ( int t = 0 ; t < chain.mtiles ; t++ ) {
        mcmcxCoord.push_back(chain.mSeeds(t,0));
    }
    for ( int t = 0 ; t < chain.mtiles ; t++ ) {
        mcmcyCoord.push_back(chain.mSeeds(t,1));
    }
    
    JtDhatJ += expectedIBD;
}
bool EEMS2::output_current_state(int chainIndex ) const {
    Chain chain = chains[chainIndex];
    ofstream out; bool error = false;
    out.open((params.mcmcpath + "/lastqtiles.txt").c_str(),ofstream::out);
    if (!out.is_open()) { error = true; return(error); }
    out << chain.qtiles << endl;
    out.close( );
    out.open((params.mcmcpath + "/lastmtiles.txt").c_str(),ofstream::out);
    if (!out.is_open()) { error = true; return(error); }
    out << chain.mtiles << endl;
    out.close( );
    out.open((params.mcmcpath + "/lastthetas.txt").c_str(),ofstream::out);
    if (!out.is_open()) { error = true; return(error); }
    out << fixed << setprecision(6) << chain.df << endl;
    out.close( );
    out.open((params.mcmcpath + "/lastdfpars.txt").c_str(),ofstream::out);
    if (!out.is_open()) { error = true; return(error); }
    out << fixed << setprecision(6) << params.dfmin << " " << params.dfmax << endl;
    out.close( );
    out.open((params.mcmcpath + "/lastmhyper.txt").c_str(),ofstream::out);
    if (!out.is_open()) { error = true; return(error); }
    out << fixed << setprecision(6) << chain.mrateMu << " " << chain.mrateS2 << endl;
    out.close( );
    out.open((params.mcmcpath + "/lastqhyper.txt").c_str(),ofstream::out);
    if (!out.is_open()) { error = true; return(error); }
    out << fixed << setprecision(6) << chain.qrateMu << " " << chain.qrateS2 << endl;
    out.close( );
    out.open((params.mcmcpath + "/lastpilogl.txt").c_str(),ofstream::out);
    if (!out.is_open()) { error = true; return(error); }
    out << fixed << setprecision(6) << chain.pi << " " << chain.ll << endl;
    out.close( );
    out.open((params.mcmcpath + "/lastmeffct.txt").c_str(),ofstream::out);
    if (!out.is_open()) { error = true; return(error); }
    out << fixed << setprecision(6) << chain.mEffcts << endl;
    out.close( );
    out.open((params.mcmcpath + "/lastmseeds.txt").c_str(),ofstream::out);
    if (!out.is_open()) { error = true; return(error); }
    out << fixed << setprecision(6) << chain.mSeeds << endl;
    out.close( );
    out.open((params.mcmcpath + "/lastqeffct.txt").c_str(),ofstream::out);
    if (!out.is_open()) { error = true; return(error); }
    out << fixed << setprecision(6) << chain.qEffcts << endl;
    out.close( );
    out.open((params.mcmcpath + "/lastqseeds.txt").c_str(),ofstream::out);
    if (!out.is_open()) { error = true; return(error); }
    out << fixed << setprecision(6) << chain.qSeeds << endl;
    out.close( );
    
    return(error);
}
bool EEMS2::output_results(const MCMC &mcmc, int chainIndex) const {
    Chain chain = chains[chainIndex];
    ofstream out; bool error = false;
    MatrixXd oDemes = MatrixXd::Zero(o,3);
    oDemes << graph.get_the_obsrv_demes(),cvec;
    out.open((params.mcmcpath + "/rdistoDemes.txt").c_str(),ofstream::out);
    if (!out.is_open()) { return false; }
    out << oDemes << endl;
    out.close( );
    out.open((params.mcmcpath + "/rdistJtDobsJ.txt").c_str(),ofstream::out);
    if (!out.is_open()) { return false; }
    out << observedIBD.array() / cMatrix.array() << endl;
    out.close( );
    out.open((params.mcmcpath + "/rdistJtDhatJ.txt").c_str(),ofstream::out);
    if (!out.is_open()) { return false; }
    double niters = mcmc.num_iters_to_save();
    out << JtDhatJ/niters << endl;
    out.close( );
    out.open((params.mcmcpath + "/mcmcqtiles.txt").c_str(),ofstream::out);
    if (!out.is_open()) { return false; }
    out << fixed << setprecision(14) << mcmcqtiles << endl;
    out.close( );
    out.open((params.mcmcpath + "/mcmcmtiles.txt").c_str(),ofstream::out);
    if (!out.is_open()) { return false; }
    out << fixed << setprecision(14) << mcmcmtiles << endl;
    out.close( );
    out.open((params.mcmcpath + "/mcmcthetas.txt").c_str(),ofstream::out);
    if (!out.is_open()) { return false; }
    out << fixed << setprecision(14) << mcmcthetas << endl;
    out.close( );
    out.open((params.mcmcpath + "/mcmcqhyper.txt").c_str(),ofstream::out);
    if (!out.is_open()) { return false; }
    out << fixed << setprecision(14) << mcmcqhyper << endl;
    out.close( );
    out.open((params.mcmcpath + "/mcmcmhyper.txt").c_str(),ofstream::out);
    if (!out.is_open()) { return false; }
    out << fixed << setprecision(14) << mcmcmhyper << endl;
    out.close( );
    out.open((params.mcmcpath + "/mcmcpilogl.txt").c_str(),ofstream::out);
    if (!out.is_open()) { return false; }
    out << fixed << setprecision(14) << mcmcpilogl << endl;
    out.close( );
    error = dlmcell(params.mcmcpath + "/mcmcmrates.txt",mcmcmtiles,mcmcmRates); if (error) { return(error); }
    error = dlmcell(params.mcmcpath + "/mcmcxcoord.txt",mcmcmtiles,mcmcxCoord); if (error) { return(error); }
    error = dlmcell(params.mcmcpath + "/mcmcycoord.txt",mcmcmtiles,mcmcyCoord); if (error) { return(error); }
    error = dlmcell(params.mcmcpath + "/mcmcqrates.txt",mcmcqtiles,mcmcqRates); if (error) { return(error); }
    error = dlmcell(params.mcmcpath + "/mcmcwcoord.txt",mcmcqtiles,mcmcwCoord); if (error) { return(error); }
    error = dlmcell(params.mcmcpath + "/mcmczcoord.txt",mcmcqtiles,mcmczCoord); if (error) { return(error); }
    error = output_current_state(chainIndex); if (error) { return(error); }
    out.open((params.mcmcpath + "/eemsrun.txt").c_str(),ofstream::out);
    if (!out.is_open( )) { return false; }
    out << "Input parameter values:" << endl << params << endl
    << "Acceptance proportions:" << endl << mcmc << endl
    << "Final log prior: " << chain.pi << endl
    << "Final log llike: " << chain.ll << endl;
    out.close( );
    cerr << "Final log prior: " << chain.pi << endl
    << "Final log llike: " << chain.ll << endl;
    return(error);
}

void EEMS2::check_ll_computation(int chainIndex ) const {
    Chain chain = chains[chainIndex];
    double pi0 = eval_prior(chain.mSeeds,chain.mEffcts,chain.mrateMu,chain.mrateS2,
                            chain.qSeeds,chain.qEffcts,chain.qrateMu,chain.qrateS2,
                            chain.df);
    double ll0 = eems2_likelihood(chain.mSeeds, chain.mEffcts, chain.mrateMu, chain.qSeeds, chain.qEffcts, chain.qrateMu, chain.df, true);
    if ((abs(chain.pi-pi0)/abs(pi0)>1e-12)||
        (abs(chain.ll-ll0)/abs(ll0)>1e-12)) {
        cerr << "[EEMS2::testing]   |ll0-ll|/|ll0| = " << abs(chain.ll - ll0)/abs(ll0) << endl;
        cerr << "[EEMS2::testing]   |pi0-pi|/|pi0| = " << abs(chain.pi - pi0)/abs(pi0) << endl;
        exit(1);
    }
}
double EEMS2::eval_prior(const MatrixXd &mSeeds, const VectorXd &mEffcts, const double mrateMu, const double mrateS2,
                         const MatrixXd &qSeeds, const VectorXd &qEffcts, const double qrateMu, const double qrateS2,
                         const double df) const {
    // Important: Do not use any of the current parameter values in this function,
    // i.e., those that start with nowXXX
    
    bool inrange = true;
    int qtiles = qEffcts.size();
    int mtiles = mEffcts.size();
    for ( int i = 0 ; i < qtiles ; i++ ) {
        if (!habitat.in_point(qSeeds(i,0),qSeeds(i,1))) { inrange = false; }
    }
    for ( int i = 0 ; i < mtiles ; i++ ) {
        if (!habitat.in_point(mSeeds(i,0),mSeeds(i,1))) { inrange = false; }
    }
    
    if (qEffcts.minCoeff() < params.qEffctLowerBound || qEffcts.maxCoeff() > params.qEffctUpperBound) { inrange = false; }
    if (mEffcts.minCoeff() < params.mEffctLowerBound || mEffcts.maxCoeff() > params.mEffctUpperBound) { inrange = false; }

    if (mrateMu>params.mrateMuUpperBound || mrateMu < params.qrateMuLowerBound) { inrange = false; }
    if (qrateMu>params.qrateMuUpperBound || qrateMu < params.qrateMuLowerBound ) { inrange = false; }
    
    if (df<params.dfmin || df>params.dfmax) { inrange = false; }
   
    if (!inrange) { return (-Inf); }
    
    double logpi =
    + dnegbinln(mtiles,params.mnegBiSize,params.mnegBiProb)
    + dnegbinln(qtiles,params.qnegBiSize,params.qnegBiProb)
    + dinvgamln(mrateS2,params.mrateShape_2,params.mrateScale_2)
    + dinvgamln(qrateS2,params.qrateShape_2,params.qrateScale_2);
    for (int i = 0 ; i < qtiles ; i++) {
        logpi += dtrnormln(qEffcts(i),0.0,qrateS2,params.qEffctLowerBound, params.qEffctUpperBound);
    }
    for (int i = 0 ; i < mtiles ; i++) {
        logpi += dtrnormln(mEffcts(i),0.0,mrateS2,params.mEffctLowerBound, params.mEffctUpperBound);
    }
    return (logpi);
}

void EEMS2::calculateIntegral(MatrixXd &eigenvals, MatrixXd &eigenvecs, const VectorXd &q, MatrixXd &integral, double bnd) const {
    
    // weights for the gaussian quadrature
    VectorXd weights(50);
    // abisca for the gaussian quadrature
    VectorXd x(50);
    
    getWeights(weights, x);
    
    // change of variable ("u substitution")
    weights = weights*(params.genomeSize/bnd)*(100/(2*bnd));
    x = x/(2*bnd/100);
    
    MatrixXd coalp = MatrixXd::Zero(o,o);
    MatrixXd P(d,d);
    integral.setZero();
    
    // x.size() is 50 but we're going to 25 to speed up the algorithm as the
    // magnnitude of the remaining 25 weights are neglible.
    for (int t = 0; t < 25; t++){
        P = eigenvecs.topRows(o) * ( ((VectorXd)((eigenvals.array() * x[t]).exp())).asDiagonal() * eigenvecs.transpose());
        coalp = P * q.asDiagonal() * P.transpose();
        integral += coalp*weights(t);
    }
    
}

double EEMS2::eems2_likelihood(const MatrixXd &mSeeds, const VectorXd &mEffcts, const double mrateMu,
                               const MatrixXd &qSeeds, const VectorXd &qEffcts, const double qrateMu,
                               const double df, const bool ismUpdate) const {
    
    // Important: Do not use any of the current parameter values in this function,
    // i.e., those that start with nowXXX
    
    // mSeeds, mEffcts and mrateMu define the migration Voronoi tessellation
    // qSeeds, qEffcts define the diversity Voronoi tessellation
    VectorXi mColors, qColors;
    // For every deme in the graph -- which diversity tile does the deme fall into?
    graph.index_closest_to_deme(qSeeds,qColors);
    // For every deme in the graph -- which migration tile does the deme fall into?
    graph.index_closest_to_deme(mSeeds,mColors);
    
    VectorXd q = VectorXd::Zero(d);
    // Transform the log10 diversity parameters into diversity rates on the original scale
    for ( int alpha = 0 ; alpha < d ; alpha++ ) {
        double log10q_alpha = qEffcts(qColors(alpha)) + qrateMu;
        q(alpha) = pow(10.0,log10q_alpha);
    }

    
    // perform costly eigen-decompositon only if updating migration rates
    if (ismUpdate){
        
        MatrixXd M = MatrixXd::Zero(d,d);
        int alpha, beta;
        // Transform the log10 migration parameters into migration rates on the original scale
        for ( int edge = 0 ; edge < graph.get_num_edges() ; edge++ ) {
            graph.get_edge(edge,alpha,beta);
            double log10m_alpha = mEffcts(mColors(alpha)) + mrateMu;
            double log10m_beta = mEffcts(mColors(beta)) + mrateMu;
            M(alpha,beta) = 0.5 * pow(10.0,log10m_alpha) + 0.5 * pow(10.0,log10m_beta);
            M(beta,alpha) = M(alpha,beta);
        }
        
        // Make M into a rate matrix
        M.diagonal() = -1* M.rowwise().sum();
        
        SelfAdjointEigenSolver<MatrixXd> es;
        es.compute(M);
        eigenvals = es.eigenvalues();
        eigenvecs = es.eigenvectors();
    }
    
    MatrixXd lowerExpectedIBD = MatrixXd::Zero(o, o);
    calculateIntegral(eigenvals, eigenvecs, q, lowerExpectedIBD, params.lowerBound);
    
    
    if (isfinite(params.upperBound)){
        MatrixXd upperExpectedIBD = MatrixXd::Zero(o, o);
        calculateIntegral(eigenvals, eigenvecs, q, upperExpectedIBD, params.upperBound);
        expectedIBD = lowerExpectedIBD - upperExpectedIBD;
    }
    else{
        expectedIBD = lowerExpectedIBD;
    }
    
    
    double phi = pow(10.0, df);
    double logll = negbiln(expectedIBD, observedIBD, cvec, cClasses, phi);
    
    return (logll);
}
double EEMS2::getMigrationRate(const int edge, int chainIndex) const {
    Chain chain = chains[chainIndex];
    int nEdges = graph.get_num_edges();
    assert((edge>=0) && (edge<nEdges));
    int alpha, beta;
    graph.get_edge(edge,alpha,beta);
    /*
     It will be highly inefficient to compute the mapping mColors every time we want to
     look up the migration rate of an edge (Therefore, nowmColors is updated after any
     update that can change the mapping of vertices/demes to mVoronoi tiles)
     VectorXi mColors;
     graph.index_closest_to_deme(nowmSeeds,mColors);
     */
    double m_alpha = pow(10, chain.mEffcts(chain.mColors(alpha)) + chain.mrateMu);
    double m_beta = pow(10, chain.mEffcts(chain.mColors(beta)) + chain.mrateMu);
    return (0.5 * m_alpha + 0.5 * m_beta);
    
}
double EEMS2::getCoalescenceRate(const int alpha, int chainIndex) const {
    Chain chain = chains[chainIndex];
    int nDemes = graph.get_num_total_demes();
    assert((alpha>=0) && (alpha<nDemes));
    /*
     It will be highly inefficient to compute the mapping qColors every time we want to
     look up the coalescence rate of a deme (Therefore, nowqColors is updated after any
     update that can change the mapping of vertices/demes to qVoronoi tiles)
     VectorXi qColors;
     graph.index_closest_to_deme(nowqSeeds,qColors);
     */
    double q_alpha = pow(10, chain.qEffcts(chain.qColors(alpha)) + chain.qrateMu);
    return (q_alpha);
}

// this is a test function
void EEMS2::writePopSizes(int chainIndex) const{
    Chain chain = chains[chainIndex];
    ofstream out;
    VectorXi mColors, qColors;
    graph.index_closest_to_deme(chain.qSeeds,qColors);
    graph.index_closest_to_deme(chain.mSeeds,mColors);
    VectorXd q = VectorXd::Zero(d);
    for ( int alpha = 0 ; alpha < d ; alpha++ ) {
        double log10q_alpha = chain.qEffcts(qColors(alpha)) + chain.qrateMu;
        q(alpha) = pow(10.0,log10q_alpha);
    }
    
    out.open((params.mcmcpath + "/coalescentrates.txt").c_str(), ios::out | ios::app);
    out << fixed << setprecision(14) << q.transpose() << endl;
    out.close( );
    
}

void EEMS2::printMigrationAndCoalescenceRates(int chainIndex ) const {
    
    int nDemes = graph.get_num_total_demes();
    cout << "Here is the coalescence (actually, diversity) rate of each deme" << endl;
    for ( int alpha = 0 ; alpha<nDemes ; alpha++ ) {
        cout << "  deme = " << alpha << ", q rate = " << getCoalescenceRate(alpha, chainIndex) << endl;
    }
    int nEdges = graph.get_num_edges();
    cout << "Here is the migration rate of each edge" << endl;
    for ( int edge = 0 ; edge<nEdges ; edge++ ) {
        cout << "  edge = " << edge << ", m rate = " << getMigrationRate(edge, chainIndex) << endl;
    }
    int alpha, beta;
    cout << "Here is the migration rate of each edge, this time edges as pairs of demes" << endl;
    for ( int edge = 0 ; edge<nEdges ; edge++ ) {
        graph.get_edge(edge,alpha,beta);
        cout << "  edge = (" << alpha << "," << beta << "), m rate = " << getMigrationRate(edge, chainIndex) << endl;
    }
    
}

int EEMS2:getnChains() const {
    return(nChains);
}

int EEMS2:getll(int chainIndex) const {
    Chain chain = chains[chainIndex];
    return(chain.ll);
}

double EEMS2:getTemperature(int chainIndex) const {
    Chain chain = chains[chainIndex];
    return(chain.Temperature);
}

void EEMS2:transferChain(int donorIndex, int receiptIndex){
    
    Chain donorChain = chains[donorIndex];
    Chain receiptChain = chains[receiptIndex];
    
    receiptChain.qtiles = donorChain.qtiles;
    receiptChain.mtiles = donorChain.mtiles;
    receiptChain.df = donorChain.df;
    receiptChain.pi = donorChain.pi;
    receiptChain.ll = donorChain.ll;
    receiptChain.mrateMu = donorChain.mrateMu;
    receiptChain.qrateMu = donorChain.qrateMu;
    receiptChain.qEffcts = donorChain.qEffcts;
    receiptChain.qSeeds = donorChain.qSeeds;
    receiptChain.mSeeds = donorChain.mSeeds;
    receiptChain.mrateS2 = donorChain.mrateS2;
    receiptChain.qrateS2 = donorChain.qrateS2;
    
}