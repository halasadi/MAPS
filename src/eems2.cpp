
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

void EEMS2::rnorm_effects(const double mu, const double rateS2, const double upperBound, VectorXd &Effcts) {
    for (int i = 0 ; i < Effcts.rows() ; i++ ) {
        Effcts(i) = draw.rtrnorm(mu,rateS2,upperBound);
    }
}

void EEMS2::initialize_sims( ) {
    cerr << "[Sims::initialize]" << endl;
    MatrixXd Sims = readMatrixXd(params.datapath + ".diffs");
    cout << Sims.rows() << endl;
    cout << Sims.cols() << endl;
    if ((Sims.rows()!=n)||(Sims.cols()!=n)) {
        cerr << "  Error reading similarities matrix " << params.datapath + ".diffs" << endl
        << "  Expect a " << n << "x" << n << " matrix of pairwise similarities" << endl; exit(1);
    }
    cerr << "  Loaded similarities matrix from " << params.datapath + ".diffs" << endl;
    
    cerr << "[Sims::initialize] Done." << endl << endl;
    
    // the number of comparisons for each deme.
    
    cMatrix = MatrixXd::Zero(o,o);
    cvec = VectorXd::Zero(o);
    totalSharingM = MatrixXd::Zero(o, o);
    int demei;
    int demej;
    // all (n choose 2) comparisons at the individual level
    for ( int i = 0 ; i < n ; i ++ ) {
        for (int j = (i+1); j < n; j++){
            demei = graph.get_deme_of_indiv(i);
            demej = graph.get_deme_of_indiv(j);
            cMatrix(demei, demej) += 1;
            cMatrix(demej, demei) += 1;
            totalSharingM(demei, demej) += Sims(i,j);
            totalSharingM(demej, demei) = totalSharingM(demei, demej);
        }
    }
    for ( int i = 0 ; i < n ; i ++ ) {
        cvec(graph.get_deme_of_indiv(i)) += 1;
    }
    
}
void EEMS2::initialize_state( ) {
    cerr << "[EEMS2::initialize_state]" << endl;
    nowdf = n;
    // Initialize the two Voronoi tessellations
    nowqtiles = draw.rnegbin(2*o,0.5); // o is the number of observed demes
    nowmtiles = draw.rnegbin(2*o,0.5);
    cerr << "  EEMS starts with " << nowqtiles << " qtiles and " << nowmtiles << " mtiles" << endl;
    // Draw the Voronoi centers Coord uniformly within the habitat
    nowqSeeds = MatrixXd::Zero(nowqtiles,2); randpoint_in_habitat(nowqSeeds);
    nowmSeeds = MatrixXd::Zero(nowmtiles,2); randpoint_in_habitat(nowmSeeds);
    nowmrateS2 = draw.rinvgam(0.5,0.5);
    nowqrateS2 = draw.rinvgam(0.5,0.5);
    
    // Assign migration rates to the Voronoi tiles
    double lowerBound = -10;
    nowmrateMu = lowerBound + draw.runif() * (params.mrateMuUpperBound - lowerBound);
    nowqrateMu = lowerBound + draw.runif() * (params.qrateMuUpperBound - lowerBound);
    
    // Assign rates to the Voronoi tiles
    nowqEffcts = VectorXd::Zero(nowqtiles); rnorm_effects(nowqrateMu, nowqrateS2, params.qEffctUpperBound, nowqEffcts);
    nowmEffcts = VectorXd::Zero(nowmtiles); rnorm_effects(nowmrateMu, nowmrateS2, params.mEffctUpperBound, nowmEffcts);
    // Initialize the mapping of demes to qVoronoi tiles
    graph.index_closest_to_deme(nowqSeeds,nowqColors);
    // Initialize the mapping of demes to mVoronoi tiles
    graph.index_closest_to_deme(nowmSeeds,nowmColors);
    cerr << "[EEMS2::initialize_state] Done." << endl << endl;
}
void EEMS2::load_final_state( ) {
    cerr << "[EEMS2::load_final_state]" << endl;
    MatrixXd tempi; bool error = false;
    tempi = readMatrixXd(params.prevpath + "/lastqtiles.txt");
    if ((tempi.rows()!=1) || (tempi.cols()!=1)) { error = true; }
    nowqtiles = tempi(0,0);
    tempi = readMatrixXd(params.prevpath + "/lastmtiles.txt");
    if ((tempi.rows()!=1) || (tempi.cols()!=1)) { error = true; }
    nowmtiles = tempi(0,0);
    cerr << "  EEMS starts with " << nowqtiles << " qtiles and " << nowmtiles << " mtiles" << endl;
    tempi = readMatrixXd(params.prevpath + "/lastthetas.txt");
    if ((tempi.rows()!=1) || (tempi.cols()!=2)) { error = true; }
    nowdf = tempi(0,0);
    tempi = readMatrixXd(params.prevpath + "/lastdfpars.txt");
    if ((tempi.rows()!=1) || (tempi.cols()!=2)) { error = true; }
    params.dfmin = tempi(0,0);
    params.dfmax = tempi(0,1);
    tempi = readMatrixXd(params.prevpath + "/lastqhyper.txt");
    if ((tempi.rows()!=1) || (tempi.cols()!=1)) { error = true; }
    nowqrateS2 = tempi(0,0);
    tempi = readMatrixXd(params.prevpath + "/lastmhyper.txt");
    if ((tempi.rows()!=1) || (tempi.cols()!=2)) { error = true; }
    nowmrateMu = tempi(0,0);
    nowmrateS2 = tempi(0,1);
    tempi = readMatrixXd(params.prevpath + "/lastqhyper.txt");
    if ((tempi.rows()!=1) || (tempi.cols()!=2)) { error = true; }
    nowqrateMu = tempi(0,0);
    nowqrateS2 = tempi(0,1);
    tempi = readMatrixXd(params.prevpath + "/lastqeffct.txt");
    if ((tempi.rows()!=nowqtiles) || (tempi.cols()!=1)) { error = true; }
    nowqEffcts = tempi.col(0);
    tempi = readMatrixXd(params.prevpath + "/lastmeffct.txt");
    if ((tempi.rows()!=nowmtiles) || (tempi.cols()!=1)) { error = true; }
    nowmEffcts = tempi.col(0);
    nowqSeeds = readMatrixXd(params.prevpath + "/lastqseeds.txt");
    if ((nowqSeeds.rows()!=nowqtiles) || (nowqSeeds.cols()!=2)) { error = true; }
    nowmSeeds = readMatrixXd(params.prevpath + "/lastmseeds.txt");
    if ((nowmSeeds.rows()!=nowmtiles) || (nowmSeeds.cols()!=2)) { error = true; }
    if (error) {
        cerr << "  Error loading MCMC state from " << params.prevpath << endl; exit(1);
    }
    // Initialize the mapping of demes to qVoronoi tiles
    graph.index_closest_to_deme(nowqSeeds,nowqColors);
    // Initialize the mapping of demes to mVoronoi tiles
    graph.index_closest_to_deme(nowmSeeds,nowmColors);
    cerr << "[EEMS2::load_final_state] Done." << endl << endl;
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
    nowpi = eval_prior(nowmSeeds,nowmEffcts,nowmrateMu,nowmrateS2,
                       nowqSeeds,nowqEffcts,nowqrateMu,nowqrateS2,
                       nowdf);
    nowll = eems2_likelihood(nowmSeeds, nowmEffcts, nowmrateMu, nowqSeeds, nowqEffcts, nowqrateMu, nowdf);
    cerr << "Input parameters: " << endl << params << endl
    << "Initial log prior: " << nowpi << endl
    << "Initial log llike: " << nowll << endl << endl;
    if ((nowpi==-Inf) || (nowpi==Inf) || (nowll==-Inf) || (nowll==Inf)) { error = true; }
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
        if (u2 < 0.5) {
            move = M_MEAN_RATE_UPDATE;
        } else{
            move = Q_MEAN_RATE_UPDATE;
        }
    }
    return(move);
}

double EEMS2::eval_proposal_rate_one_qtile(Proposal &proposal) const {
    return(eems2_likelihood(nowmSeeds, nowmEffcts, nowmrateMu, nowqSeeds, proposal.newqEffcts, nowqrateMu, nowdf));
}
double EEMS2::eval_proposal_move_one_qtile(Proposal &proposal) const {
    return(eems2_likelihood(nowmSeeds, nowmEffcts, nowmrateMu, proposal.newqSeeds, nowqEffcts, nowqrateMu, nowdf));
}
double EEMS2::eval_birthdeath_qVoronoi(Proposal &proposal) const {
    return(eems2_likelihood(nowmSeeds, nowmEffcts, nowmrateMu, proposal.newqSeeds, proposal.newqEffcts, nowqrateMu, nowdf));
}
double EEMS2::eval_proposal_rate_one_mtile(Proposal &proposal) const {
    return(eems2_likelihood(nowmSeeds, proposal.newmEffcts, nowmrateMu, nowqSeeds, nowqEffcts, nowqrateMu, nowdf));
    
}
double EEMS2::eval_proposal_overall_mrate(Proposal &proposal) const {
    return(eems2_likelihood(nowmSeeds, nowmEffcts, proposal.newmrateMu, nowqSeeds, nowqEffcts, nowqrateMu, nowdf));
}
double EEMS2::eval_proposal_overall_qrate(Proposal &proposal) const {
    return(eems2_likelihood(nowmSeeds, nowmEffcts, nowmrateMu, nowqSeeds, nowqEffcts, proposal.newqrateMu, nowdf));
}
// Propose to move one tile in the migration Voronoi tessellation
double EEMS2::eval_proposal_move_one_mtile(Proposal &proposal) const {
    return(eems2_likelihood(proposal.newmSeeds, nowmEffcts, nowmrateMu, nowqSeeds, nowqEffcts, nowqrateMu, nowdf));
}
double EEMS2::eval_birthdeath_mVoronoi(Proposal &proposal) const {
    return(eems2_likelihood(proposal.newmSeeds, proposal.newmEffcts, nowmrateMu, nowqSeeds, nowqEffcts, nowqrateMu, nowdf));
}
void EEMS2::propose_df(Proposal &proposal,const MCMC &mcmc) {
    proposal.move = DF_UPDATE;
    proposal.newpi = -Inf;
    proposal.newll = -Inf;
    // EEMS is initialized with df = nIndiv
    // Keep df = nIndiv for the first mcmc.numBurnIter/2 iterations
    // This should make it easier to move in the parameter space
    // since the likelihood is proportional to 0.5 * pdf * ll_atfixdf
    if (mcmc.currIter > (mcmc.numBurnIter/2)) {
        double newdf = draw.rnorm(nowdf,params.dfProposalS2);
        if ( (newdf>params.dfmin) && (newdf<params.dfmax) ) {
            proposal.newdf = newdf;
            proposal.newpi = eval_prior(nowmSeeds,nowmEffcts,nowmrateMu,nowmrateS2,
                                        nowqSeeds,nowqEffcts,nowqrateMu,nowqrateS2,
                                        newdf);
            proposal.newll = eems2_likelihood(nowmSeeds, nowmEffcts, nowmrateMu, nowqSeeds, nowqEffcts, nowqrateMu, newdf);
        }
    }
}

void EEMS2::propose_rate_one_qtile(Proposal &proposal) {
    // Choose a tile at random to update
    int qtile = draw.runif_int(0,nowqtiles-1);
    // Make a random-walk proposal, i.e., add small offset to current value
    double curqEffct = nowqEffcts(qtile);
    double newqEffct = draw.rnorm(curqEffct,params.qEffctProposalS2);
    proposal.move = Q_VORONOI_RATE_UPDATE;
    proposal.newqEffcts = nowqEffcts;
    proposal.newqEffcts(qtile) = newqEffct;
    // The prior distribution on the tile effects is truncated normal
    // So first check whether the proposed value is in range
    // Then update the prior and evaluate the new likelihood
    if ( newqEffct >= -100 && newqEffct <= params.qEffctUpperBound) {
        proposal.newpi = eval_prior(nowmSeeds,nowmEffcts,nowmrateMu,nowmrateS2,
                                    nowqSeeds,proposal.newqEffcts,nowqrateMu,nowqrateS2,
                                    nowdf);
        proposal.newll = eval_proposal_rate_one_qtile(proposal);
    } else {
        proposal.newpi = -Inf;
        proposal.newll = -Inf;
    }
}
void EEMS2::propose_rate_one_mtile(Proposal &proposal) {
    // Choose a tile at random to update
    int mtile = draw.runif_int(0,nowmtiles-1);
    // Make a random-walk proposal, i.e., add small offset to current value
    double curmEffct = nowmEffcts(mtile);
    double newmEffct = draw.rnorm(curmEffct,params.mEffctProposalS2);
    proposal.move = M_VORONOI_RATE_UPDATE;
    proposal.newmEffcts = nowmEffcts;
    proposal.newmEffcts(mtile) = newmEffct;
    // The prior distribution on the tile effects is truncated normal
    // So first check whether the proposed value is in range
    // Then update the prior and evaluate the new likelihood
    if ( newmEffct >= -100 && newmEffct <= params.mEffctUpperBound) {
        proposal.newpi = eval_prior(nowmSeeds,proposal.newmEffcts,nowmrateMu,nowmrateS2,
                                    nowqSeeds,nowqEffcts,nowqrateMu,nowqrateS2,
                                    nowdf);
        proposal.newll = eval_proposal_rate_one_mtile(proposal);
    } else {
        proposal.newpi = -Inf;
        proposal.newll = -Inf;
    }
}
void EEMS2::propose_overall_mrate(Proposal &proposal) {
    // Make a random-walk Metropolis-Hastings proposal
    double newmrateMu = draw.rnorm(nowmrateMu,params.mrateMuProposalS2);
    proposal.move = M_MEAN_RATE_UPDATE;
    proposal.newmrateMu = newmrateMu;
    // If the proposed value is in range, the prior probability does not change
    // as the prior distribution on mrateMu is uniform
    // Otherwise, setting the prior and the likelihood to -Inf forces a rejection
    if ( newmrateMu <= params.mrateMuUpperBound && newmrateMu >= -100) {
        proposal.newpi = nowpi;
        proposal.newll = eval_proposal_overall_mrate(proposal);
    } else {
        proposal.newpi = -Inf;
        proposal.newll = -Inf;
    }
}

void EEMS2::propose_overall_qrate(Proposal &proposal) {
    // Make a random-walk Metropolis-Hastings proposal
    double newqrateMu = draw.rnorm(nowqrateMu,params.qrateMuProposalS2);
    proposal.move = Q_MEAN_RATE_UPDATE;
    proposal.newqrateMu = newqrateMu;
    // If the proposed value is in range, the prior probability does not change
    // as the prior distribution on qrateMu is uniform
    // Otherwise, setting the prior and the likelihood to -Inf forces a rejection
    if ( newqrateMu <= params.qrateMuUpperBound && newqrateMu >= -100) {
        proposal.newpi = nowpi;
        proposal.newll = eval_proposal_overall_qrate(proposal);
    } else {
        proposal.newpi = -Inf;
        proposal.newll = -Inf;
    }
}
void EEMS2::propose_move_one_qtile(Proposal &proposal) {
    // Choose a tile at random to move
    int qtile = draw.runif_int(0,nowqtiles-1);
    // Make a random-walk proposal, i.e., add small offset to current value
    // In this case, there are actually two values -- longitude and latitude
    double newqSeedx = draw.rnorm(nowqSeeds(qtile,0),params.qSeedsProposalS2x);
    double newqSeedy = draw.rnorm(nowqSeeds(qtile,1),params.qSeedsProposalS2y);
    proposal.move = Q_VORONOI_POINT_MOVE;
    proposal.newqSeeds = nowqSeeds;
    proposal.newqSeeds(qtile,0) = newqSeedx;
    proposal.newqSeeds(qtile,1) = newqSeedy;
    if (habitat.in_point(newqSeedx,newqSeedy)) {
        proposal.newpi = nowpi;
        proposal.newll = eval_proposal_move_one_qtile(proposal);
    } else {
        proposal.newpi = -Inf;
        proposal.newll = -Inf;
    }
}
void EEMS2::propose_move_one_mtile(Proposal &proposal) {
    // Choose a tile at random to move
    int mtile = draw.runif_int(0,nowmtiles-1);
    // Make a random-walk proposal, i.e., add small offset to current value
    // In this case, there are actually two values -- longitude and latitude
    double newmSeedx = draw.rnorm(nowmSeeds(mtile,0),params.mSeedsProposalS2x);
    double newmSeedy = draw.rnorm(nowmSeeds(mtile,1),params.mSeedsProposalS2y);
    proposal.move = M_VORONOI_POINT_MOVE;
    proposal.newmSeeds = nowmSeeds;
    proposal.newmSeeds(mtile,0) = newmSeedx;
    proposal.newmSeeds(mtile,1) = newmSeedy;
    if (habitat.in_point(newmSeedx,newmSeedy)) {
        proposal.newpi = nowpi;
        proposal.newll = eval_proposal_move_one_mtile(proposal);
    } else {
        proposal.newpi = -Inf;
        proposal.newll = -Inf;
    }
}
void EEMS2::propose_birthdeath_qVoronoi(Proposal &proposal) {
    int newqtiles = nowqtiles,r;
    double u = draw.runif();
    double pBirth = 0.5;
    double pDeath = 0.5;
    proposal.newqEffcts = nowqEffcts;
    proposal.newqSeeds = nowqSeeds;
    // If there is exactly one tile, rule out a death proposal
    if ((nowqtiles==1) || (u<0.5)) { // Propose birth
        if (nowqtiles==1) { pBirth = 1.0; }
        newqtiles++;
        MatrixXd newqSeed = MatrixXd::Zero(1,2);
        randpoint_in_habitat(newqSeed);
        pairwise_distance(nowqSeeds,newqSeed).col(0).minCoeff(&r);
        // The new tile is assigned a rate by perturbing the current rate at the new seed
        double nowqEffct = nowqEffcts(r);
        double newqEffct = draw.rtrnorm(nowqEffct,params.qEffctProposalS2,params.qEffctUpperBound);
        insertRow(proposal.newqSeeds,newqSeed.row(0));
        insertElem(proposal.newqEffcts,newqEffct);
        // Compute log(proposal ratio) and log(prior ratio)
        //proposal.newratioln = log(pDeath/pBirth)
        //  - dtrnormln(newqEffct,nowqEffct,params.qEffctProposalS2,params.qEffctHalfInterval);
        proposal.newratioln = log(pDeath/pBirth)
        - dtrnormln(newqEffct,nowqEffct,params.qEffctProposalS2,params.qEffctUpperBound);
        
        proposal.newpi = eval_prior(nowmSeeds,nowmEffcts,nowmrateMu,nowmrateS2,
                                    proposal.newqSeeds,proposal.newqEffcts,nowqrateMu,nowqrateS2,
                                    nowdf)
        // CHECK
                        + log((nowqtiles+params.negBiSize)/(newqtiles/params.negBiProb));
    } else {                      // Propose death
        if (nowqtiles==2) { pBirth = 1.0; }
        newqtiles--;
        int qtileToRemove = draw.runif_int(0,newqtiles);
        MatrixXd oldqSeed = nowqSeeds.row(qtileToRemove);
        removeRow(proposal.newqSeeds,qtileToRemove);
        removeElem(proposal.newqEffcts,qtileToRemove);
        pairwise_distance(proposal.newqSeeds,oldqSeed).col(0).minCoeff(&r);
        double nowqEffct = proposal.newqEffcts(r);
        double oldqEffct = nowqEffcts(qtileToRemove);
        // Compute log(prior ratio) and log(proposal ratio)
        //proposal.newratioln = log(pBirth/pDeath)
        //  + dtrnormln(oldqEffct,nowqEffct,params.qEffctProposalS2,params.qEffctHalfInterval);
        
        proposal.newratioln = log(pBirth/pDeath)
        + dtrnormln(oldqEffct,nowqEffct,params.qEffctProposalS2,params.qEffctUpperBound);
        
        proposal.newpi = eval_prior(nowmSeeds,nowmEffcts,nowmrateMu,nowmrateS2,
                                    proposal.newqSeeds,proposal.newqEffcts,nowqrateMu,nowqrateS2,
                                    nowdf)
        // CHECK
                         + log((nowqtiles/params.negBiProb)/(newqtiles+params.negBiSize));
    }
    proposal.move = Q_VORONOI_BIRTH_DEATH;
    proposal.newqtiles = newqtiles;
    proposal.newll = eval_birthdeath_qVoronoi(proposal);
}
void EEMS2::propose_birthdeath_mVoronoi(Proposal &proposal) {
    int newmtiles = nowmtiles,r;
    double u = draw.runif();
    double pBirth = 0.5;
    double pDeath = 0.5;
    proposal.newmEffcts = nowmEffcts;
    proposal.newmSeeds = nowmSeeds;
    if ((nowmtiles==1) || (u<0.5)) { // Propose birth
        if (nowmtiles==1) { pBirth = 1.0; }
        newmtiles++;
        MatrixXd newmSeed = MatrixXd::Zero(1,2);
        randpoint_in_habitat(newmSeed);
        pairwise_distance(nowmSeeds,newmSeed).col(0).minCoeff(&r);
        double nowmEffct = nowmEffcts(r);
        double newmEffct = draw.rtrnorm(nowmEffct,params.mEffctProposalS2,params.mEffctUpperBound);
        insertRow(proposal.newmSeeds,newmSeed.row(0));
        insertElem(proposal.newmEffcts,newmEffct);
        // Compute log(prior ratio) and log(proposal ratio)
        //proposal.newratioln = log(pDeath/pBirth)
        //  - dtrnormln(newmEffct,nowmEffct,params.mEffctProposalS2,params.mEffctHalfInterval);
        proposal.newratioln = log(pDeath/pBirth)
        - dtrnormln(newmEffct,nowmEffct,params.mEffctProposalS2,params.mEffctUpperBound);
        
        proposal.newpi = eval_prior(proposal.newmSeeds,proposal.newmEffcts,nowmrateMu,nowmrateS2,
                                    nowqSeeds,nowqEffcts,nowqrateMu,nowqrateS2,
                                    nowdf)
        // CHECK
                        + log((nowmtiles+params.negBiSize)/(newmtiles/params.negBiProb));
    } else {                      // Propose death
        if (nowmtiles==2) { pBirth = 1.0; }
        newmtiles--;
        int mtileToRemove = draw.runif_int(0,newmtiles);
        MatrixXd oldmSeed = nowmSeeds.row(mtileToRemove);
        removeRow(proposal.newmSeeds,mtileToRemove);
        removeElem(proposal.newmEffcts,mtileToRemove);
        pairwise_distance(proposal.newmSeeds,oldmSeed).col(0).minCoeff(&r);
        double nowmEffct = proposal.newmEffcts(r);
        double oldmEffct = nowmEffcts(mtileToRemove);
        // Compute log(prior ratio) and log(proposal ratio)
        //proposal.newratioln = log(pBirth/pDeath)
        //  + dtrnormln(oldmEffct,nowmEffct,params.mEffctProposalS2,params.mEffctHalfInterval);
        proposal.newratioln = log(pBirth/pDeath)
        + dtrnormln(oldmEffct,nowmEffct,params.mEffctProposalS2,params.mEffctUpperBound);
        
        proposal.newpi = eval_prior(proposal.newmSeeds,proposal.newmEffcts,nowmrateMu,nowmrateS2,
                                    nowqSeeds,nowqEffcts,nowqrateMu,nowqrateS2,
                                    nowdf)
        // CHECK
                        + log((nowmtiles/params.negBiProb)/(newmtiles+params.negBiSize));
    }
    proposal.move = M_VORONOI_BIRTH_DEATH;
    proposal.newmtiles = newmtiles;
    proposal.newll = eval_birthdeath_mVoronoi(proposal);
}
void EEMS2::update_hyperparams( ) {
    double SSq = nowqEffcts.squaredNorm();
    double SSm = nowmEffcts.squaredNorm();
    
    nowqrateS2 = draw.rinvgam(params.qrateShape_2 + 0.5 * nowqtiles, params.qrateScale_2 + 0.5 * SSq);
    nowmrateS2 = draw.rinvgam(params.mrateShape_2 + 0.5 * nowmtiles, params.mrateScale_2 + 0.5 * SSm);
    nowpi = eval_prior(nowmSeeds,nowmEffcts,nowmrateMu,nowmrateS2,
                       nowqSeeds,nowqEffcts,nowqrateMu,nowqrateS2,
                       nowdf);
}
bool EEMS2::accept_proposal(Proposal &proposal) {
    double u = draw.runif( );
    // The proposal cannot be accepted because the prior is 0
    // This can happen if the proposed value falls outside the parameter's support
    if ( proposal.newpi == -Inf ) {
        proposal.newpi = nowpi;
        proposal.newll = nowll;
        return false;
    }
    double ratioln = proposal.newpi - nowpi + proposal.newll - nowll;
    // If the proposal is either birth or death, add the log(proposal ratio)
    if (proposal.move==Q_VORONOI_BIRTH_DEATH ||
        proposal.move==M_VORONOI_BIRTH_DEATH) {
        ratioln += proposal.newratioln;
    }
    if ( log(u) < min(0.0,ratioln) ) {
        switch (proposal.move) {
            case Q_VORONOI_RATE_UPDATE:
                nowqEffcts = proposal.newqEffcts;
                break;
            case Q_VORONOI_POINT_MOVE:
                nowqSeeds = proposal.newqSeeds;
                // Update the mapping of demes to qVoronoi tiles
                graph.index_closest_to_deme(nowqSeeds,nowqColors);
                break;
            case Q_VORONOI_BIRTH_DEATH:
                nowqSeeds = proposal.newqSeeds;
                nowqEffcts = proposal.newqEffcts;
                nowqtiles = proposal.newqtiles;
                // Update the mapping of demes to qVoronoi tiles
                graph.index_closest_to_deme(nowqSeeds,nowqColors);
                break;
            case M_VORONOI_RATE_UPDATE:
                nowmEffcts = proposal.newmEffcts;
                break;
            case M_MEAN_RATE_UPDATE:
                nowmrateMu = proposal.newmrateMu;
                break;
            case M_VORONOI_POINT_MOVE:
                nowmSeeds = proposal.newmSeeds;
                // Update the mapping of demes to mVoronoi tiles
                graph.index_closest_to_deme(nowmSeeds,nowmColors);
                break;
            case M_VORONOI_BIRTH_DEATH:
                nowmSeeds = proposal.newmSeeds;
                nowmEffcts = proposal.newmEffcts;
                nowmtiles = proposal.newmtiles;
                // Update the mapping of demes to mVoronoi tiles
                graph.index_closest_to_deme(nowmSeeds,nowmColors);
                break;
            case Q_MEAN_RATE_UPDATE:
                nowqrateMu = proposal.newqrateMu;
                break;
            case DF_UPDATE:
                nowdf = proposal.newdf;
                break;
            default:
                cerr << "[RJMCMC] Unknown move type" << endl;
                exit(1);
        }
        nowpi = proposal.newpi;
        nowll = proposal.newll;
        return true;
    } else {
        proposal.newpi = nowpi;
        proposal.newll = nowll;
        return false;
    }
}
///////////////////////////////////////////
void EEMS2::print_iteration(const MCMC &mcmc) const {
    cerr << " Ending iteration " << mcmc.currIter
    << " with acceptance proportions:" << endl << mcmc
    << " and effective degrees of freedom = " << nowdf << endl
    << "         number of qVoronoi tiles = " << nowqtiles << endl
    << "         number of mVoronoi tiles = " << nowmtiles << endl
    << "          Log prior = " << nowpi << endl
    << "          Log llike = " << nowll << endl;
}
void EEMS2::save_iteration(const MCMC &mcmc) {
    int iter = mcmc.to_save_iteration( );
    mcmcqhyper(iter,0) = nowqrateMu;
    mcmcqhyper(iter,1) = nowqrateS2;
    mcmcmhyper(iter,0) = nowmrateMu;
    mcmcmhyper(iter,1) = nowmrateS2;
    mcmcpilogl(iter,0) = nowpi;
    mcmcpilogl(iter,1) = nowll;
    mcmcqtiles(iter) = nowqtiles;
    mcmcmtiles(iter) = nowmtiles;
    mcmcthetas(iter) = nowdf;
    for ( int t = 0 ; t < nowqtiles ; t++ ) {
        mcmcqRates.push_back(pow(10.0, nowqEffcts(t)));
    }
    for ( int t = 0 ; t < nowqtiles ; t++ ) {
        mcmcwCoord.push_back(nowqSeeds(t,0));
    }
    for ( int t = 0 ; t < nowqtiles ; t++ ) {
        mcmczCoord.push_back(nowqSeeds(t,1));
    }
    for ( int t = 0 ; t < nowmtiles ; t++ ) {
        mcmcmRates.push_back(pow(10.0, nowmEffcts(t)));
    }
    for ( int t = 0 ; t < nowmtiles ; t++ ) {
        mcmcxCoord.push_back(nowmSeeds(t,0));
    }
    for ( int t = 0 ; t < nowmtiles ; t++ ) {
        mcmcyCoord.push_back(nowmSeeds(t,1));
    }
}
bool EEMS2::output_current_state( ) const {
    ofstream out; bool error = false;
    out.open((params.mcmcpath + "/lastqtiles.txt").c_str(),ofstream::out);
    if (!out.is_open()) { error = true; return(error); }
    out << nowqtiles << endl;
    out.close( );
    out.open((params.mcmcpath + "/lastmtiles.txt").c_str(),ofstream::out);
    if (!out.is_open()) { error = true; return(error); }
    out << nowmtiles << endl;
    out.close( );
    /*out.open((params.mcmcpath + "/lastthetas.txt").c_str(),ofstream::out);
    if (!out.is_open()) { error = true; return(error); }
    out << fixed << setprecision(6) << nowdf << endl;
    out.close( );
    */
    /*
    out.open((params.mcmcpath + "/lastdfpars.txt").c_str(),ofstream::out);
    if (!out.is_open()) { error = true; return(error); }
    out << fixed << setprecision(6) << params.dfmin << " " << params.dfmax << endl;
    out.close( );
    */
    out.open((params.mcmcpath + "/lastmhyper.txt").c_str(),ofstream::out);
    if (!out.is_open()) { error = true; return(error); }
    out << fixed << setprecision(6) << nowmrateMu << " " << nowmrateS2 << endl;
    out.close( );
    out.open((params.mcmcpath + "/lastqhyper.txt").c_str(),ofstream::out);
    if (!out.is_open()) { error = true; return(error); }
    out << fixed << setprecision(6) << nowqrateMu << " " << nowqrateS2 << endl;
    out.close( );
    out.open((params.mcmcpath + "/lastpilogl.txt").c_str(),ofstream::out);
    if (!out.is_open()) { error = true; return(error); }
    out << fixed << setprecision(6) << nowpi << " " << nowll << endl;
    out.close( );
    out.open((params.mcmcpath + "/lastmeffct.txt").c_str(),ofstream::out);
    if (!out.is_open()) { error = true; return(error); }
    out << fixed << setprecision(6) << nowmEffcts << endl;
    out.close( );
    out.open((params.mcmcpath + "/lastmseeds.txt").c_str(),ofstream::out);
    if (!out.is_open()) { error = true; return(error); }
    out << fixed << setprecision(6) << nowmSeeds << endl;
    out.close( );
    out.open((params.mcmcpath + "/lastqeffct.txt").c_str(),ofstream::out);
    if (!out.is_open()) { error = true; return(error); }
    out << fixed << setprecision(6) << nowqEffcts << endl;
    out.close( );
    out.open((params.mcmcpath + "/lastqseeds.txt").c_str(),ofstream::out);
    if (!out.is_open()) { error = true; return(error); }
    out << fixed << setprecision(6) << nowqSeeds << endl;
    out.close( );
    return(error);
}
bool EEMS2::output_results(const MCMC &mcmc) const {
    ofstream out; bool error = false;
    MatrixXd oDemes = MatrixXd::Zero(o,3);
    oDemes << graph.get_the_obsrv_demes(),cvec;
    out.open((params.mcmcpath + "/rdistoDemes.txt").c_str(),ofstream::out);
    if (!out.is_open()) { return false; }
    out << oDemes << endl;
    out.close( );
    /*out.open((params.mcmcpath + "/rdistJtDobsJ.txt").c_str(),ofstream::out);
     if (!out.is_open()) { return false; }
     MatrixXd Pairs = cvec*cvec.transpose(); Pairs -= cvec.asDiagonal();
     out << JtDobsJ.cwiseQuotient(Pairs) << endl;
     out.close( );
     out.open((params.mcmcpath + "/rdistJtDhatJ.txt").c_str(),ofstream::out);
     if (!out.is_open()) { return false; }
     int niters = mcmc.num_iters_to_save();
     out << JtDhatJ/niters << endl;
     out.close( );
     */
    out.open((params.mcmcpath + "/mcmcqtiles.txt").c_str(),ofstream::out);
    if (!out.is_open()) { return false; }
    out << mcmcqtiles << endl;
    out.close( );
    out.open((params.mcmcpath + "/mcmcmtiles.txt").c_str(),ofstream::out);
    if (!out.is_open()) { return false; }
    out << mcmcmtiles << endl;
    out.close( );
    out.open((params.mcmcpath + "/mcmcthetas.txt").c_str(),ofstream::out);
    if (!out.is_open()) { return false; }
    out << fixed << setprecision(6) << mcmcthetas << endl;
    out.close( );
    out.open((params.mcmcpath + "/mcmcqhyper.txt").c_str(),ofstream::out);
    if (!out.is_open()) { return false; }
    out << fixed << setprecision(6) << mcmcqhyper << endl;
    out.close( );
    out.open((params.mcmcpath + "/mcmcmhyper.txt").c_str(),ofstream::out);
    if (!out.is_open()) { return false; }
    out << fixed << setprecision(6) << mcmcmhyper << endl;
    out.close( );
    out.open((params.mcmcpath + "/mcmcpilogl.txt").c_str(),ofstream::out);
    if (!out.is_open()) { return false; }
    out << fixed << setprecision(6) << mcmcpilogl << endl;
    out.close( );
    error = dlmcell(params.mcmcpath + "/mcmcmrates.txt",mcmcmtiles,mcmcmRates); if (error) { return(error); }
    error = dlmcell(params.mcmcpath + "/mcmcxcoord.txt",mcmcmtiles,mcmcxCoord); if (error) { return(error); }
    error = dlmcell(params.mcmcpath + "/mcmcycoord.txt",mcmcmtiles,mcmcyCoord); if (error) { return(error); }
    error = dlmcell(params.mcmcpath + "/mcmcqrates.txt",mcmcqtiles,mcmcqRates); if (error) { return(error); }
    error = dlmcell(params.mcmcpath + "/mcmcwcoord.txt",mcmcqtiles,mcmcwCoord); if (error) { return(error); }
    error = dlmcell(params.mcmcpath + "/mcmczcoord.txt",mcmcqtiles,mcmczCoord); if (error) { return(error); }
    error = output_current_state( ); if (error) { return(error); }
    out.open((params.mcmcpath + "/eemsrun.txt").c_str(),ofstream::out);
    if (!out.is_open( )) { return false; }
    out << "Input parameter values:" << endl << params << endl
    << "Acceptance proportions:" << endl << mcmc << endl
    << "Final log prior: " << nowpi << endl
    << "Final log llike: " << nowll << endl;
    out.close( );
    cerr << "Final log prior: " << nowpi << endl
    << "Final log llike: " << nowll << endl;
    return(error);
}

void EEMS2::check_ll_computation( ) const {
    double pi0 = eval_prior(nowmSeeds,nowmEffcts,nowmrateMu,nowmrateS2,
                            nowqSeeds,nowqEffcts,nowqrateMu,nowqrateS2,
                            nowdf);
    double ll0 = eems2_likelihood(nowmSeeds, nowmEffcts, nowmrateMu, nowqSeeds, nowqEffcts, nowqrateMu, nowdf);
    if ((abs(nowpi-pi0)/abs(pi0)>1e-12)||
        (abs(nowll-ll0)/abs(ll0)>1e-12)) {
        cerr << "[EEMS2::testing]   |ll0-ll|/|ll0| = " << abs(nowll - ll0)/abs(ll0) << endl;
        cerr << "[EEMS2::testing]   |pi0-pi|/|pi0| = " << abs(nowpi - pi0)/abs(pi0) << endl;
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
    //if (qEffcts.cwiseAbs().minCoeff()>params.qEffctHalfInterval) { inrange = false; }
    //if (mEffcts.cwiseAbs().minCoeff()>params.mEffctHalfInterval) { inrange = false; }
    if (qEffcts.minCoeff() < -100 || qEffcts.maxCoeff() > params.qEffctUpperBound) { inrange = false; }
    if (mEffcts.minCoeff() < -100 || mEffcts.maxCoeff() > params.mEffctUpperBound) { inrange = false; }
    if (mrateMu>params.mrateMuUpperBound || mrateMu < -100) { inrange = false; }
    if (qrateMu>params.qrateMuUpperBound || qrateMu < -100) { inrange = false; }
    if ((df<params.dfmin) || (df>params.dfmax)) { inrange = false; }
    if (!inrange) { return (-Inf); }
    
    // remove df when done
    double logpi =
    + dnegbinln(mtiles,params.negBiSize,params.negBiProb)
    + dnegbinln(qtiles,params.negBiSize,params.negBiProb)
    + dinvgamln(mrateS2,params.mrateShape_2,params.mrateScale_2)
    + dinvgamln(qrateS2,params.qrateShape_2,params.qrateScale_2);
    for (int i = 0 ; i < qtiles ; i++) {
        //logpi += dtrnormln(qEffcts(i),0.0,qrateS2,params.qEffctHalfInterval);
        logpi += dtrnormln(qEffcts(i),qrateMu,qrateS2,params.qEffctUpperBound);
    }
    for (int i = 0 ; i < mtiles ; i++) {
        //logpi += dtrnormln(mEffcts(i),0.0,mrateS2,params.mEffctHalfInterval);
        logpi += dtrnormln(mEffcts(i),mrateMu,mrateS2,params.mEffctUpperBound);
    }
    return (logpi);
}

void EEMS2::calculateIntegral(const MatrixXd &M, const VectorXd &W, MatrixXd &lambda, double L, double r) const {
    
    // weights for the gaussian quadrature
    VectorXd w(30);
    // abisca for the gaussian quadrature
    VectorXd x(30);
    
    computeWeights(w, x, r, L);
    
    SelfAdjointEigenSolver<MatrixXd> es;
    es.compute(M);
    VectorXd eigenvalues = es.eigenvalues();
    MatrixXd V = es.eigenvectors();
    MatrixXd Dt(d, d);
    MatrixXd p = MatrixXd::Zero(o,o);
    MatrixXd P(d,d);
    lambda.setZero();
    
    // To Do: clean up this code
    for (int t = 0; t < x.size(); t++){
        Dt = ((VectorXd)((eigenvalues.array() * x[t]).exp())).asDiagonal();
        P = V*Dt*V.transpose();
        for (int i = 0; i < o; i++){
            for (int j = i; j < o; j++){
                // make P.row(i) a column vector since W.array() is a column vector
                // To do: Can vectorize this loop so p = ((MatrixXd) W'.array() * P.array()) * P^T
                // where W' = [W; W; W; ... W]
                p(i,j) = (W.array() * P.row(i).transpose().array() * P.row(j).transpose().array()).sum();
     
                //double temp = 0;
                //for (int k = 0; k < d; k++){
                //    temp += P(i,k)*P(j,k)*W(k);
                //}
                //p(i,j) = temp;
                lambda(i,j) += p(i,j)*w(t);
                lambda(j,i) = lambda(i,j);
            }
        }
    }
    
    // To do: use params.genomeSize
    // 3e9 is genomeSize, should read it in
    lambda = lambda*(3e9);
}

double EEMS2::eems2_likelihood(const MatrixXd &mSeeds, const VectorXd &mEffcts, const double mrateMu,
                               const MatrixXd &qSeeds, const VectorXd &qEffcts, const double qrateMu,
                               const double df) const {
    
    // Important: Do not use any of the current parameter values in this function,
    // i.e., those that start with nowXXX
    
    // mSeeds, mEffcts and mrateMu define the migration Voronoi tessellation
    // qSeeds, qEffcts define the diversity Voronoi tessellation
    VectorXi mColors, qColors;
    // For every deme in the graph -- which migration tile does the deme fall into?
    graph.index_closest_to_deme(mSeeds,mColors);
    // For every deme in the graph -- which diversity tile does the deme fall into?
    graph.index_closest_to_deme(qSeeds,qColors);
    VectorXd W = VectorXd::Zero(d);
    // Transform the log10 diversity parameters into diversity rates on the original scale
    for ( int alpha = 0 ; alpha < d ; alpha++ ) {
        //double log10q_alpha = qEffcts(qColors(alpha)) + qrateMu;
        //W(alpha) = pow(10.0,log10q_alpha);
        double q_alpha = pow(10.0, qEffcts(qColors(alpha)));
        W(alpha) = q_alpha;
    }
    
    // To Do: make this a sparse matrix to save space
    MatrixXd M = MatrixXd::Zero(d,d);
    int alpha, beta;
    // Transform the log10 migration parameters into migration rates on the original scale
    for ( int edge = 0 ; edge < graph.get_num_edges() ; edge++ ) {
        graph.get_edge(edge,alpha,beta);
        //double log10m_alpha = mEffcts(mColors(alpha)) + mrateMu;
        //double log10m_beta = mEffcts(mColors(beta)) + mrateMu;
        //M(alpha,beta) = 0.5 * pow(10.0,log10m_alpha) + 0.5 * pow(10.0,log10m_beta);
        double m_alpha = mEffcts(mColors(alpha));
        double m_beta = mEffcts(mColors(beta));
        M(alpha,beta) = 0.5 * pow(10.0, m_alpha) + 0.5 * pow(10.0, m_beta);
        M(beta,alpha) = M(alpha,beta);
    }
    M.diagonal() = -1* M.rowwise().sum();
    // FOR TESTING ONLY
    

    /*
    M.setZero();
    W.setZero();
    M(0,1) = M(1,0) = 0.099962;
    M(0,2) = M(2,0) = 0.099962;
    M(1,2) = M(2,1) = 0.099962;
    M(1,3) = M(3,1) = 0.099962;
    M(2,3) = M(3,2) = 0.099962;
    M.diagonal() = -1* M.rowwise().sum();

    W(0) = 0.00005;
    W(1) = 0.00005;
    W(2) = 0.00005;
    W(3) = 0.00005;
    */

    double r = 1e-8;
    MatrixXd lambda(o,o);
    double cutOff = 4e6;
    //clock_t begin_time = clock();
    calculateIntegral(M, W, lambda, cutOff, r);
    //cout << float( clock () - begin_time ) /  CLOCKS_PER_SEC << "\n\n" << endl;
    
    //cout << "OBSERVED:\n\n\n" << totalSharingM.array() / cMatrix.array() << endl;
    
    /*cout << "Migration:\n" << M << endl;
    cout << "Coalescent rates:\n" << W << endl;
    cout << "Average IBD:\n\n " << lambda << endl;
     */
 
    //cout << lambda << endl;
    //cout << "EXPECTED:\n\n\n" << lambda << endl;

    double logll = poisln(lambda, totalSharingM, cMatrix);
    //double logll = -1;
    //cout << logll << endl;
    if (logll != logll){
        cout << "trouble with ll" << endl;
        throw std::exception();
    }
    return (logll);
}
double EEMS2::getMigrationRate(const int edge) const {
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
    double m_alpha = pow(10, nowmEffcts(nowmColors(alpha)));
    double m_beta = pow(10, nowmEffcts(nowmColors(beta)));
    return (0.5 * m_alpha + 0.5 * m_beta);
    
}
double EEMS2::getCoalescenceRate(const int alpha) const {
    int nDemes = graph.get_num_total_demes();
    assert((alpha>=0) && (alpha<nDemes));
    /*
     It will be highly inefficient to compute the mapping qColors every time we want to
     look up the coalescence rate of a deme (Therefore, nowqColors is updated after any
     update that can change the mapping of vertices/demes to qVoronoi tiles)
     VectorXi qColors;
     graph.index_closest_to_deme(nowqSeeds,qColors);  
     */
    double q_alpha = pow(10, nowqEffcts(nowqColors(alpha)));
    return (q_alpha);
}

void EEMS2::printMigrationAndCoalescenceRates( ) const {
    
    int nDemes = graph.get_num_total_demes();
    cout << "Here is the coalescence (actually, diversity) rate of each deme" << endl;
    for ( int alpha = 0 ; alpha<nDemes ; alpha++ ) {
        cout << "  deme = " << alpha << ", q rate = " << getCoalescenceRate(alpha) << endl;
    }
    int nEdges = graph.get_num_edges();
    cout << "Here is the migration rate of each edge" << endl;
    for ( int edge = 0 ; edge<nEdges ; edge++ ) {
        cout << "  edge = " << edge << ", m rate = " << getMigrationRate(edge) << endl;
    }
    int alpha, beta;
    cout << "Here is the migration rate of each edge, this time edges as pairs of demes" << endl;
    for ( int edge = 0 ; edge<nEdges ; edge++ ) {
        graph.get_edge(edge,alpha,beta);
        cout << "  edge = (" << alpha << "," << beta << "), m rate = " << getMigrationRate(edge) << endl;
    }

}
