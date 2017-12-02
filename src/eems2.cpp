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
string EEMS2::olderpath( ) const { return params.olderpath; }
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

void EEMS2::rnorm_effects(const double rateS2, VectorXd &Effcts) {
    for (int i = 0 ; i < Effcts.rows() ; i++ ) {
        Effcts(i) = draw.rnorm(0.0,rateS2);
    }
}

void EEMS2::initialize_sims( ) {
    cerr << "[Sims::initialize]" << endl;
    MatrixXd sims = readMatrixXd(params.datapath + ".sims");
    Sims = sims.cast <int> ();
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

void EEMS2::initialize_state(const MCMC &mcmc) {
    cerr << "[EEMS2::initialize_state]" << endl;
    
    MatrixXd observed_means = observedIBD.array() / cMatrix.array();
    neffective = MatrixXd::Zero(o,o);
    double var;
    VectorXi indiv2deme = graph.get_indiv2deme();
    
    
    cerr << "[EEMS2::bootstrapping the data to compute variance]" << endl;
    vector<double> x;
    vector<double> y;
    double epsilon_mean = 1e-8;
    double epsilon_var = 1e-10;
    for (int alpha = 0; alpha < o; alpha++){
        for (int beta = alpha; beta < o; beta++){
            var = get_bootstrap_var(Sims, cvec, indiv2deme, 1000, alpha, beta);
            if (observed_means(alpha,beta) < epsilon_mean || var < epsilon_var){
                neffective(alpha,alpha) = -1;
            } else{
                neffective(alpha, beta) = observed_means(alpha,beta) / var;
                x.push_back(cMatrix(alpha,beta));
                y.push_back(neffective(alpha,beta));
            }
         
            neffective(beta,alpha) = neffective(alpha,beta);
        }
    }
    // now fill in missing data
    double sx = 0.0, sy = 0.0, sxx = 0.0, sxy = 0.0;
    int n = x.size();
    for (int i = 0; i < n; ++i)
    {
        sx += x[i];
        sy += y[i];
        sxx += x[i]*x[i];
        sxy += x[i]*y[i];
    }
    double delta = n*sxx - sx*sx;
    double slope = (n*sxy - sx*sy)/delta;
    // intercept is usually around zero
    double intercept = (sxx*sy - sx*sxy)/delta;
    
    for (int alpha = 0; alpha < o; alpha++){
        for (int beta = alpha; beta < o; beta++){
            if (neffective(alpha,beta) <= 0){
                neffective(alpha,beta) = slope * (float) cMatrix(alpha,beta);
                neffective(beta,alpha) = neffective(alpha,beta);
            }
            
            if (neffective(alpha,beta) < 0){
                cerr << "ERROR, effective number of samples = 0 ";
                exit(1);
            }
            
        }
    }
    
    // Initialize the two Voronoi tessellations
    nowqtiles = draw.rnegbin(2*o,0.5); // o is the number of observed demes
    nowmtiles = draw.rnegbin(2*o,0.5);
    cerr << "  EEMS starts with " << nowqtiles << " qtiles and " << nowmtiles << " mtiles" << endl;
    // Draw the Voronoi centers Coord uniformly within the habitat
    
    
    nowqtiles = o;
    nowqSeeds = MatrixXd::Zero(nowqtiles,2); //randpoint_in_habitat(nowqSeeds);
    nowqSeeds = graph.get_the_obsrv_demes();
    
    nowmSeeds = MatrixXd::Zero(nowmtiles,2); randpoint_in_habitat(nowmSeeds);
    
    
    // CHANGE THESE
    nowmrateS = -1;
    nowqrateS = -1;
    
    double abs_qrateS2 = exp(2 * nowqrateS);
    double abs_mrateS2 = exp(2 * nowmrateS);
    
    int niters = mcmc.num_iters_to_save();
    mRates = MatrixXd::Zero(niters, d);
    qRates = MatrixXd::Zero(niters, d);
    
    log10_old_mMeanRates = VectorXd::Zero(d);
    log10_old_qMeanRates = VectorXd::Zero(d);
    
    if (!params.olderpath.empty()){
        cout << "Previous Run " << params.olderpath << endl;
        load_older_rates();
    }
    
    nowmrateMu = -2 + draw.runif() * 4;
    nowqrateMu = -2 + draw.runif() * 4;
    
    // Assign rates to the Voronoi tiles
    nowqEffcts = VectorXd::Zero(nowqtiles); rnorm_effects(abs_qrateS2,nowqEffcts);
    nowmEffcts = VectorXd::Zero(nowmtiles); rnorm_effects(abs_mrateS2,nowmEffcts);
    
    // Initialize the mapping of demes to qVoronoi tiles
    graph.index_closest_to_deme(nowqSeeds,nowqColors);
    // Initialize the mapping of demes to mVoronoi tiles
    graph.index_closest_to_deme(nowmSeeds,nowmColors);
    cerr << "[EEMS2::initialize_state] Done." << endl << endl;
}

void EEMS2::store_rates(const MCMC &mcmc) {
    
    int iter = mcmc.to_save_iteration( );

    VectorXi mColors, qColors;
    // For every deme in the graph -- which diversity tile does the deme fall into?
    graph.index_closest_to_deme(nowqSeeds,qColors);
    // For every deme in the graph -- which migration tile does the deme fall into?
    graph.index_closest_to_deme(nowmSeeds,mColors);
    
    for ( int alpha = 0 ; alpha < d ; alpha++ ) {
        
        // coalescent rates
        double log10q_alpha = nowqEffcts(qColors(alpha)) + nowqrateMu + log10_old_qMeanRates(alpha);
        qRates(iter, alpha) = pow(10.0, log10q_alpha);
        
        // migration rates
        double log10m_alpha = nowmEffcts(mColors(alpha)) + nowmrateMu + log10_old_mMeanRates(alpha);
        mRates(iter, alpha) = pow(10.0, log10m_alpha);
    }
   
}

bool EEMS2::write_rates() {
    
    ofstream out; bool error = false;
    out.open((params.mcmcpath + "/qRates.txt").c_str(),ofstream::out);
    if (!out.is_open()) { error = true; return(error); }
    out << fixed << setprecision(6) << qRates << endl;
    out.close( );
    
    out.open((params.mcmcpath + "/mRates.txt").c_str(),ofstream::out);
    if (!out.is_open()) { error = true; return(error); }
    out << fixed << setprecision(6) << mRates << endl;
    out.close( );
    
    return(error);
    
}

void EEMS2::load_older_rates( ){
    
    cerr << "[EEMS2::load_older_rates]" << endl;
    MatrixXd older_rates; bool error = false;
    
    older_rates = readMatrixXd(params.olderpath + "/mRates.txt");
    if ((older_rates.rows()<1) || (older_rates.cols()<1)) { error = true; }
    for (int i = 0; i < older_rates.cols(); i++){
        log10_old_mMeanRates(i) = log10(older_rates.col(i).sum()/older_rates.col(i).size());
    }
    
    
    older_rates = readMatrixXd(params.olderpath + "/qRates.txt");
    if ((older_rates.rows()<1) || (older_rates.cols()<1)) { error = true; }
    for (int i = 0; i < older_rates.cols(); i++){
        log10_old_qMeanRates(i) = log10(older_rates.col(i).sum()/older_rates.col(i).size());
    }
    
    
    if (error) {
        cerr << "  Error loading older means from " << params.olderpath << endl; exit(1);
    }
    cerr << "[EEMS::load_older_rates] Done." << endl << endl;
    
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
    
    /*
    tempi = readMatrixXd(params.prevpath + "/lastthetas.txt");
    if ((tempi.rows()!=1) || (tempi.cols()!=1)) { error = true; }
    nowdf = tempi(0,0);
    tempi = readMatrixXd(params.prevpath + "/lastdfpars.txt");
    if ((tempi.rows()!=1) || (tempi.cols()!=2)) { error = true; }
    params.dfmin = tempi(0,0);
    params.dfmax = tempi(0,1);
     */
    tempi = readMatrixXd(params.prevpath + "/lastmhyper.txt");
    if ((tempi.rows()!=1) || (tempi.cols()!=2)) { error = true; }
    nowmrateMu = tempi(0,0);
    nowmrateS = tempi(0,1);
    tempi = readMatrixXd(params.prevpath + "/lastqhyper.txt");
    if ((tempi.rows()!=1) || (tempi.cols()!=2)) { error = true; }
    nowqrateMu = tempi(0,0);
    nowqrateS = tempi(0,1);
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
    // Initialize the mapping of demes to qVoronoi tiles
    graph.index_closest_to_deme(nowmSeeds,nowmColors);
    // Initialize the mapping of demes to mVoronoi tiles
    graph.index_closest_to_deme(nowqSeeds,nowqColors);
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
    nowpi = eval_prior(nowmSeeds,nowmEffcts,nowmrateMu,nowmrateS,
                       nowqSeeds,nowqEffcts,nowqrateMu,nowqrateS);
    nowll = eems2_likelihood(nowmSeeds, nowmEffcts, nowmrateMu, nowqSeeds, nowqEffcts, nowqrateMu, true);
    cerr << "Input parameters: " << endl << params << endl
    << "Initial log prior: " << nowpi << endl
    << "Initial log llike: " << nowll << endl << endl;
    if ((nowpi==-Inf) || (nowpi==Inf) || (nowll==-Inf) || (nowll==Inf)) { error = true; }
    return(error);
}
MoveType EEMS2::choose_move_type( ) {
    double u1 = draw.runif( );
    double u2 = draw.runif( );
    double u3 = draw.runif( );
    // There are 4 types of proposals:
    // * birth/death (with equal probability)
    // * move a tile (chosen uniformly at random)
    // * update the rate of a tile (chosen uniformly at random)
    // * update the mean migration rate or the mean coalescent rate (with equal probability)
    
    MoveType move = UNKNOWN_MOVE_TYPE;
    
    if (u3 < 0.25 && params.olderpath.empty()){
        if (u2 < 0.5) {
            move = M_MEAN_RATE_UPDATE;
        } else  {
            move = Q_MEAN_RATE_UPDATE;
        }
        
        return(move);
    }

    if (u1 < 0.33) {
        // Propose birth/death to update the Voronoi tessellation of the effective diversity,
        // with probability params.qVoronoiPr (which is 0.05 by default). Otherwise,
        // propose birth/death to update the Voronoi tessellation of the effective migration.
        if (u2 < params.qVoronoiPr) {
            move = Q_VORONOI_BIRTH_DEATH;
        } else {
            move = M_VORONOI_BIRTH_DEATH;
        }
    } else if (u1 < 0.66) {
        if (u2 < params.qVoronoiPr) {
            move = Q_VORONOI_POINT_MOVE;
        } else {
            move = M_VORONOI_POINT_MOVE;
        }
    } else {
        if (u2 < params.qVoronoiPr) {
            move = Q_VORONOI_RATE_UPDATE;
        } else {
            move = M_VORONOI_RATE_UPDATE;
        }
    }
    return(move);
}

double EEMS2::eval_proposal_rate_one_qtile(Proposal &proposal) const {
    return(eems2_likelihood(nowmSeeds, nowmEffcts, nowmrateMu, nowqSeeds, proposal.newqEffcts, nowqrateMu, false));
}
double EEMS2::eval_proposal_move_one_qtile(Proposal &proposal) const {
    return(eems2_likelihood(nowmSeeds, nowmEffcts, nowmrateMu, proposal.newqSeeds, nowqEffcts, nowqrateMu, false));
}
double EEMS2::eval_birthdeath_qVoronoi(Proposal &proposal) const {
    return(eems2_likelihood(nowmSeeds, nowmEffcts, nowmrateMu, proposal.newqSeeds, proposal.newqEffcts, nowqrateMu, false));
}
double EEMS2::eval_proposal_rate_one_mtile(Proposal &proposal) const {
    return(eems2_likelihood(nowmSeeds, proposal.newmEffcts, nowmrateMu, nowqSeeds, nowqEffcts, nowqrateMu, true));
    
}
double EEMS2::eval_proposal_overall_mrate(Proposal &proposal) const {
    return(eems2_likelihood(nowmSeeds, nowmEffcts, proposal.newmrateMu, nowqSeeds, nowqEffcts, nowqrateMu, true));
}
double EEMS2::eval_proposal_overall_qrate(Proposal &proposal) const {
    return(eems2_likelihood(nowmSeeds, nowmEffcts, nowmrateMu, nowqSeeds, nowqEffcts, proposal.newqrateMu, false));
}
// Propose to move one tile in the migration Voronoi tessellation
double EEMS2::eval_proposal_move_one_mtile(Proposal &proposal) const {
    return(eems2_likelihood(proposal.newmSeeds, nowmEffcts, nowmrateMu, nowqSeeds, nowqEffcts, nowqrateMu, true));
}
double EEMS2::eval_birthdeath_mVoronoi(Proposal &proposal) const {
    return(eems2_likelihood(proposal.newmSeeds, proposal.newmEffcts, nowmrateMu, nowqSeeds, nowqEffcts, nowqrateMu, true));
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
    proposal.newpi = eval_prior(nowmSeeds,nowmEffcts,nowmrateMu,nowmrateS,
                                nowqSeeds,proposal.newqEffcts,nowqrateMu,nowqrateS);
    proposal.newll = eval_proposal_rate_one_qtile(proposal);
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
    proposal.newpi = eval_prior(nowmSeeds,proposal.newmEffcts,nowmrateMu,nowmrateS,
                                nowqSeeds,nowqEffcts,nowqrateMu,nowqrateS);
    proposal.newll = eval_proposal_rate_one_mtile(proposal);
}
void EEMS2::propose_overall_mrate(Proposal &proposal) {
    // Make a random-walk Metropolis-Hastings proposal
    double newmrateMu = draw.rnorm(nowmrateMu,params.mrateMuProposalS2);
    proposal.move = M_MEAN_RATE_UPDATE;
    proposal.newmrateMu = newmrateMu;
    // If the proposed value is in range, the prior probability does not change
    // as the prior distribution on mrateMu is uniform
    // Otherwise, setting the prior and the likelihood to -Inf forces a rejection
    if (newmrateMu <= params.mrateMuUpperBound && newmrateMu >= params.mrateMuLowerBound) {
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
    if (newqrateMu <= params.qrateMuUpperBound && newqrateMu >= params.qrateMuLowerBound) {
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
        double newqEffct = draw.rnorm(nowqEffct,params.qEffctProposalS2);
        insertRow(proposal.newqSeeds,newqSeed.row(0));
        insertElem(proposal.newqEffcts,newqEffct);
        // Compute log(proposal ratio) and log(prior ratio)
        proposal.newratioln = log(pDeath/pBirth)
        - dnormln(newqEffct,nowqEffct,params.qEffctProposalS2);
        
        proposal.newpi = eval_prior(nowmSeeds,nowmEffcts,nowmrateMu,nowmrateS,
                                    proposal.newqSeeds,proposal.newqEffcts,nowqrateMu,nowqrateS)
        + log((nowqtiles+params.qnegBiSize)/(newqtiles/params.qnegBiProb));
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
        proposal.newratioln = log(pBirth/pDeath)
        + dnormln(oldqEffct,nowqEffct,params.qEffctProposalS2);
        
        proposal.newpi = eval_prior(nowmSeeds,nowmEffcts,nowmrateMu,nowmrateS,
                                    proposal.newqSeeds,proposal.newqEffcts,nowqrateMu,nowqrateS)
        + log((nowqtiles/params.qnegBiProb)/(newqtiles+params.qnegBiSize));
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
        double newmEffct = draw.rnorm(nowmEffct,params.mEffctProposalS2);
        insertRow(proposal.newmSeeds,newmSeed.row(0));
        insertElem(proposal.newmEffcts,newmEffct);
        // Compute log(prior ratio) and log(proposal ratio)
        proposal.newratioln = log(pDeath/pBirth)
        - dnormln(newmEffct,nowmEffct,params.mEffctProposalS2);
        
        proposal.newpi = eval_prior(proposal.newmSeeds,proposal.newmEffcts,nowmrateMu,nowmrateS,
                                    nowqSeeds,nowqEffcts,nowqrateMu,nowqrateS)
        + log((nowmtiles+params.mnegBiSize)/(newmtiles/params.mnegBiProb));
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
        proposal.newratioln = log(pBirth/pDeath)
        + dnormln(oldmEffct,nowmEffct,params.mEffctProposalS2);
        
        proposal.newpi = eval_prior(proposal.newmSeeds,proposal.newmEffcts,nowmrateMu,nowmrateS,
                                    nowqSeeds,nowqEffcts,nowqrateMu,nowqrateS)
        + log((nowmtiles/params.mnegBiProb)/(newmtiles+params.mnegBiSize));
    }
    proposal.move = M_VORONOI_BIRTH_DEATH;
    proposal.newmtiles = newmtiles;
    proposal.newll = eval_birthdeath_mVoronoi(proposal);
}
void EEMS2::update_hyperparams( ) {
    double SSq = nowqEffcts.squaredNorm();
    double SSm = nowmEffcts.squaredNorm();
    
    double abs_qrateS2 = draw.rtrinvgam(nowqtiles, 0.5 * SSq, 0.001, 1);
    double abs_mrateS2 = draw.rtrinvgam(nowmtiles, 0.5 * SSm, 0.001, exp(1));
    
    nowqrateS = log(abs_qrateS2)/2;
    nowmrateS = log(abs_mrateS2)/2;
    
    nowpi = eval_prior(nowmSeeds,nowmEffcts,nowmrateMu,nowmrateS,
                       nowqSeeds,nowqEffcts,nowqrateMu,nowqrateS);
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
    << "         number of qVoronoi tiles = " << nowqtiles << endl
    << "         number of mVoronoi tiles = " << nowmtiles << endl
    << "          Log prior = " << nowpi << setprecision(8) << endl
    << "          Log llike = " << nowll << setprecision(8) << endl;
}
void EEMS2::save_iteration(const MCMC &mcmc) {
    int iter = mcmc.to_save_iteration( );
    mcmcqhyper(iter,0) = nowqrateMu;
    mcmcqhyper(iter,1) = nowqrateS;
    mcmcmhyper(iter,0) = nowmrateMu;
    mcmcmhyper(iter,1) = nowmrateS;
    mcmcpilogl(iter,0) = nowpi;
    mcmcpilogl(iter,1) = nowll;
    mcmcqtiles(iter) = nowqtiles;
    mcmcmtiles(iter) = nowmtiles;
    //mcmcthetas(iter) = nowdf;
    
    for ( int t = 0 ; t < nowqtiles ; t++ ) {
        mcmcqRates.push_back(pow(10.0,nowqEffcts(t) + nowqrateMu));
    }
    for ( int t = 0 ; t < nowqtiles ; t++ ) {
        mcmcwCoord.push_back(nowqSeeds(t,0));
    }
    for ( int t = 0 ; t < nowqtiles ; t++ ) {
        mcmczCoord.push_back(nowqSeeds(t,1));
    }
    for ( int t = 0 ; t < nowmtiles ; t++ ) {
        mcmcmRates.push_back(pow(10.0,nowmEffcts(t) + nowmrateMu));
    }
    for ( int t = 0 ; t < nowmtiles ; t++ ) {
        mcmcxCoord.push_back(nowmSeeds(t,0));
    }
    for ( int t = 0 ; t < nowmtiles ; t++ ) {
        mcmcyCoord.push_back(nowmSeeds(t,1));
    }
    
    JtDhatJ += expectedIBD;
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
    /*
    out.open((params.mcmcpath + "/lastthetas.txt").c_str(),ofstream::out);
    if (!out.is_open()) { error = true; return(error); }
    out << fixed << setprecision(6) << nowdf << endl;
    out.close( );
    out.open((params.mcmcpath + "/lastdfpars.txt").c_str(),ofstream::out);
    if (!out.is_open()) { error = true; return(error); }
    out << fixed << setprecision(6) << params.dfmin << " " << params.dfmax << endl;
    out.close( );
     */
    out.open((params.mcmcpath + "/lastmhyper.txt").c_str(),ofstream::out);
    if (!out.is_open()) { error = true; return(error); }
    out << fixed << setprecision(6) << nowmrateMu << " " << nowmrateS << endl;
    out.close( );
    out.open((params.mcmcpath + "/lastqhyper.txt").c_str(),ofstream::out);
    if (!out.is_open()) { error = true; return(error); }
    out << fixed << setprecision(6) << nowqrateMu << " " << nowqrateS << endl;
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
    double pi0 = eval_prior(nowmSeeds,nowmEffcts,nowmrateMu,nowmrateS,
                            nowqSeeds,nowqEffcts,nowqrateMu,nowqrateS);
    double ll0 = eems2_likelihood(nowmSeeds, nowmEffcts, nowmrateMu, nowqSeeds, nowqEffcts, nowqrateMu, true);
    if ((abs(nowpi-pi0)/abs(pi0)>1e-12)||
        (abs(nowll-ll0)/abs(ll0)>1e-12)) {
        cerr << "[EEMS2::testing]   |ll0-ll|/|ll0| = " << abs(nowll - ll0)/abs(ll0) << endl;
        cerr << "[EEMS2::testing]   |pi0-pi|/|pi0| = " << abs(nowpi - pi0)/abs(pi0) << endl;
        exit(1);
    }
}
double EEMS2::eval_prior(const MatrixXd &mSeeds, const VectorXd &mEffcts, const double mrateMu, const double mrateS,
                         const MatrixXd &qSeeds, const VectorXd &qEffcts, const double qrateMu, const double qrateS) const {
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
    
    if (mrateMu>params.mrateMuUpperBound || mrateMu < params.qrateMuLowerBound) { inrange = false; }
    if (qrateMu>params.qrateMuUpperBound || qrateMu < params.qrateMuLowerBound ) { inrange = false; }
    if (qrateS < -6.9077 || qrateS > 0) { inrange = false;}
    if (mrateS < -6.9077 || mrateS > 1) { inrange = false;}
    if (!inrange) { return (-Inf); }
    
    double abs_qrateS2 = exp(2 * qrateS);
    double abs_mrateS2 = exp(2 * mrateS);
    
    double logpi =
    + dnegbinln(mtiles,params.mnegBiSize,params.mnegBiProb)
    + dnegbinln(qtiles,params.qnegBiSize,params.qnegBiProb);
    for (int i = 0 ; i < qtiles ; i++) {
        logpi += dnormln(qEffcts(i),0.0,abs_qrateS2);
    }
    for (int i = 0 ; i < mtiles ; i++) {
        logpi += dnormln(mEffcts(i),0.0,abs_mrateS2);
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
    
    MatrixXd coalp;
    MatrixXd P;
    //VectorXd qp = qMeanRate * VectorXd::Ones(d);
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
                               const bool ismUpdate) const {
    
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
        double log10q_alpha = qEffcts(qColors(alpha)) + qrateMu + log10_old_qMeanRates(alpha);
        q(alpha) = pow(10.0,log10q_alpha);
    }

    int alpha, beta;
    
    // perform costly eigen-decompositon only if updating migration rates
    if (ismUpdate){
        
        MatrixXd M = MatrixXd::Zero(d,d);
        //int alpha, beta;
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
    
    double logll = poisln(expectedIBD, observedIBD, neffective, cMatrix);
 
    if (logll != logll){
        cerr << "trouble with ll" << endl;
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
    double m_alpha = pow(10, nowmEffcts(nowmColors(alpha)) + nowmrateMu);
    double m_beta = pow(10, nowmEffcts(nowmColors(beta)) + nowmrateMu);
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
    double q_alpha = pow(10, nowqEffcts(nowqColors(alpha)) + nowqrateMu);
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
