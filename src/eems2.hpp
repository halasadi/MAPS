#pragma once

#include "util.hpp"
#include "mcmc.hpp"
#include "draw.hpp"

#include "graph.hpp"
#include "habitat.hpp"

#ifndef EEMS2_H
#define EEMS2_H

/*
 An updated set of parameter values
 The type of move is necessary in order to know which parameters have a new proposed value;
 the rest of the parameters won't be set to their current values (to avoid unnecessary copying)
 For example, if move = M_VORONOI_BIRTH_DEATH,
 then newmtiles, newmSeeds, nowmColors, newmEffcts (and of course, newpi, newll, newratioln) would have changed
 if move = M_VORONOI_POINT_MOVE,
 then newmSeeds, nowmColors (and of course, newpi, newll) would have changed
 The ratioln is the proposal ratio for birth/death proposal.
 For the usual Metropolis-Hastings updates, the acceptance probability is
 alpha = (prior ratio) * (likelihood ratio)
 For the birth/deatch RJ-MCMC updates, the acceptance probability is
 alpha = (proposal ratio) * (prior ratio) * (likelihood ratio)
 See Green, "Reversible jump Markov chain Monte Carlo computation and Bayesian model determination"
 */
struct Proposal {
    MoveType move; // the type of proposal/update
    int newqtiles; // number of m and q tiles, respectively
    int newmtiles;
    
    double newmrateS, newqrateS;
    double newpi; // log prior
    double newll; // log likelihood
    double newratioln; // RJ-MCMC proposal ratio, on the log scale
    double newmrateMu; // overall (mean) migration rate,
    double newqrateMu;
    
    VectorXd newqEffcts; // the diversity rate of each q tile
    VectorXd newmEffcts; // the migration rate of each m tile, relative to the ovarall mrateMu
    MatrixXd newqSeeds;  // the location of each q tile within the habitat
    MatrixXd newmSeeds;  // the location of each m tile within the habitat
    
};

class EEMS2 {
public:
    
    EEMS2(const Params &params);
    ~EEMS2( );
    
    void initialize_state(const MCMC &mcmc);
    void load_final_state( );
    void load_rates( );

    bool start_eems(const MCMC &mcmc);
    double eval_prior(const MatrixXd &mSeeds, const VectorXd &mEffcts, const double mrateMu, const double mrateS,
                      const MatrixXd &qSeeds, const VectorXd &qEffcts, const double qrateMu, const double qrateS) const;
    double eems2_likelihood(MatrixXd newmSeeds, MatrixXd newqSeeds, VectorXd newmEffcts,
                            VectorXd newqEffcts, double newmrateMu, double newdf, bool ismUpdate,
                            double nowmrateS, double nowqrateS) const;
    
    void calculateIntegral(MatrixXd &eigenvals, MatrixXd &eigenvecs, const VectorXd &q,
                           MatrixXd &integral, double bnd) const;
    
    MoveType choose_move_type(const MCMC &mcmc);
    // These functions change the within demes component:
    double eval_proposal_rate_one_qtile(Proposal &proposal) const;
    double eval_proposal_move_one_qtile(Proposal &proposal) const;
    double eval_birthdeath_qVoronoi(Proposal &proposal) const;
    // These functions change the between demes component:
    double eval_proposal_rate_one_mtile(Proposal &proposal) const;
    double eval_proposal_overall_mrate(Proposal &proposal) const;
    double eval_proposal_overall_qrate(Proposal &proposal) const;
    double eval_proposal_move_one_mtile(Proposal &proposal) const;
    double eval_birthdeath_mVoronoi(Proposal &proposal) const;
    
    // Random-walk Metropolis-Hastings proposals:
    void propose_omegaq(Proposal &proposal);
    void propose_omegam(Proposal &proposal);
    void propose_rate_one_qtile(Proposal &proposal);
    void propose_rate_one_mtile(Proposal &proposal);
    void propose_overall_mrate(Proposal &proposal);
    void propose_overall_qrate(Proposal &proposal);
    void propose_move_one_qtile(Proposal &proposal);
    void propose_move_one_mtile(Proposal &proposal);
    void propose_birthdeath_qVoronoi(Proposal &proposal);
    void propose_birthdeath_mVoronoi(Proposal &proposal);
    bool accept_proposal(Proposal &proposal, const MCMC &mcmc);
    
    void print_iteration(const MCMC &mcmc) const;
    void save_iteration(const MCMC &mcmc);
    bool output_results(const MCMC &mcmc) const;
    bool output_current_state() const;
    void check_ll_computation() const;
    string datapath() const;
    string mcmcpath() const;
    string prevpath() const;
    string olderpath() const;
    string gridpath() const;
    void store_rates(const MCMC &mcmc);
    bool write_rates();
    
private:
    
    Draw draw; // Random number generator
    Graph graph;
    Params params;
    Habitat habitat;
    
    int o; // number of observed demes
    int d; // total number of demes
    int n; // number of samples
    MatrixXd observedIBD; // observed means (for number of IBD blocks)
    MatrixXd cMatrix; // number of pairwise observations between observed populations
    VectorXd cvec; // c is the vector of counts
    VectorXd cClasses; // cClasses is a vector of count of the number of 0's, number of 1's, etc. For likelihood;
    double maxCnt; // the maximum number of IBD segments shared (over all pairs)
    MatrixXd neffective;
    
    MatrixXd mRates;
    MatrixXd qRates;
    MatrixXi Sims;
    
    VectorXd log10_old_mMeanRates;
    VectorXd log10_old_qMeanRates;
    
    MatrixXd JtDhatJ;
    mutable MatrixXd expectedIBD;
    mutable MatrixXd eigenvals;
    mutable MatrixXd eigenvecs;
    
    // The current set of parameter values:
    int nowmtiles, nowqtiles; // number of m and q tiles, respectively
    MatrixXd nowmSeeds; VectorXd nowmEffcts; double nowmrateMu; // parameters to describe the m Voronoi tessellation
    MatrixXd nowqSeeds; VectorXd nowqEffcts;                    // parameters to describe the q Voronoi tessellation
    double nowqrateS, nowmrateS; // two hyperparameters -- the standard deviation of nowqEffcts and nowmEffcts, respectively
    //double nowsigma2, nowpi, nowll, nowdf; // variance scale, log prior, log likelihood, degrees of freedom
    double nowqrateMu, nowpi, nowll; // variance scale, log prior, log likelihood, degrees of freedom
    
    VectorXi nowqColors; // mapping that indicates which q tiles each vertex/deme falls into
    VectorXi nowmColors; // mapping that indicates which m tiles each vertex/deme falls into
    
    // Variables to store the results in:
    // Fixed size:
    
    MatrixXd mcmcmhyper;
    MatrixXd mcmcqhyper;
    MatrixXd mcmcthetas;
    MatrixXd mcmcpilogl;
    VectorXd mcmcmtiles;
    VectorXd mcmcqtiles;
    // Variable length:
    vector<double> mcmcmRates;
    vector<double> mcmcqRates;
    vector<double> mcmcxCoord;
    vector<double> mcmcyCoord;
    vector<double> mcmcwCoord;
    vector<double> mcmczCoord;
    
    void initialize_sims();
    void randpoint_in_habitat(MatrixXd &Seeds);
    void rnorm_effects(const double HalfInterval, const double rateS2, VectorXd &Effcts);
    
    double eems2_likelihood(const MatrixXd &mSeeds, const VectorXd &mEffcts, const double mrateMu,
                            const MatrixXd &qSeeds, const VectorXd &qEffcts,
                            const double qrateMu, const bool ismUpdate,
                            const double nowmrateS, const double nowqrateS) const;
};

#endif
