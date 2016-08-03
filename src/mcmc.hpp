#pragma once

#include "util.hpp"

#ifndef MCMC_H
#define MCMC_H

enum MoveType {
    CHAIN_SWAP,
    Q_VORONOI_RATE_UPDATE,
    Q_VORONOI_POINT_MOVE,
    Q_VORONOI_BIRTH_DEATH,
    Q_MEAN_RATE_UPDATE,
    M_VORONOI_RATE_UPDATE,
    M_VORONOI_POINT_MOVE,
    M_VORONOI_BIRTH_DEATH,
    M_MEAN_RATE_UPDATE,
    DF_UPDATE,
    UNKNOWN_MOVE_TYPE
};

class MCMC {
public:
    
    MCMC(const Params &params);
    ~MCMC();
    
    int currIter;
    int numMCMCIter;
    int numBurnIter;
    int numThinIter;
    bool isbaseChain;
    bool finished;
    
    void end_iteration( );
    void add_to_okay_moves(const int type);
    void add_to_total_moves(const int type);
    int num_iters_to_save( ) const;
    int to_store_iteration() const;
    int to_write_iteration( ) const;
    void restart(const Params &params, bool isbaseChain);
    
    friend ostream& operator<<(ostream& out, const MCMC& mcmc);
    
private:
    
    int numTypes;
    vector<double> okayMoves;
    vector<double> totalMoves;
    
};

#endif