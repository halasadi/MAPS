
#include "mcmc.hpp"

MCMC::MCMC(const Params &params, vector<double> temperatures_in) {
    numMCMCIter = params.numMCMCIter;
    numBurnIter = params.numBurnIter;
    numThinIter = params.numThinIter;
    currIter = 0;
    numTypes = 10;
    finished = false;
    temperatures = temperatures_in;
    nChains = temperatures.size();
    chain_no = nChains-1;
    okayMoves = vector<double>(numTypes,0);
    totalMoves = vector<double>(numTypes,0);
}
MCMC::~MCMC( ) { }
int MCMC::num_iters_to_save( ) const {
    int a = (numMCMCIter - numBurnIter) / (numThinIter + 1);
    return (a);
}
int MCMC::to_write_iteration( ) const {
    if (currIter>numBurnIter & temperatures[chain_no] == 1) {
        int a = (currIter - numBurnIter) / (numThinIter + 1);
        int b = (currIter - numBurnIter) % (numThinIter + 1);
        if (b==0) { return (a-1); }
    }
    return (-1);
}


bool MCMC::to_print_iteration( ) const {
    if ((currIter % 1000) == 0){
        return(true);
    } else{
        return(false);
    }
}

int MCMC::to_save_iteration( ) const {
    if (currIter > (numBurnIter/2)){
        return(1);
    }
    return(-1);
}

void MCMC::restart(const Params &params, int chain_in){
    chain_no = chain_in;
    if (temperatures[chain_no] == 1){
        numMCMCIter = params.numMCMCIter;
    } else{
        numMCMCIter = numBurnIter;
    }
    currIter = 0;
    finished = false;
    okayMoves.clear();
    okayMoves = vector<double>(numTypes,0);
    totalMoves.clear();
    totalMoves = vector<double>(numTypes,0);
}

ostream& operator<<(ostream& out, const MCMC& mcmc) {
    for ( int i = 0 ; i < mcmc.numTypes ; i++ ) {
        double a = mcmc.okayMoves.at(i);
        double A = mcmc.totalMoves.at(i);
        out << setprecision(2) << "\t(" << (int)a << "/" << (int)A << ") = " << 100.0*(a/A) << "% for type ";
        switch (i) {
            case Q_VORONOI_RATE_UPDATE:
                out << "\"qTileRate\"" << endl;
                break;
            case Q_VORONOI_POINT_MOVE:
                out << "\"qTileMove\"" << endl;
                break;
            case Q_VORONOI_BIRTH_DEATH:
                out << "\"qBirthDeath\"" << endl;
                break;
            case M_VORONOI_RATE_UPDATE:
                out << "\"mTileRate\"" << endl;
                break;
            case M_MEAN_RATE_UPDATE:
                out << "\"mMeanRate\"" << endl;
                break;
            case M_VORONOI_POINT_MOVE:
                out << "\"mTileMove\"" << endl;
                break;
            case M_VORONOI_BIRTH_DEATH:
                out << "\"mBirthDeath\"" << endl;
                break;
            case Q_MEAN_RATE_UPDATE:
                out << "\"qMeanRate\"" << endl;
                break;
            case DF_UPDATE:
                out << "\"d.f.\"" << endl;
                break;
            case CHAIN_SWAP:
                out << "\"chainSwap\"" << endl;
                break;
            default:
                cerr << "[RJMCMC] Unknown move type" << endl;
                exit(1);
        }
    }
    return out;
}
void MCMC::end_iteration( ) {
    if (++currIter==numMCMCIter) {
        finished = true;
    }
    if (!(currIter%1000)) {
        cerr << "Iteration " << currIter << " of " << numMCMCIter << "..." << endl;
    }
}
void MCMC::add_to_okay_moves(const int type) {
    okayMoves[type]++;
}
void MCMC::add_to_total_moves(const int type) {
    totalMoves[type]++;
}