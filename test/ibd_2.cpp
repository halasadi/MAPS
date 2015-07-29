#include "ibd.hpp"
#include <math.h>
#include <fstream>
#include <ctime>
#include <vector>


const int nrow = 15;
const int ncol = 15;
// total number of nodes
const int ndemes = nrow*ncol;

// number of states in the markov chain
// eqn = no of deme pairs + coalescent state
const int nstates = (int) (ndemes*(ndemes+1))/2 + 1;

void swap(int &i, int &j){
    int temp = i;
    i = j;
    j = temp;
}

// define the format you want, you only need one instance of this...
const static IOFormat CSVFormat(FullPrecision, DontAlignCols, ", ", "\n");

// a way to go from index in the Q matrix to (i,j)
// excluding the coalescent state
int lookup[nstates-1][2];

int revLookup(double i, double j){
    // subtract the one at the end because matrices in C++ start with 0
    // sum_{k=0}^{i-1) n-k  + sum_{k=i}^{k=j} 1  - 1
 
    // we restrict i <= j
    if (i > j){
        swap(i,j);
    }
    double index = 0;
    if (i >= 1){
        index += 0.5*i*(2.0*ndemes+1-i);
    }
    index += (j-i);
    return ((int) index);
}

struct node {
    vector<int> neighbors;
    int label;
};


double max(double a, double b){
    if (a > b){
        return(a);
    }
    return(b);
}

void padm(MatrixXd &H, MatrixXd &E){
    // REQUIRES: H and E to be both MatrixXd of size NxN.
    // MODIFIES: E
    // EFFECTS: calculates the exponential of H and stores it into E.
    // padm.m from expokit package translated to C++ by Hussein Al-Asadi
    
    
    const int N = H.rows();
    if (N != H.cols() || N != E.rows() || N != E.cols()){
        cout << "ERROR: DIMENSIONS OF INPUT MATRICES ARE INCONSISTENT!" << endl;
    }
    
    // recommended (6,6)-degree rational Pade approximation
    const int p = 6;
    
    // pade coefficients
    VectorXd c(p+1);
    c[0] = 1.0;
    for (int k = 1; k <= p; k++){
        c(k) = c(k-1)*((p+1.0-k)/(k*(2.0*p+1.0-k)));
    }
    
    // L-Infinity norm as defined by norm(H, 'inf') in Matlab
    double s = H.cwiseAbs().rowwise().sum().maxCoeff();
    
    if (s > 0.5){
        s = max(0, floor(log(s)/log(2)) + 2);
        H = pow(2, -s)*H;
    }
    
    // Horner evaluation of the irreducible fraction
    MatrixXd I = MatrixXd::Identity(N,N);
    MatrixXd H2(N,N);
    H2 = H*H;
    MatrixXd Q = c(p)*I;
    MatrixXd P = c(p-1)*I;
    
    int odd = 1;
    for (int k = p-1; k > 0; k--){
        if (odd == 1){
            Q = Q*H2 + c(k-1)*I;
        } else{
            P = P*H2 + c(k-1)*I;
        }
        odd = 1 - odd;
    }
    if (odd == 1){
        Q = Q*H;
        Q = Q-P;
        E = -1*(I + 2*Q.lu().solve(P));
    }
    else{
        P = P*H;
        Q = Q - P;
        E = I + 2*Q.lu().solve(P);
    }
    
    // Squaring
    // loop floor(s) times
    for (int k = 0; k < floor(s); k++){
        E = E*E;
    }
}


void writeToCSVfile(string name, MatrixXd matrix)
{
    ofstream file(name.c_str());
    file << matrix.format(CSVFormat);
}

void calculateProduct(VectorXd &z, VectorXd &q, const MatrixXd &M, const VectorXd &W, vector<node> nodes) {
    node demei;
    node demej;
    
    // sum is to keep track of the row sum. So then we can fill in the diagonals
    double sum;
    
    // sweeping across the entries of the vector z where z = A*q
    int index = 0;
    
    // going from index to (i,j)
    int state;
    
    for (int i = 0; i < ndemes; i++){
        for (int j = i; j < ndemes; j++){
            demei = nodes[i];
            demej = nodes[j];
            sum = 0.0;
            
            // let i move and fix j since only one step transitions are allowed. Need to look up the migration rate from i to the neighbor of i
            for (int k = 0; k < demei.neighbors.size(); k++){
                sum += M(i,demei.neighbors[k]);
                // we're on the "indexth" row of A and looking at the "stateth" entry of this row
                state = revLookup(demei.neighbors[k], j);
                z[index] += M(i,demei.neighbors[k])*q(state);
            }
            
            // let j move and fix i since only one step transitions are allowed. Need to look up the migration rate from j to the neightbor of j
            for (int k = 0; k < demej.neighbors.size() ; k++){
                sum += M(j,demej.neighbors[k]);
                // we're on the "indexth" row of A and looking at the "stateth" entry of this row
                state = revLookup(i, demej.neighbors[k]);
                z[index] += M(j,demej.neighbors[k])*q(state);
            }
            
            // both i and j coalescece
            if (i == j){
                sum += W(i);
                z[index] += W(i) * q(nstates-1);
            }
            
            // the diagonal
            z[index] -= sum*q(index);
            
            index += 1;
        }
    }
    z[nstates-1] = 0;
}


void krylovProj(MatrixXd &H, MatrixXd &Q, const MatrixXd &M, const VectorXd &W, vector<node> nodes, const int m) {
    // REQUIRES: dimKrylov < nstates
    // MODIFIES: H, Q
    // EFFECTS: krylov projection of the rate matrix,e.g. if  A is rate matrix then
    // finds the decomposition A=Q'HQ
    
    // set up storage
    H.setZero();
    Q.setZero();
    VectorXd z(nstates);
    VectorXd q(nstates);
    
    // initialize first kyrlov basis
    Q(nstates-1, 0) = 1;
    
    // Arnoldi iteration
    for (int k = 1; k < m; k++){
        q = Q.col(k-1);
        z.setZero();
        calculateProduct(z, q, M,W, nodes);
        for (int i = 0; i < k; i++){
            H(i, k-1) = Q.col(i).dot(z);
            z = z - H(i, k-1) * Q.col(i);
        }
        
        H(k, k-1) = z.norm();
        if (H(k,k-1) == 0){
            return;
        }
        Q.col(k) = z / H(k,k-1);
    }
}

// For testing purposes only
// order goes like this for 2 deme model:
// (1,1) -> 0
// (1,2) -> 1
// (2,2) -> 2
// C     -> 3
// Lookup and revLookup is just to go from (i,j) <-> k
void makeFullMatrix(vector<node> &nodes, const MatrixXd &M, const MatrixXd &W, MatrixXd &Q){
    int index;
    node demei;
    node demej;
    int neighbor;
    Q.setZero();
    // i < (nstates-1) because the last row of Q is all zeros
    for (int i = 0; i < (nstates-1); i++){
        demei = nodes[lookup[i][0]];
        demej = nodes[lookup[i][1]];
        
        // fix deme i and look at all the possble demes lineage i can go to (fix lineage j).
        for (int k = 0; k < demei.neighbors.size(); k++){
            neighbor = demei.neighbors[k];
            index = revLookup(neighbor, demej.label);
            Q(i, index) += M(neighbor, demei.label);
        }
        
        for (int k = 0; k < demej.neighbors.size(); k++){
            neighbor = demej.neighbors[k];
            index = revLookup(demei.label, neighbor);
            Q(i, index) += M(neighbor, demej.label);
        }
        
        if (demei.label == demej.label){
            Q(i, nstates-1) += W(demei.label);
        }
        
        Q(i,i) = -1*(Q.row(i).sum());
        
    }
}

void computeWeights(VectorXd &w, VectorXd &x, double r, double L){
    // REQUIRES: w and x vectors of length 30, r is recombination rate, and L (in base pairs) of cutoff.
    // MODIFIES: w and x
    // EFFECTS: x will contain the x-values telling you where to evaluate P(T_mrca = x); w will contains the weights
    // This function allows user to compute an integral by computing \sum_i P(T_mrca = x_i) * w_i
    
    x << 0.118440697736960550688, 0.3973475034735802657556, 0.8365549141880933313119, 1.437175158191620443607,
    2.200789508440616292336, 3.129448303166859096349, 4.225699164493802071261, 5.492626704368934083587,
    6.933903364122364597039, 8.553853192793023779194, 10.35753137020864105106, 12.35082332811269876439,
    14.54056869943518703492, 16.93471724415800802837, 19.54252664684054185266, 22.37481610233449499411,
    25.44429563058376261798, 28.76600031447167014762, 32.35787326932856805551, 36.24156497875364752439,
    40.44355691460364227197, 44.99678841355200250088, 49.94309754094208987181, 55.33704611950810443499,
    61.25224904369593075136, 67.79260716731075303985, 75.11420274687672563149, 83.47405073153149030595,
    93.36359463048878316735, 106.0462505962874034422;
    
    w << 0.02093564741472521761, 0.09585049298017654367, 0.18833296435057945936, 0.23281944819987904471,
    0.2060782293528492151, 0.138528960450616358, 0.07293919110208096649, 0.030605607903988887905,
    0.010333948458420042431, 0.002821608083735993584, 6.2402663742264620427E-4, 1.1168849922460852198E-4,
    1.6129719270580565631E-5, 1.87044426274856472768E-6, 1.72995513372709914535E-7, 1.26506996496773906645E-8,
    7.2352574135703022224E-10, 3.19320138447436406004E-11, 1.069761647687436460972E-12, 2.66597906070505518515E-14,
    4.82019019925788439097E-16, 6.12740480626441608041E-18, 5.26125812567892365789E-20, 2.89562589607893296815E-22,
    9.51695437836864011982E-25, 1.69046847745875738033E-27, 1.39738002075239812243E-30, 4.20697826929603166432E-34,
    2.89826026866498969507E-38, 1.411587124593531584E-43;
    
    // integral_0^{\inf} 2rte^(-2trL)f(t) = (1/L^2) \integral_0^{\inf} f(u/2rL) ue^(-u)
    w = w*(1/(L*2.0*r*L));
    x = x/(2*r*L);
}


void calculateIntegralKrylov(const MatrixXd &M, const MatrixXd &W, MatrixXd &lambda, double L, double r, vector<node> &nodes, const int m){
    MatrixXd Q(nstates, m);
    MatrixXd H(m, m);
    krylovProj(H, Q, M, W, nodes, m);
    
    // weights for the gaussian quadrature
    VectorXd w(30);
    // abisca for the gaussian quadrature
    VectorXd x(30);
    
    computeWeights(w, x, r, L);
    
    // where to store the matrix exponential
    MatrixXd E(m,m);
    
    // to get the last column
    VectorXd l = VectorXd::Zero(nstates);
    l[nstates-1] = 1.0;
    
    // storing the probabilities
    MatrixXd P(nstates, 30);
    
    MatrixXd Ht(m,m);
    
    for (int i = 0; i < x.size() ; i++){
        Ht = H*x[i];
        padm(Ht, E);
        P.col(i) = (Q*E)*(Q.transpose()*l);
    }
    
    VectorXd p(30);
    int state;
    for (int i = 0; i < ndemes; i++){
        for (int j = i; j < ndemes; j++)
        {
            state = revLookup(i,j);
            // estimate the probability mass function
            p(0) = 0;
            p.tail(29) = (P.row(state).tail(29)- P.row(state).head(29)).array()/(x.tail(29)-x.head(29)).transpose().array();
            // compute the integral
            // 3e9 is genome size
            lambda(i,j) = (3e9)*(w.dot(p));
            //if (lambda(i,j) < 0){
            //    throw std::exception();
            //}
            lambda(j,i) = lambda(i,j);
        }
    }
   
    
}


void calculateIntegral(const MatrixXd &M, const MatrixXd &W, MatrixXd &lambda, double L, double r, vector<node> &nodes){
    VectorXd w(30);
    VectorXd x(30);
    computeWeights(w, x, r, L);

    MatrixXd E(nstates, nstates);
    MatrixXd A(nstates, nstates);
    A.setZero();
    makeFullMatrix(nodes, M, W, A);
    
    //cout << "A: \n" << A << endl;
    
    // to get the last column
    VectorXd l = VectorXd::Zero(nstates);
    l[nstates-1] = 1.0;
    
    MatrixXd At(nstates, nstates);
    MatrixXd P(nstates, 30);
    for (int i = 0; i < x.size() ; i++){
        At = A*x[i];
        padm(At, E);
        P.col(i) = E*l;
    }
    
    VectorXd p(30);
    int state;
    for (int i = 0; i < ndemes; i++){
        for (int j = i; j < ndemes; j++)
        {
            state = revLookup(i,j);
            p(0) = 0;
            p.tail(29) = (P.row(state).tail(29)- P.row(state).head(29)).array()/(x.tail(29)-x.head(29)).transpose().array();
            lambda(i,j) = (3e9)*(w.dot(p));
            lambda(j,i) = lambda(i,j);
        }
    }
    
    
}


void populate_nodes(vector<node> &nodes, MatrixXi &DemePairs){
    // MODIFIES: nodes
    nodes.resize(ndemes);
    
    for (int i = 0; i < ndemes; i++){
        //std::cout << "neighbors of deme " << i <<  " are:";
        //for (unsigned ii=0; ii< nodes[i].neighbors.size(); ii++)
        //    std::cout << ' ' << nodes[i].neighbors[ii];
        //std::cout << '\n';
        nodes[i].neighbors.clear();
    }
    
    for (int i = 0; i < DemePairs.rows(); i++){
        int alpha = DemePairs(i,0);
        int beta = DemePairs(i, 1);
        nodes[alpha].neighbors.push_back(beta);
        nodes[beta].neighbors.push_back(alpha);
        nodes[alpha].label = alpha;
        nodes[beta].label = beta;
    }
}


double poisln(const MatrixXd &Lambda, const MatrixXd &lambda, const MatrixXd &cMatrix){
    double ll = 0;
    int n = Lambda.rows();
    for (int i = 0; i < n; i++){
        for (int j = i; j < n; j++){
            cout << "i: " << i << ", j: " << j << endl;
            ll += lambda(i,j)*log(Lambda(i,j))-cMatrix(i,j)*Lambda(i,j);
        }
    }
    cout << "logl: " << ll << endl;
    return(ll);
}



void get_edge(int edge, int &alpha, int &beta, MatrixXi &DemePairs)
{
    alpha = DemePairs(edge,0); beta = DemePairs(edge,1);
}

void make_edges(int nrow, int ncol, MatrixXi &DemePairs){
    int edge = 0;
    for (int i = 0; i < ndemes; i++){
        
        // if node not on bottom
        if ((i+1) <= ncol*(nrow-1)){
            DemePairs(edge,0) = i;
            DemePairs(edge,1) = i+ncol;
            edge += 1;
        }
        
        // if node not on the right edge
        if (((i+1) % ncol) != 0){
            DemePairs(edge,0) = i;
            DemePairs(edge,1) = i+1;
            edge += 1;
        }
        
        // if node not on the right edge AND not on bottom
        if  ((((i+1) % ncol) != 0)  & ((i+1) <= ncol*(nrow-1)) ){
            DemePairs(edge, 0) = i;
            DemePairs(edge, 1) = i+ncol + 1;
            edge += 1;
        }
    }
}

int main()
{
    // populate lookup array
    int ind = 0;
    for (int i = 0; i < ndemes; i++){
        for (int j=i; j < ndemes; j++){
            lookup[ind][0] = i;
            lookup[ind][1] = j;
            ind += 1;
        }
    }
    
    VectorXd mrates(ndemes);
    VectorXd W(ndemes);
    VectorXd ones = VectorXd::Ones(ndemes);
    W.setRandom(ndemes);
    mrates.setRandom(ndemes);
    W = (0.00001/2) * (W + ones);
    mrates = (0.1/2) * (mrates + ones);
    
    //cout << W << endl;
    //cout << mrates << endl;
    
    // Must agree with ndemes above
    // set up the graph here
    const int nedges = (ncol-1)*nrow + (nrow-1)*ncol + (ncol-1)*(nrow-1);
    MatrixXi DemePairs = MatrixXd::Zero(nedges, 2).cast <int> ();
    
    make_edges(nrow, ncol, DemePairs);
    
    /*
    // edge 0
    DemePairs(0,0) = 0;
    DemePairs(0,1) = 1;
    DemePairs(1,0) = 0;
    DemePairs(1,1) = 2;
    DemePairs(2,0) = 1;
    DemePairs(2,1) = 2;
    DemePairs(3,0) = 1;
    DemePairs(3,1) = 3;
    DemePairs(4,0) = 2;
    DemePairs(4,1) = 4;
    DemePairs(5,0) = 2;
    DemePairs(5,1) = 5;
    DemePairs(6,0) = 2;
    DemePairs(6,1) = 3;
    DemePairs(7,0) = 3;
    DemePairs(7,1) = 5;
    DemePairs(8,0) = 4;
    DemePairs(8,1) = 5;
    */
     
    vector<node> nodes;
    
    populate_nodes(nodes, DemePairs);
    
    MatrixXd M = MatrixXd::Zero(ndemes,ndemes);
    int alpha, beta;
    for ( int edge = 0 ; edge < nedges ; edge++ ) {
        get_edge(edge,alpha,beta, DemePairs);
        double m_alpha = mrates(alpha);
        double m_beta = mrates(beta);
        M(alpha,beta) = 0.5 * m_alpha + 0.5 * m_beta;
        M(beta,alpha) = M(alpha,beta);
    }
    
    /*
    M(0,1) = M(1,0) = 0.01;
    M(0,2) = M(2,0) = 0.01;
    M(1,2) = M(2,1) = 0.01;
    M(1,3) = M(3,1) = 0.01;
    M(2,3) = M(3,2) = 0.01;
    M(2,4) = M(4,2) = 0.01;
    M(2,5) = M(5,2) = 0.01;
    M(3,5) = M(5,3) = 0.01;
    M(4,5) = M(5,4) = 0.01;
    
    W(0) = 0.00000000001;
    W(1) = 0.00000000002;
    W(2) = 0.00000000003;
    W(3) = 0.00000000004;
    W(4) = 0.00000000005;
    W(5) = 0.00000000006;
    */

    
    double r = 1e-8;
    double L = 4e6;
    
    
    //cout << "W:\n" << W << endl;
    //cout << "M:\n" << M << endl;
    
    // dimension of krylov subpsace
    int m;
    m = 100;
    MatrixXd lambda = MatrixXd::Zero(ndemes, ndemes);
    calculateIntegralKrylov(M, W, lambda, L, r, nodes, m);
    cout << "almost exact lambda: \n" << lambda.row(0) << endl;
    
    lambda.setZero();
    m = 15;
    calculateIntegralKrylov(M, W, lambda, L, r, nodes, m);
 
    cout << "approx lambda:\n" << lambda.row(0) << endl;
    
}
