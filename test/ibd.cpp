#include "ibd.hpp"
#include <math.h>
#include <fstream>
#include <ctime>


// number of nodes row-wise
const int nrow = 20;

// number of nodes column-wise
const int ncol = 20;

// total number of nodes
const int ndemes = ncol * nrow;

// number of states in the markov chain
// eqn = no of deme pairs + coalescent state
const int nstates = (int) (ndemes*(ndemes+1))/2 + 1;

// dimension of Krylov subspace
const int m = 10;

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
// if any of these items are -1 then => not defined
// e.g. a node in the graph located on the left corner does
// not have a left neighbor so llabel=-1
    
    
// each node is labeled by an integer from 0 to (g-1) where
// g is the number of nodes in the graph. Labeled from left to right
// for example in a graph of 8 nodes: 0-1-2-3
//                                    | | | |
//                                    4-5-6-7

    int label;
    
//  coalescent rate
    double q;
    
// following Desi's style where each deme has a migration parameter
// and migration rate between connecting deme A and deme B = mA+mB/2
//  migration rate:
    double m;
    
// left label, right label, top label, bottom label
    int llabel;
    int rlabel;
    int tlabel;
    int blabel;
    
};

void makeGraph (node nodes[], VectorXd &mrates, VectorXd &qrates, const int nrow, const int ncol){
    // for now, just constant migration rates and coalescent rates
    // and simple rectangular grid

    for (int i = 0, length = ndemes; i < length; i++){
        node &myNode = nodes[i];
        
        myNode.label = i;
        myNode.q = qrates(i);
        myNode.m = mrates(i);
        
        myNode.llabel = i-1;
        myNode.rlabel = i+1;
        myNode.tlabel = i-ncol;
        myNode.blabel = i+ncol;
        
        
        // node on the left edge
        if (!(i % ncol)){
            myNode.llabel = -1;
        }
        // node on the right edge
        if (!( (i+1) % ncol )){
            myNode.rlabel = -1;
        }
        // node on the top
        if (i < ncol){
            myNode.tlabel = -1;
        }
        // node on the bottom
        if ((i+1) > ncol*(nrow-1)){
            myNode.blabel = -1;
        }
        
    }
    
    
}

void calculateProduct(node nodes[], VectorXd &z, VectorXd &q){
    
    double sum;
    int index = 0;
    int state;
    double mi;
    double mj;
    node demei;
    node demej;
    
    // i < (nstates-1) because the last row of Q is all zeros
    for (int i =0; i < ndemes; i++){
        for (int j = i; j < ndemes; j++){
            demei = nodes[i];
            demej = nodes[j];
        
            mi = demei.m;
            mj = demej.m;
            sum = 0.0;
        
        // look at first lineage. e.g. (i,j) then look at movements of lineage i
            if (demei.llabel != -1){
                state = revLookup(demei.llabel, demej.label);
                sum += 0.5*(nodes[demei.llabel].m + mi);
                z[index] += 0.5*(nodes[demei.llabel].m + mi) * q(state);
            }
            if (demei.rlabel != -1){
                state = revLookup(demei.rlabel, demej.label);
                sum += 0.5*(nodes[demei.rlabel].m + mi);
                z[index] += 0.5*(nodes[demei.rlabel].m + mi) * q(state);
            }
            if (demei.tlabel != -1){
                state = revLookup(demei.tlabel, demej.label);
                sum += 0.5*(nodes[demei.tlabel].m + mi);
                z[index] += 0.5*(nodes[demei.tlabel].m + mi) * q(state);

            }
            if (demei.blabel != -1){
                state = revLookup(demei.blabel, demej.label);
                sum += 0.5*(nodes[demei.blabel].m + mi);
                z[index] += 0.5*(nodes[demei.blabel].m + mi) * q(state);
            
            }
        
            // look at second lineage
            if (demej.llabel != -1){
                state = revLookup(demei.label, demej.llabel);
                sum += 0.5*(nodes[demej.llabel].m + mj);
                z[index] += 0.5*(nodes[demej.llabel].m + mj) * q(state);
            
            }
            if (demej.rlabel != -1){
                state = revLookup(demei.label, demej.rlabel);
                sum += 0.5*(nodes[demej.rlabel].m + mj);
                z[index] += 0.5*(nodes[demej.rlabel].m + mj) * q(state);
            
            }
            if (demej.tlabel != -1){
                state = revLookup(demei.label, demej.tlabel);
                sum += 0.5*(nodes[demej.tlabel].m + mj);
                z[index] += 0.5*(nodes[demej.tlabel].m + mj) * q(state);
            
            }
            if (demej.blabel != -1){
                state = revLookup(demei.label, demej.blabel);
                sum += 0.5*(nodes[demej.blabel].m + mj);
                z[index] += 0.5*(nodes[demej.blabel].m + mj) * q(state);
            
            }
        
            // if both in the same deme, they can coalesce
            if (demei.label == demej.label){
                // matrix indices start with 0
                sum += demei.q;
                z[index] += demei.q * q(nstates-1);
            }
        
            z[index] -= sum*q(index);
            index +=1;

        }
    }
    
    // representing dot product with zero vector
    z[nstates-1] = 0;
    
}

void krylovProj(node nodes[], MatrixXd &H, MatrixXd &Q, const int m){
    H.setZero();
    VectorXd z(nstates);
    VectorXd q(nstates);
    Q.col(0).setZero();
    Q(nstates-1, 0) = 1;
    for (int k = 1; k < m; k++){
        q = Q.col(k-1);
        z.setZero();
        calculateProduct(nodes, z, q);
        for (int i = 0; i < k; i++){
            H(i, k-1) = Q.col(i).dot(z);
            z = z - H(i, k-1) * Q.col(i);
        }
        
        // need to compare this to Matlab code
        //if (k < nstates){
            H(k, k-1) = z.norm();
            if (H(k,k-1) == 0){
                    return;
            }
            Q.col(k) = z / H(k,k-1);
        //}
        
    }
    
}

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

// For testing purposes only
void makeFullMatrix(node nodes[], MatrixXd &Q){
  // Requires: Q be NxN
  int index;
  double mi;
  double mj;
  node demei;
  node demej;
    Q.setZero();
    // i < (nstates-1) because the last row of Q is all zeros
  for (int i = 0; i < (nstates-1); i++){
    demei = nodes[lookup[i][0]];
    demej = nodes[lookup[i][1]];
        
    mi = demei.m;
    mj = demej.m;
    
    // look at first lineage. e.g. (i,j) then look at movements of lineage i
    if (demei.llabel != -1){
      index = revLookup(demei.llabel, demej.label);
      Q(i, index) += 0.5*(nodes[demei.llabel].m + mi);
    }
    if (demei.rlabel != -1){
      index = revLookup(demei.rlabel, demej.label);
      Q(i, index) += 0.5*(nodes[demei.rlabel].m + mi);
    }
    if (demei.tlabel != -1){
      index = revLookup(demei.tlabel, demej.label);
      Q(i, index) += 0.5*(nodes[demei.tlabel].m + mi);
    }
    if (demei.blabel != -1){
      index = revLookup(demei.blabel, demej.label);
      Q(i, index) += 0.5*(nodes[demei.blabel].m + mi);
    }
        
        // look at second lineage
    if (demej.llabel != -1){
      index = revLookup(demei.label, demej.llabel);
      Q(i, index) += 0.5*(nodes[demej.llabel].m + mj);
    }
    if (demej.rlabel != -1){
      index = revLookup(demei.label, demej.rlabel);
      Q(i, index) += 0.5*(nodes[demej.rlabel].m + mj);
    }
    if (demej.tlabel != -1){
      index = revLookup(demei.label, demej.tlabel);
      Q(i, index) += 0.5*(nodes[demej.tlabel].m + mj);
      
    }
    if (demej.blabel != -1){
      index = revLookup(demei.label, demej.blabel);
      Q(i, index) += 0.5*(nodes[demej.blabel].m + mj);
    }
    
    // if both in the same deme, they can coalesce
    if (demei.label == demej.label){
      // matrix indices start with 0
      Q(i,nstates-1) += demei.q;
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

void calculateIntegralKrylov(node nodes[], double r, double L, VectorXd &approximate){
    MatrixXd Q(nstates, m);
    MatrixXd H(m, m);
    krylovProj(nodes, H, Q, m);
    
    // abisca for the gaussian quadrature
    VectorXd w(30);
    VectorXd x(30);
    computeWeights(w, x, r, L);
    // where to store the matrix exponential
    MatrixXd E(m,m);
    
    // to get the last column
    VectorXd l = VectorXd::Zero(nstates);
    l[nstates-1] = 1.0;
    
    // storing the probabilities
    MatrixXd Papprox(nstates, 30);
    
    MatrixXd Ht(m,m);
    for (int i = 0; i < x.size() ; i++){
        Ht = H*x[i];
        padm(Ht, E);
        Papprox.col(i) = (Q*E)*(Q.transpose()*l);
    }
    
    int cnt = 0;
    VectorXd papprox(30);
    int state;
    for (int i = 0; i < ndemes; i++){
        for (int j = i; j < ndemes; j++)
        {
            state = revLookup(i,j);
            // estimate the probability mass function
            papprox(0) = 0;
            papprox.tail(29) = (Papprox.row(state).tail(29)- Papprox.row(state).head(29)).array()/(x.tail(29)-x.head(29)).transpose().array();
            approximate[cnt] = w.dot(papprox);
            cnt += 1;
        }
    }
    
}

void calculateIntegral(node nodes[], double r, double L, VectorXd &exact){
    VectorXd w(30);
    VectorXd x(30);
    computeWeights(w, x, r, L);

    MatrixXd E(nstates, nstates);
    MatrixXd A(nstates, nstates);
    makeFullMatrix(nodes, A);
    
    // to get the last column
    VectorXd l = VectorXd::Zero(nstates);
    l[nstates-1] = 1.0;
    
    MatrixXd At(nstates, nstates);
    MatrixXd Pexact(nstates, 30);
    for (int i = 0; i < x.size() ; i++){
        At = A*x[i];
        padm(At, E);
        Pexact.col(i) = E*l;
    }
    
    int cnt = 0;
    VectorXd pexact(30);
    int state;
    for (int i = 0; i < ndemes; i++){
        for (int j = i; j < ndemes; j++)
        {
            state = revLookup(i,j);
            pexact(0) = 0;
            pexact.tail(29) = (Pexact.row(state).tail(29)- Pexact.row(state).head(29)).array()/(x.tail(29)-x.head(29)).transpose().array();
            exact[cnt] = w.dot(pexact);
            cnt += 1;
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
    double r = 1e-8;
    double L = 2e6;


    // array of nodes
    node nodes[ndemes];
    
    const int nreps = 1;
    for (int i = 0; i < nreps; i++){
        VectorXd mrates(ndemes);
        VectorXd qrates(ndemes);
        
        //mrates << 0.00532767, 0.00218959, 0.000470446, 0.00678865;
        //qrates << 7.82637e-09, 0.000131538, 0.000755605,  0.00045865;
        
        VectorXd ones = VectorXd::Ones(ndemes);
        qrates.setRandom(ndemes);
        mrates.setRandom(ndemes);
        qrates = (0.001/2) * (qrates + ones);
        mrates = (0.01/2) * (mrates + ones);

        
        // initalize graph
        makeGraph(nodes, mrates, qrates, nrow, ncol);
        
 
        VectorXd approximate(nstates-1);
        //VectorXd exact(nstates-1);
        //clock_t begin_time = clock();
        //calculateIntegral(nodes, r, L, exact);
        //cout << "For Full: " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;;
        clock_t begin_time = clock();
        calculateIntegralKrylov(nodes, r, L, approximate);
        cout << "For Krylov: " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;;
        //cout << "error: " << (approximate-exact).norm() << "\n\n\n" << endl;

        
        
         }
}
