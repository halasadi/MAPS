#include "ibd.hpp"
#include <math.h>

// number of nodes row-wise
const int nrow = 5;

// number of nodes column-wise
const int ncol = 5;

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

// For testing purposes only
void test_FullMatrixExp(node nodes[], MatrixXd &Q, MatrixXd &E){
  // Requires: Q and E be NxN
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
      Q(i, index) += 0.5*(nodes[demei.llabel].m + mj);
    }
    if (demei.rlabel != -1){
      index = revLookup(demei.rlabel, demej.label);
      Q(i, index) += 0.5*(nodes[demei.rlabel].m + mj);
    }
    if (demei.tlabel != -1){
      index = revLookup(demei.tlabel, demej.label);
      Q(i, index) += 0.5*(nodes[demei.tlabel].m + mj);
    }
    if (demei.blabel != -1){
      index = revLookup(demei.blabel, demej.label);
      Q(i, index) += 0.5*(nodes[demei.blabel].m + mj);
    }
        
        // look at second lineage
    if (demej.llabel != -1){
      index = revLookup(demei.label, demej.llabel);
      Q(i, index) += 0.5*(nodes[demej.llabel].m + mi);
    }
    if (demej.rlabel != -1){
      index = revLookup(demei.label, demej.rlabel);
      Q(i, index) += 0.5*(nodes[demej.rlabel].m + mi);
    }
    if (demej.tlabel != -1){
      index = revLookup(demei.label, demej.tlabel);
      Q(i, index) += 0.5*(nodes[demej.tlabel].m + mi);
      
    }
    if (demej.blabel != -1){
      index = revLookup(demei.label, demej.blabel);
      Q(i, index) += 0.5*(nodes[demej.blabel].m + mi);
    }
    
    // if both in the same deme, they can coalesce
    if (demei.label == demej.label){
      // matrix indices start with 0
      Q(i,nstates-1) += demei.q;
    }

    Q(i,i) = -1*(Q.row(i).sum());
    
  }
    
  padm(Q, E);
}

void computeWeights(VectorXd &w, VectorXd &x, double u){
  // REQUIRES: w and x vectors of length 20 and u > 0.
  // MODIFIES: w and x
  // EFFECTS: x will contain the x-values telling you where to evaluate P(T_mrca = x); w will contains the weights
  // This function allows user to compute an integral by computing \sum_i P(T_mrca = x_i) * w_i 

    x << 0.07053988969199015, 0.3721268180016133, 0.9165821024832778, 1.707306531028346,
        2.749199255309431, 4.048925313850894, 5.615174970861625, 7.459017453671069,
        9.594392869581098, 12.03880254696431, 14.81429344263073, 17.94889552051938,
        21.47878824028502, 25.45170279318691, 29.9325546317006, 35.01343424047902,
        40.83305705672853, 47.61999404734653, 55.8107957500639, 66.52441652561578;
  
    w << 0.1687468018511152, 0.2912543620060685, 0.266686102867001, 0.1660024532695055,
        0.07482606466879245, 0.02496441730928314, 0.006202550844572207, 0.001144962386476896,
        0.0001557417730278113, 1.54014408652249e-05, 1.086486366517979e-06, 5.3301209095567e-08,
        1.757981179050555e-09, 3.725502402512292e-11, 4.767529251578155e-13, 3.372844243362382e-15,
        1.155014339500436e-17, 1.539522140582338e-20, 5.28644272556909e-24, 1.656456612498994e-28;
    
    // integral f(t) e^(-2*m*t) (1 + 2*m*t) dt = integral f(u/2m) e^-u (1+u) du/2m
    w = w.cwiseProduct(VectorXd::Ones(20)+x)/(2*u);
    x = x/(2*u);

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


    // array of nodes
    node nodes[ndemes];
    const int nreps = 1;
    for (int i = 0; i < nreps; i++){
        VectorXd mrates(ndemes);
        VectorXd qrates(ndemes);
        qrates << 1e-3, 1e-3, 1e-3, 1e-3;
        mrates << 0.464932, 0.464932, 0.464932, 0.464932;
        
        //VectorXd ones = VectorXd::Ones(ndemes);
        //qrates.setRandom(ndemes);
        //mrates.setRandom(ndemes);
        //qrates = 0.001 * (qrates + ones);
        //mrates = 0.001 * (mrates + ones);
        //qrates = (qrates + ones);
        //mrates = (mrates + ones);
        // initalize graph
        makeGraph(nodes, mrates, qrates, nrow, ncol);
    
        MatrixXd Q(nstates, m);
        MatrixXd H(m, m);
        krylovProj(nodes, H, Q, m);
        MatrixXd E(m,m);
        cout << H << endl;
        
        double t = 1;
        H = H*t;
        padm(H, E);
        VectorXd Y_1(nstates);
        VectorXd Y_2(nstates);
        MatrixXd E2(nstates, nstates);
        MatrixXd Q2(nstates, nstates);
        test_FullMatrixExp(nodes, Q2, E2);
        VectorXd l = VectorXd::Zero(nstates);
        l[nstates-1] = 1.0;
        VectorXd Y(nstates);
        Y_1 = (Q*E)*(Q.transpose()*l);
        Y_2 = E2*l;
        //cout << (Y_1-Y_2).norm() << endl;
        

        //VectorXd w(20);
        //VectorXd x(20);
        //computeWeights(w, x, 2);
        
        
         }
}
