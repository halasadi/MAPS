
#include "util.hpp"

extern string dist_metric;

Params::Params( ) { }
Params::~Params( ) { }
Params::Params(const string &params_file, const long seed_from_command_line) {
    try {
        po::options_description eems_options("EEMS options from parameter file");
        eems_options.add_options()
        ("seed", po::value<long>(&seed)->default_value(seed_from_command_line), "Random seed")
        ("datapath", po::value<string>(&datapath)->required(), "Path to coord/sims/outer files")
        ("mcmcpath", po::value<string>(&mcmcpath)->required(), "Path to output directory")
        ("prevpath", po::value<string>(&prevpath)->default_value(""), "Path to previous output directory")
        ("gridpath", po::value<string>(&gridpath)->default_value(""), "Path to demes/edges/ipmap files")
        ("nIndiv", po::value<int>(&nIndiv)->required(), "nIndiv")
        ("genomeSize", po::value<double>(&genomeSize)->default_value(3000), "genomeSize")
        ("lowerBound", po::value<double>(&lowerBound)->required(), "lowerBound")
        ("upperBound", po::value<double>(&upperBound)->default_value(numeric_limits<double>::infinity()), "upperBound")
        ("nDemes", po::value<int>(&nDemes)->default_value(1), "nDemes")
        ("distance", po::value<string>(&distance)->default_value("euclidean"), "distance")
        ("numMCMCIter", po::value<int>(&numMCMCIter)->default_value(1), "numMCMCIter")
        ("numBurnIter", po::value<int>(&numBurnIter)->default_value(0), "numBurnIter")
        ("numThinIter", po::value<int>(&numThinIter)->default_value(0), "numThinIter")
        ("mSeedsProposalS2", po::value<double>(&mSeedsProposalS2)->default_value(0.01), "mSeedsProposalS2")
        ("qSeedsProposalS2", po::value<double>(&qSeedsProposalS2)->default_value(0.1), "qSeedsProposalS2")
        ("mEffctProposalS2", po::value<double>(&mEffctProposalS2)->default_value(0.1), "mEffctProposalS2")
        ("qEffctProposalS2", po::value<double>(&qEffctProposalS2)->default_value(0.001), "qEffctProposalS2")
        ("mrateMuProposalS2", po::value<double>(&mrateMuProposalS2)->default_value(0.01), "mrateMuProposalS2")
        ("qrateMuProposalS2", po::value<double>(&qrateMuProposalS2)->default_value(0.01), "qrateMuProposalS2")
        ("omegaProposalS2", po::value<double>(&omegaProposalS2)->default_value(0.1), "omegaProposalS2")
        ("qVoronoiPr", po::value<double>(&qVoronoiPr)->default_value(0.5), "qVoronoiPr")
        ("mrateShape", po::value<double>(&mrateShape_2)->default_value(0.001), "mrateShape")
        ("qrateShape", po::value<double>(&qrateShape_2)->default_value(0.001), "qrateShape")
        ("qrateScale", po::value<double>(&qrateScale_2)->default_value(1.0), "qrateScale")
        ("mrateScale", po::value<double>(&mrateScale_2)->default_value(1.0), "mrateScale")
        ("mnegBiProb", po::value<double>(&mnegBiProb)->default_value(0.67), "mnegBiProb")
        ("mnegBiSize", po::value<int>(&mnegBiSize)->default_value(10), "mnegBiSize")
        ("qnegBiProb", po::value<double>(&qnegBiProb)->default_value(0.67), "qnegBiProb")
        ("olderpath", po::value<string>(&olderpath)->default_value(""), "Path to a run with a older time period")
        ("qnegBiSize", po::value<int>(&qnegBiSize)->default_value(10), "qnegBiSize");
        ifstream instrm(params_file.c_str());
        po::variables_map vm;
        po::store(po::parse_config_file(instrm,eems_options,true),vm);
        po::notify(vm);
        instrm.close();
    } catch(exception& e) {
        cerr << "[EEMS::Params] Error parsing input parameters in " << params_file << ": " << endl;
        cerr << e.what() << endl; exit(1);
    }
    mrateShape_2 /= 2.0;
    qrateShape_2 /= 2.0;
    mrateScale_2 /= 2.0;
    qrateScale_2 /= 2.0;
    
    testing = false;
    
    mrateMuUpperBound = 10; 
    qrateMuUpperBound = 10;
    mrateMuLowerBound = -10.0;
    qrateMuLowerBound = -10.0;
    
    mEffctHalfInterval = 20;
    qEffctHalfInterval = 20;
    
    min_omegaq = -6.9077;
    max_omegaq = 0;
    min_omegam = -6.9077;
    max_omegam = 0.6937;
}
ostream& operator<<(ostream& out, const Params& params) {
    out << "               datapath = " << params.datapath << endl
    << "               mcmcpath = " << params.mcmcpath << endl
    << "               prevpath = " << params.prevpath << endl
    << "               gridpath = " << params.gridpath << endl
    << "               distance = " << params.distance << endl
    << "                 nIndiv = " << params.nIndiv << endl
    << "             genomeSize = " << params.genomeSize << endl
    << "             lowerBound = " << params.lowerBound << endl
    << "             upperBound = " << params.upperBound << endl
    << "                 nDemes = " << params.nDemes << endl
    << "                   seed = " << params.seed << endl
    << "            numMCMCIter = " << params.numMCMCIter << endl
    << "            numBurnIter = " << params.numBurnIter << endl
    << "            numThinIter = " << params.numThinIter << endl
    << "              mnegBiSize = " << params.mnegBiSize << endl
    << "              mnegBiProb = " << params.mnegBiProb << endl
    << "              qnegBiSize = " << params.qnegBiSize << endl
    << "              qnegBiProb = " << params.qnegBiProb << endl
    << "             qVoronoiPr = " << params.qVoronoiPr << endl
    << "             mrateShape = " << 2.0*params.mrateShape_2 << endl
    << "             qrateShape = " << 2.0*params.qrateShape_2 << endl
    << "             qrateScale = " << 2.0*params.qrateScale_2 << endl
    << "             mrateScale = " << 2.0*params.mrateScale_2 << endl
    << "       mSeedsProposalS2 = " << params.mSeedsProposalS2 << endl
    << "       qSeedsProposalS2 = " << params.qSeedsProposalS2 << endl
    << "       mEffctProposalS2 = " << params.mEffctProposalS2 << endl
    << "       qEffctProposalS2 = " << params.qEffctProposalS2 << endl
    << "      mrateMuProposalS2 = " << params.mrateMuProposalS2 << endl
    << "      qrateMuProposalS2 = " << params.qrateMuProposalS2 << endl;
    return out;
}
bool Params::check_input_params( ) const {
    bool error = false;
    boost::filesystem::path mcmcdir(mcmcpath.c_str());
    boost::filesystem::path prevdir(prevpath.c_str());
    cerr << "Using Boost " << BOOST_LIB_VERSION
    << " and Eigen " << EIGEN_WORLD_VERSION << "." << EIGEN_MAJOR_VERSION << "." << EIGEN_MINOR_VERSION << endl
    << "  EEMS was tested with Boost 1_57 and Eigen 3.2.4" << endl << endl;
    if (!numeric_limits<double>::has_infinity || !numeric_limits<double>::infinity()) {
        cerr << "  Infinity not supported on this platform" << endl;
        error = true;
    }
    if (!boost::filesystem::create_directory(mcmcdir) && !exists(mcmcdir)) {
        cerr << "  Failed to create output directory " << mcmcpath << endl;
        error = true;
    }
    if (!(!distance.compare("euclidean") || !distance.compare("greatcirc"))) {
        cerr << "  Choose either 'euclidean' or 'greatcirc' distance metric" << endl;
        error = true;
    }
    if (!boost::filesystem::exists(datapath + ".coord") ||
        !boost::filesystem::exists(datapath + ".sims") ||
        !boost::filesystem::exists(datapath + ".outer")) {
        cerr << "  Failed to find input files " << datapath << ".coord/sims/outer" << endl;
        error = true;
    }
    if (!gridpath.empty() &&
        (!boost::filesystem::exists(gridpath + ".demes") ||
         !boost::filesystem::exists(gridpath + ".edges") ||
         !boost::filesystem::exists(gridpath + ".ipmap"))) {
            // Path to input population grid is specified
            // but one of the required files is missing
            cerr << "  Failed to find graph files " << gridpath << ".demes/edges/ipmap" << endl;
            error = true;
        }
    if (!prevpath.empty() && !exists(prevdir)) {
        cerr << "  Failed to find directory " << prevpath << " to previous EEMS output" << endl;
        error = true;
    }
    if (!(mSeedsProposalS2>0) || !(mEffctProposalS2>0) || !(mrateMuProposalS2>0) ||
        !(qSeedsProposalS2>0) || !(qEffctProposalS2>0) || !(qrateMuProposalS2>0)) {
        cerr << "  Choose positive variance parameters for the proposal distributions:" << endl
        << "  mrateMuProposalS2 = " << mrateMuProposalS2 << ", qrateMuProposalS2 = " << qrateMuProposalS2 << endl
        << "   mSeedsProposalS2 = " << mSeedsProposalS2 << ", mEffctProposalS2 = " << mEffctProposalS2 << endl
        << "   qSeedsProposalS2 = " << qSeedsProposalS2 << ", qEffctProposalS2 = " << qEffctProposalS2 << endl;
        error = true;
    }
    
    if (genomeSize > 3.3e3){
        cerr << "  Error with genome size: " << endl
        << " genomeSize = " << genomeSize << endl;
        error = true;
    }
    
    
    if (lowerBound < 0){
        cerr << "  Error with IBD cut off: " << endl
        << " lowerBound = " << lowerBound << endl;
        error = true;
    }
    
    
    if (!(numMCMCIter>0) || !(numBurnIter>=0) || !(numThinIter>=0) || !(numMCMCIter>numBurnIter) ||
        !(numMCMCIter>numBurnIter+numThinIter)) {
        cerr << "  Error with the MCMC parameters:" << endl
        << "  numMCMCIter = " << numMCMCIter << ", numBurnIter = " << numBurnIter << ", numThinIter " << numThinIter << endl;
        error = true;
    }
    if (!(qrateShape_2>0) || !(mrateShape_2>0) ||
        !(qrateScale_2>0) || !(mrateScale_2>0)) {
        cerr << "  Error with the Inverse Gamma hyperparameters:" << endl
        << "  qrateShape = " << 2.0*qrateShape_2 << ", qrateScale = " << 2.0*qrateScale_2 << endl
        << "  mrateShape = " << 2.0*mrateShape_2 << ", mrateScale = " << 2.0*mrateScale_2 << endl;
        error = true;
    }
    if (!(mnegBiSize>0) || !( (mnegBiProb>0) && (mnegBiProb<1) )) {
        cerr << "  Error with the m Negative Binomial hyperparameters:" << endl
        << "  mnegBiSize = " << mnegBiSize << ", mnegBiProb = " << mnegBiProb << endl;
        error = true;
    }
    if (!(qnegBiSize>0) || !( (qnegBiProb>0) && (qnegBiProb<1) )) {
        cerr << "  Error with the q Negative Binomial hyperparameters:" << endl
        << "  qnegBiSize = " << qnegBiSize << ", qnegBiProb = " << qnegBiProb << endl;
        error = true;
    }
    return(error);
}
VectorXd split(const string &line) {
    istringstream in(line);
    vector<double> numbers;
    double number;
    while (!in.eof()) { in >> number;
        if (!in.fail()) { numbers.push_back(number); } else { break; }
    }
    if (in.fail()||in.bad()) { return (VectorXd::Zero(0)); }
    // Since we don't know how many numbers there are,
    // first we store the data in a std::vector of doubles,
    // and then typecast it to Eigen::VectorXd
    return (VectorXd::Map(&numbers[0],numbers.size()));
}
bool isposdef(const MatrixXd &A) {
    SelfAdjointEigenSolver<MatrixXd> eig(A,EigenvaluesOnly);
    double minval = eig.eigenvalues().minCoeff();
    return (minval>0);
}

double logdet(const MatrixXd &A) {
    return (A.selfadjointView<Lower>().ldlt().vectorD().array().log().sum());
}
double pseudologdet(const MatrixXd &A, const int rank) {
    SelfAdjointEigenSolver<MatrixXd> x(A);
    return (x.eigenvalues().reverse().array().head(rank).log().sum());
}


double get_bootstrap_var(const MatrixXi &Sims, VectorXd cvec, const VectorXi &indiv2deme, int nb, int alpha, int beta){
    
    int n_alpha = cvec(alpha);
    int n_beta = cvec(beta);
    
    VectorXi deme1_indices = VectorXi::Zero(n_alpha);
    VectorXi deme2_indices = VectorXi::Zero(n_beta);
    
    int cnt_alpha = 0;
    int cnt_beta = 0;

    
    for (int i = 0; i < indiv2deme.size(); i++){
        
        if (indiv2deme(i) == alpha){
            deme1_indices(cnt_alpha) = i;
            cnt_alpha += 1;
        }
        if (indiv2deme(i) == beta){
            deme2_indices(cnt_beta) = i;
            cnt_beta += 1;
        }
        
    }
    
    
    VectorXi deme1_subsamples = VectorXi::Zero(n_alpha);
    VectorXi deme2_subsamples = VectorXi::Zero(n_beta);
    double running_sum;
    double n;
    
    VectorXd myMeans = VectorXd::Zero(nb);
    
    for (int b = 0; b < nb; b++){
        
        for (int i = 0; i < n_alpha; i++){
            deme1_subsamples(i) = deme1_indices(rand() % n_alpha);
        }
        for (int i = 0; i < n_beta; i++){
            deme2_subsamples(i) = deme2_indices(rand() % n_beta);
        }
        
        running_sum = 0;
        n = 0;
        
        if (alpha == beta){
            for (int i = 0; i < (n_alpha-1); i++){
                for (int j = (i+1); j < n_alpha; j++){
                    running_sum += Sims(deme1_subsamples(i), deme1_subsamples(j));
                    n += 1;
                }
            }
        } else {
            for (int i = 0; i < n_alpha; i++){
                for (int j = 0; j < n_beta; j++){
                    running_sum += Sims(deme1_subsamples(i), deme2_subsamples(j));
                    n += 1;
                }
            }
        }
        
        myMeans(b) = running_sum / n;
        
    }
    
    
    double avg = myMeans.sum() / nb;
    double variance = 0;
    for (int i = 0; i < nb; i++){
        variance = variance + (myMeans(i) - avg) * (myMeans(i) - avg);
    }
    variance = variance / (nb - 1);
    
    return(variance);

}

double poisln(const MatrixXd &expectedIBD, const MatrixXd &observedIBDCnt, const MatrixXd &ceffective, const MatrixXd &cMatrix){
    double ll = 0;
    int o = expectedIBD.rows();
    double lamda;
    double xbar;
    for (int alpha = 0; alpha < o; alpha++){
        for (int beta = alpha; beta < o; beta++){
            if (expectedIBD(alpha,beta) < 1e-8){
                lamda = 1e-8;
            } else {
                lamda = expectedIBD(alpha,beta);
            }
            xbar = observedIBDCnt(alpha,beta) / cMatrix(alpha,beta);
            
            ll += ceffective(alpha,beta) * (xbar * log(lamda) - lamda);
        }
    }
    return(ll);
}

double mvgammaln(const double a, const int p) {
    double val = 0.25*log_pi*p*(p-1);
    for (int i = 0 ; i < p ; i++ ) {
        val += lgamma(a-0.5*i);
    }
    return (val);
}
/*
 Squared Euclidean distance: Taking the square root is not necessary
 because EEMS uses the distances to find closest points. For example,
 pairwise_distance(X,Y).col(0).minCoeff( &closest )
 finds the row/point in X that is closest to the first row/point in Y
 */
MatrixXd euclidean_dist(const MatrixXd &X, const MatrixXd &Y) {
    return (  X.rowwise().squaredNorm().eval().replicate(1,Y.rows())
            + Y.rowwise().squaredNorm().eval().transpose().replicate(X.rows(),1)
            - 2.0*X*Y.transpose() );
}
// Great circle distance, up to a constant of proportionality equal to 2*R
// where R is the earth's radius
MatrixXd greatcirc_dist(const MatrixXd &X, const MatrixXd &Y) {
    int nr = X.rows();
    int nc = Y.rows();
    MatrixXd lon1 = X.col(0).replicate(1,nc) * pi_180;
    MatrixXd lat1 = X.col(1).replicate(1,nc) * pi_180;
    MatrixXd lon2 = Y.col(0).transpose().replicate(nr,1) * pi_180;
    MatrixXd lat2 = Y.col(1).transpose().replicate(nr,1) * pi_180;
    MatrixXd dlon = 0.5*(lon2 - lon1);
    MatrixXd dlat = 0.5*(lat2 - lat1);
    MatrixXd a = dlat.array().sin().square().matrix() +
    (dlon.array().sin().square() * lat1.array().cos() * lat2.array().cos()).matrix();
    MatrixXd c = (a.array()<1.0).select(a.array().sqrt(),MatrixXd::Ones(nr,nc)).array().asin();
    return (c); // Instead of (2*R*c) where R = 6378137 is the Earth's radius.
}
// Compute pairwise distances between the rows of X and the rows of Y.
// Choose either Euclidean or great circle distance
MatrixXd pairwise_distance(const MatrixXd &X, const MatrixXd &Y) {
    if (!dist_metric.compare("greatcirc")) {
        return (greatcirc_dist(X,Y));
    } else {
        return (euclidean_dist(X,Y));
    }
}
// Read a matrix, with unknown dimensions, from a text file
// Return an empty matrix (0 rows, 0 columns) if there is an error
MatrixXd readMatrixXd(const string &filename) {
    ifstream instrm(filename.c_str(), ios::in);
    if (!instrm.is_open( )) { return (MatrixXd::Zero(0,0)); }
    string line; getline(instrm,line);
    boost::algorithm::trim(line);
    // Split the first line into numbers
    VectorXd row = split(line);
    // This tells us the number of columns
    int cols = row.size();
    if (!cols) { return (MatrixXd::Zero(0,0)); }
    MatrixXd mat(1,cols); mat.row(0) = row;
    while (getline(instrm,line)) {
        boost::algorithm::trim(line);
        row = split(line);
        if (row.size()==cols) {
            // Resize the matrix to add each row (done once at the start of EEMS, so okay)
            insertRow(mat,row);
        } else {
            return (MatrixXd::Zero(0,0));
        }
    }
    return(mat);
}
double trace_AxB(const MatrixXd &A, const MatrixXd &B) {
    return (A.cwiseProduct(B).sum());
}
// If there are two draws and the first has two tiles and the second -- three tiles,
// then sizes = c(2,3) and array = c(m_{1t_1},m_{1t_2},m_{2t_1},m_{2t_2},m_{2t_3})
bool dlmcell(const string &filename, const VectorXd &sizes, const vector<double> &array) {
    bool error = false;
    if (array.size()!=sizes.sum()) { error = true; return(error); }
    ofstream out(filename.c_str(),ofstream::out);
    if (!out.is_open()) { error = true; return(error); }
    vector<double>::const_iterator it = array.begin();
    for ( int i = 0 ; i < sizes.size() ; i++ ) {
        for ( int j = 0 ; j < sizes(i) ; j++ ) {
            out << fixed << setprecision(14) << *it << " "; it++;
        }
        out << endl;
    }
    out.close( );
    return(error);
}
void removeElem(VectorXd &vec, const int elemToRemove)
{
    int elems = vec.size() - 1;
    if (elemToRemove <= elems) {
        // Copy the last element into the element to remove
        if (elemToRemove < elems) {
            vec(elemToRemove) = vec(elems);
        }
        vec.conservativeResize(elems);
    }
}
void removeRow(MatrixXd &mat, const int rowToRemove)
{
    int rows = mat.rows() - 1;
    int cols = mat.cols();
    if( rowToRemove <= rows ) {
        // Copy the last row into the row to remove
        if ( rowToRemove < rows ) {
            mat.row(rowToRemove) = mat.row(rows);
        }
        mat.conservativeResize(rows,cols);
    }
}
void insertElem(VectorXd &vec, const double &elem)
{
    int elems = vec.size();
    vec.noalias() = (VectorXd(elems+1) << vec, elem).finished();
}
void insertRow(MatrixXd &mat, const VectorXd &row)
{
    int rows = mat.rows();
    int cols = mat.cols();
    mat.noalias() = (MatrixXd(rows+1,cols) << mat, row.transpose()).finished();
}
/*
 Log probability mass function of the negative binomial distribution, with
 parameters size (number of failures) and prob (probability of success)
 * up to a constant of proportionality that depends on size and prob,
 but not on the random variable k (number of successes)
 Also -- log pmf of the zero-truncated negative binomial, with the same parameters,
 as long as k>0, since truncating only changes the constant of proportionality
 Uses gamma(n) = (n-1)!
 */
double dnegbinln(const int k, const int size, const double prob) {
    double pln = -Inf;
    if ( (k>=0) && (size>0) && (prob>0) && (prob<1) ) {
        pln = lgamma(size+k) - lgamma(k+1) + k*log(prob);
    }
    return (pln);
}
/*
 Log probability density function of the inverse gamma distribution,
 with parameters shape and scale
 * Up to a constant of proportionality that depends on the shape and the scale,
 but not on the random variable x
 */
double dinvgamln(const double x, const double shape, const double scale) {
    double pln = -Inf;
    if ( (x>0.0) && (shape>0.0) && (scale>0.0) ) {
        pln = - (shape+1.0)*log(x) - scale/x;
    }
    return (pln);
}
/*
 Multivariate normal probability density function, on the log scale,
 up to a constant of proportionality
 */
double dmvnormln(const VectorXd &x, const VectorXd &mu, const MatrixXd &sigma) {
    double pln = -Inf;
    if ( isposdef(sigma) ) {
        LLT<MatrixXd> Chol(sigma);
        MatrixXd L = Chol.matrixL();
        pln = - L.diagonal().array().log().sum();
        pln -= 0.5 * (x-mu).transpose() * Chol.solve(x-mu);
    }
    return (pln);
}
/*
 Truncated normal probability density function, on the log scale,
 with support [0,bnd], including the normalizing constant
 */
double dtrnormln(const double x, const double mu, const double sigma2, const double bnd) {
    double pln = -Inf;
    if ( (sigma2>0) && (x>=-bnd) && (x<=bnd) ) {
        boost::math::normal pnorm(mu,sqrt(sigma2));
        pln = - 0.5 * log(sigma2) - 0.5 * (x-mu) * (x-mu) / sigma2
        - log(cdf(pnorm,bnd) - cdf(pnorm,-bnd));
    }
    return (pln);
}
VectorXd slice(const VectorXd &A, const VectorXi &I) {
    int elems = I.size();
    VectorXd B(elems);
    assert(I.minCoeff() >= 0);
    assert(I.maxCoeff() < A.size());
    for ( int i = 0 ; i < elems ; i++ ) {
        B(i) = A(I(i));
    }
    return (B);
}
MatrixXd slice(const MatrixXd &A, const VectorXi &R, const VectorXi &C) {
    int rows = R.size();
    int cols = C.size();
    MatrixXd B(rows,cols);
    assert(R.minCoeff() >= 0);
    assert(C.minCoeff() >= 0);
    assert(R.maxCoeff() < A.rows());
    assert(C.maxCoeff() < A.cols());
    for ( int i = 0 ; i < rows ; i++ ) {
        for ( int j = 0 ; j< cols ; j++ ) {
            B(i,j) = A(R(i),C(j));
        } }
    return (B);
}

double max(double a, double b){
    if (a > b){
        return(a);
    }
    return(b);
}

void getWeights(VectorXd &w, VectorXd &x){
    // REQUIRES: w and x vectors of length 50
    // MODIFIES: w and x
    // EFFECTS: x will contain the x-values telling you where to evaluate P(T_mrca = x); w will contains the weights
    // This function allows user to compute an integral by computing \sum_i P(T_mrca = x_i) * w_i
    
    x << 0.07197890982430907685, 0.2413621356214323113832, 0.50777161206496736682,
    0.87144100848215091489, 1.3327153593876555612, 1.89203857575589673578,
    2.54995389165696269159, 3.30710638809426104509, 4.1642464309382815239,
    5.1222338489683650003, 6.1820428555624776171, 7.3447677682201374241,
    8.611629605039126607, 9.9839836572644003588, 11.4633281577394397878,
    13.0513141887293622298, 14.7497570005632491357, 16.5606489462107181103,
    18.4861742778364391121, 20.5287261015344429523, 22.6909258483688038223,
    24.97564569685792480173, 27.3860344785262187279, 29.92554771997283930857,
    32.59798262998117745422, 35.40751903929353023831, 38.35876755865307164698,
    41.4568265582708885324, 44.7073500182295270902, 48.11662889629256334459,
    51.6916894678465181416, 55.4404132017820128574, 59.37168428037916420825,
    63.49557305617946339417, 67.82356688527364443624, 72.36886439711898866969,
    77.14675619634162666911, 82.17512565948067597586, 87.4751203582245020515,
    93.0720721704444912416, 98.99679073227190311995, 105.2874371482835905836,
    111.9923375735336617652, 119.1743972669017647382, 126.91841438735756732,
    135.344080011584273835, 144.6313615494928725149, 155.0771275144866916574,
    167.2505316308244871017, 182.620207348251479189;
    
    w << 0.008098150669659729617, 0.04130873125538665997, 0.09625940978218537466,
    0.1503491696588344311, 0.17934986299937562831, 0.173581394953721460356,
    0.140856655908618615123, 0.097739526671480248591, 0.0587261120368546837416,
    0.030808796814653968376, 0.014192143542695654478, 0.0057625101944473758948,
    0.0020676112211018121615, 6.566117266315217618E-4, 1.84713846428817139761E-4,
    4.6041662947613873147E-5, 1.016612790318865366769E-5, 1.9870817337963512747E-6,
    3.4344907131572547167E-7, 5.2416191882540332129E-8, 7.050822455756935318E-9,
    8.3415284815468702479E-10, 8.6573754732312158465E-11, 7.8596173841489468557E-12,
    6.2209642223163649024E-13, 4.2769613995502895646E-14, 2.5433779837237564877E-15,
    1.302075024421022199E-16, 5.7083496285476396191E-18, 2.13034598442665506488E-19,
    6.72273042229833935E-21, 1.7803851060592980154E-22, 3.9231584139344725937E-24,
    7.1232728017092643659E-26, 1.05390989562966620855E-27, 1.25438816254046964E-29,
    1.1832996034202271289E-31, 8.693916313933935085E-34, 4.8733576082223977928E-36,
    2.0332444973728568667E-38, 6.12678969991989233E-41, 1.28463695773429860451E-43,
    1.787969839412486007E-46, 1.55365992711589439292E-49, 7.761542881245283525E-53,
    1.984698518925299219E-56, 2.18270822517498151914E-60, 7.7576333601861023782E-65,
    5.1724748561078432042E-70, 1.6224693284923917835E-76;
    
}
