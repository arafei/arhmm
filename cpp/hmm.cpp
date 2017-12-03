
// Final project

#include <iostream>
#include <vector>
#include <cstdlib>
#include <iomanip>
#include <cmath>
#include <limits>
#include <climits>
#include <algorithm>
#include <iterator>
#include <fstream>
//#include <boost/tokenizer.hpp>
//#include <boost/lexical_cast.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/exponential.hpp>
#include <boost/math/distributions/poisson.hpp>
#include <boost/math/distributions/lognormal.hpp>
#include </afs/umich.edu/user/a/r/arafei/Private/biostat682/test/Eigen/Cholesky>
#include </afs/umich.edu/user/a/r/arafei/Private/biostat682/test/Eigen/QR>
#include </afs/umich.edu/user/a/r/arafei/Private/biostat682/test/Eigen/Dense>
#include </afs/umich.edu/user/a/r/arafei/Private/biostat682/test/Eigen/Core>
#include </afs/umich.edu/user/a/r/arafei/Private/biostat682/test/Eigen/SVD>
#include "Matrix6.h"
using namespace std;
#include "HMM6.h"


/////////////////////////////////////////// Read from file //////////////


/////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////// main function ///////////////

int main(int argc, char* argv[]){


	// read data from CSV file
	vector<std::vector<double> > values;
	ifstream fin(argv[1]);
	for(string line; getline(fin, line);){
		replace(line.begin(), line.end(), ',', ' ');
		istringstream in(line);
		values.push_back(
		vector<double>(std::istream_iterator<double>(in),
			std::istream_iterator<double>()));
	}

	int m, tok;
	int n=(int)values.size();
	int nX=(int)values[0].size()-1;

cout << "n= " << n << endl;
cout << "n= " << nX << endl;

	
	double tok1;
	vector<int> dist;
	cout << "Insert the number of hidden states:" << endl;
	cin >> m;
	if (m<2){
		cerr << "The number of states should be >2. Please try again." << endl;
		return -1;
	}
	cout << "Insert the distribution type for each state:" << endl;
	cout << "1=Gaussian, 2=Exponential, 3=Poisson, and 4=LogNormal" << endl;
	for(int i=0; i<m; i++){
		cin >> tok;
		if(tok!=1 && tok!=2 && tok!=3 && tok!=4){
			cerr << "Wrong input. Please try again." << endl;
			return -1;
		}
		else
		dist.push_back(tok);
	}
	    
	cout << "Insert the vector of trends:" << endl << "0=No trend, 1=Trend" << endl;
	vector<int> trend;
	for(int i=0; i<m; i++){
		cin >> tok;
		trend.push_back(tok);
    	}

    	int pp=0;
	int maxparam=0;
	vector<int> p(m,4);
	for(int i=0; i<m; i++){
		if(dist[i]==1 || dist[i]==4){
			if (trend[i]==1) p[i]=5+nX;
			else p[i]=2;
		}
	else {
		if (trend[i]==1) p[i]=4+nX;
		else p[i]=1;}
	pp = max(pp,p[i]);
	maxparam=maxparam+p[i];
	}    

	HMM6 hmm(m, n, nX, pp, maxparam);
	
	for(int i=0; i<n; i++){
		hmm.outs.push_back(values[i][0]);
		for(int j=0; j<nX; j++)
			hmm.X.data[i][j]=values[i][j+1];
	}

	hmm.dist=dist;
	hmm.trend=trend;

	for(int i=0; i<m; i++){
        cout << "Insert the initial values of parameters for the state " << i+1 << ": " << endl;
        cout << "Example: B0, B1, B2, B3, Sigma for Normal Dist" << endl;
        for(int j=0; j<p[i]; j++){
            cin >> tok1;
            hmm.theta.data[j][i] = tok1;
        }
    }

	//insert values of initial probabilities 
	cout << "insert initial probabilities: " << endl;
	double pis=0;
	for(int i=0; i<m; i++){
		cin >> hmm.pis[i];
		if(hmm.pis[i]>1 || hmm.pis[i]<0){
			cerr << "Initial probabilities shoud be between 0 and 1." << endl;
			return -1;
		}
	pis=pis+hmm.pis[i];
	}
	if(pis!=1){
		cerr << "Sum of the initial probabilities shoud be 1." << endl;
		return -1;
	}

	//insert values of transition probabilities 	
	for (int i=0; i<m; i++){
		pis=0;
		cout << "insert transition probabilities from state: " << i+1 << endl;
		for (int j=0; j<m; j++){
			cin >> hmm.trans.data[i][j];
			if(hmm.trans.data[i][j]>1 || hmm.trans.data[i][j]<0){
				cerr << "Transition probabilities shoud be between 0 and 1." << endl;
				return -1;
			}
			pis=pis+hmm.trans.data[i][j];
		}
		if(pis!=1){
			cerr << "Sum of the transition probabilities in each state shoud be 1." << endl;
			return -1;
		}
	}

	//Run HMM model
	hmm.covariate();
	hmm.trainHMM();
	hmm.viterbi();
	hmm.predict();
	
	// Prints out outputs

	// Prints out the maximum likelihood estimation of parameters
	cout << "Estimation of parameters:" << endl;
	for(int i=0; i<m; i++){
		for(int j=0; j<p[i]; j++) 
			cout <<  setprecision(4) << setw(10) << hmm.theta.data[j][i];
		cout << endl << endl;
	}

	// Prints out estimation of initial probabilities
	cout << "Estimation of initial probabilities:" << endl;
		for(int j=0; j<m; j++) 
			cout << setprecision(4) << setw(10) << hmm.pis[j];
	cout << endl << endl;

	// Prints out estimation of transition probabilities
	cout << "Estimation of transition probabilities:" << endl;
	for(int i=0; i<m; i++){
		for(int j=0; j<m; j++) 
			cout << setprecision(4) << setw(10) << hmm.trans.data[i][j];
	cout << endl << endl;
	}

	// Prints out the most likely hidden sequence
	cout << "The most likely state sequence:" << endl;
	for(int i=0; i<n; i++) 
		cout<< hmm.path[i] << " ";
	cout << endl << endl;
	
	// Prints out the fitted values of the model
	cout << "Fitted values of the model:" << endl;
	for(int i=0; i<n; i++) 
		cout<< hmm.fitted[i] << " ";
	cout << endl << endl;
	
	// Prints out Goodness of fit measures
	cout << setw(15) << "#Parameters" << setw(15) << "#Iterations" << setw(15) << "Log-likelihood" << setw(15) << "R-square" << setw(15) << "adj.R-Square" << setw(15) << "AIC" << setw(15) << "BIC" << endl;
	cout << setprecision(4) << setw(15) << hmm.gof[0] << setw(15) << hmm.gof[1] << setw(15) << hmm.gof[2] << setw(15) << hmm.gof[3] << setw(15) << hmm.gof[4] << setw(15) << hmm.gof[5] << setw(15) << hmm.gof[6] << endl << endl;

	return 0;
}