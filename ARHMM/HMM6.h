#ifndef __HMM_6_H
#define __HMM_6_H
#include "Matrix6.h"

using namespace Eigen;
using namespace boost::math;

const int r=12;
const int maxit=10000;
const double convthrd=1e-5;

class HMM6{
public:

bool tr;
int nStates, nTimes, nX, nCov, nParams, numit, maxParams;
double logLik;
vector<int> dist, trend, path;
vector<double> outs, pis, fitted, gof;
Matrix6<double> theta, X, trans, emis, alphas, betas, cov, deltas;

HMM6(int states, int times, int nX, int params, int maxp) : nStates(states), nTimes(times), nCov(nX), nParams(params), maxParams(maxp),
theta(params, states, 0), X(times, nX, 0), trans(states, states, 0), emis(states, times, 0),
alphas(states, times, 0), betas(states, times, 0), cov(times, nX+4, 0),
deltas(states, times, 0){
	pis.resize(nStates);
	path.resize(nTimes);
	fitted.resize(nTimes);
	path.resize(nTimes);
	gof.resize(7);
}


void covariate();
void density();
void forwardBackward();
void trainHMM();
void viterbi();
void predict();
};
#endif

void HMM6::covariate(){
	bool tr=false;
	for(int i=0; i<nStates; i++){
		if(trend[i]==1){
			tr=true;
			break;
		}
	}
	if(tr){
		for(int t=0; t<nTimes; t++){
			cov.data[t][0]=1.0;
			cov.data[t][1]=t+1.0;
			cov.data[t][2]=cos(2.0*M_PI*(double)((t+1.0)/(double)r));
			cov.data[t][3]=sin(2.0*M_PI*(double)((t+1.0)/(double)r));
			for(int j=0; j<nCov; j++){
				cov.data[t][j+4]=X.data[t][j];			
			}
		}
	}
}

void HMM6::density(){
	double mu;
	for(int i=0; i<nStates; i++){
	// Nomal distribution
		if(dist[i]==1){
			if(trend[i]==1){
				for(int t=0; t<nTimes; t++){
					mu=0;
					for(int j=0; j<nCov+4; j++)
						mu = mu + theta.data[j][i] * cov.data[t][j];
					emis.data[i][t]=pdf(normal_distribution<>(mu, sqrt(theta.data[nCov+4][i])), outs[t]);
				}
			}else{
				for(int t=0; t<nTimes; t++)
					emis.data[i][t]=pdf(normal_distribution<>(theta.data[0][i], theta.data[1][i]), outs[t]);
			}
		}
	// Exponential distribution
		else if(dist[i]==2){
			if(trend[i]==1){
				for(int t=0; t<nTimes; t++){
					mu=0;
					for(int j=0; j<nCov+4; j++)
						mu = mu + theta.data[j][i] * cov.data[t][j];
					emis.data[i][t]=pdf(exponential_distribution<>(mu), outs[t]);
				}
			}else{
				for(int t=0; t<nTimes; t++)
					emis.data[i][t]=pdf(exponential_distribution<>(theta.data[0][i]), outs[t]);
			}
		}
	// Poisson distribution
		else if(dist[i]==3){
			if(trend[i]==1){
				for(int t=0; t<nTimes; t++){
					mu=0;
					for(int j=0; j<nCov+4; j++)
						mu = mu + theta.data[j][i] * cov.data[t][j];
					emis.data[i][t]=pdf(poisson_distribution<>(mu), outs[t]);
				}
			}else{
				for(int t=0; t<nTimes; t++)
					emis.data[i][t]=pdf(poisson_distribution<>(theta.data[0][i]), outs[t]);
			}
		}
	// Log-normal distribution
		else if(dist[i]==4){
			if(trend[i]==1){
				for(int t=0; t<nTimes; t++){
					mu=0;
					for(int j=0; j<nCov+4; j++)
						mu = mu + theta.data[j][i] * cov.data[t][j];
					emis.data[i][t]=pdf(lognormal_distribution<>(mu, sqrt(theta.data[nCov+4][i])), outs[t]);
				}
			}else{
				for(int t=0; t<nTimes; t++)
					emis.data[i][t]=pdf(normal_distribution<>(theta.data[0][i], theta.data[1][i]), outs[t]);
			}
		}
	}
}

void HMM6::forwardBackward(){

    for(int i=0; i<nStates; i++){
        for (int j=0; j<nTimes; j++){
		alphas.data[i][j]=0; betas.data[i][j]=0;
	}
    }
    int scale =1;
    vector<double> rr(nTimes,0);
    double ss =0;
    for(int j=0; j<nStates; j++){ 
		alphas.data[j][0]=pis[j]*emis.data[j][0]; 
		ss=ss+alphas.data[j][0];
	}
    if(scale==1){
        rr[0]=ceil(log10(ss));
        for(int j=0; j<nStates; j++) 
			alphas.data[j][0]=alphas.data[j][0]/pow(10,rr[0]);
    }
    for (int j=0; j<(nTimes-1); j++){
        ss=0;
        for(int i=0; i<nStates; i++){
            for (int h=0; h<nStates; h++){
                alphas.data[i][j+1]=alphas.data[i][j+1] + alphas.data[h][j]*trans.data[h][i];
            }
            alphas.data[i][j+1]=alphas.data[i][j+1]*emis.data[i][j+1];
            ss=ss+alphas.data[i][j+1];
        }
        rr[j+1]=ceil(log10(ss));
        if(scale==1){
            for(int ii=0; ii<nStates; ii++) 
				alphas.data[ii][j+1]=alphas.data[ii][j+1]/pow(10,rr[j+1]);
        }
    }
    double sf=0;
    for(int j=0; j<nTimes; j++) sf=sf+rr[j];
    ss=0;
    for(int j=0; j<nStates; j++) 
		ss = ss+alphas.data[j][nTimes-1];
		logLik = sf*log(10)+log(ss);
    for(int j=0; j<nStates; j++) 
		betas.data[j][nTimes-1]=1;
    ss=0;
    for (int j=0; j<nStates; j++) 
		ss=ss+betas.data[j][nTimes-1];
    rr[nTimes-1]=ceil(log10(nStates));
    
    for (int j=0; j<nStates; j++) 
		betas.data[j][nTimes-1]=betas.data[j][nTimes-1]/pow(10, rr[nTimes-1]);
    
    for (int j=nTimes-2; j>=0; j--){
        ss=0;
        for(int h=0; h<nStates; h++){
            for (int i=0; i<nStates; i++){
                betas.data[h][j]=betas.data[h][j] + trans.data[h][i]*emis.data[i][j+1]*betas.data[i][j+1];
            }
            ss=ss+betas.data[h][j];
        }
        rr[j]=ceil(log10(ss));
        for(int ii=0; ii<nStates; ii++) 
			betas.data[ii][j]=betas.data[ii][j]/pow(10,rr[j]);
    }
}

void HMM6::trainHMM(){
	//Matrix6<double> initheta=theta;
	//int maxparams=theta.rowNums();
	vector<vector<vector<double> > > tau3(nStates, vector<vector<double> >(nStates,vector<double>(nTimes-1,0)));
	Matrix6<double> tau2(nStates, nTimes, 0);
	numit=0;
	bool converged=false;
	Matrix6<double> tauslice(nStates, nStates);
	double oldL =LLONG_MIN;
	vector<double> oldpis;
	Matrix6<double> oldtheta(nParams, nStates), oldtrans(nStates, nStates);
	double normf, sum;
	while((numit < maxit) && !converged){
		//keep the ininitial values of parameters
		oldtheta=theta;
		oldpis=pis;
		oldtrans=trans;
		//E step///////////////////////////////////////////////////////

		//precalculate conditional likelihood values
		density();

		//run forward-backward to compute a and b
		forwardBackward();

		//logLik.push_back(newL);

		//calculate tau_{hij}
		for(int j=0; j<(nTimes-1); j++){
			for(int h=0; h<nStates; h++){
				for(int i=0; i<nStates; i++){
					tau3[h][i][j]=alphas.data[h][j]*trans.data[h][i]*emis.data[i][j+1]*betas.data[i][j+1];
				}
			}
			normf=0;
			for(int h=0; h<nStates; h++){
				for(int i=0; i<nStates; i++){
					tauslice.data[h][i]=tau3[h][i][j];
					normf=normf+tauslice.data[h][i];
				}
			}
			for(int h=0; h<nStates; h++){
				for(int i=0; i<nStates; i++){
					tau3[h][i][j]=tau3[h][i][j]/normf;
				}
			}
		}
		//calculate tau_{i1}
		normf=0;
		for(int i=0; i<nStates; i++){
			tau2.data[i][0]=pis[i]*emis.data[i][0];
			normf=normf+tau2.data[i][0];
		}
		if(normf != 0){
			for(int i=0; i<nStates; i++)
				tau2.data[i][0]=tau2.data[i][0]/normf;
		}else{
			cerr << "normalization factor is zero" << endl;
			exit(EXIT_FAILURE);
		}		

		//calculate tau_{ij}
		for(int i=0; i<nStates; i++){
			for(int j=1; j<nTimes; j++){
				sum=0;
				for(int h=0; h<nStates; h++)
					sum=sum+tau3[h][i][j-1];
				tau2.data[i][j]=sum;
			}
		}

		////////// M Step /////////////////////////////////////

		// update priors
		for(int i=0; i<nStates; i++)
			pis[i]=tau2.data[i][0];

		// update transition probabilities
		for(int h=0; h<nStates; h++){
			for(int i=0; i<nStates; i++){
				trans.data[h][i]=0;
				sum=0;
				for(int j=0; j<(nTimes-1); j++){
					trans.data[h][i]=trans.data[h][i]+tau3[h][i][j];
					sum=sum+tau2.data[h][j];
				}
				trans.data[h][i]=trans.data[h][i]/sum;
			}
		}
		// normalize transition probabilities
		for(int h=0; h<nStates; h++){
			sum=0;
			for(int j=0; j<nStates; j++)
				sum=sum+trans.data[h][j];
			for(int i=0; i<nStates; i++)
				trans.data[h][i]=trans.data[h][i]/sum;
		}
		//update of conditional distribution parameters

		//Matrix6<double> oldtheta=theta;
		Matrix6<double> mu(nStates, nTimes, 0);
		
		//update seasonality and trend parameters of all states
		//for(int i=0; i<nStates; i++)
			//theta[i].clear();

		for(int i=0; i<nStates; i++){
		vector<double> param;
			if(dist[i]==1){
				if(trend[i]==1){
					MatrixXd M = MatrixXd::Zero(nCov+4, nCov+4);
					VectorXd V = VectorXd::Zero(nCov+4);
					VectorXd W(nCov+4);
					for (int j = 0; j<nTimes; j++) {
						for(int k=0; k<(nCov+4); k++)
							W(k) = cov.data[j][k];
						M = M + tau2.data[i][j] * (W*W.transpose());
						V = V + tau2.data[i][j] * outs[j] * W;
					}
					VectorXd tmpParam = M.colPivHouseholderQr().solve(V);
					for(int k=0; k<(nCov+4); k++)
						param.push_back( tmpParam(k) );

					// update variances
					// construct means from trend/seasonality parameters
					for (int j = 0; j<nTimes; j++) {
						for(int k=0; k<(nCov+4); k++)
							W(k) = cov.data[j][k];
						mu.data[i][j] = W.transpose()*tmpParam;
					}					
					double sigma=0.0;
					for(int j=0; j<nTimes; j++)
						sigma=sigma+tau2.data[i][j]*(outs[j]-mu.data[i][j])*(outs[j]-mu.data[i][j]);
					sum=0;
					for(int j=0; j<nTimes; j++)
						sum=sum+tau2.data[i][j];
					sigma=sigma/sum;
					param.push_back(sigma);
				}
				else{
					//no trend/seasonality modeling
					double mu=0;
					double sigma=0;
					for(int j=0; j<nTimes; j++){
						mu=mu+tau2.data[i][j]*outs[j];
						sigma=sigma+tau2.data[i][j]*outs[j]*outs[j];
					}
					normf=0;
					for(int j=0; j<nTimes; j++)
						normf=normf+tau2.data[i][j];
					mu=mu/normf;
					sigma=sigma/normf;
					param.push_back(mu);
					param.push_back(sigma-mu*mu);
				}			
			}
			
			// Exponential density
			if(dist[i]==2){
				// with trend/seasonality modeling
				if(trend[i]==1){
					MatrixXd M = MatrixXd::Zero(nCov+4, nCov+4);
					VectorXd V = VectorXd::Zero(nCov+4);
					VectorXd W(nCov+4);
					for (int j = 0; j<nTimes; j++) {
						for(int k=0; k<(nCov+4); k++)
							W(k) = cov.data[j][k];
						M = M + tau2.data[i][j] * (W*W.transpose());
						V = V + tau2.data[i][j] * outs[j] * W;
					}
					VectorXd tmpParam = M.colPivHouseholderQr().solve(V);
					for(int k=0; k<(nCov+4); k++)
						param.push_back( tmpParam(k) );

					// update variances
					// construct means from trend/seasonality parameters
					for(int j = 0; j<nTimes; j++) {
						for(int k=0; k<(nCov+4); k++)
							W(k) = cov.data[j][k];
						mu.data[i][j] = W.transpose()*tmpParam;
					}					
				}
				else{
					//no trend/seasonality modeling
					double mu=0;
					for(int j=0; j<nTimes; j++){
						mu=mu+tau2.data[i][j]*outs[j];
					}
					normf=0;
					for(int j=0; j<nTimes; j++)
						normf=normf+tau2.data[i][j];
					mu=mu/normf;
					param.push_back(mu);
				}			
			}

			// Poisson density
			if(dist[i]==3){
				// with trend/seasonality modeling
				if(trend[i]==1){
					MatrixXd M = MatrixXd::Zero(nCov+4, nCov+4);
					VectorXd V = VectorXd::Zero(nCov+4);
					VectorXd W(nCov+4);
					for (int j = 0; j<nTimes; j++) {
						for(int k=0; k<(nCov+4); k++)
							W(k) = cov.data[j][k];
						M = M + tau2.data[i][j] * (W*W.transpose());
						V = V + tau2.data[i][j] * outs[j] * W;
					}
					VectorXd tmpParam = M.colPivHouseholderQr().solve(V);
					for(int k=0; k<(nCov+4); k++)
						param.push_back( tmpParam(k) );


					// update variances
					// construct means from trend/seasonality parameters
					for (int j = 0; j<nTimes; j++) {
						for(int k=0; k<(nCov+4); k++)
							W << cov.data[j][k];
						mu.data[i][j] = W.transpose()*tmpParam;
					}					
				}
				else{
					//no trend/seasonality modeling
					double mu=0;
					for(int j=0; j<nTimes; j++){
						mu=mu+tau2.data[i][j]*outs[j];
					}
					normf=0;
					for(int j=0; j<nTimes; j++)
						normf=normf+tau2.data[i][j];
					mu=mu/normf;
					param.push_back(mu);
				}			
			}
			// log-normal
			if(dist[i]==4){
				if(trend[i]==1){
					MatrixXd M = MatrixXd::Zero(nCov+4, nCov+4);
					VectorXd V = VectorXd::Zero(nCov+4);
					VectorXd W(nCov+4);
					for (int j = 0; j<nTimes; j++) {
						for(int k=0; k<(nCov+4); k++)
							W(k) = cov.data[j][k];
						M = M + tau2.data[i][j] * (W*W.transpose());
						V = V + tau2.data[i][j] * outs[j] * W;
					}
					VectorXd tmpParam = M.colPivHouseholderQr().solve(V);
					for(int k=0; k<(nCov+4); k++)
						param.push_back( tmpParam(k) );

					// update variances
					// construct means from trend/seasonality parameters
					for (int j = 0; j<nTimes; j++) {
						for(int k=0; k<(nCov+4); k++)
							W(k) = cov.data[j][k];
						mu.data[i][j] = W.transpose()*tmpParam;
					}					
					double sigma=0.0;
					for(int j=0; j<nTimes; j++)
						sigma=sigma+tau2.data[i][j]*(outs[j]-mu.data[i][j])*(outs[j]-mu.data[i][j]);
					sum=0;
					for(int j=0; j<nTimes; j++)
						sum=sum+tau2.data[i][j];
					sigma=sigma/sum;
					param.push_back(sigma);
				}
				else{
					//no trend/seasonality modeling
					// replace zero data with .01
					vector<double> yp=outs;
					for(int t=0; t<nTimes; t++){
						if(yp[t]==0)
							yp[t]=yp[t]+0.01;
					}
					double mu=0;
					double sigma=0;
					for(int j=0; j<nTimes; j++){
						mu=mu+tau2.data[i][j]*log(yp[j])/sum;
						sigma=sigma+tau2.data[i][j]*log(yp[j])*log(yp[j]);
					}
					normf=0;
					for(int j=0; j<nTimes; j++)
						normf=normf+tau2.data[i][j];
					mu=mu/normf;
					sigma=sigma/normf;
					param.push_back(mu);
					param.push_back(sigma-mu*mu);
				}			
			}
			
			//store updated density parameters for this state
			for(unsigned int j=0; j<param.size(); j++)
				theta.data[j][i]=param[j];
		}
		//check for convergence
		converged=(fabs(oldL-logLik)<convthrd) ? true : false;
		if(logLik < oldL){
			cerr << endl << endl << "WARNING: likelihood decreased after " << numit << " iterations." << endl << "Please try again using new initial vaules." << endl;
			exit(EXIT_FAILURE);
		}
		oldL=logLik;
		numit=numit+1;
	}
	logLik=oldL;
	theta=oldtheta;
	pis=oldpis;
	trans=oldtrans;
}

void HMM6::viterbi(){

	//density();
	
	// initialize delta, psi
	for(int i=0; i<nStates; i++)
		deltas.data[i][0]=pis[i]*emis.data[i][0];

	// recursion
	vector<double> pr(nStates);
	double maxp, sum;
	Matrix6<int> psi(nStates, nTimes, 0);
	for(int j=1; j<nTimes; j++){
		for(int i=0; i<nStates; i++){
			maxp=0;
			for(int h=0; h<nStates; h++){
				pr[h]=deltas.data[h][j-1]*trans.data[h][i];
				maxp=max(maxp, pr[h]);
			}
			deltas.data[i][j]=maxp*emis.data[i][j];
			vector<int> idx;
			for(int h=0; h<(int)pr.size(); h++){
				if(pr[h]>=maxp)
					idx.push_back(h);
			}
			// in case we find more than one maximum value, use first one
			psi.data[i][j]=idx[0];
		}
		// scale delta values to keep them from growing too large or too small
		sum=0;
		for(int i=0; i<nStates; i++)
			sum=sum+deltas.data[i][j];
		double scalef=ceil(log10(sum));
		for(int i=0; i<nStates; i++)
			deltas.data[i][j]=deltas.data[i][j]/pow(10.0, scalef);
	}
	// backtracking
	// initialize
	maxp=0;
	vector<int> S;
	for(int h=0; h<nStates; h++)
		maxp=max(maxp, deltas.data[h][nTimes-1]);
	for(int h=0; h<nStates; h++){
		if(deltas.data[h][nTimes-1]>=maxp)
			S.push_back(h);
	}
	// break ties
	path.push_back(S[0]);
	
	// recurse
	for(int j=1; j<nTimes; ++j){
		path[j]=psi.data[S[0]][j];
	}
}

//////////////////////////////////////////// Prediction fuction //////////

void HMM6::predict(){
	for(int t=0; t<nTimes; t++){
		for(int i=0; i<nStates; i++){
			if(path[t]==i){
				if(trend[i]==1)
					for(int k=0; k<(nCov+4); k++)
						fitted[t] = fitted[t] + theta.data[k][i] * cov.data[t][k];
				else
					fitted[t] = theta.data[0][i];
			}
		}
	}

	// goodness of fit
	double SSE=(outs[0]-fitted[0])*(outs[0]-fitted[0]);
	double ybar=outs[0];
	double SSY=0;
	for(int i=1; i<nTimes; i++){
		SSY=SSY+(double)(i)*(outs[i]-ybar)*(outs[i]-ybar)/(double)(i+1);
		ybar=ybar+(outs[i]-ybar)/(double)(i+1);
		SSE=SSE+(outs[i]-fitted[i])*(outs[i]-fitted[i]);
	}
	int maxpar=maxParams+(nStates-1)+nStates*(nStates-1);
	gof[0]=maxpar;
	gof[1]=numit;
	gof[2]=logLik;
	gof[3]=1.0-SSE/SSY;
	gof[4]=1.0-(1.0-gof[3])*(double)(nTimes-1.0)/(nTimes-maxpar-1.0);
	gof[5]=2*maxpar-2.0*logLik;	
	gof[6]=maxpar*log(nTimes)-2.0*logLik;

}