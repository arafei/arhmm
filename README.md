# arhmm

Extension of a Hidden Markov Model with time-dependent parameters: a package in C++/R

Author: Rafei, Ali

Summary:
The present program is designed for the final project of the course Biostat615, and is aimed to fit a hidden Markov model, where the parameter of the conditional distribution function varies over time. Such model has wide applications in the area of public health especially in early detection epidemics of vector-host infections. Below we have briefly described how to compile and implement the program efficiently. The program utilizes an EM-algorithm to compute the maximum likelihood estimation of parameters, and then, using a Viterbi algorithm, the most likely states sequence of the model is decoded based on the updated parameters. Eventually, some measures of goodness-of-fit including R-square, Adj. R-square, AIC, and BIC are computed to evaluate the fitness of model. We have also provided a real example of 8 years of monthly incidence rates of Malaria data, and showed how the model can apply on such data and detect the epidemic states of the disease. 
Version: 0.0.1
Date: 12-16-2016
Language: The program is written in C++.
Requirements: Boost, Eigen libraries


Introduction:
Hidden Markov Model (HMM) is known as a powerful statistical tool among clustering techniques, whose wide applications have been demonstrated in many research areas including automated speech recognition, electrocardiographic signal analysis, epileptic seizure frequency analysis, DNA sequence analysis, and early detection of outbreaks. Generally, an HMM consists of two sequences, observed sequence {Yt} and hidden sequence {St}, where the observed sequence are conditionally independent, and the state sequence follow a Markov property (Fig. 1). 
 
Fig. 1. Graphical representation of the dependence structure of two sequences in a hidden Markov
In the case of first-order Markov chain of the state sequence with k states, the initial and transition probabilities can be denoted as follow:
pi = P (St=i)   i = 1, 2, … , k;   t = 1, 2,…, n
aij = P (St=j|St-1=i)   i, j = 1, 2, … , k;   t = 1, 2,…, n
Depending on the nature of the observed sequence, the conditional distribution of Yt | St=st can be defined discrete or continuous. For example, in case of count data, we might assume and fit a Poisson distribution, which can be defined: 
 
Three essential questions arise in the context of HMM and the interest is to respond to each based on the observed sequence at hand:
	What’s the probability of observing a certain sequence of data?
	For a given set of observed sequence, how could we estimate the model parameters as well as the initial and transition probabilities?
	Based on the observed sequence, how could we decode the most likely hidden state sequence?
A short answer to each of the above questions would be: (1) based on the conditional rules of probability, the answer to the first question is straightforward. However, a recursive forward-backward algorithm can helps us reduce the volume of computations. (2) EM-algorithm provides an efficient iterative framework to estimate the model parameters. And finally, (3) the algorithm proposed by Viterbi could give us the most optimal solution for the hidden state sequence in terms of the observed sequence.
Aim and Scope:
In some application of HMM, the model parameters are not be constant and may vary over time. For example, in the disease surveillance practice, the epidemic state of a particular infectious disease may be assumed as a hidden state, and the aim is to detect the epidemic state based on a historical time series of the disease incidence data. In such scenarios, the distribution parameter may depends on some external covariates such as the weather conditions at each time period t, or seasonal fluctuation of the disease incidence itself. As a result, we might assume that the parameters can be associated with a number of covariates as well as cyclic terms characterizing the seasonal trends of the disease through a GLM. So at each state of the epidemic or non-epidemic, we can define:
Epidemic:  
Non-Epidemic:  
Where g1 and g2 are the link functions in GLM in the epidemic and non-epidemic states respectively, and X is the matrix of covariates including the cyclic terms   and   where r represents the cyclic period.
The present package has been developed to apply and fit a k-state hidden Markov Model (first-order) on a sequence of observed data, based on a given link function and a given conditional distribution form at each stage. The program will employ an efficient EM-algorithm to estimate the model parameters as well as the initial and transition distribution and the Viterbi algorithm to decode the most likely hidden state sequence of the model. As the disease incidence data typically exist in form of rate or counts, we will limit the scope of conditional distribution to Gaussian, Poisson, and Exponential distribution and the link functions to identity (assuming no link function), Poisson (log) and Logit functions. In addition, the package will be able to evaluate the goodness-of-fit in terms of the some criteria including likelihood-ratio, Akaike Information Criterion (AIC), and Bayesian Information Criterion (BIC), enabling user to explore the best choices for the number of hidden state (k), the assumed conditional distribution, and the link function. Finally, we will assess the performance of the developed package by applying on the real monthly incidence data of Malaria from 2004 to 2012, and will explore the best fitted model in terms of goodness-of-fit. 
Setting up:
The present package is written in C++, and consists of three files, MAIN.cpp, HMM615.h, and Matrix615.h. In order to run the program, you should place them in the same folder. This package uses two external library, BOOST, and EIGEN, so before running, you have to install both package in the main directory ~/include/. Please note that the Eigen library does not require to be installed. You just need to copy to the mentioned path. To compile the program, please go to the directory where the program files copied and then, write the following codes:

g++ -Wall MAIN.cpp -o hmm -I ~/include


If you have already installed the requirements correctly, the program will be compiled and run.
Inputs:
The program has several inputs fit the HMM model that are indicated respectively below: 
	The sequence of observation, which can be either counts or continuous. Examples are the time series of monthly counts of disease incidence or monthly incidence rate of a particular infection. 
To insert the observed sequence, you need to type or copy and paste the values preceded by the name of main command (here hmm for example). Below you can see an example of disease incidence rates for 20 consecutive months:

hmm 3.652672609 8.55293358 11.49462427 12.19306069 16.43053215 15.77670443 11.40707474 5.450159481 1.999244287 1.249250106 1.155633561 1.488427383 2.508935895 5.100489899 9.839811476 9.837792014 7.79769001 8.054552774 6.72421453 5.924722219 


	After calling the program and inserting the sequence of observations, you will be asked to insert the number of hidden states (m). The number of states should be integer and >1. Remember higher values of m increases the number of parameter exponentially and complicates the model and estimation of parameter. For instance, below you can see m is inserted 2:

hmm 3.652672609 8.55293358 11.49462427 12.19306069 16.43053215 15.77670443 11.40707474 5.450159481 1.999244287 1.249250106 1.155633561 1.488427383 2.508935895 5.100489899 9.839811476 9.837792014 7.79769001 8.054552774 6.72421453 5.924722219 

Insert the number of hidden states:
2


	In the next step, we need to specify the type of distribution in each state. The available distributions for this packages with corresponding code are listed below:
1=Gaussian
2=Exponential
3=Poisson, and 
4=Log-Normal
Here Poisson distribution is the only available distribution for counts data. Please note that you are not allowed to use Poisson distribution for continuous observations. The number of distribution parameters varies by type of distribution and depending on assuming trends model in each state. In the below example, a mixture of Gaussian-Gaussian has been assumed for m=2:

hmm 3.652672609 8.55293358 11.49462427 12.19306069 16.43053215 15.77670443 11.40707474 5.450159481 1.999244287 1.249250106 1.155633561 1.488427383 2.508935895 5.100489899 9.839811476 9.837792014 7.79769001 8.054552774 6.72421453 5.924722219 

Insert the number of hidden states:
2

Insert the distribution type for each state:
1=Gaussian, 2=Exponential, 3=Poisson, and 4=LogNormal
1 1


	Now you will be asked to determine if you want to apply a periodic model for the location parameter of the conditional distribution at each state of HMM. 1 represents trend and 0 represents non-trend. You are allowed to specify trends parameters for one states and without trends in the second state. Example of output of trends status for m=2 can be seen below. That is, we assume to apply a periodic regression model for each state of HMM:


hmm 3.652672609 8.55293358 11.49462427 12.19306069 16.43053215 15.77670443 11.40707474 5.450159481 1.999244287 1.249250106 1.155633561 1.488427383 2.508935895 5.100489899 9.839811476 9.837792014 7.79769001 8.054552774 6.72421453 5.924722219 

Insert the number of hidden states:
2

Insert the distribution type for each state:
1=Gaussian, 2=Exponential, 3=Poisson, and 4=LogNormal
1 1

Insert the vector of trends:
0=No trend, 1=Trend
1 1



	In the next three steps, we are asked to insert the initial values of parameters at each state of the model. The first input relates to the conditional distribution parameters. Please note that if you have specified trends in a state, the minimum number of parameters would be 4 (ß_0,ß_1,ß_2,ß_3). In our example with mixture distribution of Normal-Normal and trends at each state, we need 5 parameters at each state, four for trends and one for the distribution variance. For more information please refer to the theories in the first section.
Then, the initial values of the initial probabilities should be inserted. A number of m probabilities should be entered. Please note that the sum of initial probabilities should be 1. The final input is the initial values of transition probabilities. For each state, a number of probabilities is m. Again the sum of the probabilities at each state should be equal to 1. Please you can see an example of initial values of parameters.
Please note that if the iterative process of EM-algorithm doesn’t converge, the program returns error and asks the user to change new initial values of parameters.


hmm 3.652672609 8.55293358 11.49462427 12.19306069 16.43053215 15.77670443 11.40707474 5.450159481 1.999244287 1.249250106 1.155633561 1.488427383 2.508935895 5.100489899 9.839811476 9.837792014 7.79769001 8.054552774 6.72421453 5.924722219 

Insert the number of hidden states:
2

Insert the distribution type for each state:
1=Gaussian, 2=Exponential, 3=Poisson, and 4=LogNormal
1 1

Insert the vector of trends:
0=No trend, 1=Trend
1 1

Insert the initial values of parameters for the state 1:
Example: B0, B1, B2, B3, Sigma for Normal Dist
3  0.0   -1  0.7    1
Insert the initial values of parameters for the state 2:
Example: B0, B1, B2, B3, Sigma for Normal Dist
8 -0.1   -3  1.0    4

insert initial probabilities:
0 1

insert transition probabilities from state: 1
0.9 0.1
insert transition probabilities from state: 2
0.2 0.8


Outputs:
The output of the program can be summarized as following:
	The maximum likelihood estimation of parameters of the conditional distributions at each state of the model using EM-algorithm.
	The updated values of initial probabilities and transition probabilities using EM-algorithm.
	The most likely state sequence of the model using Viterbi algorithm.
	The fitted values of the estimated model
	Some measures of goodness-of-fit including R2, Adj-R2, AIC and BIC, enabling us to compare different models in terms of goodness of fit.

Estimation of parameters:
     1.426     0.213    -2.499      1.24   0.02673

     8.507   -0.2559    -7.176     3.919     1.391

Estimation of initial probabilities:
         0         1

Estimation of transition probabilities:
    0.8688    0.1312

    0.1794    0.8206

The most likely state sequence:
0 1 1 1 1 1 1 1 1 0 0 0 0 0 1 1 1 0 0 0

Fitted values of the model:
0.09514 7.801 11.66 14.47 15.4 14.15 10.97 6.654 2.286 1.234 0.9857 1.484 2.652 4.233 8.588 11.39 12.33 7.76 7.018 5.863

    #Parameters    #Iterations     Log-likelihood       R-square   adj.R-Square              AIC                 BIC
                      13                      8                  -21.87           0.8801             0.8002          69.74             82.68

Example:
We applied and tested the model on the monthly incidence rate (per 100,000) of Malaria, reported from 1996 to 2004 in a high-risk area in Iran. Malaria is known as an infection with intensive seasonal effects.
Below you can see the final estimate of model as well as the time series of observations and fitted values predicted by model.

Estimation of parameters:
     3.053  -0.01209    -1.639    0.5957     0.259

     7.465  -0.06704    -3.671      1.54     5.639

Estimation of initial probabilities:
         0         1

Estimation of transition probabilities:
    0.9017   0.09831

    0.2028    0.7972

The most likely state sequence:
0 1 1 1 1 1 1 1 1 0 0 0 0 0 1 1 1 1 1 1 1 0 0 0 0 0 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0

Fitted values of the model:
1.919 6.83 8.805 10.37 11.08 10.73 9.405 7.43 5.321 1.596 1.202 1.268 1.773 2.58 8 9.562 10.27 9.929 8.6 6.626 4.517 1.451 1.057 1.123 1.628 2.434 7.196 8.758 9.47 9.125 7.796 5.821 3.712 1.306 0.9118 0.978 1.483 2.289 3.177 3.905 4.275 4.184 6.991 5.017 1.913 1.161 0.7668 0.833 1.338 2.144 3.032 3.76 4.13 4.039 3.51 2.679 1.768 0.4072 0.6218 0.6879 1.193 1.999 2.887 3.615 3.985 3.894 3.365 2.534 1.623 0.8708 0.4767 0.5429 1.048 1.854 2.742 3.47 3.839 3.749 3.22 2.389 1.478 0.7258 0.3317 0.3978 0.9033 1.709 2.597 3.325 3.694 3.604 3.773 1.799 1.333 0.5807 0.1866

    #Parameters    #Iterations     Log-likelihood       R-square   adj.R-Square            AIC            BIC
                       13                   25                  -146.2           0.7682               0.731        318.4        351.6




The final estimates of the model are:
Non-Epidemic:  
Epidemic:  
The results of model goodness-of-fit show high degrees of fitness.

 
Figure 1 : Monthly incidence rate of Malaria vs fitted values by model. Red segments are the periods detected as epidemic by HMM

Discussion:
We can develop this model further in the future by including auto-regressive terms as well as other external covariates and also assuming GLM with specific link function to improve the fitness of model. One other interesting capability that can be added to the model in the future is inference and interval estimates for parameters through bootstrapping techniques.
