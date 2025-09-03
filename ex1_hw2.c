#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include <string>

using namespace std;

//Determination of the ln(L(σ)) function for rayleigh distribution:
double log_likelihood_function(double sum, double sum_square, double product, double sigma, double N)
{
    double log_likelihood = log(product) - 2 * N * log(sigma) - 1 / (2 * pow(sigma,2)) * sum_square;

    return log_likelihood;
}

void ex1_hw2()
{   
    //Determination of the observations and the size N:
    double observations[] = {0.96, 1.12, 0.85, 1.02, 1.58, 1.86, 0.79, 0.82, 0.45, 1.52};
    double N = sizeof(observations) / sizeof(observations[0]);

    //Sum of the observations:
    double sum = 0.0;
    for (int i = 0; i < N; i++) { sum += sum + observations[i]; }

    //Sum square of the observations:
    double sum_square = 0.0;
    for (int i = 0; i < N; i++) { sum_square += pow(observations[i], 2); }

    //Product of the observations:
    double product = 1.0;
    for (int i = 0; i < N; i++) { product *= observations[i]; }
    
    //The sigma from the maximization of the log likelihood:
    double sigma_hat = pow(sum_square / (2 * N), 0.5);

    //The Hessian Matrix:
    double Hessian = 4 * N / pow(sigma_hat, 2);

    //The covariance matrix (variance) for the sigma hat:
    double Variance = 1 / Hessian;

    //The error of the estimator for sigma:
    double error_of_sigma_hat = pow(Variance, 0.5);

    std::cout << "SIGMA FOR MAXIMUM LOG LIKELIHOOD     =  " << sigma_hat          << " " << std::endl;
    std::cout << "HESSIAN MATRIX                       =  " << Hessian            << " " << std::endl;
    std::cout << "VARIANCE OF THE SIGMA                =  " << Variance           << " " << std::endl;
    std::cout << "ERROR OF THE SIGMA                   =  " << error_of_sigma_hat << " " << std::endl;
    
    //The maximum log-likelihood:
    double maximum_likelihood = log_likelihood_function(sum, sum_square, product, sigma_hat, N);

    std::cout << "THE MAXIMUM LOG LIKELIHOOD IS        = " << maximum_likelihood  << " " << std::endl;
    
    //Create vectors to store sigma and log likelihood values:
    double min_sigma_value = 0.5;
    double max_sigma_value = 1.4;
    int noOfSigmaValues = 501;

    double sigma_values[noOfSigmaValues];
    double log_likelihood_values[noOfSigmaValues];
    double parabolic_estimation_values[noOfSigmaValues];

    //Calculate log likelihood for different sigma values:
    for (int i = 0; i < noOfSigmaValues; i++)
    {   
        double sigma = min_sigma_value + i * (max_sigma_value - min_sigma_value) / noOfSigmaValues; 
        double log_likelihood = log_likelihood_function(sum, sum_square, product, sigma, N);

        //Parabolic approach of log likelihood:
        double parabolic_estimation = maximum_likelihood - 0.5 / Variance * pow(sigma - sigma_hat, 2);

        sigma_values[i] = sigma;
        log_likelihood_values[i] = log_likelihood; 
        parabolic_estimation_values[i] = parabolic_estimation;

    }

    // Create the log likelihood plot in terms of sigma:
    TCanvas *c1 = new TCanvas("c1", "Likelihood function", 800, 600);
    TGraph *graph1 = new TGraph(noOfSigmaValues, sigma_values, log_likelihood_values);
    
    //Plot the ln(L(sigma_hat))-0.5 line to find the error:
    double maximum_likelihood_array[noOfSigmaValues];
    for (int i = 0; i < noOfSigmaValues; i++) { maximum_likelihood_array[i] = maximum_likelihood - 0.5; }
    TGraph *graph2 = new TGraph(noOfSigmaValues, sigma_values, maximum_likelihood_array);

    //Plot the parabolic estimation:
    TGraph *graph3 = new TGraph(noOfSigmaValues, sigma_values, parabolic_estimation_values);

    // Set graph title and axis labels
    graph1->SetTitle("Log-Likelihood Function; #sigma; ln[L(#sigma)]");

    // Draw the graphs:
    graph1->Draw("AC");
    graph1->SetLineWidth(2);
    graph1->SetLineColor(kGreen);
    graph2->Draw("SAME");
    graph2->SetLineWidth(2);
    graph2->SetLineColor(kRed);
    graph3->Draw("SAME");
    graph3->SetLineWidth(2);
    graph3->SetLineColor(kBlue);

    // Create legend:
    TLegend *legend = new TLegend(0.15, 0.75, 0.35, 0.85);
    legend->AddEntry(graph1, "Log-Likelihood Function", "l");
    legend->AddEntry(graph2, "ln[L(sigma)]-1/2", "l");
    legend->AddEntry(graph3, "Parabolic approach", "l");
    legend->Draw();

    // Display the canvas
    c1->Draw();

}