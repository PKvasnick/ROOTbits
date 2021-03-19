/*********************************************************************/
/* Bootstrap demo                                                    */
/*                                                                   */
/* To run from ROOT prompt, do                                       */
/* root [1] .L corr_bootstrap.C++                                    */
/* root [2] corr_bootstrap()                                         */
/*********************************************************************/

#include <TRandomGen.h>
#include <TMath.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <iostream>

/***********************************************************************/
/* Helper functions                                                    */
/***********************************************************************/

/** Calculate dot producct of double arrays x and y. */
double dot_product(int n, double x[], double y[])
{
	double dot = 0;
	for(int i = 0; i<n; ++i) dot += x[i]*y[i];
	return dot;
}

/** Calculate Pearson ocrrelation coefficient of x and y. */
double correlation(int n, double x[], double y[])
{
	double mu_x = TMath::Mean(x, x+n);
	double mu_y = TMath::Mean(y, y+n);
	double std_x = TMath::StdDev(x, x+n);
	double std_y = TMath::StdDev(y, y+n);
	double dot = 0;
	for(int i=0; i<n; ++i) dot+=(x[i]-mu_x)*(y[i]-mu_y);
	return dot/std_x/std_y/(n-1);
}

/***********************************************************************/
/* Main function: Data simulation and calculation of bootstrap CI      */
/***********************************************************************/
void corr_bootstrap() {

/*** Part 1. Data simulation
 * Start by simulating data with a given vector x and correlation coefficient.
 */
	// First create some array for x:
	const int arraysize = 20;
	double *x = new double[20];
	for(int i=0; i<arraysize; ++i) x[i] = i;

	const double rho0 = 0.7;

	// Create y perpendicular to x:
	// Start with initial guess - a random Gaussian vector,
	// Then project onto space perpendicular to x and adjust variance. 

	auto generator= new TRandomMT64(); 

	double *y = new double[arraysize];

	double mu_x = TMath::Mean(x, x + arraysize);
	double std_x = TMath::RMS(x, x + arraysize);

	for(int i=0; i<arraysize; ++i)
	    y[i] = generator->Gaus(0,std_x);
	double mu_y = TMath::Mean(y, y + arraysize);
	double std_y = TMath::RMS(y, y + arraysize);

	for(int i=0; i<arraysize; ++i)
	    y[i] = (y[i]-mu_y)*std_x/std_y;

	mu_y = TMath::Mean(y, y + arraysize);
	std_y = TMath::RMS(y, y + arraysize);

	// Project out the part parallel to x:

	double *normed_x = new double[arraysize];
	double norm_x = 0;
	for(int i=0; i<arraysize; ++i) {
		norm_x += (x[i]-mu_x)*(x[i]-mu_x);
		normed_x[i] = x[i]-mu_x;
	}
	norm_x = TMath::Sqrt(norm_x);
	for (int i = 0; i<arraysize; ++i)
		normed_x[i] = normed_x[i]/norm_x;

	double y_dot_normed_x = 0;
	for(int i=0; i<arraysize; ++i)
		y_dot_normed_x += (y[i]-mu_y)*normed_x[i];

	for(int i=0; i<arraysize; ++i)
		y[i] -= y_dot_normed_x * normed_x[i];

	mu_y = TMath::Mean(y, y+arraysize);
	std_y = TMath::StdDev(y, y+arraysize);

	for(int i=0; i<arraysize; ++i)
		y[i] *= mu_y + std_x/std_y;

	// Test: dot product of x and y should be 0
	double x_dot_y = 0;
	for(int i=0; i<arraysize; ++i)
		x_dot_y += normed_x[i]*(y[i]-mu_y);

	std::cout << "x.y\t" << x_dot_y << std::endl;

	mu_y = TMath::Mean(x, x+arraysize);
	std_y = TMath::StdDev(y, y+arraysize);

	std::cout << "x: Mean\t" << mu_x << "\tStd\t" << std_x << endl;
	std::cout << "y: Mean\t" << mu_y << "\tStd\t" << std_y << endl;

	double *z = new double[arraysize];
	double rho0_comp = TMath::Sqrt(1-rho0*rho0);
	for(int i=0; i<arraysize; ++i) {
		z[i] = rho0 * x[i] + rho0_comp * y[i];
	}

	double rho_rep = correlation(arraysize, x, z);
	std::cout << "Correlation: required:\t" << rho0 << "\treproduced\t"
		<< rho_rep << std::endl;

	// Plot the z vs. y
	TCanvas *cxz = new TCanvas("cxz", "Source data");
	TGraph *g = new TGraph(arraysize, x, z);
	g->Draw("AP*");
	g->SetMarkerColor(kRed);
	g->GetYaxis()->SetTitle("z");
	g->GetXaxis()->SetTitle("x");
	g->SetTitle("Simulated data");

/*** Part 2. Calculate bootstrap CI 	
 * We want statistical properties of the correlation coefficient estimate .
 * Clearly, no assumptions based on (x,y) having bivariate gaussian distribution 
 * are plausible. So we use the (x,y) data to form replicas by randomly sampling 
 * (x_i, y_i) pairs and computing correlation coefficients for each replica. 
 */
	const int n_samples = 1000;
	TH1F *hRep = new TH1F("hRep", "Rho bootstrap distribution", 100, -1, 1);
	double *xs = new double[arraysize];
	double *zs = new double[arraysize];

	for (int i_sample = 0; i_sample < n_samples; ++i_sample) {
		// Generate replica
		for(int i=0; i<arraysize; ++i) {
			int index = generator->Integer(arraysize);
			xs[i] = x[index];
			zs[i] = z[index];
		}
		double rho = correlation(arraysize, xs, zs);
		hRep->Fill(rho);
	}
	TCanvas *cRep = new TCanvas("cRep", "Bootstrap demo");
	hRep->Draw();
	hRep->GetYaxis()->SetTitle("cases");
	hRep->GetXaxis()->SetTitle("rho:");

	// Now we can estimate 95% confidence interval.
	// Note that this is NOT an exact confidence interval, though it is 
	// better than what you can get from textbook formulas. 

	double probs[2] = {0.025, 0.975};
	double quantiles[2];

	hRep->GetQuantiles(2,quantiles, probs);
	std::cout << "Approximate 95% CI for rho:\t" 
		<< quantiles[0] << " - " << quantiles[1] << std::endl;

	double height[2] = {0, 0};

	// ... and make a plot.

	TGraph *gq = new TGraph(2,quantiles, height);
	gq -> Draw("L");
	gq -> SetLineColor(kRed);
	cRep->Update();

	delete [] x;
	delete [] y;
	delete [] normed_x;
	delete [] z;
	delete [] xs;
	delete [] zs;
	
}
