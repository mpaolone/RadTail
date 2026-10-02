#include "intExact.h"

void calcRadTail(){
	
	//CHANGE THESE INPUTS TO MATCH YOUR EXPERIMET:
	double Z = 6;
	double A = 12;
	double be = 400.0; //beam energy in MeV
	double epmin = 200.0; //minimum scattered momentum, in MeV
	double ang = 90.0; // central angle of spectrometer, in deg
	TString rosfile = "inputFF2.dat"; //name of input FF2 file
	double Tb = 0.005; //total radiation thickness before primary scattering
	double Ta = 0.005; //total radiation thickness after primary scattering
	double Mt = 12.000; //Molar mass of target material
	double dth = 0.04: //this the 1/2 size of the spectrometer-theta acceptance in rads.  i.e. if your spectrometer acceptance after cuts is +/- 0.04 radians, then dth=0.04
	double dph = 0.02; //same as dth but for spectrometer-phi acceptance.
	double dp = 0.01; //this is the 1/2 width of your spectrometer-dp bin in percentage of 1.  i.e. if your data point covers a +/- 5% dp, then dp=0.05
	
	//for plotting values
	double epmax = 300.0; //maximum falue of scattered electron energy (for plot)
	double step = 2.0; // step-size of points in momentum.
	
	
	
	//set up integrators;
	intExact Ae(be,epmin,ang,Z,A);
	Ae.LoadFF(rosfile);
	Ae.addMatt(Tb,Ta,Mt);
	Ae.SetAccept(dth, dph, dp);  //this is the 1/2 theta-spectrometer solid angle (rad), 1/2 phi-spectrometer solid angle (rad), and the 1/2 dp-spectrometer bin-width:  dp=0.2 would mean a dp-spectromter cut of +/- 20%.
	ROOT::Math::GaussLegendreIntegrator ig;
	TGraph *gr = new TGraph(); //plot of cross-section versus scattered electron momentum.
	double npoint = (epmax - epmin)/step + 1.0;
	
	for(int k = 0; k < npoint; k++){
		double Ep = epmax - (double)k*step;
		Ae.SetEp(Ep);
		
		ROOT::Math::Functor1D f1d(Ae);
		ig.SetFunction(f1d);
		ig.SetRelTolerance(0.000001); //tolerance vale in integration
		
		float npint1 = 4000; //number of points to use in integration
		ig.SetNumberPoints(npint1);
		double fa = Ae.finiteAcc();
		
		double sigExt = ig.Integral(-1.0,1.0); //does the exact integration from cos(theta) of -1 to 1
		double sigb = Ae.sigb()*conv; //conv is converting units to ub/GeV-sr
		double sigtot = fa*Ae.Fsoft()*(sigExt*conv + sigb);
		gr->SetPoint(k, Ep, sigtot);
	}
	
	gr->SetMarkerStyle(20);
	gr->Draw("AP");
	gr->GetYaxis()->SetTitle("#mub/GeV/sr");
	gr->GetXaxis()->SetTitle("Scattered Electron Energy (MeV)");
}