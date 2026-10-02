//Written by Michael Paolone, NMSU, 2021

#include "Math/Functor.h"
#include "Math/WrappedTF1.h"
#include "Math/GaussIntegrator.h"
#include "Math/GaussLegendreIntegrator.h"
#include "Math/IntegratorMultiDim.h"
#include "TLorentzVector.h"
#include <TROOT.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TMath.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>

const float alpha =   1.0/137.03604;
const float pi = 3.1415926535;
const float HBARC = 197.32e-13;
const float Ma = 931.478;
const float Mp = 938.272;
const float Me = 0.511;
const float conv = HBARC*HBARC*1.0e36; //from cm/MeV-sr to ub/GeV-sr

class intExact{
private:
	double Es; //Incomming electron energy, MeV
	double Ep; //Outgoing electron energy, MeV
	double th; //central angle of scattering, degrees
	double Z; //atomic number Z of target material
	double A; //atomic number A of target material
	TGraph *FF; //Form factor squared plot
	double retVal_; //return value
	double bttot; //the total summed bt (radiation thickness with "b" factor)
	double tb; //radiation thickness before primary scattering (unitless:  thickness(cm)*density(g/cm^3)/radlength(g/cm^2))
	double ta; //radiation thickness before primary scattering (unitless:  thickness(cm)*density(g/cm^3)/radlength(g/cm^2))
	double tr; //internal radiator radiation thickness
	double q2min; //minimum Q2 in FF input.  Auto set when loading extrnal form factor input
	double q2max; //maximum Q2 in FF input. Auto set when loading extrnal form factor input
	double Q2;
	double eta; //factor used in calculating "b" factor
	double br; // the "b" factor used in radiation thickness "bt"
	double Mt; // mass of the target at the vertex of elastic scattering
	double dth; //the half theta(spectrometer) window for integration in radians. (usually the +/- solid angle cut on the spectrometer)
	double dph; //the half phi(spectrometer) window for integration in radians. (usually the +/- solid angle cut on the spectrometer)
	double dp; //the dp(spectrometer) bin size for the calculated point. Units in percentage of 1:   i.e. 100% would be a value of dp = 1,  20% would be dp = 0.2
public:
	intExact(double Es_, double Ep_, double th_, double Z_, double A_):Es(Es_),Ep(Ep_),th(th_),Z(Z_),A(A_){
		FF = new TGraph();
		retVal_ = 0;
		Q2 = 4.0*Ep*Es*pow(sin(th*pi/180.0/2.0),2);
		q2max = 0.0;
		q2min = 0.0;
		bttot = 0.0;
		eta = log(1440.0*pow(Z,-2.0/3.0))/log(183.0*pow(Z,-1.0/3.0));
		br = 4.0/3.0*(1.0 + 1.0/9.0*(Z+1)/(Z+eta)/log(183.0*pow(Z,-1.0/3.0)));
		tb = 0.0;
		ta = 0.0;
		tr = alpha/pi*(log(Q2/Me/Me) - 1.0)/br;
		Mt = A*Ma;
		dth = 0.04*0.7;  //generic from CSR analysis, can be set later.
		dph = 0.02*0.7;  //generic from CSR analysis, can be set later. 
		dp = 0.035/3.0;  //generic from CSR analysis, can be set later. 
		if(th_ < 20.0){
			dp = 0.035/11.0;
		};
	}
	void SetEp(double newEp){
		Ep = newEp;
		Q2 = 4.0*Ep*Es*pow(sin(th*pi/180.0/2.0),2);
		tr = alpha/pi*(log(Q2/Me/Me) - 1.0)/br;
	}
	void SetTh(double newTh){
		th = newTh;
		Q2 = 4.0*Ep*Es*pow(sin(th*pi/180.0/2.0),2);
		tr = alpha/pi*(log(Q2/Me/Me) - 1.0)/br;
	}
	void SetAccept(double dthi, double dphi, double xdp){
		//number of points in XS calc of data
		if(xdp > 0.99){
			dp = 0.035/xdp;
		}else{
			//exact dp bin size
			dp = xdp;
		}
		dth = dthi;
		dph = dphi;
	}
	void LoadFF(TString infile){
		FF = new TGraph(infile);
		FF->Set(FF->GetN() -1);
		q2min = TMath::MinElement(FF->GetN(),FF->GetX());
		q2max = TMath::MaxElement(FF->GetN(),FF->GetX());
	}
	//NOTE:  tbi and tai can be the sum radiation length of mixed materials, but MolM MUST BE the molar mass of the primary scattering target.
	void addMatt(double tbi, double tai, double MolM){
		tb = tbi;
		ta = tai;
		bttot = br*(tb + ta);
		Mt = Ma*MolM;
	}
	//double calcRet(const double &cosk) const{
	double calcRet(const double &cosk){
		
		double Mt = A*Ma;
		double pmag = sqrt(Ep*Ep - Me*Me);
		double smag = sqrt(Es*Es - Me*Me);
		double cosScat = cos(th*pi/180.0);
		
		double u0 = Es + Mt - Ep; //A27
		double sdotp = Es*Ep - smag*pmag*cosScat; //A26
		double u2 = 2.0*Me*Me + Mt*Mt -2.0*sdotp + 2.0*Mt*(Es -Ep); //A28
		double umag = sqrt(u0*u0 - u2);  //A28
		double w = 0.5*(u2 - Mt*Mt)/(u0 - umag*cosk); // A25
		double q2 = 2.0*Me*Me - 2.0*sdotp - 2.0*w*(Es - Ep) + 2.0*w*umag*cosk; //A30
		double cosp = (smag*cosScat - pmag)/umag; //A35
		double coss = (smag - pmag*cosScat)/umag; //A36
		double a = w*(Ep - pmag*cosp*cosk); //A31
		double ap = w*(Es - smag*coss*cosk); //A32
		double bp = - w*pmag*sqrt(1.0 - cosp*cosp)*sqrt(1.0 - cosk*cosk); //A33
		double x = sqrt(a*a -bp*bp);
		double y = sqrt(ap*ap - bp*bp);
		double v = 1.0/(ap - a);
		
		double TERM1 = pow(alpha,3)/(2.0*pi)*(Ep/Es);
		double TERM2 = 2.0*Mt*w/q2/q2/(u0 - umag*cosk);
		double TERM3 = W2bar(q2)*(
			-a*Me*Me/x/x/x*(2.0*Es*(Ep+w) + q2/2.0) - ap*Me*Me/y/y/y*(2.0*Ep*(Es -w) + q2/2.0)
				-2.0 + 2.0*v*(1.0/x - 1.0/y)*(Me*Me*(sdotp - w*w) + sdotp*(2.0*Es*Ep - sdotp + w*(Es -Ep)))
					+ 1.0/x*(2.0*(Es*Ep + Es*w + Ep*Ep) + q2/2.0 - sdotp - Me*Me)
						- 1.0/y*(2.0*(Es*Ep - Ep*w + Es*Es) + q2/2.0 - sdotp - Me*Me));
		
		double TERM4 = 0.0; //for heavy nuclie C12 and above
		if(Z < 6){
		
			//TERM4 = W1bar(q2)*((a/x/x/x + ap/y/y/y)*Me*Me*(2.0*Me*Me + q2) + 4.0 + 4.0*v*(1.0/x - 1.0/y)*sdotp*(sdotp - 2.0*Me*Me)
			//	+(1.0/x + 1.0/y)*(2.0*sdotp + 2.0*Me*Me - q2));
			TERM4 = 0.0; //NOT CORRECT FOR Z < 6:  Need to add W1bar(Q2) function and use form above!
		}
		return TERM1*TERM2*(TERM3 + TERM4);
	}
	//double Phi(double x) const{
	double Phi(double x){
		return 1.0 - x + 3.0/4.0*x*x;
	}
	//double FBar(double q2) const{
	double FBar(double q2){
		return 1.0 + 0.5772*bttot + 2.0*alpha/pi*(-14.0/9.0 + 13.0/12.0*log(-q2/Me/Me)) 
			- alpha/(2.0*pi)*pow(log(Es/Ep),2) 
			+ alpha/pi*(1.0/6.0*pi*pi - TMath::DiLog(Phi(pow(cos(th*pi/180.0/2.0),2))));
	}
	//double W2bar(double q2) const{
	double W2bar(double q2){
		double Q2loc = sqrt(-q2)/197.0;
		if(Q2loc < q2min){
			cout << "W2bar Q2 out of rang for FF!!  " << Q2loc << " < " << q2min << endl;
			return 0.0;
		}
		if(Q2loc > q2max) Q2loc = q2max; 
		return FBar(q2)*FF->Eval(Q2loc)*Z*Z;
	}
	double XSel(double Esl){
		double eprm = recoil(Esl);
		double Q2rec = 4.0*Esl*eprm*pow(sin(th*pi/180.0/2.0),2);
		double Q2loc = sqrt(Q2rec)/197.0;
		if(Q2loc < q2min || Q2loc > q2max){
			cout << "XSel Q2 out of rang for FF!!  " << Q2loc << endl;
			return 0.0;
		}
		
		double mott = pow(2.0*alpha*eprm*Z*cos(th*pi/180.0/2.0)/Q2rec,2);
		return FF->Eval(Q2loc)*mott;
	}
	double XSel(double Esl, double thi){
		double eprm = recoil(Esl, thi);
		double Q2rec = 4.0*Esl*eprm*pow(sin(thi*pi/180.0/2.0),2);
		double Q2loc = sqrt(Q2rec)/197.0;
		if(Q2loc < q2min || Q2loc > q2max){
			cout << "XSel Q2 out of rang for FF!!  " << Q2loc << endl;
			return 0.0;
		}
		
		double mott = pow(2.0*alpha*eprm*Z*cos(thi*pi/180.0/2.0)/Q2rec,2);
		return FF->Eval(Q2loc)*mott;
	}
	double sigb(){
		double sin2 = pow(sin(th*pi/180.0/2.0),2);
		double ws = Es - Ep/(1.0 - (2.0*Ep/Mt)*sin2);
		double wp = Es/(1.0 + (2.0*Es/Mt)*sin2) - Ep;
		double zeta = pi*Me/(2.0*alpha)*(tb + ta)/(Z + eta)/log(183.0*pow(Z,-1.0/3.0));
		double vs = ws/Es;
		double vp = wp/(Ep + wp);
		double TERM1 = (Mt + 2.0*(Es - ws)*sin2)/(Mt - 2.0*Ep*sin2);
		double Q2a = 4.0*(Es-ws)*recoil(Es-ws)*sin2;
		double TERM2 = FBar(-Q2a)*XSel(Es - ws)*(br*tb/ws*Phi(vs) + zeta/(2.0*ws*ws));
		double TERM3 = FBar(-Q2)*XSel(Es)*(br*ta/wp*Phi(vp) + zeta/(2.0*wp*wp));
		return TERM1*TERM2 + TERM3;
	}
	double sigb(double thi, double Epi){
		double sin2 = pow(sin(thi*pi/180.0/2.0),2);
		double ws = Es - Epi/(1.0 - (2.0*Epi/Mt)*sin2);
		double wp = Es/(1.0 + (2.0*Es/Mt)*sin2) - Epi;
		double zeta = pi*Me/(2.0*alpha)*(tb + ta)/(Z + eta)/log(183.0*pow(Z,-1.0/3.0));
		double vs = ws/Es;
		double vp = wp/(Epi + wp);
		double Q2a = 4.0*(Es-ws)*recoil(Es-ws, thi)*sin2;
		double Q2b = 4.0*(Es)*recoil(Es, thi)*sin2;
		double TERM1 = (Mt + 2.0*(Es - ws)*sin2)/(Mt - 2.0*Epi*sin2);
		double TERM2 = FBar(-Q2a)*XSel(Es - ws, thi)*(br*tb/ws*Phi(vs) + zeta/(2.0*ws*ws));
		double TERM3 = FBar(-Q2b)*XSel(Es, thi)*(br*ta/wp*Phi(vp) + zeta/(2.0*wp*wp));
		return TERM1*TERM2 + TERM3;
	}
	double Fsoft(){
		double sin2 = pow(sin(th*pi/180.0/2.0),2);
		double ws = Es - Ep/(1.0 - (2.0*Ep/Mt)*sin2);
		double wp = Es/(1.0 + (2.0*Es/Mt)*sin2) - Ep;
		return pow(ws/Es,br*(tb + tr))*pow(wp/(Ep + wp),br*(ta + tr));
	}
	double Fsoft(double thi, double Epi){
		double sin2 = pow(sin(thi*pi/180.0/2.0),2);
		double ws = Es - Epi/(1.0 - (2.0*Epi/Mt)*sin2);
		double wp = Es/(1.0 + (2.0*Es/Mt)*sin2) - Epi;
		double Q2a = 4.0*Epi*Es*sin2;
		double tra = alpha/pi*(log(Q2a/Me/Me) - 1.0)/br;
		return pow(ws/Es,br*(tb + tra))*pow(wp/(Epi + wp),br*(ta + tra));
	}
	double recoil(double Esl){
		return Esl/(1.0 + 2.0*Esl*pow(sin(th*pi/180.0/2.0),2)/Mt);
	}
	double recoil(double Esl, double thi){
		return Esl/(1.0 + 2.0*Esl*pow(sin(thi*pi/180.0/2.0),2)/Mt);
	}
	
	double sigPeak(){
		double sin2 = pow(sin(th*pi/180.0/2.0),2);
		double ws = Es - Ep/(1.0 - (2.0*Ep/Mt)*sin2);
		double wp = Es/(1.0 + (2.0*Es/Mt)*sin2) - Ep;
		double vs = ws/Es;
		double vp = wp/(Ep + wp);
		double TERM1 = (Mt + 2.0*(Es - ws)*sin2)/(Mt - 2.0*Ep*sin2);
		double Q2a = 4.0*(Es-ws)*recoil(Es-ws)*sin2;
		double TERM2 = FBar(-Q2a)*XSel(Es - ws)*(br*tr*Phi(vs)/ws);
		double TERM3 = FBar(-Q2)*XSel(Es)*(br*tr*Phi(vp)/wp);
		
		return TERM1*TERM2 + TERM3;
	}
	
	double sigPeak(double thi, double Epi){
		double sin2 = pow(sin(thi*pi/180.0/2.0),2);
		double ws = Es - Epi/(1.0 - (2.0*Epi/Mt)*sin2);
		double wp = Es/(1.0 + (2.0*Es/Mt)*sin2) - Epi;
		double vs = ws/Es;
		double vp = wp/(Epi + wp);
		double Q2a = 4.0*(Es-ws)*recoil(Es-ws, thi)*sin2;
		double Q2b = 4.0*(Es)*recoil(Es, thi)*sin2;
		double TERM1 = (Mt + 2.0*(Es - ws)*sin2)/(Mt - 2.0*Epi*sin2);
		double TERM2 = FBar(-Q2a)*XSel(Es - ws,thi)*(br*tr*Phi(vs)/ws);
		double TERM3 = FBar(-Q2b)*XSel(Es, thi)*(br*tr*Phi(vp)/wp);
		
		return TERM1*TERM2 + TERM3;
	}
	
	
	//double operator()(const double &x) const {
	//for cosx integration
	
	double finiteAcc(){
		double np = 10; // number of points to integrate over for finite acceptance correction in each dimmension. Since there are 3 dimensions, this is 10*10*10 = 10000 calculated points that integrated.
		double dthstep = 2.0*dth/(np -1.0);
		double dphstep = 2.0*dph/(np -1.0);
		double dpstep = 2.0*dp/(np -1.0);
		
		double sum = 0.0;
		
		for(int i = 0; i < np; i++){
			float dthi = -dth + (double)i*dthstep;
			for(int j = 0; j < np; j++){
				float dphi = -dph + (double)j*dphstep;
				double sin2 = pow(sin((th*pi/180.0 + dthi)/2),2) + pow(sin((dphi)*pi/180.0/2),2) - 2.0*pow(sin((th*pi/180.0 + dthi)/2),2)* pow(sin((dphi)*pi/180.0/2),2);
				double thi = 2.0*asin(sqrt(sin2))*180.0/pi;
				for(int k = 0; k < np; k++){
					float dpi = -dp + (double)k*dpstep;
					double Epi = Ep*(1.0 + dpi);
					if(Epi > recoil(Es)){
						//do nothing, not physical
					}else
						sum += Fsoft(thi,Epi)*(sigPeak(thi,Epi) + sigb(thi,Epi));
				}
			}
		}
		
		sum *= 1.0/pow(np,3);
		double cent = Fsoft()*(sigPeak() + sigb());
		return sum/cent;
	}
	
	double GetW(){
		return sqrt(Mp*Mp -Q2 + 2.0*(Mp*(Es - Ep)));
	}
	
	double operator()(const double &x){
		//const double retval = calcRet(x);
		double retval = calcRet(x);
		return retval;
	}
	//for finite accept integration
	
};

