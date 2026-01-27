//********************************************************************************
//  TridentEventGenerator
//
//  version 1.0 by Wolfgang Altmannshofer (01/18/2019)
//  version 2.0 modification by Joseph Bateman (01/15/2025)
//
//********************************************************************************

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <random>
#include <ctime>
#include <iomanip>
#include <sstream>
#include <iterator>
#include <vector>
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TKey.h"

#include "form_factors.h"

using namespace std;
using std::fstream;



//*****************************************************************
// global constants and input parameters
//*****************************************************************
bool debug_mode = false;
//*****************************************************************
double pi = 3.1415926535897;
double Twopi = 2*pi;
//*****************************************************************
double GF = 0.0000116637; // Fermi constant
double aem = 1/137.036; // em coupling
double sW2 = 0.23119; // sin^2 of the weak mixing angle
//*****************************************************************
double me = 0.000511; // electron mass
double me2 = me*me;
double mmu = 0.105658; // muon mass
double mmu2 = mmu*mmu;
double mtau = 1.778; // tau mass
double mtau2 = mtau*mtau;

//*****************************************************************
// Define masses of various nuclei
//*****************************************************************
double MCarbon = 12.01*0.9315;
double MNitrogen = 14.01*0.9315;
double MOxygen = 16.0*0.9315;
double MSilicon = 28.09*0.9315;
double MArgon = 39.95*0.9315;
double MIron = 55.85*0.9315;
//*****************************************************************
double Mproton = 0.938272;
double Mneutron = 0.939565;
//*****************************************************************
// magnetic moments of proton and neutron
double mup = 2.7928;
double mun = -1.913;

//*****************************************************************
// form factor parameters for proton and neutron
//*****************************************************************
double bpM1 = 12.31;
double bpM2 = 25.57;
double bpM3 = 30.61;
double apM1 = 1.09;
double bpE1 = 11.12;
double bpE2 = 15.16;
double bpE3 = 21.25;
double apE1 = -0.19;
double ffA = 1.68;
double ffB = 3.63;
double bnM1 = 21.3;
double bnM2 = 77;
double bnM3 = 238;
double anM1 = 8.28;

//*****************************************************************
// nuclear parameters, effective couplings, lepton masses ...
//*****************************************************************
double A, Z, M, M2;
double GV, GA, GVSM, GASM;
double m3, m4;
double InterpolationList[192][2]; // array used for interpolation of form factor grid

//*****************************************************************
// variables for reading the neutrino flux distribution file
//*****************************************************************
vector<double> distribution_list;
vector<double> Enumin_list;
vector<double> Enumax_list;
vector<double> probability_list;

int length_probability_list;
int bin;
string flux_file;

TH2D Evsz_hist;

//*****************************************************************
// variables for getting an interaction vertex
//*****************************************************************

// cm and ns
double xmin = 0;
double xmax = 256.0; 
double ymin =  -116.3; 
double ymax = 116.3 ;
double zmin = 0;
double zmax = 1036.8;

double radius_decay_pipe = 35.0; // cm

// BNB coordinates origin in MicroBooNE coordinates (no rotation required)
double BNB_x_offset  = 128.175; // cm
double BNB_y_offset  = -0.93; // cm
double BNB_z_offset  = -46336.3525; // cm

// time offsets
double MicroBooNEGlobalTimeOffset = 3125.0; //ns
double MicroBooNERandomTimeOffset = 1600.0; // ns


//*****************************************************************
// flags
//*****************************************************************
int zeroweight;
int anti;
int readerror = 0;
int flagkinematics = 0;
int errorcounter = 0;
int eventcounter = 0;
int reweightcounter = 0;
int output_format;

//*****************************************************************
// variables for data output
//*****************************************************************
double *data;
int PDG1, PDG2, PDG3, PDG4;


//*****************************************************************
// Variables entered in the users interface
//*****************************************************************
int process;
string material;
string energy_type;
string model;
string command;
string filename_in, filename_out;
double Enu, GP, MZP;
double GVtot, GAtot; 
int Nevents;


//*****************************************************************
// Variables used for random point generation
//*****************************************************************
double eps1, g;
double x1max, x1min;
double x0, x1, x2, x3, x4, x5, x6, x7, x8;
double u1max, u2max, u3max, u4max, u5max, u6max, u7max;
double u1min, u2min, u3min, u4min, u5min, u6min, u7min;
double u1, u2, u3, u4, u5, u6, u7;
double D3max, D3min, D3, D32;
double D4max, D4min, D4, D42, D3D4;
double yplus, yminus;
double phi3, phi4, P2, P5, WC, WC2, WCqc, WB, WB2, qc;
double P0C, PC, q0C, qC, eps1C, eps4C, p4C, P4; 
double cos4, cosPq, sin4, sinPq, cosP1q, cosPp1, sinP1q, sinPp1, cosPhi, sinPhi;
double p1P, q2, qp1, qp3, qp4, p1p3, p1p4, p3P, p4P, Pq, qP, p2P, qp2, p1p2, p2p3, p2p4, p3p4;


//*****************************************************************
// Variables used to define events
//*****************************************************************
double Theta, sintheta, costheta;
double Eq, qx, qy, qz;
double Enuin, pnuinx, pnuiny, pnuinz;
double Enuout, pnuoutx, pnuoutz, pnuouty; 
double Elminus, plminusx, plminusy, plminusz; 
double Elplus, plplusz, plplusy, plplusx;
double Eprime, Pprimez, Pprimey, Pprimex;
double vtx_x, vtx_y, vtx_z, vtx_t;
double prod_x, prod_y, prod_z, prod_t;

double event[7][4];


//*****************************************************************
// Variables for weight and crosssection calculations
//*****************************************************************
double J, PLP, Laa, deltaSigma, weight, maxweight, averageweight;
double crosssectionresult;
double deltacrosssectionresult;


//*****************************************************************
// functions
//*****************************************************************
void SetTridentProcess();
void SetNuclearParameters();
void GenerateEvent(bool);
void GenerateEvents();
void GenerateRandomPoint();
void DetermineWeight();
void FindMaxWeight();
void ComputeCrossSection();
void ReadDistribution();
void LoadFluxFromROOT(string, int);
void GenerateVertex();
void WriteEventFile(string);
void WriteEventHepMC3(string);
void WriteEventHepevt(string);
void WriteXsecTempFile(string);
double SquaredMatrixElementPLP(double, double);
double SquaredMatrixElementPLPanti(double, double);
double SquaredMatrixElementL(double, double);
double SquaredMatrixElementLanti(double, double);
double Interpolation(double, int);
double InverseInterpolation(double, int);
double GMp(double);
double GEp(double);
double GMn(double);
double GEn(double);
double GD(double);
double Pauli(double);
vector<double>  SampleEvsZhist(double, string);
 
//***************************************************
    bool is_nan(double x) { return x != x; }  
//***************************************************

//***************************************************
    std::random_device entropySource;
    std::mt19937 generator(entropySource());
    std::uniform_real_distribution<double> realdistribution(0.0,1.0);
//***************************************************

//*******************************************************
// Interpolation functions for the form factor grids
//*******************************************************

double Interpolation(double x, int r){

    double x1, y1, x2, y2, slope, interp;

    for (int i = 0; i < r-1; i++){
        if (InterpolationList[i][0] < x && x < InterpolationList[i+1][0]){
            // Set the x and y values of the two points
            x1 = InterpolationList[i][0];
            x2 = InterpolationList[i+1][0];
            y1 = InterpolationList[i][1];
            y2 = InterpolationList[i+1][1];
        }
    }
    
    // Calculate the slope of the line between the two points
    slope = (y2-y1)/(x2-x1);
    // Calculate the y value for the given x value
    interp = slope*(x-x1)+y1;

    return interp;
}

double InverseInterpolation(double y, int r){

    double x1, y1, x2, y2, slope, invinterp;

    for (int i = 0; i < r-1; i++){
        if (InterpolationList[i][1] > y && y > InterpolationList[i+1][1]){
            // Set the x and y values of the two points
            x1 = InterpolationList[i][0];
            x2 = InterpolationList[i+1][0];
            y1 = InterpolationList[i][1];
            y2 = InterpolationList[i+1][1];
        }
    }
    							    
    // Calculate the slope of the line between the two points
    slope = (x2-x1)/(y2-y1);   
    // Calculate the y value for the given x value
    invinterp = slope*(y-y1)+x1;

    return invinterp;
}



//***********************************************
// proton and neutron form factors
//***********************************************

// dipole form factor
double GD(double q2){double ff;
  
  ff = 1/(1+q2/0.71)/(1+q2/0.71);
  
  return ff;}

  
// magnetic form factor of the proton (from 0812.3539)
double GMp(double q2){double ff;
  
  ff = mup*(1+apM1*q2/4/Mproton/Mproton)/(1+bpM1*q2/4/Mproton/Mproton
             +bpM2*q2*q2/4/4/Mproton/Mproton/Mproton/Mproton
	     +bpM3*q2*q2*q2/4/4/4/Mproton/Mproton/Mproton/Mproton/Mproton/Mproton);
  
  return ff;}

  
// electric form factor of the proton (from 0812.3539)
double GEp(double q2){double ff;
  
  ff = (1+apE1*q2/4/Mproton/Mproton)/(1+bpE1*q2/4/Mproton/Mproton
             +bpE2*q2*q2/4/4/Mproton/Mproton/Mproton/Mproton
	     +bpE3*q2*q2*q2/4/4/4/Mproton/Mproton/Mproton/Mproton/Mproton/Mproton);
  
  return ff;}
  
  
// magnetic form factor of the neutron (from 0812.3539)
double GMn(double q2){double ff;
  
  ff = mun*(1+anM1*q2/4/Mneutron/Mneutron)/(1+bnM1*q2/4/Mneutron/Mneutron
             +bnM2*q2*q2/4/4/Mneutron/Mneutron/Mneutron/Mneutron
	     +bnM3*q2*q2*q2/4/4/4/Mneutron/Mneutron/Mneutron/Mneutron/Mneutron/Mneutron);
  
  return ff;}
  
  
// electric form factor of the neutron (from 0812.3539)
double GEn(double q2){double ff;
  
  ff = ffA*q2/4/Mneutron/Mneutron/(1+ffB*q2/4/Mneutron/Mneutron)*GD(q2);
  
  return ff;}


  
//***********************************************
// Pauli blocking factor for nucleons scattering inside the nucleus (from Lovseth/Radomski)
//***********************************************

double Pauli(double p){double blockingfactor;
  
  if (p > 2*0.235){blockingfactor = 1;}
  else {blockingfactor = 1.5*p/2/0.235 - 0.5*p*p*p/8/0.235/0.235/0.235;}
  
  return blockingfactor;}
  

//***********************************************
// define the trident processes
//***********************************************

void SetTridentProcess(){

     if(process == 1){ // nu_e -> nu_e e+ e-
         m3 = me; m4 = me;
         GVSM = 0.5 + 2*sW2; GASM = -0.5;
	     anti = 0;
         PDG1 = 12; PDG2 = 12; PDG3 = -11; PDG4 = 11;}
    
    else if(process == 2){ // nu_e -> nu_e mu+ mu-
         m3 = mmu; m4 = mmu;
         GVSM = -0.5 + 2*sW2; GASM = 0.5;
	     anti = 0;
         PDG1 = 12; PDG2 = 12; PDG3 = -13; PDG4 = 13;}
         
    else if(process == 3){ // nu_e -> nu_mu mu+ e-
         m3 = mmu; m4 = me; 
         GVSM = 1.0; GASM = -1.0;
         anti = 0; 
         PDG1 = 12; PDG2 = 14; PDG3 = -13; PDG4 = 11;}
         
    else if(process == 4){ // anti-nu_e -> anti-nu_e e+ e-
         m3 = me; m4 = me; 
         GVSM = 0.5 + 2*sW2; GASM = -0.5;
         anti = 1;
         PDG1 = -12; PDG2 = -12; PDG3 = -11; PDG4 = 11;}
         
    else if(process == 5){ // anti-nu_e -> anti-nu_e mu+ mu-
         m3 = mmu; m4 = mmu; 
         GVSM = -0.5 + 2*sW2; GASM = 0.5;
         anti = 1;
         PDG1 = -12; PDG2 = -12; PDG3 = -13; PDG4 = 13;}
    
    else if(process == 6){ // anti-nu_e -> anti-nu_mu e+ mu-
         m3 = me; m4 = mmu; 
         GVSM = 1.0; GASM = -1.0;
         anti = 1;
         PDG1 = -12; PDG2 = -14; PDG3 = -11; PDG4 = 13;}
    
    else if(process == 7){ // nu_mu -> nu_mu e+ e-
         m3 = me; m4 = me;
         GVSM = -0.5 + 2*sW2; GASM = 0.5;
	     anti = 0;
         PDG1 = 14; PDG2 = 14; PDG3 = -11; PDG4 = 11;}
         
    else if(process == 8){ // nu_mu -> nu_mu mu+ mu-
         m3 = mmu; m4 = mmu;
         GVSM = 0.5 + 2*sW2; GASM = -0.5;
         anti = 0;  
         PDG1 = 14; PDG2 = 14; PDG3 = -13; PDG4 = 13;}
         
    else if(process == 9){ // nu_mu -> nu_e e+ mu-
         m3 = me; m4 = mmu; 
         GVSM = 1.0; GASM = -1.0;
         anti = 0;
         PDG1 = 14; PDG2 = 12; PDG3 = -11; PDG4 = 13;}
         
    else if(process == 10){ // anti-nu_mu -> anti-nu_mu e+ e-
         m3 = me; m4 = me; 
         GVSM = -0.5 + 2*sW2; GASM = 0.5;
         anti = 1;  
         PDG1 = -14; PDG2 = -14; PDG3 = -11; PDG4 = 11;}

    else if(process == 11){ // anti-nu_mu -> anti-nu_mu mu+ mu-
         m3 = mmu; m4 = mmu; 
         GVSM = 0.5 + 2*sW2; GASM = -0.5;
         anti = 1;
         PDG1 = -14; PDG2 = -14; PDG3 = -13; PDG4 = 13;}

    else if(process == 12){ // anti-nu_mu -> anti-nu_e mu+ e-
         m3 = mmu; m4 = me; 
         GVSM = 1.0; GASM = -1.0; 
         anti = 1;
         PDG1 = -14; PDG2 = -12; PDG3 = -13; PDG4 = 11;}

    return;
  
}


//***********************************************
// set the nuclear parameters
//***********************************************

void SetNuclearParameters(){
            
    if (material.compare("Ar") == 0){
        A = 40;
        Z = 18;
        M = MArgon;
	    M2 = M*M;
	    for (int jj = 0; jj < 192; jj++){
            for (int kk = 0; kk < 2; kk++){
                InterpolationList[jj][kk] = InterpolationListAr[jj][kk];}}}
            
    else if (material.compare("Fe") == 0){
        A = 56;
        Z = 26;
        M = MIron;
	    M2 = M*M;
	    for (int jj = 0; jj < 192; jj++){
            for (int kk = 0; kk < 2; kk++){
                InterpolationList[jj][kk] = InterpolationListFe[jj][kk];}}}
            
    else if (material.compare("proton") == 0){
        A = 1;
        Z = 1;
        M = Mproton;
	    M2 = M*M;
	    for (int jj = 0; jj < 192; jj++){
            for (int kk = 0; kk < 2; kk++){
                InterpolationList[jj][kk] = InterpolationListProton[jj][kk];}}}
            
    else if (material.compare("neutron") == 0){
        A = 1;
        Z = 1; // "charge" of the neutron is taken care by the form factor
        M = Mneutron;
	    M2 = M*M;
	    for (int jj = 0; jj < 192; jj++){
            for (int kk = 0; kk < 2; kk++){
                InterpolationList[jj][kk] = InterpolationListNeutron[jj][kk];}}}
            
    return;

}
    
    
//*******************************************************
// main program
//*******************************************************

int main(){
    
        // Define the trident process
        std::cout << "\n";
	    std::cout << "Select the trident process  [enter 1 - 12]\n\n";
	    std::cout << " [1] nu_e -> nu_e e+ e-                  [7] nu_mu -> nu_mu e+ e- \n";
	    std::cout << " [2] nu_e -> nu_e mu+ mu-                [8] nu_mu -> nu_mu mu+ mu- \n";
	    std::cout << " [3] nu_e -> nu_mu mu+ e-                [9] nu_mu -> nu_e e+ mu- \n";
	    std::cout << " [4] anti-nu_e -> anti-nu_e e+ e-        [10] anti-nu_mu -> anti-nu_mu e+ e- \n";
	    std::cout << " [5] anti-nu_e -> anti-nu_e mu+ mu-      [11] anti-nu_mu -> anti-nu_mu mu+ mu- \n";
	    std::cout << " [6] anti-nu_e -> anti-nu_mu e+ mu-      [12] anti-nu_mu -> anti-nu_e mu+ e- \n\n";
        std::cin >> process;
	    if(process != 1 && process != 2 && process != 3 && process != 4 && process != 5 && 
	       process != 6 && process != 7 && process != 8 && process != 9 && process != 10 && 
	       process != 11 && process != 12){
        std::cout << "\n Invalid selection \n";
	    return 0;}
	   
	SetTridentProcess();
 	
 	
	// Define the target material
	    std::cout << "\n";
        std::cout << "Select the target material  [enter Ar, Fe, proton, neutron] \n\n";
	    std::cout << " [Ar] Argon      [Fe] Iron \n\n";
	    std::cout << " [proton]  Proton inside an ideal Fermi gas  \n";  
	    std::cout << " [neutron] Neutron inside an ideal Fermi gas \n\n";
        std::cin >> material;
	    if(material.compare("Ar") != 0 && material.compare("Fe") != 0 &&
	       material.compare("proton") != 0 && material.compare("neutron") != 0 ){
	    std::cout << "\n Invalid selection \n";
	    return 0;}

	SetNuclearParameters();
	
	
	// Define the incoming neutrino spectrum
	    std::cout << "\n";
        std::cout << "Are you using a fixed neutrino energy or an energy distribution? \n\n";
	    std::cout << "[1] fixed neutrino energy \n\n";
        std::cout << "[2] Use the uBooNE BNB flux \n\n";
        std::cout << "[3] Load flux from a ROOT file \n\n";
        std::cin >> energy_type;
	    if(energy_type.compare("1") != 0 && energy_type.compare("2") != 0 ){
	    std::cout << "\n Invalid selection \n";
	    return 0;}

        // Use the selected process to determine neutrino flavor and particle/antiparticle
        

        if(energy_type.compare("1") == 0){
	    Enu = 0.0;
	    std::cout << "\n";
        std::cout << "Enter the energy of the neutrino beam (in GeV) \n\n";
        std::cin >> Enu;}

        if(energy_type.compare("2") == 0){
            string is_antinu;
            string flavor;
            flux_file = "/exp/uboone/app/users/jbateman/workdir/DarkNews/Trident/data/flux/bnb/MCC9_FluxHist_volTPCActive_w2D_hists.root";
            
            LoadFluxFromROOT(flux_file, PDG1);}
        
        if(energy_type.compare("3") == 0){
            string is_antinu;
            string flavor;
            std::cout << "\n";
            std::cout << "Enter the path to the ROOT flux file \n\n";
            std::cin >> flux_file;
            
            LoadFluxFromROOT(flux_file, PDG1);}
        
        
        // Define the particle physics model
        std::cout << "\n";
        std::cout << "Select the model \n\n";
	    std::cout << "[4F]      Model independent 4 Fermi operators \n";
	    std::cout << "[SM]      Standard Model \n";
	    std::cout << "[LmuLtau] Standard Model + Z' gauge boson based on L_mu - L_tau \n\n";
        std::cin >> model;
    	if(model.compare("SM") != 0 && model.compare("4F") != 0 && model.compare("LmuLtau") != 0 ){
	    std::cout << "\n Invalid selection \n";
	    return 0;}
	   	
        if(model.compare("SM") == 0){
        GAtot=0.0; GVtot=0.0; GP=0.0; MZP=1000.0;}
	
	    if(model.compare("4F") == 0){
        GP=0.0; MZP=1000.0;
	    std::cout << "\n";
        std::cout << "Enter the value of the vector coupling g_V \n\n";
        std::cin >> GVtot; 
	    std::cout << "\n";
        std::cout << "Enter the value of the axial-vector coupling g_A \n\n";
        std::cin >> GAtot;}
        
        if(model.compare("LmuLtau") == 0){
	    GAtot=0.0; GVtot=0.0;
        std::cout << "\n";
        std::cout << "Enter the value of the L_mu - L_tau gauge coupling g' \n\n";
        std::cin >> GP; 
	    std::cout << "\n";
        std::cout << "Enter the value of the Z' mass (in GeV) \n\n";
        std::cin >> MZP;}
        
        
        // Cross section or event generation
        std::cout << "\n";
        std::cout << "You can compute the trident [CrossSection] or [GenerateEvents] \n\n";
        std::cin >> command;
	    if(command.compare("CrossSection") != 0 && command.compare("GenerateEvents") != 0 ){
	    std::cout << "\n Invalid selection \n";
	    return 0;}
        
        
        if(command.compare("GenerateEvents") == 0){
	    std::cout << "\n";
        std::cout << "Enter the number of events to be generated \n\n";
        std::cin >> Nevents; 
	    ::data = new double[28*Nevents];
	    std::cout << "\n";
        std::cout << "Enter the name of the output file \n\n";
        std::cin >> filename_out;
        std::cout << "Choose the output format: \n\n";
        std::cout << "[1] Original TEG format \n";
        std::cout << "[2] Hepevt format \n";
        std::cout << "[3] HepMC3 format \n";

        
        std::cin >> output_format;
        if(output_format != 1 && output_format != 2 && output_format != 3 ){
        std::cout << "\n Invalid selection \n";
        return 0;
        }
    }
        
        
        // Compute cross section
        if(command.compare("CrossSection") == 0){
	    ComputeCrossSection();
        std::cout << "\n\n";
	    std::cout << "The trident cross section is  ( " << crosssectionresult << " +- " << deltacrosssectionresult << " ) fb  \n";
        std::cout << "(uncertainty is the statistical uncertainty of the numerical phase space integration) \n\n";
        std::cout << "A temporary file with the cross section values has been create/added to: temp_xsec_file.txt \nPlease clear this file between different processes/runs. \n\n";
        WriteXsecTempFile("temp_xsec_file.txt");}
        
        // Generate Events
        if(command.compare("GenerateEvents") == 0){
        FindMaxWeight();
	    GenerateEvents();

        if (output_format == 1){
            filename_out.append(".teg");
            WriteEventFile(filename_out.c_str());}
        else if (output_format == 3){
            filename_out.append(".hepmc3");
            WriteEventHepMC3(filename_out.c_str());}
        else if (output_format == 2){
            filename_out.append(".hepevt");
            WriteEventHepevt(filename_out.c_str());}        
    }
	   
	      
	return 0;
}

//********************************************************************
// Generate a random point in phase space
//********************************************************************

void GenerateEvent(bool get_vertex){
	
    if (energy_type.compare("1") == 0){
      eps1 = Enu;
    }
    else if (energy_type.compare("1") != 0){  
      bin = realdistribution(generator)*length_probability_list;
      eps1 = realdistribution(generator)*(Enumax_list[bin]-Enumin_list[bin])+Enumin_list[bin];       
    }

    // Generate a random point if the neutrino energy is high enough.
    // Otherwise set flag to return a zero weight event.
   
    if ((1-(m3+m4)*(m3+m4)/2/(eps1*eps1)*(1+eps1/M))*(1-(m3+m4)*(m3+m4)/2/(eps1*eps1)*(1+eps1/M)) - (m3+m4)*(m3+m4)*(m3+m4)*(m3+m4)/4/(eps1*eps1*eps1*eps1)*(1+2*eps1/M) < 0) {
        zeroweight=1;
	for (int jj = 0; jj < 4; jj++){
        for (int kk = 0; kk < 4; kk++){
        event[jj][kk] = 0.0;}}
        //std::cout << "The neutrino energy is too small to create the charged leptons";
        return;}
        
    else {zeroweight=0; GenerateRandomPoint();}
    
    vector<double> prod(3);

    // Rotate around z-axis by random angle BEFORE rotating neutrino direction
    Theta = Twopi*realdistribution(generator);
    sintheta = sin(Theta);
    costheta = cos(Theta);

    double pnuoutx_updated = pnuoutx*costheta+pnuouty*sintheta;
    double pnuouty_updated = -pnuoutx*sintheta+pnuouty*costheta;
    pnuoutx = pnuoutx_updated;
    pnuouty = pnuouty_updated;
    

    double plminusx_updated = plminusx*costheta+plminusy*sintheta;
    double plminusy_updated = -plminusx*sintheta+plminusy*costheta;
    plminusx = plminusx_updated;
    plminusy = plminusy_updated;
    
    double plplusx_updated = plplusx*costheta+plplusy*sintheta;
    double plplusy_updated = -plplusx*sintheta+plplusy*costheta;
    plplusx = plplusx_updated;
    plplusy = plplusy_updated;

    double Pprimex_updated = Pprimex*costheta+Pprimey*sintheta;
    double Pprimey_updated = -Pprimex*sintheta+Pprimey*costheta;
    Pprimex = Pprimex_updated;
    Pprimey = Pprimey_updated;

    if (get_vertex){
        // Sample incoming neutrino energy vs production vertex z from flux histogram
        prod = SampleEvsZhist(Enuin, flux_file);
        
        // Generate a random vertex in the active volume
        vtx_x = xmin + (xmax - xmin) * realdistribution(generator);
        vtx_y = ymin + (ymax - ymin)  * realdistribution(generator);
        vtx_z = zmin + (zmax - zmin) * realdistribution(generator);

        // // Temporarily, use the cryostat volume instead, ignoring end caps
        // double r_t = 191.61; // cm
        // double z_t = 1086.49; // cm

        // double r_polar = r_t * sqrt(realdistribution(generator));;
        // double phi = Twopi * realdistribution(generator);
        // vtx_x = r_polar * cos(phi) + (xmax - xmin)/2; // uB coordinate origin is at the at the collection plane in x (0, 256.35 cm)
        // vtx_y = r_polar * sin(phi); // uB coordinate origin is at the center of the cryostat in y (-116.3 cm to 116.3 cm)
        // vtx_z = z_t * realdistribution(generator); // uB coordinate origin is at the front face of the cryostat in z (0 to 1036.8 cm)
        
        vtx_t = MicroBooNEGlobalTimeOffset + MicroBooNERandomTimeOffset * realdistribution(generator);
    
        // Get the initial neutrino direction from production point to interaction vertex
        // Then rotate all final state momenta into this frame

        vector<double> dir(3);
        dir[0] = vtx_x - prod[0];
        dir[1] = vtx_y - prod[1];
        dir[2] = vtx_z - prod[2];

        double norm = sqrt(dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]);
        dir = {dir[0]/norm, dir[1]/norm, dir[2]/norm};

        // ***********************************
        // Rotate momenta into neutrino direction frame
        // ***********************************

        vector<double> V = {0.0, 0.0, 1.0};

        vector<double> u = {
            V[1]*dir[2] - V[2]*dir[1],
            V[2]*dir[0] - V[0]*dir[2],
            V[0]*dir[1] - V[1]*dir[0]
        };

        double normFactor = sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
        double cos_phi = V[0]*dir[0] + V[1]*dir[1] + V[2]*dir[2];
        double sin_phi = normFactor;  // Use normFactor directly as sin_phi
        
        if (normFactor < 1.0e-10) {
            // Direction is already aligned with z-axis, no rotation needed
            normFactor = 1.0;
            sin_phi = 0.0;
        }
        
        u = {u[0]/normFactor, u[1]/normFactor, u[2]/normFactor};

        // Construct rodrigues rotation matrix
        double Su[3][3] = { {0, -u[2], u[1]},
                            {u[2], 0, -u[0]},
                            {-u[1], u[0], 0}};

        double uuT[3][3] = {{u[0]*u[0], u[0]*u[1], u[0]*u[2]},
                            {u[1]*u[0], u[1]*u[1], u[1]*u[2]},
                            {u[2]*u[0], u[2]*u[1], u[2]*u[2]}};

        double I_minus_uuT[3][3] = {{1 - uuT[0][0], -uuT[0][1], -uuT[0][2]},
                                    {-uuT[1][0], 1 - uuT[1][1], -uuT[1][2]},
                                    {-uuT[2][0], -uuT[2][1], 1 - uuT[2][2]}};
        

        double R[3][3] = {
            {uuT[0][0] + I_minus_uuT[0][0]* cos_phi + Su[0][0]* sin_phi, uuT[0][1] + I_minus_uuT[0][1]* cos_phi + Su[0][1]* sin_phi, uuT[0][2] + I_minus_uuT[0][2]* cos_phi + Su[0][2]* sin_phi},
            {uuT[1][0] + I_minus_uuT[1][0]* cos_phi + Su[1][0]* sin_phi, uuT[1][1] + I_minus_uuT[1][1]* cos_phi + Su[1][1]* sin_phi, uuT[1][2] + I_minus_uuT[1][2]* cos_phi + Su[1][2]* sin_phi},
            {uuT[2][0] + I_minus_uuT[2][0]* cos_phi + Su[2][0]* sin_phi, uuT[2][1] + I_minus_uuT[2][1]* cos_phi + Su[2][1]* sin_phi, uuT[2][2] + I_minus_uuT[2][2]* cos_phi + Su[2][2]* sin_phi}
        };      
        
        double pnuinx_new = R[0][0]*pnuinx + R[0][1]*pnuiny + R[0][2]*pnuinz;
        double pnuiny_new = R[1][0]*pnuinx + R[1][1]*pnuiny + R[1][2]*pnuinz;
        double pnuinz_new = R[2][0]*pnuinx + R[2][1]*pnuiny + R[2][2]*pnuinz;

        pnuinx = pnuinx_new;
        pnuiny = pnuiny_new;
        pnuinz = pnuinz_new;

        double pnuoutx_new = R[0][0]*pnuoutx + R[0][1]*pnuouty + R[0][2]*pnuoutz;
        double pnuouty_new = R[1][0]*pnuoutx + R[1][1]*pnuouty + R[1][2]*pnuoutz;
        double pnuoutz_new = R[2][0]*pnuoutx + R[2][1]*pnuouty + R[2][2]*pnuoutz;

        pnuoutx = pnuoutx_new;
        pnuouty = pnuouty_new;
        pnuoutz = pnuoutz_new;

        double plminusx_new = R[0][0]*plminusx + R[0][1]*plminusy + R[0][2]*plminusz;
        double plminusy_new = R[1][0]*plminusx + R[1][1]*plminusy + R[1][2]*plminusz;
        double plminusz_new = R[2][0]*plminusx + R[2][1]*plminusy + R[2][2]*plminusz;
        
        plminusx = plminusx_new;
        plminusy = plminusy_new;
        plminusz = plminusz_new;

        double plplusx_new = R[0][0]*plplusx + R[0][1]*plplusy + R[0][2]*plplusz;
        double plplusy_new = R[1][0]*plplusx + R[1][1]*plplusy + R[1][2]*plplusz;
        double plplusz_new = R[2][0]*plplusx + R[2][1]*plplusy + R[2][2]*plplusz;

        plplusx = plplusx_new;
        plplusy = plplusy_new;
        plplusz = plplusz_new;
        
        double Pprimex_new = R[0][0]*Pprimex + R[0][1]*Pprimey + R[0][2]*Pprimez;
        double Pprimey_new = R[1][0]*Pprimex + R[1][1]*Pprimey + R[1][2]*Pprimez;
        double Pprimez_new = R[2][0]*Pprimex + R[2][1]*Pprimey + R[2][2]*Pprimez;

        Pprimex = Pprimex_new;
        Pprimey = Pprimey_new;
        Pprimez = Pprimez_new;
    }
    
    // Define the event
    event[0][0] = 0.;
    event[0][1] = prod[0];
    event[0][2] = prod[1];
    event[0][3] = prod[2];
    
    event[1][0] = Enuin;
    event[1][1] = pnuinx;
    event[1][2] = pnuiny;
    event[1][3] = pnuinz;

    event[2][0] = vtx_t;
    event[2][1] = vtx_x;
    event[2][2] = vtx_y;
    event[2][3] = vtx_z;
    
    event[3][0] = Enuout;
    event[3][1] = pnuoutx;
    event[3][2] = pnuouty;
    event[3][3] = pnuoutz;
    
    event[4][0] = Elminus;
    event[4][1] = plminusx;
    event[4][2] = plminusy;
    event[4][3] = plminusz;
    
    event[5][0] = Elplus;
    event[5][1] = plplusx;
    event[5][2] = plplusy;
    event[5][3] = plplusz;

    event[6][0] = Eprime;
    event[6][1] = Pprimex;
    event[6][2] = Pprimey;
    event[6][3] = Pprimez;
    
    return;
}

//*************************************************************************

void GenerateRandomPoint(){
 
    g = eps1/M;
    x0 = M*eps1; // (B.1) of Lovseth, Radomski
    
    x1max = 2*eps1*eps1/(1+2*g)*(1-(m3+m4)*(m3+m4)/2/(eps1*eps1)*(1+g)+sqrt((1-(m3+m4)*(m3+m4)/2/(eps1*eps1)*(1+g))*(1-(m3+m4)*(m3+m4)/2/(eps1*eps1)*(1+g))-(m3+m4)*(m3+m4)*(m3+m4)*(m3+m4)/4/(eps1*eps1*eps1*eps1)*(1+2*g)));
    x1min = (m3+m4)*(m3+m4)*(m3+m4)*(m3+m4)/x1max/(1+2*g); // (B.4) of Lovseth, Radomski
    
    u1max = Interpolation(x1min, 192);
    // put hard cut off for the form factor
    // at some point the interpolation function becomes ureliable
    if (x1max > 1.0){u1min = Interpolation(1.0, 192);}
    else{u1min = Interpolation(x1max, 192);}
    if (x1min > 1.0){u1min = Interpolation(x1min, 192);}

    u1 = u1min + (u1max-u1min)*realdistribution(generator);
    x1 = InverseInterpolation(u1, 192);
    
    
    u2max = (1+2*g)*(x1max-x1)*(x1-x1min)/((m3+m4)*(m3+m4)+x1*(1+g)+2*eps1*sqrt((x1+(x1*x1)/4/M2)));
    u2min = 0.0; // (B.5) of Lovseth, Radomski
    
    u2 = u2min + (u2max-u2min)*realdistribution(generator);
    x2 = (u2+x1+(m3+m4)*(m3+m4))/2;
    
    
    WB = sqrt(u2+(m3+m4)*(m3+m4));
    WB2 = WB*WB;
    yplus = (WB2-m3*m3-m4*m4)+sqrt(WB2*WB2+(m3*m3-m4*m4)*(m3*m3-m4*m4)-2*m3*m3*WB2-2*m4*m4*WB2);
    yminus = 4*m3*m3*m4*m4/yplus; // (B.9b) of Lovseth, Radomski
    
    D3max = m3*m3+x1*m4*m4/(WB2)+x2*yplus/(WB2);
    D3min = m3*m3+x1*m4*m4/(WB2)+x2*yminus/(WB2); // (B.9b) of Lovseth, Radomski
    
    u3max = log(D3max);
    u3min = log(D3min); // (B.9a) of Lovseth, Radomski
    u3 = u3min + (u3max-u3min)*realdistribution(generator);

    x3 = (exp(u3)-x1)/2; // (B.8) of Lovseth, Radomski
    D3 = x1+2*x3;
    D32 = D3*D3;
    
    u4max = ((D3max-D3)*(D3-D3min)*WB2/x1/(u2+2*m3*(m3+m4)+2*x3*(x2-x1)/x1+2*x2*sqrt(x3*x3+m3*m3*x1)/x1))*((D3max-D3)*(D3-D3min)*WB2/x1/(u2+2*m3*(m3+m4)+2*x3*(x2-x1)/x1+2*x2*sqrt(x3*x3+m3*m3*x1)/x1));
    //  (B.10) of Lovseth, Radomski (square is according to Brown et al.)
    u4min = 0.0;
    u4 = u4min + (u4max-u4min)*realdistribution(generator);

    x5 = x2-x3-(m4*m4-m3*m3+x1+sqrt(u4))/2; // (B.8) of Lovseth, Radomski
    
    u5max = Twopi; // (B.11) of Lovseth, Radomski
    u5min = 0.0;
    u5 = u5min + (u5max-u5min)*realdistribution(generator);

    phi3 = u5;
    
    P2 = M2*(u2max-u2)/2*(x2+eps1*x1/2/M+eps1*sqrt(x1+x1*x1/4/M2))/(x2*x2); // (B.13) of Lovseth, Radomski
    P5 = x1/2*(sqrt(u4max)-sqrt(u4))*((sqrt(u4max)-sqrt(u4))/2+2*x2*sqrt(x3*x3+m3*m3*x1)/x1)/(x2*x2);
    
    WC = sqrt(sqrt(u4)+m4*m4); // from (24) of Lovseth, Radomski
    WC2 = WC*WC;
    qc = 1/WC*sqrt((x2-x3)*(x2-x3)-2*x1*x5+m3*m3*x1); // (B.16) of Lovseth, Radomski
    WCqc = sqrt((x2-x3)*(x2-x3)-2*x1*x5+m3*m3*x1);
    
    x7 = (x0*x1*x5+x0*x2*x3-x1*x2*x5/2)/(x2*x2)-sqrt(P2*P5)*cos(u5); // (B.12) of Lovseth, Radomski
    
    u7max = Twopi; // (B.17) of Lovseth, Radomski
    u7min = 0.0;
    u7 = u7min + (u7max-u7min)*realdistribution(generator);

    phi4 = u7;
    
    D4max = m4*m4*(WC2+2*x5-m3*m3)/WC2+(1-m4*m4/WC2)/2*(WC2+x1+2*x5-m3*m3+sqrt((WC2+2*x5+x1-m3*m3)*(WC2+2*x5+x1-m3*m3)-4*x1*(2*x5-m3*m3))); // (B.18) of Lovseth, Radomski
    D4min = m4*m4*(WC2+2*x5-m3*m3)/WC2+2*(1-m4*m4/WC2)*x1*(2*x5-m3*m3)/((WC2+x1+2*x5-m3*m3)+sqrt((WC2+x1+2*x5-m3*m3)*(WC2+x1+2*x5-m3*m3)-4*x1*(2*x5-m3*m3))); // (B.19) of Lovseth, Radomski

    u6max = log(D4max); // (B.17) of Lovseth, Radomski
    u6min = log(D4min);
    u6 = u6min + (u6max-u6min)*realdistribution(generator);

    
    x4 = (-x1+exp(u6))/2; // (B.15) of Lovseth, Radomski
    D4 = x1+2*x4;
    D42 = D4*D4;
    D3D4 = D3*D4;
    
    P4 = (D4max-D4)*(D4-D4min)/4/(qc*qc); // below (B.21) of Lovseth, Radomski
    
    x6 = ((WC2+m4*m4)/2*(x2*(x2-x3)-x1*x5)+x4*(x2-x5)*(x2-x1-x3)-x4*x2*WC2)/WCqc/WCqc-x2*sqrt(P4*P5)*cos(phi4)/WCqc; // (B.21) of Lovseth, Radomski
    
    // x8 result in (B.21) of Lovseth, Radomski is missing a term according to Brown et al. and is not used in the following
    
    // below (B.20) in Lovseth, Radomski
    P0C = (x0-0.5*x1-x7)/WC; // note the typo ...x7/WC) in Lovseth, Radomski
    PC = sqrt(P0C*P0C-M2);
    q0C = (x2-x1-x3)/WC;
    qC = sqrt(q0C*q0C+x1);
    eps1C = (x2-x5)/WC;
    eps4C = (WC2+m4*m4)/2/WC;
    p4C = sqrt(u4)/2/WC;
    
    cos4 = (eps4C*q0C-x4)/p4C/qC;
    sin4 = sqrt(1-cos4*cos4);
    
    flagkinematics = 0;
      
    if(cos4*cos4>1.01){
      flagkinematics = 1; 
//      std::cout << "Warning: event with extreme kinematics. Numerical precision is insufficient. \n";
//      std::cout << "cos4^2-1 = " << (eps4C*q0C-x4)/p4C/qC*(eps4C*q0C-x4)/p4C/qC-1 << "\n\n";
    }
    
    if(1-cos4*cos4<1.0e-16){
      sin4=0.0; // adhoc fix
      cos4=1.0; // adhoc fix
    }
      
    cosPq = (P0C*q0C+x1/2)/PC/qC;
    sinPq = sqrt(1-cosPq*cosPq);
    
    if(1-cosPq*cosPq<1.0e-16){
      flagkinematics = 1;
//      std::cout << "Warning: event with extreme kinematics. Numerical precision is insufficient. \n";
//      std::cout << "cosPq^2-1 = " << (P0C*q0C+x1/2)/PC/qC*(P0C*q0C+x1/2)/PC/qC-1 << "\n\n";
    }
    
    cosP1q = (eps1C*q0C-x2)/eps1C/qC;
    sinP1q = sqrt(1-cosP1q*cosP1q);
    
    if(1-cosP1q*cosP1q<1.0e-16){
      flagkinematics = 1;
//      std::cout << "Warning: event with extreme kinematics. Numerical precision is insufficient. \n";
//      std::cout << "cosP1q^2-1 = " << (eps1C*q0C-x2)/eps1C/qC*(eps1C*q0C-x2)/eps1C/qC-1 << "\n\n";
    }
   
    cosPp1 = (P0C*eps1C-x0)/eps1C/PC;
    sinPp1 = sqrt(1-cosPp1*cosPp1);
    
    if(1-cosPp1*cosPp1<1.0e-16){
      flagkinematics = 1;
//      std::cout << "Warning: event with extreme kinematics. Numerical precision is insufficient. \n";
//      std::cout << "cosPp1^2-1 = " << (P0C*eps1C-x0)/eps1C/PC*(P0C*eps1C-x0)/eps1C/PC-1 << "\n\n";
    }
   
    cosPhi = (cosPp1-cosPq*cosP1q)/(sinPq*sinP1q);
    sinPhi = sqrt(1-cosPhi*cosPhi);
    
    if(cosPhi>1.01){
      flagkinematics = 1;
//      std::cout << "Warning: event with extreme kinematics. Numerical precision is insufficient. \n";
//      std::cout << "cosPhi^2-1 = " << (cosPp1-cosPq*cosP1q)/(sinPq*sinP1q)*(cosPp1-cosPq*cosP1q)/(sinPq*sinP1q)-1 << "\n\n";
    }
      if(1-cosPhi*cosPhi<1.0e-16){
      sinPhi=0.0; // adhoc fix
      cosPhi=1.0; // adhoc fix
    }

    x8 = P0C*eps4C-PC*p4C*(cos4*cosPq+sin4*sinPq*(cos(phi4)*cosPhi-sin(phi4)*sinPhi)); // (B .20) in Lovseth, Radomski

    // (B.1) and (B.2) of Lovseth, Radomski
    p1P = -x0;
    q2 = x1;
    qp1 = -x2;
    qp3 = -x3;
    qp4 = -x4;
    p1p3 = -x5;
    p1p4 = -x6;
    p3P = -x7;
    p4P = -x8;
    Pq = x1*0.5;
    qP = Pq;
    p2P = x1*0.5-x0+x7+x8;
    qp2 = x1-x2+x3+x4;
    p1p2 = -x2+x5+x6;
    p2p3 = (m3*m3-m4*m4+x1)*0.5-x2+x4+x6;
    p2p4 = (m4*m4-m3*m3+x1)*0.5-x2+x3+x5;
    p3p4 = (m3*m3+m4*m4-x1)*0.5+x2-x3-x4-x5-x6;
    
   
    // Define the 4 momenta of the events
    
    Enuin = eps1;
    pnuinx = 0.0;
    pnuiny = 0.0;
    pnuinz = eps1;
    
    // note that signs are such that the metric of Lovseth, Radomski is taken into account
    
    Eq = -qP/M;
    qz = Eq+qp1/eps1;
    qy = sqrt(q2-2*Eq*qp1/eps1-qp1*qp1/eps1/eps1);
    qx = 0.0;
    
    Elplus = -p3P/M;
    plplusz = Elplus+p1p3/eps1; 
    plplusy = (Eq*Elplus+qp3-qz*plplusz)/qy;
    plplusx = sqrt(-2*Elplus*p1p3/eps1-p1p3*p1p3/eps1/eps1-plplusy*plplusy-m3*m3);

    Enuout = -p2P/M; 
    pnuoutz = Enuout+p1p2/eps1;
    pnuouty = (Eq*Enuout+qp2-qz*pnuoutz)/qy;
    pnuoutx = 1/plplusx*(p2p3-p1p2*p1p3/eps1/eps1-p1p2*Elplus/eps1-p1p3*Enuout/eps1-pnuouty*plplusy);
    
    Elminus = -p4P/M;
    plminusz = Elminus+p1p4/eps1;
    plminusy = (Eq*Elminus+qp4-qz*plminusz)/qy;
    plminusx = -plplusx-pnuoutx;
        
    // special case where the anti-muon has no momentum in x direction
    if (-2*Elplus*p1p3/eps1-p1p3*p1p3/eps1/eps1-plplusy*plplusy < 0){
       flagkinematics = 1;
//       std::cout << "Warning: event with extreme kinematics. Numerical precision is insufficient. \n";
//       std::cout << "Elplus = " << Elplus << "\n";
//       std::cout << "plplusz = " << plplusz << "\n";
//       std::cout << "plplusy = " << plplusy << "\n";
//       std::cout << "plplusx = " << plplusx << "\n\n";
    }
    if (-2*Elplus*p1p3/eps1-p1p3*p1p3/eps1/eps1-plplusy*plplusy-m3*m3 < 0){
        plplusx = 0.0;
	    plminusx = sqrt(-2*Elminus*p1p4/eps1-p1p4*p1p4/eps1/eps1-plminusy*plminusy-m4*m4);
	    pnuoutx = -plminusx;
    }
    
    // special case where both muon and anti-muon have no momentum in x direction
    if (-2*Elplus*p1p3/eps1-p1p3*p1p3/eps1/eps1-plplusy*plplusy-m3*m3 < 0 && -2*Elminus*p1p4/eps1-p1p4*p1p4/eps1/eps1-plminusy*plminusy < 0){
       flagkinematics = 1;
//       std::cout << "Warning: event with extreme kinematics. Numerical precision is insufficient. \n";
//       std::cout << "Elminus = " << Elminus << "\n";
//       std::cout << "plminusz = " << plminusz << "\n";
//       std::cout << "plminusy = " << plminusy << "\n";
//       std::cout << "plminusx = " << plminusx << "\n\n";
    }
    if (-2*Elplus*p1p3/eps1-p1p3*p1p3/eps1/eps1-plplusy*plplusy-m3*m3 < 0 && -2*Elminus*p1p4/eps1-p1p4*p1p4/eps1/eps1-plminusy*plminusy-m4*m4 < 0){
        plplusx = 0.0;
	    plminusx = 0.0;
	    pnuoutx = 0.0;
    }
   
    
    // special case where the momentum transfer to the nucleon has no y component    
    if (q2-2*Eq*qp1/eps1-qp1*qp1/eps1/eps1 < -q2){
       flagkinematics = 1;
//       std::cout << "Warning: event with extreme kinematics. Numerical precision is insufficient. \n";
//       std::cout << "Eq = " << Eq << "\n";
//       std::cout << "qz = " << qz << "\n";
//       std::cout << "qy = " << qy << "\n\n";
    }
    if (q2-2*Eq*qp1/eps1-qp1*qp1/eps1/eps1 < 0){
        qy = 0.0;
	
	    plplusy = sqrt(-2*Elplus*p1p3/eps1-p1p3*p1p3/eps1/eps1-m3*m3);
        plplusx = 0.0;
	
	    pnuouty = (Elplus*Enuout+p2p3-plplusz*pnuoutz)/plplusy;
        pnuoutx = sqrt(-2*Enuout*p1p2/eps1-p1p2*p1p2/eps1/eps1-pnuouty*pnuouty);
	
 	    plminusy = (Elplus*Elminus+p3p4-plplusz*plminusz)/plplusy;
        plminusx = -plplusx-pnuoutx;
    }
    
    
    // special case where the momentum transfer to the nucleon has no y component and the outgoing neutrino has no momentum in x direction   
    if (q2-2*Eq*qp1/eps1-qp1*qp1/eps1/eps1 < 0 && -2*Enuout*p1p2/eps1-p1p2*p1p2/eps1/eps1 < 0){
       flagkinematics = 1;
//       std::cout << "Warning: event with extreme kinematics. Numerical precision is insufficient. \n";
//       std::cout << "Enuout = " << Enuout << "\n";
//       std::cout << "pnuoutz = " << pnuoutz << "\n";
//       std::cout << "pnuouty = " << pnuouty << "\n";
//       std::cout << "pnuoutx = " << pnuoutx << "\n\n";
    }
    if (q2-2*Eq*qp1/eps1-qp1*qp1/eps1/eps1 < 0 && -2*Enuout*p1p2/eps1-p1p2*p1p2/eps1/eps1-pnuouty*pnuouty < 0){
        qy = 0.0;
	
	    plplusy = sqrt(-2*Elplus*p1p3/eps1-p1p3*p1p3/eps1/eps1-m3*m3);
        plplusx = 0.0;
	
	    pnuouty = (Elplus*Enuout+p2p3-plplusz*pnuoutz)/plplusy;
        pnuoutx = 0.0;
	
	    plminusy = (Elplus*Elminus+p3p4-plplusz*plminusz)/plplusy;
        plminusx = 0.0;
    }
    
    
    // special case where the momentum transfer to the nucleon and the momentum of the anti-muon have no y component
    if (q2-2*Eq*qp1/eps1-qp1*qp1/eps1/eps1 < 0 && -2*Elplus*p1p3/eps1-p1p3*p1p3/eps1/eps1 < 0){
       flagkinematics = 1;
//       std::cout << "Warning: event with extreme kinematics. Numerical precision is insufficient. \n";
//       std::cout << "Elplus = " << Elplus << "\n";
//       std::cout << "plplusz = " << plplusz << "\n";
//       std::cout << "plplusy = " << plplusy << "\n";
//       std::cout << "plplusx = " << plplusx << "\n\n";
    }
    if (q2-2*Eq*qp1/eps1-qp1*qp1/eps1/eps1 < 0 && -2*Elplus*p1p3/eps1-p1p3*p1p3/eps1/eps1-m3*m3 < 0){
        qy = 0.0;
	
	    plplusy = 0.0;
        plplusx = 0.0;
	
	    plminusy = sqrt(-2*Elminus*p1p4/eps1-p1p4*p1p4/eps1/eps1-m4*m4);
        plminusx = 0.0;
	
        pnuouty = -plminusy;
        pnuoutx = 0.0;
    }
    
    
    // special case where the momentum transfer to the nucleon and the momenta of the muon and anti-muon have no y component
    if (q2-2*Eq*qp1/eps1-qp1*qp1/eps1/eps1 < 0 && -2*Elplus*p1p3/eps1-p1p3*p1p3/eps1/eps1-m3*m3 < 0 && -2*Elminus*p1p4/eps1-p1p4*p1p4/eps1/eps1 < 0){
       flagkinematics = 1;
//       std::cout << "Warning: event with extreme kinematics. Numerical precision is insufficient. \n";
//       std::cout << "Elminus = " << Elminus << "\n";
//       std::cout << "plminusz = " << plminusz << "\n";
//       std::cout << "plminusy = " << plminusy << "\n";
//       std::cout << "plminusx = " << plminusx << "\n\n";
    }
    if (q2-2*Eq*qp1/eps1-qp1*qp1/eps1/eps1 < 0 && -2*Elplus*p1p3/eps1-p1p3*p1p3/eps1/eps1-m3*m3 < 0 && -2*Elminus*p1p4/eps1-p1p4*p1p4/eps1/eps1-m4*m4 < 0){
        qy = 0.0;
	
	    plplusy = 0.0;
        plplusx = 0.0;
	
	    plminusy = 0.0;
        plminusx = 0.0;
	
        pnuouty = 0.0;
        pnuoutx = 0.0;
    }
    
    
    // nuclear recoil 4-momentum
    
        Eprime = M + Enuin - Enuout - Elplus - Elminus;
        Pprimex = - pnuoutx - plplusx - plminusx;
        Pprimey = - pnuouty - plplusy - plminusy;
        Pprimez = pnuinz - pnuoutz - plplusz - plminusz;
    
    return;
}


//*****************************************************************
// squared matrix elements for neutrino and anti-neutrino tridents
//*****************************************************************

double SquaredMatrixElementPLP(double GFV, double GFA){

    double GFA2, GFV2, GFAGFV, result;
    
    GFA2 = GFA*GFA;
    GFV2 = GFV*GFV;
    GFAGFV = GFA*GFV;

    result = 4*((64*GFA2*m3*m4*p1p2*p3P*p3P)/D32 - 
                (64*GFV2*m3*m4*p1p2*p3P*p3P)/D32 + 
                (64*GFA2*p1p4*p2p3*p3P*p3P)/D32 + 
                (128*GFAGFV*p1p4*p2p3*p3P*p3P)/D32 + 
                (64*GFV2*p1p4*p2p3*p3P*p3P)/D32 + 
                (64*GFA2*p1p3*p2p4*p3P*p3P)/D32 - 
                (128*GFAGFV*p1p3*p2p4*p3P*p3P)/D32 + 
                (64*GFV2*p1p3*p2p4*p3P*p3P)/D32 - 
                (128*GFA2*m3*m4*p1p2*p3P*p4P)/D3D4 + 
                (128*GFV2*m3*m4*p1p2*p3P*p4P)/D3D4 - 
                (128*GFA2*p1p4*p2p3*p3P*p4P)/D3D4 - 
                (256*GFAGFV*p1p4*p2p3*p3P*p4P)/D3D4 - 
                (128*GFV2*p1p4*p2p3*p3P*p4P)/D3D4 - 
                (128*GFA2*p1p3*p2p4*p3P*p4P)/D3D4 + 
                (256*GFAGFV*p1p3*p2p4*p3P*p4P)/D3D4 - 
                (128*GFV2*p1p3*p2p4*p3P*p4P)/D3D4 + 
                (64*GFA2*m3*m4*p1p2*p4P*p4P)/D42 - 
                (64*GFV2*m3*m4*p1p2*p4P*p4P)/D42 + 
                (64*GFA2*p1p4*p2p3*p4P*p4P)/D42 + 
                (128*GFAGFV*p1p4*p2p3*p4P*p4P)/D42 + 
                (64*GFV2*p1p4*p2p3*p4P*p4P)/D42 + 
                (64*GFA2*p1p3*p2p4*p4P*p4P)/D42 - 
                (128*GFAGFV*p1p3*p2p4*p4P*p4P)/D42 + 
                (64*GFV2*p1p3*p2p4*p4P*p4P)/D42 - 
                (16*GFA2*M2*m3*m4*p1p2*q2)/D32 - 
                (16*GFA2*M2*m3*m4*p1p2*q2)/D42 + 
                (32*GFA2*M2*m3*m4*p1p2*q2)/D3D4 + 
                (16*GFV2*M2*m3*m4*p1p2*q2)/D32 + 
                (16*GFV2*M2*m3*m4*p1p2*q2)/D42 - 
                (32*GFV2*M2*m3*m4*p1p2*q2)/D3D4 + 
                (64*GFA2*m3*m4*p1P*p2P*q2)/D3D4 - 
                (64*GFV2*m3*m4*p1P*p2P*q2)/D3D4 - 
                (16*GFA2*M2*p1p4*p2p3*q2)/D32 - 
                (16*GFA2*M2*p1p4*p2p3*q2)/D42 - 
                (32*GFA2*M2*p1p4*p2p3*q2)/D3D4 - 
                (32*GFAGFV*M2*p1p4*p2p3*q2)/D32 - 
                (32*GFAGFV*M2*p1p4*p2p3*q2)/D42 + 
                (64*GFAGFV*M2*p1p4*p2p3*q2)/D3D4 - 
                (16*GFV2*M2*p1p4*p2p3*q2)/D32 - 
                (16*GFV2*M2*p1p4*p2p3*q2)/D42 - 
                (32*GFV2*M2*p1p4*p2p3*q2)/D3D4 - 
                (16*GFA2*M2*p1p3*p2p4*q2)/D32 - 
                (16*GFA2*M2*p1p3*p2p4*q2)/D42 - 
                (32*GFA2*M2*p1p3*p2p4*q2)/D3D4 + 
                (32*GFAGFV*M2*p1p3*p2p4*q2)/D32 + 
                (32*GFAGFV*M2*p1p3*p2p4*q2)/D42 - 
                (64*GFAGFV*M2*p1p3*p2p4*q2)/D3D4 - 
                (16*GFV2*M2*p1p3*p2p4*q2)/D32 - 
                (16*GFV2*M2*p1p3*p2p4*q2)/D42 - 
                (32*GFV2*M2*p1p3*p2p4*q2)/D3D4 - 
                (32*GFA2*p1p4*p2P*p3P*q2)/D32 - 
                (32*GFA2*p1p4*p2P*p3P*q2)/D3D4 - 
                (64*GFAGFV*p1p4*p2P*p3P*q2)/D32 + 
                (64*GFAGFV*p1p4*p2P*p3P*q2)/D3D4 - 
                (32*GFV2*p1p4*p2P*p3P*q2)/D32 - 
                (32*GFV2*p1p4*p2P*p3P*q2)/D3D4 - 
                (32*GFA2*p1P*p2p4*p3P*q2)/D32 - 
                (32*GFA2*p1P*p2p4*p3P*q2)/D3D4 + 
                (64*GFAGFV*p1P*p2p4*p3P*q2)/D32 - 
                (64*GFAGFV*p1P*p2p4*p3P*q2)/D3D4 - 
                (32*GFV2*p1P*p2p4*p3P*q2)/D32 - 
                (32*GFV2*p1P*p2p4*p3P*q2)/D3D4 + 
                (64*GFA2*M2*p1p2*p3p4*q2)/D3D4 + 
                (64*GFV2*M2*p1p2*p3p4*q2)/D3D4 + 
                (64*GFA2*p1P*p2P*p3p4*q2)/D3D4 + 
                (64*GFV2*p1P*p2P*p3p4*q2)/D3D4 - 
                (32*GFA2*p1p3*p2P*p4P*q2)/D42 - 
                (32*GFA2*p1p3*p2P*p4P*q2)/D3D4 + 
                (64*GFAGFV*p1p3*p2P*p4P*q2)/D42 - 
                (64*GFAGFV*p1p3*p2P*p4P*q2)/D3D4 - 
                (32*GFV2*p1p3*p2P*p4P*q2)/D42 -
                (32*GFV2*p1p3*p2P*p4P*q2)/D3D4 - 
                (32*GFA2*p1P*p2p3*p4P*q2)/D42 - 
                (32*GFA2*p1P*p2p3*p4P*q2)/D3D4 - 
                (64*GFAGFV*p1P*p2p3*p4P*q2)/D42 + 
                (64*GFAGFV*p1P*p2p3*p4P*q2)/D3D4 - 
                (32*GFV2*p1P*p2p3*p4P*q2)/D42 - 
                (32*GFV2*p1P*p2p3*p4P*q2)/D3D4 + 
                (64*GFA2*p1p2*p3P*p4P*q2)/D3D4 + 
                (64*GFV2*p1p2*p3P*p4P*q2)/D3D4 - 
                (64*GFA2*m3*m4*p1p2*p3P*qP)/D32 + 
                (64*GFA2*m3*m4*p1p2*p3P*qP)/D3D4 + 
                (64*GFV2*m3*m4*p1p2*p3P*qP)/D32 - 
                (64*GFV2*m3*m4*p1p2*p3P*qP)/D3D4 - 
                (64*GFA2*p1p4*p2p3*p3P*qP)/D32 + 
                (64*GFA2*p1p4*p2p3*p3P*qP)/D3D4 - 
                (128*GFAGFV*p1p4*p2p3*p3P*qP)/D32 + 
                (128*GFAGFV*p1p4*p2p3*p3P*qP)/D3D4 - 
                (64*GFV2*p1p4*p2p3*p3P*qP)/D32 + 
                (64*GFV2*p1p4*p2p3*p3P*qP)/D3D4 - 
                (64*GFA2*p1p3*p2p4*p3P*qP)/D32 + 
                (64*GFA2*p1p3*p2p4*p3P*qP)/D3D4 + 
                (128*GFAGFV*p1p3*p2p4*p3P*qP)/D32 - 
                (128*GFAGFV*p1p3*p2p4*p3P*qP)/D3D4 - 
                (64*GFV2*p1p3*p2p4*p3P*qP)/D32 + 
                (64*GFV2*p1p3*p2p4*p3P*qP)/D3D4 - 
                (64*GFA2*m3*m4*p1p2*p4P*qP)/D42 + 
                (64*GFA2*m3*m4*p1p2*p4P*qP)/D3D4 + 
                (64*GFV2*m3*m4*p1p2*p4P*qP)/D42 - 
                (64*GFV2*m3*m4*p1p2*p4P*qP)/D3D4 - 
                (64*GFA2*p1p4*p2p3*p4P*qP)/D42 + 
                (64*GFA2*p1p4*p2p3*p4P*qP)/D3D4 - 
                (128*GFAGFV*p1p4*p2p3*p4P*qP)/D42 + 
                (128*GFAGFV*p1p4*p2p3*p4P*qP)/D3D4 - 
                (64*GFV2*p1p4*p2p3*p4P*qP)/D42 + 
                (64*GFV2*p1p4*p2p3*p4P*qP)/D3D4 - 
                (64*GFA2*p1p3*p2p4*p4P*qP)/D42 + 
                (64*GFA2*p1p3*p2p4*p4P*qP)/D3D4 + 
                (128*GFAGFV*p1p3*p2p4*p4P*qP)/D42 - 
                (128*GFAGFV*p1p3*p2p4*p4P*qP)/D3D4 - 
                (64*GFV2*p1p3*p2p4*p4P*qP)/D42 + 
                (64*GFV2*p1p3*p2p4*p4P*qP)/D3D4 - 
                (64*GFA2*p1p4*p2p3*qP*qP)/D3D4 - 
                (64*GFV2*p1p4*p2p3*qP*qP)/D3D4 - 
                (64*GFA2*p1p3*p2p4*qP*qP)/D3D4 - 
                (64*GFV2*p1p3*p2p4*qP*qP)/D3D4 + 
                (64*GFA2*p1p2*p3p4*qP*qP)/D3D4 + 
                (64*GFV2*p1p2*p3p4*qP*qP)/D3D4 - 
                (64*GFA2*p2p4*p3P*p3P*qp1)/D32 + 
                (128*GFAGFV*p2p4*p3P*p3P*qp1)/D32 - 
                (64*GFV2*p2p4*p3P*p3P*qp1)/D32 + 
                (64*GFA2*p2p3*p3P*p4P*qp1)/D3D4 + 
                (128*GFAGFV*p2p3*p3P*p4P*qp1)/D3D4 + 
                (64*GFV2*p2p3*p3P*p4P*qp1)/D3D4 + 
                (64*GFA2*p2p4*p3P*p4P*qp1)/D3D4 - 
                (128*GFAGFV*p2p4*p3P*p4P*qp1)/D3D4 + 
                (64*GFV2*p2p4*p3P*p4P*qp1)/D3D4 - 
                (64*GFA2*p2p3*p4P*p4P*qp1)/D42 - 
                (128*GFAGFV*p2p3*p4P*p4P*qp1)/D42 - 
                (64*GFV2*p2p3*p4P*p4P*qp1)/D42 - 
                (64*GFA2*m3*m4*p2P*qP*qp1)/D3D4 + 
                (64*GFV2*m3*m4*p2P*qP*qp1)/D3D4 + 
                (64*GFA2*p2p4*p3P*qP*qp1)/D32 - 
                (128*GFAGFV*p2p4*p3P*qP*qp1)/D32 + 
                (128*GFAGFV*p2p4*p3P*qP*qp1)/D3D4 + 
                (64*GFV2*p2p4*p3P*qP*qp1)/D32 - 
                (64*GFA2*p2P*p3p4*qP*qp1)/D3D4 - 
                (64*GFV2*p2P*p3p4*qP*qp1)/D3D4 + 
                (64*GFA2*p2p3*p4P*qP*qp1)/D42 + 
                (128*GFAGFV*p2p3*p4P*qP*qp1)/D42 - 
                (128*GFAGFV*p2p3*p4P*qP*qp1)/D3D4 + 
                (64*GFV2*p2p3*p4P*qP*qp1)/D42 - 
                (64*GFA2*p1p4*p3P*p3P*qp2)/D32 - 
                (128*GFAGFV*p1p4*p3P*p3P*qp2)/D32 - 
                (64*GFV2*p1p4*p3P*p3P*qp2)/D32 + 
                (64*GFA2*p1p3*p3P*p4P*qp2)/D3D4 - 
                (128*GFAGFV*p1p3*p3P*p4P*qp2)/D3D4 + 
                (64*GFV2*p1p3*p3P*p4P*qp2)/D3D4 + 
                (64*GFA2*p1p4*p3P*p4P*qp2)/D3D4 + 
                (128*GFAGFV*p1p4*p3P*p4P*qp2)/D3D4 + 
                (64*GFV2*p1p4*p3P*p4P*qp2)/D3D4 - 
                (64*GFA2*p1p3*p4P*p4P*qp2)/D42 + 
                (128*GFAGFV*p1p3*p4P*p4P*qp2)/D42 - 
                (64*GFV2*p1p3*p4P*p4P*qp2)/D42 - 
                (64*GFA2*m3*m4*p1P*qP*qp2)/D3D4 + 
                (64*GFV2*m3*m4*p1P*qP*qp2)/D3D4 + 
                (64*GFA2*p1p4*p3P*qP*qp2)/D32 + 
                (128*GFAGFV*p1p4*p3P*qP*qp2)/D32 - 
                (128*GFAGFV*p1p4*p3P*qP*qp2)/D3D4 + 
                (64*GFV2*p1p4*p3P*qP*qp2)/D32 - 
                (64*GFA2*p1P*p3p4*qP*qp2)/D3D4 - 
                (64*GFV2*p1P*p3p4*qP*qp2)/D3D4 + 
                (64*GFA2*p1p3*p4P*qP*qp2)/D42 - 
                (128*GFAGFV*p1p3*p4P*qP*qp2)/D42 + 
                (128*GFAGFV*p1p3*p4P*qP*qp2)/D3D4 + 
                (64*GFV2*p1p3*p4P*qP*qp2)/D42 - 
                (64*GFA2*M2*m3*m4*qp1*qp2)/D3D4 + 
                (64*GFV2*M2*m3*m4*qp1*qp2)/D3D4 - 
                (64*GFA2*M2*p3p4*qp1*qp2)/D3D4 - 
                (64*GFV2*M2*p3p4*qp1*qp2)/D3D4 - 
                (128*GFA2*p3P*p4P*qp1*qp2)/D3D4 - 
                (128*GFV2*p3P*p4P*qp1*qp2)/D3D4 + 
                (64*GFA2*p1p4*p2P*p3P*qp3)/D32 + 
                (128*GFAGFV*p1p4*p2P*p3P*qp3)/D32 + 
                (64*GFV2*p1p4*p2P*p3P*qp3)/D32 + 
                (64*GFA2*p1P*p2p4*p3P*qp3)/D32 - 
                (128*GFAGFV*p1P*p2p4*p3P*qp3)/D32 + 
                (64*GFV2*p1P*p2p4*p3P*qp3)/D32 - 
                (64*GFA2*p1p4*p2P*p4P*qp3)/D3D4 - 
                (128*GFAGFV*p1p4*p2P*p4P*qp3)/D3D4 - 
                (64*GFV2*p1p4*p2P*p4P*qp3)/D3D4 - 
                (64*GFA2*p1P*p2p4*p4P*qp3)/D3D4 + 
                (128*GFAGFV*p1P*p2p4*p4P*qp3)/D3D4 - 
                (64*GFV2*p1P*p2p4*p4P*qp3)/D3D4 + 
                (64*GFA2*p1p4*p2P*qP*qp3)/D3D4 + 
                (64*GFV2*p1p4*p2P*qP*qp3)/D3D4 + 
                (64*GFA2*p1P*p2p4*qP*qp3)/D3D4 + 
                (64*GFV2*p1P*p2p4*qP*qp3)/D3D4 - 
                (64*GFA2*p1p2*p4P*qP*qp3)/D3D4 - 
                (64*GFV2*p1p2*p4P*qP*qp3)/D3D4 + 
                (32*GFA2*M2*p2p4*qp1*qp3)/D32 + 
                (32*GFA2*M2*p2p4*qp1*qp3)/D3D4 - 
                (64*GFAGFV*M2*p2p4*qp1*qp3)/D32 + 
                (64*GFAGFV*M2*p2p4*qp1*qp3)/D3D4 + 
                (32*GFV2*M2*p2p4*qp1*qp3)/D32 + 
                (32*GFV2*M2*p2p4*qp1*qp3)/D3D4 + 
                (64*GFA2*p2P*p4P*qp1*qp3)/D3D4 + 
                (128*GFAGFV*p2P*p4P*qp1*qp3)/D3D4 + 
                (64*GFV2*p2P*p4P*qp1*qp3)/D3D4 + 
                (32*GFA2*M2*p1p4*qp2*qp3)/D32 + 
                (32*GFA2*M2*p1p4*qp2*qp3)/D3D4 + 
                (64*GFAGFV*M2*p1p4*qp2*qp3)/D32 - 
                (64*GFAGFV*M2*p1p4*qp2*qp3)/D3D4 + 
                (32*GFV2*M2*p1p4*qp2*qp3)/D32 + 
                (32*GFV2*M2*p1p4*qp2*qp3)/D3D4 + 
                (64*GFA2*p1P*p4P*qp2*qp3)/D3D4 - 
                (128*GFAGFV*p1P*p4P*qp2*qp3)/D3D4 + 
                (64*GFV2*p1P*p4P*qp2*qp3)/D3D4 - 
                (64*GFA2*p1p3*p2P*p3P*qp4)/D3D4 + 
                (128*GFAGFV*p1p3*p2P*p3P*qp4)/D3D4 - 
                (64*GFV2*p1p3*p2P*p3P*qp4)/D3D4 - 
                (64*GFA2*p1P*p2p3*p3P*qp4)/D3D4 - 
                (128*GFAGFV*p1P*p2p3*p3P*qp4)/D3D4 - 
                (64*GFV2*p1P*p2p3*p3P*qp4)/D3D4 + 
                (64*GFA2*p1p3*p2P*p4P*qp4)/D42 - 
                (128*GFAGFV*p1p3*p2P*p4P*qp4)/D42 + 
                (64*GFV2*p1p3*p2P*p4P*qp4)/D42 + 
                (64*GFA2*p1P*p2p3*p4P*qp4)/D42 + 
                (128*GFAGFV*p1P*p2p3*p4P*qp4)/D42 + 
                (64*GFV2*p1P*p2p3*p4P*qp4)/D42 + 
                (64*GFA2*p1p3*p2P*qP*qp4)/D3D4 + 
                (64*GFV2*p1p3*p2P*qP*qp4)/D3D4 + 
                (64*GFA2*p1P*p2p3*qP*qp4)/D3D4 + 
                (64*GFV2*p1P*p2p3*qP*qp4)/D3D4 - 
                (64*GFA2*p1p2*p3P*qP*qp4)/D3D4 - 
                (64*GFV2*p1p2*p3P*qP*qp4)/D3D4 + 
                (32*GFA2*M2*p2p3*qp1*qp4)/D42 + 
                (32*GFA2*M2*p2p3*qp1*qp4)/D3D4 + 
                (64*GFAGFV*M2*p2p3*qp1*qp4)/D42 - 
                (64*GFAGFV*M2*p2p3*qp1*qp4)/D3D4 + 
                (32*GFV2*M2*p2p3*qp1*qp4)/D42 + 
                (32*GFV2*M2*p2p3*qp1*qp4)/D3D4 + 
                (64*GFA2*p2P*p3P*qp1*qp4)/D3D4 - 
                (128*GFAGFV*p2P*p3P*qp1*qp4)/D3D4 + 
                (64*GFV2*p2P*p3P*qp1*qp4)/D3D4 + 
                (32*GFA2*M2*p1p3*qp2*qp4)/D42 + 
                (32*GFA2*M2*p1p3*qp2*qp4)/D3D4 - 
                (64*GFAGFV*M2*p1p3*qp2*qp4)/D42 + 
                (64*GFAGFV*M2*p1p3*qp2*qp4)/D3D4 + 
                (32*GFV2*M2*p1p3*qp2*qp4)/D42 + 
                (32*GFV2*M2*p1p3*qp2*qp4)/D3D4 + 
                (64*GFA2*p1P*p3P*qp2*qp4)/D3D4 + 
                (128*GFAGFV*p1P*p3P*qp2*qp4)/D3D4 + 
                (64*GFV2*p1P*p3P*qp2*qp4)/D3D4 - 
                (64*GFA2*M2*p1p2*qp3*qp4)/D3D4 - 
                (64*GFV2*M2*p1p2*qp3*qp4)/D3D4 - 
                (128*GFA2*p1P*p2P*qp3*qp4)/D3D4 - 
                (128*GFV2*p1P*p2P*qp3*qp4)/D3D4);

    return result;
}


double SquaredMatrixElementPLPanti(double GFV, double GFA){

    double GFA2, GFV2, GFAGFV, result;
    
    GFA2 = GFA*GFA;
    GFV2 = GFV*GFV;
    GFAGFV = GFA*GFV;

    result = 4*((64*GFA2*m3*m4*p1p2*p3P*p3P)/D32 - 
                (64*GFV2*m3*m4*p1p2*p3P*p3P)/D32 + 
                (64*GFA2*p1p4*p2p3*p3P*p3P)/D32 - 
                (128*GFAGFV*p1p4*p2p3*p3P*p3P)/D32 + 
                (64*GFV2*p1p4*p2p3*p3P*p3P)/D32 + 
                (64*GFA2*p1p3*p2p4*p3P*p3P)/D32 + 
                (128*GFAGFV*p1p3*p2p4*p3P*p3P)/D32 + 
                (64*GFV2*p1p3*p2p4*p3P*p3P)/D32 - 
                (128*GFA2*m3*m4*p1p2*p3P*p4P)/D3D4 + 
                (128*GFV2*m3*m4*p1p2*p3P*p4P)/D3D4 - 
                (128*GFA2*p1p4*p2p3*p3P*p4P)/D3D4 + 
                (256*GFAGFV*p1p4*p2p3*p3P*p4P)/D3D4 - 
                (128*GFV2*p1p4*p2p3*p3P*p4P)/D3D4 - 
                (128*GFA2*p1p3*p2p4*p3P*p4P)/D3D4 - 
                (256*GFAGFV*p1p3*p2p4*p3P*p4P)/D3D4 - 
                (128*GFV2*p1p3*p2p4*p3P*p4P)/D3D4 + 
                (64*GFA2*m3*m4*p1p2*p4P*p4P)/D42 - 
                (64*GFV2*m3*m4*p1p2*p4P*p4P)/D42 + 
                (64*GFA2*p1p4*p2p3*p4P*p4P)/D42 - 
                (128*GFAGFV*p1p4*p2p3*p4P*p4P)/D42 + 
                (64*GFV2*p1p4*p2p3*p4P*p4P)/D42 + 
                (64*GFA2*p1p3*p2p4*p4P*p4P)/D42 + 
                (128*GFAGFV*p1p3*p2p4*p4P*p4P)/D42 + 
                (64*GFV2*p1p3*p2p4*p4P*p4P)/D42 - 
                (16*GFA2*M2*m3*m4*p1p2*q2)/D32 - 
                (16*GFA2*M2*m3*m4*p1p2*q2)/D42 + 
                (32*GFA2*M2*m3*m4*p1p2*q2)/D3D4 + 
                (16*GFV2*M2*m3*m4*p1p2*q2)/D32 + 
                (16*GFV2*M2*m3*m4*p1p2*q2)/D42 - 
                (32*GFV2*M2*m3*m4*p1p2*q2)/D3D4 + 
                (64*GFA2*m3*m4*p1P*p2P*q2)/D3D4 - 
                (64*GFV2*m3*m4*p1P*p2P*q2)/D3D4 - 
                (16*GFA2*M2*p1p4*p2p3*q2)/D32 - 
                (16*GFA2*M2*p1p4*p2p3*q2)/D42 - 
                (32*GFA2*M2*p1p4*p2p3*q2)/D3D4 + 
                (32*GFAGFV*M2*p1p4*p2p3*q2)/D32 + 
                (32*GFAGFV*M2*p1p4*p2p3*q2)/D42 - 
                (64*GFAGFV*M2*p1p4*p2p3*q2)/D3D4 - 
                (16*GFV2*M2*p1p4*p2p3*q2)/D32 - 
                (16*GFV2*M2*p1p4*p2p3*q2)/D42 - 
                (32*GFV2*M2*p1p4*p2p3*q2)/D3D4 - 
                (16*GFA2*M2*p1p3*p2p4*q2)/D32 - 
                (16*GFA2*M2*p1p3*p2p4*q2)/D42 - 
                (32*GFA2*M2*p1p3*p2p4*q2)/D3D4 - 
                (32*GFAGFV*M2*p1p3*p2p4*q2)/D32 - 
                (32*GFAGFV*M2*p1p3*p2p4*q2)/D42 + 
                (64*GFAGFV*M2*p1p3*p2p4*q2)/D3D4 - 
                (16*GFV2*M2*p1p3*p2p4*q2)/D32 - 
                (16*GFV2*M2*p1p3*p2p4*q2)/D42 - 
                (32*GFV2*M2*p1p3*p2p4*q2)/D3D4 - 
                (32*GFA2*p1p4*p2P*p3P*q2)/D32 - 
                (32*GFA2*p1p4*p2P*p3P*q2)/D3D4 + 
                (64*GFAGFV*p1p4*p2P*p3P*q2)/D32 - 
                (64*GFAGFV*p1p4*p2P*p3P*q2)/D3D4 - 
                (32*GFV2*p1p4*p2P*p3P*q2)/D32 - 
                (32*GFV2*p1p4*p2P*p3P*q2)/D3D4 - 
                (32*GFA2*p1P*p2p4*p3P*q2)/D32 - 
                (32*GFA2*p1P*p2p4*p3P*q2)/D3D4 - 
                (64*GFAGFV*p1P*p2p4*p3P*q2)/D32 + 
                (64*GFAGFV*p1P*p2p4*p3P*q2)/D3D4 - 
                (32*GFV2*p1P*p2p4*p3P*q2)/D32 - 
                (32*GFV2*p1P*p2p4*p3P*q2)/D3D4 + 
                (64*GFA2*M2*p1p2*p3p4*q2)/D3D4 + 
                (64*GFV2*M2*p1p2*p3p4*q2)/D3D4 + 
                (64*GFA2*p1P*p2P*p3p4*q2)/D3D4 + 
                (64*GFV2*p1P*p2P*p3p4*q2)/D3D4 - 
                (32*GFA2*p1p3*p2P*p4P*q2)/D42 - 
                (32*GFA2*p1p3*p2P*p4P*q2)/D3D4 - 
                (64*GFAGFV*p1p3*p2P*p4P*q2)/D42 + 
                (64*GFAGFV*p1p3*p2P*p4P*q2)/D3D4 - 
                (32*GFV2*p1p3*p2P*p4P*q2)/D42 - 
                (32*GFV2*p1p3*p2P*p4P*q2)/D3D4 - 
                (32*GFA2*p1P*p2p3*p4P*q2)/D42 - 
                (32*GFA2*p1P*p2p3*p4P*q2)/D3D4 + 
                (64*GFAGFV*p1P*p2p3*p4P*q2)/D42 - 
                (64*GFAGFV*p1P*p2p3*p4P*q2)/D3D4 - 
                (32*GFV2*p1P*p2p3*p4P*q2)/D42 - 
                (32*GFV2*p1P*p2p3*p4P*q2)/D3D4 + 
                (64*GFA2*p1p2*p3P*p4P*q2)/D3D4 + 
                (64*GFV2*p1p2*p3P*p4P*q2)/D3D4 - 
                (64*GFA2*m3*m4*p1p2*p3P*qP)/D32 + 
                (64*GFA2*m3*m4*p1p2*p3P*qP)/D3D4 + 
                (64*GFV2*m3*m4*p1p2*p3P*qP)/D32 - 
                (64*GFV2*m3*m4*p1p2*p3P*qP)/D3D4 - 
                (64*GFA2*p1p4*p2p3*p3P*qP)/D32 + 
                (64*GFA2*p1p4*p2p3*p3P*qP)/D3D4 + 
                (128*GFAGFV*p1p4*p2p3*p3P*qP)/D32 - 
                (128*GFAGFV*p1p4*p2p3*p3P*qP)/D3D4 - 
                (64*GFV2*p1p4*p2p3*p3P*qP)/D32 + 
                (64*GFV2*p1p4*p2p3*p3P*qP)/D3D4 - 
                (64*GFA2*p1p3*p2p4*p3P*qP)/D32 + 
                (64*GFA2*p1p3*p2p4*p3P*qP)/D3D4 - 
                (128*GFAGFV*p1p3*p2p4*p3P*qP)/D32 + 
                (128*GFAGFV*p1p3*p2p4*p3P*qP)/D3D4 - 
                (64*GFV2*p1p3*p2p4*p3P*qP)/D32 + 
                (64*GFV2*p1p3*p2p4*p3P*qP)/D3D4 - 
                (64*GFA2*m3*m4*p1p2*p4P*qP)/D42 + 
                (64*GFA2*m3*m4*p1p2*p4P*qP)/D3D4 + 
                (64*GFV2*m3*m4*p1p2*p4P*qP)/D42 - 
                (64*GFV2*m3*m4*p1p2*p4P*qP)/D3D4 - 
                (64*GFA2*p1p4*p2p3*p4P*qP)/D42 + 
                (64*GFA2*p1p4*p2p3*p4P*qP)/D3D4 + 
                (128*GFAGFV*p1p4*p2p3*p4P*qP)/D42 - 
                (128*GFAGFV*p1p4*p2p3*p4P*qP)/D3D4 - 
                (64*GFV2*p1p4*p2p3*p4P*qP)/D42 + 
                (64*GFV2*p1p4*p2p3*p4P*qP)/D3D4 - 
                (64*GFA2*p1p3*p2p4*p4P*qP)/D42 + 
                (64*GFA2*p1p3*p2p4*p4P*qP)/D3D4 - 
                (128*GFAGFV*p1p3*p2p4*p4P*qP)/D42 + 
                (128*GFAGFV*p1p3*p2p4*p4P*qP)/D3D4 - 
                (64*GFV2*p1p3*p2p4*p4P*qP)/D42 + 
                (64*GFV2*p1p3*p2p4*p4P*qP)/D3D4 - 
                (64*GFA2*p1p4*p2p3*qP*qP)/D3D4 - 
                (64*GFV2*p1p4*p2p3*qP*qP)/D3D4 - 
                (64*GFA2*p1p3*p2p4*qP*qP)/D3D4 - 
                (64*GFV2*p1p3*p2p4*qP*qP)/D3D4 + 
                (64*GFA2*p1p2*p3p4*qP*qP)/D3D4 + 
                (64*GFV2*p1p2*p3p4*qP*qP)/D3D4 - 
                (64*GFA2*p2p4*p3P*p3P*qp1)/D32 - 
                (128*GFAGFV*p2p4*p3P*p3P*qp1)/D32 - 
                (64*GFV2*p2p4*p3P*p3P*qp1)/D32 + 
                (64*GFA2*p2p3*p3P*p4P*qp1)/D3D4 - 
                (128*GFAGFV*p2p3*p3P*p4P*qp1)/D3D4 + 
                (64*GFV2*p2p3*p3P*p4P*qp1)/D3D4 + 
                (64*GFA2*p2p4*p3P*p4P*qp1)/D3D4 + 
                (128*GFAGFV*p2p4*p3P*p4P*qp1)/D3D4 + 
                (64*GFV2*p2p4*p3P*p4P*qp1)/D3D4 - 
                (64*GFA2*p2p3*p4P*p4P*qp1)/D42 + 
                (128*GFAGFV*p2p3*p4P*p4P*qp1)/D42 - 
                (64*GFV2*p2p3*p4P*p4P*qp1)/D42 - 
                (64*GFA2*m3*m4*p2P*qP*qp1)/D3D4 + 
                (64*GFV2*m3*m4*p2P*qP*qp1)/D3D4 + 
                (64*GFA2*p2p4*p3P*qP*qp1)/D32 + 
                (128*GFAGFV*p2p4*p3P*qP*qp1)/D32 - 
                (128*GFAGFV*p2p4*p3P*qP*qp1)/D3D4 + 
                (64*GFV2*p2p4*p3P*qP*qp1)/D32 - 
                (64*GFA2*p2P*p3p4*qP*qp1)/D3D4 - 
                (64*GFV2*p2P*p3p4*qP*qp1)/D3D4 + 
                (64*GFA2*p2p3*p4P*qP*qp1)/D42 - 
                (128*GFAGFV*p2p3*p4P*qP*qp1)/D42 + 
                (128*GFAGFV*p2p3*p4P*qP*qp1)/D3D4 + 
                (64*GFV2*p2p3*p4P*qP*qp1)/D42 - 
                (64*GFA2*p1p4*p3P*p3P*qp2)/D32 + 
                (128*GFAGFV*p1p4*p3P*p3P*qp2)/D32 - 
                (64*GFV2*p1p4*p3P*p3P*qp2)/D32 + 
                (64*GFA2*p1p3*p3P*p4P*qp2)/D3D4 + 
                (128*GFAGFV*p1p3*p3P*p4P*qp2)/D3D4 + 
                (64*GFV2*p1p3*p3P*p4P*qp2)/D3D4 + 
                (64*GFA2*p1p4*p3P*p4P*qp2)/D3D4 - 
                (128*GFAGFV*p1p4*p3P*p4P*qp2)/D3D4 + 
                (64*GFV2*p1p4*p3P*p4P*qp2)/D3D4 - 
                (64*GFA2*p1p3*p4P*p4P*qp2)/D42 - 
                (128*GFAGFV*p1p3*p4P*p4P*qp2)/D42 - 
                (64*GFV2*p1p3*p4P*p4P*qp2)/D42 - 
                (64*GFA2*m3*m4*p1P*qP*qp2)/D3D4 + 
                (64*GFV2*m3*m4*p1P*qP*qp2)/D3D4 + 
                (64*GFA2*p1p4*p3P*qP*qp2)/D32 - 
                (128*GFAGFV*p1p4*p3P*qP*qp2)/D32 + 
                (128*GFAGFV*p1p4*p3P*qP*qp2)/D3D4 + 
                (64*GFV2*p1p4*p3P*qP*qp2)/D32 - 
                (64*GFA2*p1P*p3p4*qP*qp2)/D3D4 - 
                (64*GFV2*p1P*p3p4*qP*qp2)/D3D4 + 
                (64*GFA2*p1p3*p4P*qP*qp2)/D42 + 
                (128*GFAGFV*p1p3*p4P*qP*qp2)/D42 - 
                (128*GFAGFV*p1p3*p4P*qP*qp2)/D3D4 + 
                (64*GFV2*p1p3*p4P*qP*qp2)/D42 - 
                (64*GFA2*M2*m3*m4*qp1*qp2)/D3D4 + 
                (64*GFV2*M2*m3*m4*qp1*qp2)/D3D4 - 
                (64*GFA2*M2*p3p4*qp1*qp2)/D3D4 - 
                (64*GFV2*M2*p3p4*qp1*qp2)/D3D4 - 
                (128*GFA2*p3P*p4P*qp1*qp2)/D3D4 - 
                (128*GFV2*p3P*p4P*qp1*qp2)/D3D4 + 
                (64*GFA2*p1p4*p2P*p3P*qp3)/D32 - 
                (128*GFAGFV*p1p4*p2P*p3P*qp3)/D32 + 
                (64*GFV2*p1p4*p2P*p3P*qp3)/D32 + 
                (64*GFA2*p1P*p2p4*p3P*qp3)/D32 + 
                (128*GFAGFV*p1P*p2p4*p3P*qp3)/D32 + 
                (64*GFV2*p1P*p2p4*p3P*qp3)/D32 - 
                (64*GFA2*p1p4*p2P*p4P*qp3)/D3D4 + 
                (128*GFAGFV*p1p4*p2P*p4P*qp3)/D3D4 - 
                (64*GFV2*p1p4*p2P*p4P*qp3)/D3D4 - 
                (64*GFA2*p1P*p2p4*p4P*qp3)/D3D4 - 
                (128*GFAGFV*p1P*p2p4*p4P*qp3)/D3D4 - 
                (64*GFV2*p1P*p2p4*p4P*qp3)/D3D4 + 
                (64*GFA2*p1p4*p2P*qP*qp3)/D3D4 + 
                (64*GFV2*p1p4*p2P*qP*qp3)/D3D4 + 
                (64*GFA2*p1P*p2p4*qP*qp3)/D3D4 + 
                (64*GFV2*p1P*p2p4*qP*qp3)/D3D4 - 
                (64*GFA2*p1p2*p4P*qP*qp3)/D3D4 - 
                (64*GFV2*p1p2*p4P*qP*qp3)/D3D4 + 
                (32*GFA2*M2*p2p4*qp1*qp3)/D32 + 
                (32*GFA2*M2*p2p4*qp1*qp3)/D3D4 + 
                (64*GFAGFV*M2*p2p4*qp1*qp3)/D32 - 
                (64*GFAGFV*M2*p2p4*qp1*qp3)/D3D4 + 
                (32*GFV2*M2*p2p4*qp1*qp3)/D32 + 
                (32*GFV2*M2*p2p4*qp1*qp3)/D3D4 + 
                (64*GFA2*p2P*p4P*qp1*qp3)/D3D4 - 
                (128*GFAGFV*p2P*p4P*qp1*qp3)/D3D4 + 
                (64*GFV2*p2P*p4P*qp1*qp3)/D3D4 + 
                (32*GFA2*M2*p1p4*qp2*qp3)/D32 + 
                (32*GFA2*M2*p1p4*qp2*qp3)/D3D4 - 
                (64*GFAGFV*M2*p1p4*qp2*qp3)/D32 + 
                (64*GFAGFV*M2*p1p4*qp2*qp3)/D3D4 + 
                (32*GFV2*M2*p1p4*qp2*qp3)/D32 + 
                (32*GFV2*M2*p1p4*qp2*qp3)/D3D4 + 
                (64*GFA2*p1P*p4P*qp2*qp3)/D3D4 + 
                (128*GFAGFV*p1P*p4P*qp2*qp3)/D3D4 + 
                (64*GFV2*p1P*p4P*qp2*qp3)/D3D4 - 
                (64*GFA2*p1p3*p2P*p3P*qp4)/D3D4 - 
                (128*GFAGFV*p1p3*p2P*p3P*qp4)/D3D4 - 
                (64*GFV2*p1p3*p2P*p3P*qp4)/D3D4 - 
                (64*GFA2*p1P*p2p3*p3P*qp4)/D3D4 + 
                (128*GFAGFV*p1P*p2p3*p3P*qp4)/D3D4 - 
                (64*GFV2*p1P*p2p3*p3P*qp4)/D3D4 + 
                (64*GFA2*p1p3*p2P*p4P*qp4)/D42 + 
                (128*GFAGFV*p1p3*p2P*p4P*qp4)/D42 + 
                (64*GFV2*p1p3*p2P*p4P*qp4)/D42 + 
                (64*GFA2*p1P*p2p3*p4P*qp4)/D42 - 
                (128*GFAGFV*p1P*p2p3*p4P*qp4)/D42 + 
                (64*GFV2*p1P*p2p3*p4P*qp4)/D42 + 
                (64*GFA2*p1p3*p2P*qP*qp4)/D3D4 + 
                (64*GFV2*p1p3*p2P*qP*qp4)/D3D4 + 
                (64*GFA2*p1P*p2p3*qP*qp4)/D3D4 + 
                (64*GFV2*p1P*p2p3*qP*qp4)/D3D4 - 
                (64*GFA2*p1p2*p3P*qP*qp4)/D3D4 - 
                (64*GFV2*p1p2*p3P*qP*qp4)/D3D4 + 
                (32*GFA2*M2*p2p3*qp1*qp4)/D42 + 
                (32*GFA2*M2*p2p3*qp1*qp4)/D3D4 - 
                (64*GFAGFV*M2*p2p3*qp1*qp4)/D42 + 
                (64*GFAGFV*M2*p2p3*qp1*qp4)/D3D4 + 
                (32*GFV2*M2*p2p3*qp1*qp4)/D42 + 
                (32*GFV2*M2*p2p3*qp1*qp4)/D3D4 + 
                (64*GFA2*p2P*p3P*qp1*qp4)/D3D4 + 
                (128*GFAGFV*p2P*p3P*qp1*qp4)/D3D4 + 
                (64*GFV2*p2P*p3P*qp1*qp4)/D3D4 + 
                (32*GFA2*M2*p1p3*qp2*qp4)/D42 + 
                (32*GFA2*M2*p1p3*qp2*qp4)/D3D4 + 
                (64*GFAGFV*M2*p1p3*qp2*qp4)/D42 - 
                (64*GFAGFV*M2*p1p3*qp2*qp4)/D3D4 + 
                (32*GFV2*M2*p1p3*qp2*qp4)/D42 + 
                (32*GFV2*M2*p1p3*qp2*qp4)/D3D4 + 
                (64*GFA2*p1P*p3P*qp2*qp4)/D3D4 - 
                (128*GFAGFV*p1P*p3P*qp2*qp4)/D3D4 + 
                (64*GFV2*p1P*p3P*qp2*qp4)/D3D4 - 
                (64*GFA2*M2*p1p2*qp3*qp4)/D3D4 - 
                (64*GFV2*M2*p1p2*qp3*qp4)/D3D4 - 
                (128*GFA2*p1P*p2P*qp3*qp4)/D3D4 - 
                (128*GFV2*p1P*p2P*qp3*qp4)/D3D4);

    return result;
}



double SquaredMatrixElementL(double GFV, double GFA){

    double GFA2, GFV2, GFAGFV, result;
    
    GFA2 = GFA*GFA;
    GFV2 = GFV*GFV;
    GFAGFV = GFA*GFV;
    
    result = 4*((-64*GFA2*m3*m3*m3*m4*p1p2)/D32 + 
                (64*GFV2*m3*m3*m3*m4*p1p2)/D32 - 
                (64*GFA2*m3*m4*m4*m4*p1p2)/D42 + 
                (64*GFV2*m3*m4*m4*m4*p1p2)/D42 - 
                (64*GFA2*m3*m3*p1p4*p2p3)/D32 - 
                (128*GFAGFV*m3*m3*p1p4*p2p3)/D32 - 
                (64*GFV2*m3*m3*p1p4*p2p3)/D32 - 
                (64*GFA2*m4*m4*p1p4*p2p3)/D42 - 
                (128*GFAGFV*m4*m4*p1p4*p2p3)/D42 - 
                (64*GFV2*m4*m4*p1p4*p2p3)/D42 - 
                (64*GFA2*m3*m3*p1p3*p2p4)/D32 + 
                (128*GFAGFV*m3*m3*p1p3*p2p4)/D32 - 
                (64*GFV2*m3*m3*p1p3*p2p4)/D32 - 
                (64*GFA2*m4*m4*p1p3*p2p4)/D42 + 
                (128*GFAGFV*m4*m4*p1p3*p2p4)/D42 - 
                (64*GFV2*m4*m4*p1p3*p2p4)/D42 - 
                (128*GFA2*m3*m4*p1p2*p3p4)/D3D4 + 
                (128*GFV2*m3*m4*p1p2*p3p4)/D3D4 - 
                (128*GFA2*p1p4*p2p3*p3p4)/D3D4 - 
                (256*GFAGFV*p1p4*p2p3*p3p4)/D3D4 - 
                (128*GFV2*p1p4*p2p3*p3p4)/D3D4 - 
                (128*GFA2*p1p3*p2p4*p3p4)/D3D4 + 
                (256*GFAGFV*p1p3*p2p4*p3p4)/D3D4 - 
                (128*GFV2*p1p3*p2p4*p3p4)/D3D4 + 
                (64*GFA2*m3*m4*p1p2*q2)/D32 + 
                (64*GFA2*m3*m4*p1p2*q2)/D42 - 
                (64*GFA2*m3*m4*p1p2*q2)/D3D4 - 
                (64*GFV2*m3*m4*p1p2*q2)/D32 - 
                (64*GFV2*m3*m4*p1p2*q2)/D42 + 
                (64*GFV2*m3*m4*p1p2*q2)/D3D4 + 
                (32*GFA2*p1p4*p2p3*q2)/D32 + 
                (32*GFA2*p1p4*p2p3*q2)/D42 + 
                (64*GFAGFV*p1p4*p2p3*q2)/D32 + 
                (64*GFAGFV*p1p4*p2p3*q2)/D42 - 
                (128*GFAGFV*p1p4*p2p3*q2)/D3D4 + 
                (32*GFV2*p1p4*p2p3*q2)/D32 + 
                (32*GFV2*p1p4*p2p3*q2)/D42 + 
                (32*GFA2*p1p3*p2p4*q2)/D32 + 
                (32*GFA2*p1p3*p2p4*q2)/D42 - 
                (64*GFAGFV*p1p3*p2p4*q2)/D32 - 
                (64*GFAGFV*p1p3*p2p4*q2)/D42 + 
                (128*GFAGFV*p1p3*p2p4*q2)/D3D4 + 
                (32*GFV2*p1p3*p2p4*q2)/D32 + 
                (32*GFV2*p1p3*p2p4*q2)/D42 - 
                (64*GFA2*p1p2*p3p4*q2)/D3D4 - 
                (64*GFV2*p1p2*p3p4*q2)/D3D4 + 
                (64*GFA2*m4*m4*p2p3*qp1)/D42 + 
                (128*GFAGFV*m4*m4*p2p3*qp1)/D42 + 
                (64*GFV2*m4*m4*p2p3*qp1)/D42 + 
                (64*GFA2*m3*m3*p2p4*qp1)/D32 - 
                (128*GFAGFV*m3*m3*p2p4*qp1)/D32 + 
                (64*GFV2*m3*m3*p2p4*qp1)/D32 + 
                (64*GFA2*p2p3*p3p4*qp1)/D3D4 + 
                (128*GFAGFV*p2p3*p3p4*qp1)/D3D4 + 
                (64*GFV2*p2p3*p3p4*qp1)/D3D4 + 
                (64*GFA2*p2p4*p3p4*qp1)/D3D4 - 
                (128*GFAGFV*p2p4*p3p4*qp1)/D3D4 + 
                (64*GFV2*p2p4*p3p4*qp1)/D3D4 + 
                (64*GFA2*m4*m4*p1p3*qp2)/D42 - 
                (128*GFAGFV*m4*m4*p1p3*qp2)/D42 + 
                (64*GFV2*m4*m4*p1p3*qp2)/D42 + 
                (64*GFA2*m3*m3*p1p4*qp2)/D32 + 
                (128*GFAGFV*m3*m3*p1p4*qp2)/D32 + 
                (64*GFV2*m3*m3*p1p4*qp2)/D32 + 
                (64*GFA2*p1p3*p3p4*qp2)/D3D4 - 
                (128*GFAGFV*p1p3*p3p4*qp2)/D3D4 + 
                (64*GFV2*p1p3*p3p4*qp2)/D3D4 + 
                (64*GFA2*p1p4*p3p4*qp2)/D3D4 + 
                (128*GFAGFV*p1p4*p3p4*qp2)/D3D4 + 
                (64*GFV2*p1p4*p3p4*qp2)/D3D4 + 
                (128*GFA2*m3*m4*qp1*qp2)/D3D4 - 
                (128*GFV2*m3*m4*qp1*qp2)/D3D4 - 
                (64*GFA2*m3*m4*p1p2*qp3)/D32 + 
                (64*GFA2*m3*m4*p1p2*qp3)/D3D4 + 
                (64*GFV2*m3*m4*p1p2*qp3)/D32 - 
                (64*GFV2*m3*m4*p1p2*qp3)/D3D4 + 
                (64*GFA2*p1p4*p2p3*qp3)/D3D4 + 
                (128*GFAGFV*p1p4*p2p3*qp3)/D3D4 + 
                (64*GFV2*p1p4*p2p3*qp3)/D3D4 + 
                (64*GFA2*p1p3*p2p4*qp3)/D3D4 - 
                (128*GFAGFV*p1p3*p2p4*qp3)/D3D4 + 
                (64*GFV2*p1p3*p2p4*qp3)/D3D4 - 
                (128*GFA2*p1p4*p2p4*qp3)/D3D4 - 
                (128*GFV2*p1p4*p2p4*qp3)/D3D4 - 
                (64*GFA2*p2p4*qp1*qp3)/D32 + 
                (128*GFAGFV*p2p4*qp1*qp3)/D32 - 
                (64*GFV2*p2p4*qp1*qp3)/D32 - 
                (64*GFA2*p1p4*qp2*qp3)/D32 - 
                (128*GFAGFV*p1p4*qp2*qp3)/D32 - 
                (64*GFV2*p1p4*qp2*qp3)/D32 - 
                (64*GFA2*m3*m4*p1p2*qp4)/D42 + 
                (64*GFA2*m3*m4*p1p2*qp4)/D3D4 + 
                (64*GFV2*m3*m4*p1p2*qp4)/D42 - 
                (64*GFV2*m3*m4*p1p2*qp4)/D3D4 - 
                (128*GFA2*p1p3*p2p3*qp4)/D3D4 - 
                (128*GFV2*p1p3*p2p3*qp4)/D3D4 + 
                (64*GFA2*p1p4*p2p3*qp4)/D3D4 + 
                (128*GFAGFV*p1p4*p2p3*qp4)/D3D4 + 
                (64*GFV2*p1p4*p2p3*qp4)/D3D4 + 
                (64*GFA2*p1p3*p2p4*qp4)/D3D4 - 
                (128*GFAGFV*p1p3*p2p4*qp4)/D3D4 + 
                (64*GFV2*p1p3*p2p4*qp4)/D3D4 - 
                (64*GFA2*p2p3*qp1*qp4)/D42 - 
                (128*GFAGFV*p2p3*qp1*qp4)/D42 - 
                (64*GFV2*p2p3*qp1*qp4)/D42 - 
                (64*GFA2*p1p3*qp2*qp4)/D42 + 
                (128*GFAGFV*p1p3*qp2*qp4)/D42 - 
                (64*GFV2*p1p3*qp2*qp4)/D42);
    
    return result;}
    
    
    
double SquaredMatrixElementLanti(double GFV, double GFA){

    double GFA2, GFV2, GFAGFV, result;
    
    GFA2 = GFA*GFA;
    GFV2 = GFV*GFV;
    GFAGFV = GFA*GFV;
    
    result = 4*((-64*GFA2*m3*m3*m3*m4*p1p2)/D32 + 
                (64*GFV2*m3*m3*m3*m4*p1p2)/D32 - 
                (64*GFA2*m3*m4*m4*m4*p1p2)/D42 + 
                (64*GFV2*m3*m4*m4*m4*p1p2)/D42 - 
                (64*GFA2*m3*m3*p1p4*p2p3)/D32 - 
                (128*GFAGFV*m3*m3*p1p4*p2p3)/D32 - 
                (64*GFV2*m3*m3*p1p4*p2p3)/D32 - 
                (64*GFA2*m4*m4*p1p4*p2p3)/D42 - 
                (128*GFAGFV*m4*m4*p1p4*p2p3)/D42 - 
                (64*GFV2*m4*m4*p1p4*p2p3)/D42 - 
                (64*GFA2*m3*m3*p1p3*p2p4)/D32 + 
                (128*GFAGFV*m3*m3*p1p3*p2p4)/D32 - 
                (64*GFV2*m3*m3*p1p3*p2p4)/D32 - 
                (64*GFA2*m4*m4*p1p3*p2p4)/D42 + 
                (128*GFAGFV*m4*m4*p1p3*p2p4)/D42 - 
                (64*GFV2*m4*m4*p1p3*p2p4)/D42 - 
                (128*GFA2*m3*m4*p1p2*p3p4)/D3D4 + 
                (128*GFV2*m3*m4*p1p2*p3p4)/D3D4 - 
                (128*GFA2*p1p4*p2p3*p3p4)/D3D4 - 
                (256*GFAGFV*p1p4*p2p3*p3p4)/D3D4 - 
                (128*GFV2*p1p4*p2p3*p3p4)/D3D4 - 
                (128*GFA2*p1p3*p2p4*p3p4)/D3D4 + 
                (256*GFAGFV*p1p3*p2p4*p3p4)/D3D4 - 
                (128*GFV2*p1p3*p2p4*p3p4)/D3D4 + 
                (64*GFA2*m3*m4*p1p2*q2)/D32 + 
                (64*GFA2*m3*m4*p1p2*q2)/D42 - 
                (64*GFA2*m3*m4*p1p2*q2)/D3D4 - 
                (64*GFV2*m3*m4*p1p2*q2)/D32 - 
                (64*GFV2*m3*m4*p1p2*q2)/D42 + 
                (64*GFV2*m3*m4*p1p2*q2)/D3D4 + 
                (32*GFA2*p1p4*p2p3*q2)/D32 + 
                (32*GFA2*p1p4*p2p3*q2)/D42 + 
                (64*GFAGFV*p1p4*p2p3*q2)/D32 + 
                (64*GFAGFV*p1p4*p2p3*q2)/D42 - 
                (128*GFAGFV*p1p4*p2p3*q2)/D3D4 + 
                (32*GFV2*p1p4*p2p3*q2)/D32 + 
                (32*GFV2*p1p4*p2p3*q2)/D42 + 
                (32*GFA2*p1p3*p2p4*q2)/D32 + 
                (32*GFA2*p1p3*p2p4*q2)/D42 - 
                (64*GFAGFV*p1p3*p2p4*q2)/D32 - 
                (64*GFAGFV*p1p3*p2p4*q2)/D42 + 
                (128*GFAGFV*p1p3*p2p4*q2)/D3D4 + 
                (32*GFV2*p1p3*p2p4*q2)/D32 + 
                (32*GFV2*p1p3*p2p4*q2)/D42 - 
                (64*GFA2*p1p2*p3p4*q2)/D3D4 - 
                (64*GFV2*p1p2*p3p4*q2)/D3D4 + 
                (64*GFA2*m4*m4*p2p3*qp1)/D42 + 
                (128*GFAGFV*m4*m4*p2p3*qp1)/D42 + 
                (64*GFV2*m4*m4*p2p3*qp1)/D42 + 
                (64*GFA2*m3*m3*p2p4*qp1)/D32 - 
                (128*GFAGFV*m3*m3*p2p4*qp1)/D32 + 
                (64*GFV2*m3*m3*p2p4*qp1)/D32 + 
                (64*GFA2*p2p3*p3p4*qp1)/D3D4 + 
                (128*GFAGFV*p2p3*p3p4*qp1)/D3D4 + 
                (64*GFV2*p2p3*p3p4*qp1)/D3D4 + 
                (64*GFA2*p2p4*p3p4*qp1)/D3D4 - 
                (128*GFAGFV*p2p4*p3p4*qp1)/D3D4 + 
                (64*GFV2*p2p4*p3p4*qp1)/D3D4 + 
                (64*GFA2*m4*m4*p1p3*qp2)/D42 - 
                (128*GFAGFV*m4*m4*p1p3*qp2)/D42 + 
                (64*GFV2*m4*m4*p1p3*qp2)/D42 + 
                (64*GFA2*m3*m3*p1p4*qp2)/D32 + 
                (128*GFAGFV*m3*m3*p1p4*qp2)/D32 + 
                (64*GFV2*m3*m3*p1p4*qp2)/D32 + 
                (64*GFA2*p1p3*p3p4*qp2)/D3D4 - 
                (128*GFAGFV*p1p3*p3p4*qp2)/D3D4 + 
                (64*GFV2*p1p3*p3p4*qp2)/D3D4 + 
                (64*GFA2*p1p4*p3p4*qp2)/D3D4 + 
                (128*GFAGFV*p1p4*p3p4*qp2)/D3D4 + 
                (64*GFV2*p1p4*p3p4*qp2)/D3D4 + 
                (128*GFA2*m3*m4*qp1*qp2)/D3D4 - 
                (128*GFV2*m3*m4*qp1*qp2)/D3D4 - 
                (64*GFA2*m3*m4*p1p2*qp3)/D32 + 
                (64*GFA2*m3*m4*p1p2*qp3)/D3D4 + 
                (64*GFV2*m3*m4*p1p2*qp3)/D32 - 
                (64*GFV2*m3*m4*p1p2*qp3)/D3D4 + 
                (64*GFA2*p1p4*p2p3*qp3)/D3D4 + 
                (128*GFAGFV*p1p4*p2p3*qp3)/D3D4 + 
                (64*GFV2*p1p4*p2p3*qp3)/D3D4 + 
                (64*GFA2*p1p3*p2p4*qp3)/D3D4 - 
                (128*GFAGFV*p1p3*p2p4*qp3)/D3D4 + 
                (64*GFV2*p1p3*p2p4*qp3)/D3D4 - 
                (128*GFA2*p1p4*p2p4*qp3)/D3D4 - 
                (128*GFV2*p1p4*p2p4*qp3)/D3D4 - 
                (64*GFA2*p2p4*qp1*qp3)/D32 + 
                (128*GFAGFV*p2p4*qp1*qp3)/D32 - 
                (64*GFV2*p2p4*qp1*qp3)/D32 - 
                (64*GFA2*p1p4*qp2*qp3)/D32 - 
                (128*GFAGFV*p1p4*qp2*qp3)/D32 - 
                (64*GFV2*p1p4*qp2*qp3)/D32 - 
                (64*GFA2*m3*m4*p1p2*qp4)/D42 + 
                (64*GFA2*m3*m4*p1p2*qp4)/D3D4 + 
                (64*GFV2*m3*m4*p1p2*qp4)/D42 - 
                (64*GFV2*m3*m4*p1p2*qp4)/D3D4 - 
                (128*GFA2*p1p3*p2p3*qp4)/D3D4 - 
                (128*GFV2*p1p3*p2p3*qp4)/D3D4 + 
                (64*GFA2*p1p4*p2p3*qp4)/D3D4 + 
                (128*GFAGFV*p1p4*p2p3*qp4)/D3D4 + 
                (64*GFV2*p1p4*p2p3*qp4)/D3D4 + 
                (64*GFA2*p1p3*p2p4*qp4)/D3D4 - 
                (128*GFAGFV*p1p3*p2p4*qp4)/D3D4 + 
                (64*GFV2*p1p3*p2p4*qp4)/D3D4 - 
                (64*GFA2*p2p3*qp1*qp4)/D42 - 
                (128*GFAGFV*p2p3*qp1*qp4)/D42 - 
                (64*GFV2*p2p3*qp1*qp4)/D42 - 
                (64*GFA2*p1p3*qp2*qp4)/D42 + 
                (128*GFAGFV*p1p3*qp2*qp4)/D42 - 
                (64*GFV2*p1p3*qp2*qp4)/D42);
    
    return result;}
    
//*************************************************
// Compute cross section for fixed neutrino energy
//*************************************************


void ComputeCrossSection(){
    
    double varold, intold, sigma;
    double integral = 0.0;
    double var = 0.0;
    
    errorcounter=0;
    float progress = 0.0;
    int barWidth = 50;
    
    std::cout << " \n\n";
    
    for (int ii = 1; ii < 30000001; ii++){
      bool get_vertex=false;
      GenerateEvent(get_vertex);
      DetermineWeight();
      
      intold = integral;
      varold = var;
	    
      integral = (intold*(ii-1)+weight)/ii;
      var = ((varold+intold*intold)*(ii-1)+weight*weight)/ii-integral*integral;
      sigma = sqrt(var/ii);
      
      //progress bar
      if ( ii % 500000 == 0){
      progress = ii/30000000.0;
      std::cout << "computing cross section:  [";
      int pos = barWidth * progress;
      for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
      }
      std::cout << "]  " << int(progress * 100.0) << " %\r";
      std::cout.flush();
      }
      
    }

    crosssectionresult=integral;
    deltacrosssectionresult=sigma;
    
    if (errorcounter > 0){
    std::cout << " \n\n";
    std::cout << "Out of 30,000,000 weighted events, " << errorcounter << " events with unphysical kinematics were ignored \n\n";}
    
}
  
  
//**********************************
// Determine the weight of an event
//**********************************

void DetermineWeight(){
    
    if (zeroweight == 1){weight = 0.0;}
    else {
  
    // Set vector and axialvector couplings
    // note that the sign in the propagator corresponds to the unconventional metric of Lovseth, Radomski 
    
    if (model.compare("4F") == 0){
    GV = GVtot/2;  // factor 1/2 to agree with the normalization of the effective Hamiltonian of the paper
    GA = GAtot/2;} // factor 1/2 to agree with the normalization of the effective Hamiltonian of the paper
    else {  
    GV = GVSM - GP*GP*sqrt(2)/2/GF/(2*p1p2-MZP*MZP);
    GA = GASM;}
    
    // Jacobian factor (27) of Lovseth, Radomski, without the form factor 
    // (form factor is cancelled by the corresponding term in PLP)
    J = Twopi/8/eps1/M*D3/16/(-qp1)/sqrt(u4)*D4/8/qc/WC;
    
    // Define the squared matrix element
    if (anti == 0){PLP = SquaredMatrixElementPLP(GV, GA);}
    else {PLP = SquaredMatrixElementPLPanti(GV, GA);}
    
    Laa = 0.0;
    if (material.compare("proton") == 0 || material.compare("neutron") == 0){
    if (anti == 0){Laa = SquaredMatrixElementL(GV, GA);}
    else {Laa = SquaredMatrixElementLanti(GV, GA);}
    }
	      
    // Cross section (26) of Lovseth, Radomski
    if (material.compare("proton") == 0){
      deltaSigma=(GF*GF)*aem*aem*(Z*Z)*0.5/(2*pi*2*pi*2*pi*2*pi*2*pi*2*pi)/M/eps1*
      (4*PLP+4*M2*Laa*(1+q2/4/M2)*q2/4/M2*GMp(q2)*GMp(q2)/(GEp(q2)*GEp(q2)+q2/4/M2*GMp(q2)*GMp(q2)))/
      (q2*q2)*J*(0.389379e12)*Pauli(sqrt(q2*(1+q2/4/M2)));
    }
    else if (material.compare("neutron") == 0){
      deltaSigma=(GF*GF)*aem*aem*(Z*Z)*0.5/(2*pi*2*pi*2*pi*2*pi*2*pi*2*pi)/M/eps1*
      (4*PLP+4*M2*Laa*(1+q2/4/M2)*q2/4/M2*GMn(q2)*GMn(q2)/(GEn(q2)*GEn(q2)+q2/4/M2*GMn(q2)*GMn(q2)))/
      (q2*q2)*J*(0.389379e12)*Pauli(sqrt(q2*(1+q2/4/M2)));
    }
    else {
      deltaSigma=(GF*GF)*aem*aem*(Z*Z)*0.5/(2*pi*2*pi*2*pi*2*pi*2*pi*2*pi)/M/eps1*4*PLP/(q2*q2)*J*(0.389379e12);
    }
    
       
    // Define the weight of the event
    weight = deltaSigma*((u1max-u1min)*(u2max-u2min)*(u3max-u3min)*(u4max-u4min)*(u5max-u5min)*(u6max-u6min)*(u7max-u7min));}
    
    if (energy_type.compare("1") != 0){
      weight=weight*probability_list[bin]*length_probability_list;
      }
      
    if (flagkinematics == 1){
      weight = 0.0; errorcounter = errorcounter+1; 
    }
    
    if (is_nan(weight) == 1){
      weight = 0.0; errorcounter = errorcounter+1;
    }
   
    return;
}


//**********************************
// Determine the maximum weight
//**********************************

void FindMaxWeight(){
  
    maxweight = 0.0;
    averageweight = 0.0;
    
    float progress = 0.0;
    int barWidth = 50;

    std::cout << "\n\n";
    
    for(int n=1; n<30000001; n++){
      
      //progress bar
      if (n % 500000 == 0){
      progress = n/30000000.0;
      std::cout << "finding max weight:  [";
      int pos = barWidth * progress;
      for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
      }
      std::cout << "]  " << int(progress * 100.0) << " %\r";
      std::cout.flush();
      }
    bool get_vertex = false;
    GenerateEvent(get_vertex);
    DetermineWeight();
    
    averageweight = ((n-1)*averageweight+weight)/n;
    
    if (weight > maxweight){maxweight = weight;}}
          
    std::cout << "\n\n";
    std::cout << "To determine the maximum weight, 30,000,000 weighted events were computed \n\n";
    std::cout << "max_weight =  " << maxweight*2 << "\n";
    std::cout << "average_weight =  " << averageweight << "\n";
 
    return;
 
}

//*************************************************
// Read the file containing the neutrino flux distribution
//*************************************************

void ReadDistribution(){
    
    length_probability_list = distribution_list.size()/3;
    
    double Enumin;
    double Enumax;
    double probability;
    
    for(int n=0; n < length_probability_list; n++){
      
      Enumin=distribution_list[n*3];
      Enumax=distribution_list[n*3+1];
      probability=distribution_list[n*3+2];
      
      Enumin_list.push_back(Enumin);
      Enumax_list.push_back(Enumax);
      probability_list.push_back(probability);
      
    }
    
    return;
    
}

//*************************************************
// Load flux from a ROOT file
//*************************************************

void MakeFluxFromROOT(string filename, int target_pdg){
    
    TFile* infile = TFile::Open(filename.c_str());
    if (!infile || infile->IsZombie()) {
        std::cerr << "Error: Could not open ROOT file: " << filename << std::endl;
        return;
    }
    
    // Get the flux TTree
    TTree* flux_tree = (TTree*)infile->Get("flux");
    if (!flux_tree) {
        std::cerr << "Error: Could not find 'flux' TTree in ROOT file" << std::endl;
        infile->Close();
        return;
    }
    
    // Set up branch addresses
    Double_t wgt, E;
    Int_t pdg;
    Double_t vtxx, vtxy, vtxz;
    Double_t px, py, pz;
    
    flux_tree->SetBranchAddress("wgt", &wgt);
    flux_tree->SetBranchAddress("E", &E);
    flux_tree->SetBranchAddress("pdg", &pdg);
    flux_tree->SetBranchAddress("vtxx", &vtxx);
    flux_tree->SetBranchAddress("vtxy", &vtxy);
    flux_tree->SetBranchAddress("vtxz", &vtxz);
    flux_tree->SetBranchAddress("px", &px);
    flux_tree->SetBranchAddress("py", &py);
    flux_tree->SetBranchAddress("pz", &pz);
    
    // Find energy range and create histogram
    Long64_t nentries = flux_tree->GetEntries();
    std::cout << "Found " << nentries << " entries in flux tree" << std::endl;
    
    double min_E = 1e9;
    double max_E = 0.0;
    
    // First pass: find energy range for matching PDG
    for (Long64_t i = 0; i < nentries; i++) {
        flux_tree->GetEntry(i);
        if (pdg == target_pdg) {
            if (E < min_E) min_E = E;
            if (E > max_E) max_E = E;
        }
    }
    
    if (max_E <= min_E) {
        std::cerr << "Error: No entries found for PDG code " << target_pdg << std::endl;
        infile->Close();
        return;
    }
    
    std::cout << "Energy range: " << min_E << " to " << max_E << " GeV" << std::endl;
    
    // Create histogram with appropriate binning
    int nbins = 100;
    TH1D* flux_hist = new TH1D("flux_hist", "Neutrino Flux", nbins, min_E, max_E);
    
    // Second pass: fill histogram with weights for matching PDG
    for (Long64_t i = 0; i < nentries; i++) {
        flux_tree->GetEntry(i);
        if (pdg == target_pdg) {
            flux_hist->Fill(E, wgt);
        }
    }
    
    // Convert histogram to flux distribution format
    distribution_list.clear();
    Enumin_list.clear();
    Enumax_list.clear();
    probability_list.clear();
    
    double total_content = 0.0;
    
    // Calculate total content for normalization
    for (int i = 1; i <= nbins; i++) {
        total_content += flux_hist->GetBinContent(i);
    }
    
    if (total_content <= 0) {
        std::cerr << "Error: Histogram has no positive content for PDG " << target_pdg << std::endl;
        delete flux_hist;
        infile->Close();
        return;
    }
    
    // Fill distribution list
    for (int i = 1; i <= nbins; i++) {
        double emin = flux_hist->GetBinLowEdge(i);
        double emax = flux_hist->GetBinLowEdge(i) + flux_hist->GetBinWidth(i);
        double content = flux_hist->GetBinContent(i);
        double probability = content / total_content;
        
        if (probability > 0) {  // Only add bins with non-zero probability
            distribution_list.push_back(emin);
            distribution_list.push_back(emax);
            distribution_list.push_back(probability);
        }
    }
    
    delete flux_hist;
    infile->Close();
    
    std::cout << "Successfully loaded flux from ROOT file for PDG " << target_pdg << std::endl;
    std::cout << "Created distribution with " << distribution_list.size()/3 << " energy bins\n" << std::endl;
    
    // Now read the distribution as usual
    ReadDistribution();
    
    return;
    
}

void LoadFluxFromROOT(string filename, int target_pdg){
    
    TFile* infile = TFile::Open(filename.c_str());
    if (!infile || infile->IsZombie()) {
        std::cerr << "Error: Could not open ROOT file: " << filename << std::endl;
        return;
    }

    // Turn PDG code into name string to grab the right histogram
    string pdg_name;
    if (abs(target_pdg) == 12) pdg_name = "nue";
    else if (abs(target_pdg) == 14) pdg_name = "numu";
    else if (abs(target_pdg) == 16) pdg_name = "nutau";
    else {
        std::cerr << "Error: Unsupported PDG code " << target_pdg << std::endl;
        infile->Close();
        return;
    }
    if (target_pdg < 0) pdg_name = pdg_name+"bar";

    string key; 
    // Get the histogram
    if (filename.find("MCC9") != string::npos){
        key = "hE"+pdg_name+"_cv";
    }
    else if (filename.find("g4lbne") != string::npos){
        key = pdg_name+"_flux";
    }
    else {
        std::cerr << "Error: Unrecognized ROOT file format for flux" << std::endl;
        infile->Close();
        return;
    }

    TH1D* flux_hist = (TH1D*)infile->Get(key.c_str());
    if (!flux_hist) {
        std::cerr << "Error: Could not find histogram " << key << " in ROOT file" << std::endl;
        infile->Close();
        return;
    }

    // Convert histogram to flux distribution format
    distribution_list.clear();
    Enumin_list.clear();
    Enumax_list.clear();
    probability_list.clear();

    double total_content = 0.0;
    int nbins = flux_hist->GetNbinsX();

    // Calculate total content for normalization
    for (int i = 1; i <= nbins; i++) {
        total_content += flux_hist->GetBinContent(i);
    }
    if (total_content <= 0) {
        std::cerr << "Error: Histogram has no positive content for PDG " << target_pdg << std::endl;
        infile->Close();
        return;
    }

    // Fill distribution list
    for (int i = 1; i <= nbins; i++) {
        double emin = flux_hist->GetBinLowEdge(i);
        double emax = flux_hist->GetBinLowEdge(i) + flux_hist->GetBinWidth(i);
        double content = flux_hist->GetBinContent(i);
        double probability = content / total_content;
        if (probability > 0) {  // Only add bins with non-zero probability
            distribution_list.push_back(emin);
            distribution_list.push_back(emax);
            distribution_list.push_back(probability);
        }
    }
    infile->Close();

    std::cout << "Successfully loaded flux from ROOT file for PDG " << target_pdg << std::endl;
    std::cout << "Created distribution with " << distribution_list.size()/3 << " energy bins\n" << std::endl;

    // Now read the distribution as usual
    ReadDistribution();

    return;

}
//*************************************************
// Write file with the generated events
//*************************************************

void WriteEventFile(string filename){
  
    ofstream outfile;
    outfile.open(filename.c_str(), ios_base::trunc | ios_base::out | ios_base::in);
    
    // write the header
    outfile << "This is an event file created by TEG version 1.0 \n\n\n";
    outfile << "The format of the events is similar to that of MadGraph \n\n";
    outfile << "======================================== \n";
    outfile << "PDG code, in/out state, px, py, pz, E, m \n"; 
    outfile << "======================================== \n\n";   
    outfile << "PDG codes are       11 for a electron \n";
    outfile << "                   -11 for a positron \n";
    outfile << "                    12 for a electron neutrino \n";
    outfile << "                   -12 for a electron anti-neutrino \n";
    outfile << "                    13 for a muon \n";
    outfile << "                   -13 for a anti-muon \n";
    outfile << "                    14 for a muon neutrino \n";
    outfile << "                   -14 for a muon anti-neutrino \n";
    outfile << "                    15 for a tau \n";
    outfile << "                   -15 for a anti-tau \n";
    outfile << "                    16 for a tau neutrino \n";
    outfile << "                   -16 for a tau anti-neutrino \n";
    outfile << "                  2212 for a proton \n";
    outfile << "                  2112 for a neutron \n\n";
    outfile << "in/out state is    -1 for a particle in the initial state \n";
    outfile << "                    1 for a particle in the final state \n\n";
    outfile << "px is the particle momentum in x direction \n";
    outfile << "py is the particle momentum in y direction \n";
    outfile << "pz is the particle momentum in z direction (= beam direction) \n";
    outfile << "E is the particle Energy \n";
    outfile << "m is the particle mass \n\n";
    outfile << "(momenta, energies and masses are all in GeV) \n\n\n";
    
    int n = 0;

    while(n < 28*Nevents){
            n +=4; // Skipping the neutrino production vertex

            outfile << "\n<event>\n";
            outfile << PDG1 << " -1 ";
            for(int d=1;d<4;d++){
                outfile << std::fixed << std::setprecision(8) << ::data[n+d] << " ";
            }
            outfile << std::fixed << std::setprecision(8) << ::data[n] << " ";
            outfile << std::fixed << std::setprecision(8) << "0.0\n";
            n += 4;

            if (material.compare("Ar") == 0){
            outfile << 1000180400 << "  -1 ";
                // At rest argon nucleus
                for(int d=1;d<4;d++){
                    outfile << std::fixed << std::setprecision(8) << 0. << " ";
                }
                // Argon rest mass
                outfile << std::fixed << std::setprecision(8) << MArgon << " ";
                outfile << std::fixed << std::setprecision(8) << MArgon << "\n";
            }   

            n+=4; // Skipping the vertex position 
            
            outfile << PDG2 << "  1 ";
            for(int d=1;d<4;d++){
                outfile << std::fixed << std::setprecision(8) << ::data[n+d] << " ";
            }
            outfile << std::fixed << std::setprecision(8) << ::data[n] << " ";
            outfile << std::fixed << std::setprecision(8) << "0.0\n";
            n += 4;
            
            outfile << PDG4 << "  1 ";
            for(int d=1;d<4;d++){
                outfile << std::fixed << std::setprecision(8) << ::data[n+d] << " ";
            }
            outfile << std::fixed << std::setprecision(8) << ::data[n] << " ";
            outfile << std::fixed << std::setprecision(8) << m4 << "\n";
            n += 4;
            
            outfile << PDG3 << "  1 ";
            for(int d=1;d<4;d++){
                outfile << std::fixed << std::setprecision(8) << ::data[n+d] << " ";
            }
            outfile << std::fixed << std::setprecision(8) << ::data[n] << " ";
            outfile << std::fixed << std::setprecision(8) << m3 << "\n";
            n += 4;
	    	    
            if (material.compare("proton") == 0){
            outfile << 2212 << "  1 ";
                for(int d=1;d<4;d++){
                    outfile << std::fixed << std::setprecision(8) << ::data[n+d] << " ";
                }
                outfile << std::fixed << std::setprecision(8) << ::data[n] << " ";
                outfile << std::fixed << std::setprecision(8) << M << "\n";
            }
            else if (material.compare("neutron") == 0){
            outfile << 2112 << "  1 ";
                for(int d=1;d<4;d++){
                    outfile << std::fixed << std::setprecision(8) << ::data[n+d] << " ";
                }
                outfile << std::fixed << std::setprecision(8) << ::data[n] << " ";
                outfile << std::fixed << std::setprecision(8) << M << "\n";
            }
            else if (material.compare("Ar") == 0){
            outfile << 1000180400 << "  1 ";
                for(int d=1;d<4;d++){
                    outfile << std::fixed << std::setprecision(8) << ::data[n+d] << " ";
                }
                outfile << std::fixed << std::setprecision(8) << ::data[n] << " ";
                outfile << std::fixed << std::setprecision(8) << M << "\n";
            }   
            n += 8; // Skipping the final state nucleus
            outfile << "</event>";
    }

    outfile.close();
}

void WriteEventHepMC3(string filename){
  
    ofstream outfile;
    outfile.open(filename.c_str(), ios_base::trunc | ios_base::out | ios_base::in);
    
    // write the header
    outfile << "HepMC::Version 3.02.05\n";
    outfile << "HepMC::Asciiv3-START_EVENT_LISTING\n";
    
    int n = 0;
    int event_num = 0;

    while(n < 28*Nevents){
        n+=4 ; // Skipping the neutrino production vertex

        int particle = 0;
        int vertex = 0;

        outfile << "E " << event_num << " 1 6\n";
        outfile << "U GEV CM\n";

        particle +=1;
        outfile << "P " << particle << " 0 " << std::fixed  << PDG1  << " ";
        for(int d=1;d<4;d++){
            outfile << std::fixed << std::setprecision(8) << ::data[n+d] << " ";
        }
        outfile << std::fixed << std::setprecision(8) << ::data[n] << " ";
        outfile << std::fixed << std::setprecision(8) << "0.0" << " 4\n";
        n += 4;

        int target_pdg = -1;
        double target_mass = -1.;
        if (material.compare("Ar") == 0){
            target_pdg = 1000180400;
            target_mass = MArgon;
        }   
        else if (material.compare("neutron") == 0){
            target_pdg = 2112;
            target_mass = Mneutron;
        }
        else if (material.compare("proton") == 0){
            target_pdg == 2212;
            target_mass = Mproton;
        }
        
        particle +=1;
        outfile << "P " << std::fixed << particle << std::fixed  << " 0 " << std::fixed  << target_pdg  << " ";
        // At rest argon nucleus
        for(int d=1;d<4;d++){
            outfile << std::fixed << std::setprecision(8) << 0. << " ";
        }
        // Argon rest mass
        outfile << std::fixed << std::setprecision(8) << target_mass << " ";
        outfile << std::fixed << std::setprecision(8) << target_mass << " 4\n";
        // No need to increment n here since we didn't read from data

        vertex += -1;
        outfile << "V " << std::fixed << vertex << std::fixed  << " 0 [1,2] @ " << std::fixed;
        for(int d=1;d<4;d++){
            outfile << std::fixed << std::setprecision(8) << ::data[n+d] << " ";
        }
        outfile << std::fixed << std::setprecision(8) << ::data[n] << "\n";
        n += 4;

        particle +=1;
        outfile << "P " << particle << " -1 " << std::fixed  << PDG2  << " ";
        for(int d=1;d<4;d++){
            outfile << std::fixed << std::setprecision(8) << ::data[n+d] << " ";
        }
        outfile << std::fixed << std::setprecision(8) << ::data[n] << " ";
        outfile << std::fixed << std::setprecision(8) << "0.0" <<" 1\n";
        n += 4;

        particle +=1;
        outfile << "P " << particle << " -1 " << std::fixed  << PDG4  << " ";
        for(int d=1;d<4;d++){
            outfile << std::fixed << std::setprecision(8) << ::data[n+d] << " ";
        }
        outfile << std::fixed << std::setprecision(8) << ::data[n] << " ";
        outfile << std::fixed << std::setprecision(8) << m4 << " 1\n";
        n += 4;

        particle +=1;
        outfile << "P " << particle << " -1 " << std::fixed  << PDG3  << " ";
        for(int d=1;d<4;d++){
            outfile << std::fixed << std::setprecision(8) << ::data[n+d] << " ";
        }
        outfile << std::fixed << std::setprecision(8) << ::data[n] << " ";
        outfile << std::fixed << std::setprecision(8) << m3 << " 1\n";
        n += 4;

        particle +=1;
        outfile << "P " << particle << " -1 " << std::fixed  << target_pdg  << " ";
        for(int d=1;d<4;d++){
            outfile << std::fixed << std::setprecision(8) << ::data[n+d] << " ";
        }
        outfile << std::fixed << std::setprecision(8) << ::data[n] << " ";
        outfile << std::fixed << std::setprecision(8) << target_mass << " 1\n";
        n += 4;

        event_num++;
    }
    outfile << "HepMC::Asciiv3-END_EVENT_LISTING";
    outfile.close();
}

void WriteEventHepevt(string filename){
  
    ofstream outfile;
    outfile.open(filename.c_str(), ios_base::trunc | ios_base::out | ios_base::in);
    
    int n = 0;
    int event_num = 1;

    while(n < 28*Nevents){
        
        int particle = 0;
        int vertex = 0;

        outfile << event_num << " 6\n";
        // HEPEVT format:
        // <event no.> <number of particles> 
        // <status (1 for final state, 2 for intermediate, 3 for beam particle) > <PDG ID> <parent1> <parent2> <child1> <child2> <px> <py> <pz> <E> <mass> <prod x> <prod y> <prod z> <prod t>
        
        // neutrino production vertex
        double nu_t = ::data[n];
        double nu_x = ::data[n+1];
        double nu_y = ::data[n+2];
        double nu_z = ::data[n+3];
        // Interaction vertex
        double int_t = ::data[n+8];
        double int_x = ::data[n+9];
        double int_y = ::data[n+10];
        double int_z = ::data[n+11];
        
        n +=4 ; // Skipping the neutrino production vertex as we have stored it already

        outfile << 3 << " " << PDG1 << " 0 0 3 6 ";
        for(int d=1;d<4;d++){
            outfile << std::fixed << std::setprecision(8) << ::data[n+d] << " ";
        }
        outfile << std::fixed << std::setprecision(8) << ::data[n] << " ";
        outfile << std::fixed << std::setprecision(8) << "0 " << nu_x << " " << nu_y << " " << nu_z << " " << nu_t << "\n";
        n += 8; // Skipping to interaction vertex

        int target_pdg = -1;
        double target_mass = -1.;
        if (material.compare("Ar") == 0){
            target_pdg = 1000180400;
            target_mass = MArgon;
        }   
        else if (material.compare("neutron") == 0){
            target_pdg = 2112;
            target_mass = Mneutron;
        }
        else if (material.compare("proton") == 0){
            target_pdg == 2212;
            target_mass = Mproton;
        }
        outfile << 3 << " " << target_pdg << " 0 0 3 6 ";
        // At rest argon nucleus
        for(int d=1;d<4;d++){
            outfile << std::fixed << std::setprecision(8) << 0. << " ";
        }
        // Argon rest mass
        outfile << std::fixed << std::setprecision(8) << target_mass << " ";
        outfile << std::fixed << std::setprecision(8) << target_mass << " " << int_x << " " << int_y << " " << int_z << " " << 0 << "\n";
    

        outfile << 1 << " " << PDG2 << " 1 2 0 0 ";
        for(int d=1;d<4;d++){
            outfile << std::fixed << std::setprecision(8) << ::data[n+d] << " ";
        }
        outfile << std::fixed << std::setprecision(8) << ::data[n] << " ";
        outfile << std::fixed << std::setprecision(8) << "0 " << int_x << " " << int_y << " " << int_z << " " << int_t << "\n";
        n += 4;

        outfile << 1 << " " << PDG4 << " 1 2 0 0 ";
        for(int d=1;d<4;d++){
            outfile << std::fixed << std::setprecision(8) << ::data[n+d] << " ";
        }
        outfile << std::fixed << std::setprecision(8) << ::data[n] << " ";
        outfile << std::fixed << std::setprecision(8) << m4 << " " << int_x << " " << int_y << " " << int_z << " " << int_t << "\n";
        n += 4;

        outfile << 1 << " " << PDG3 << " 1 2 0 0 ";
        for(int d=1;d<4;d++){
            outfile << std::fixed << std::setprecision(8) << ::data[n+d] << " ";
        }
        outfile << std::fixed << std::setprecision(8) << ::data[n] << " ";
        outfile << std::fixed << std::setprecision(8) << m3 << " " << int_x << " " << int_y << " " << int_z << " " << int_t << "\n";
        n += 4;

        // Write argon as an intermediate final state, as g4 doesnt know what to do with it
        outfile << 1 << " " << target_pdg << " 1 2 0 0 ";
        for(int d=1;d<4;d++){
            outfile << std::fixed << std::setprecision(8) << ::data[n+d] << " ";
        }
        outfile << std::fixed << std::setprecision(8) << ::data[n] << " ";
        outfile << std::fixed << std::setprecision(8) << target_mass << " " << int_x << " " << int_y << " " << int_z << " " << int_t << "\n";
        n += 4;
        event_num++;
    }
    outfile.close();
}

bool exists_test (const std::string& name) {
    if (FILE *file = fopen(name.c_str(), "r")) {
        fclose(file);
        return true;
    } else {
        return false;
    }   
}

void WriteXsecTempFile(string filename){
  
    ofstream outfile;
    if(exists_test(filename)){
        std::cout << "Appending to existing cross section temporary file \n";
        outfile.open(filename, std::ios_base::app); // append instead of overwrite     
    }
    else {
        outfile.open(filename.c_str(), ios_base::trunc | ios_base::out | ios_base::in);
        outfile << "# This is a temporary file created by TEG version 2.0 \n";
        outfile << "# It contains the computed cross section for given neutrino energy \n\n";
        outfile << "# Neutrino energy (GeV)    Cross section (fb)    Uncertainty (fb) \n";
    }
    // else {
    //     std::cout << "Appending to existing cross section temporary file \n";
    //     outfile.seekp(0, ios::end);
    // }
    outfile << std::fixed << std::setprecision(6) << Enu << "    ";
    outfile << std::scientific << std::setprecision(6) << crosssectionresult << "    ";
    outfile << std::scientific << std::setprecision(6) << deltacrosssectionresult << "\n";
    outfile.close();
}

//*************************************************
// Sample the E-z distribution to get the 
// production vertex and momentum direction of the 
// incoming neutrino
//*************************************************

void LoadEvsZhist(string filename){
    // Loading the relevant histogram once so we don't have to do it for every event
    std::cout << "Loading E vs Z histogram from file: " << filename << std::endl;
    TFile* infile = TFile::Open(filename.c_str());
    if (!infile || infile->IsZombie()) {
        std::cerr << "Error: Could not open ROOT file: " << filename << ". Saving empty histogram." << std::endl;
        return;
    }
    string pdg_name;
    if (abs(PDG1) == 12) pdg_name = "nue";
    else if (abs(PDG1) == 14) pdg_name = "numu";
    else if (abs(PDG1) == 16) pdg_name = "nutau";
    else {
        std::cerr << "Error: Could not find a histogram for PDG " << PDG1 << " in ROOT file. Saving empty histogram." << std::endl;
        infile->Close();
        return;
    }
    if (PDG1 < 0) pdg_name = pdg_name+"bar";
    
    string key = "hE_vs_z_"+pdg_name;
    TH2D* hist = (TH2D*)infile->Get(key.c_str());
    if (!hist) {
        std::cerr << "Error: Could not find histogram " << key << " in ROOT file. Saving empty histogram." << std::endl;
        infile->Close();
        return;
    }
    
    hist->SetDirectory(0); // Detach from file ownership
    Evsz_hist = *hist;
    infile->Close();
    std::cout << "Successfully loaded E vs Z histogram for PDG " << PDG1 << std::endl;
    return;
} 

vector<double> SampleEvsZhist(double E_nu, string filename){

    vector<double> production_vertex(3); // x,y,z
    
    TH2D* hist = &Evsz_hist;
    if (!hist) {
        std::cerr << "Error: E vs Z histogram not loaded. Returning empty vector." << std::endl;
        return production_vertex;
    }

    // Find the bin corresponding to E_nu
    int ebin = hist->GetXaxis()->FindBin(E_nu);
    // Project the histogram onto the z-axis for the given energy bin
    TH1D* zproj = hist->ProjectionY("zproj", ebin, ebin);
    // Normalize the projection to create a probability distribution
    zproj->Scale(1.0 / zproj->Integral());
    // Sample a z value from the distribution
    double rand = realdistribution(generator);
    double cumulative = 0.0;
    double sampled_z = 0.0;
    for (int bin = 1; bin <= zproj->GetNbinsX(); bin++) {
        cumulative += zproj->GetBinContent(bin);
        if (rand < cumulative) {
            sampled_z = zproj->GetBinCenter(bin);
            break;
        }
    }
    if (bin == zproj->GetNbinsX() + 1) {
        sampled_z = zproj->GetBinCenter(zproj->GetNbinsX());
    }

    // Get random x and y within the beam radius 
    double r = radius_decay_pipe * sqrt(realdistribution(generator));
    double theta = 2 * pi * realdistribution(generator);
    double sampled_x = r * cos(theta);
    double sampled_y = r * sin(theta);

    production_vertex[0] = sampled_x + BNB_x_offset;
    production_vertex[1] = sampled_y + BNB_y_offset;
    production_vertex[2] = sampled_z + BNB_z_offset;

    if (debug_mode){
        // write production vertex to a temp file for debugging
        ofstream outfile;
        string temp_filename = "production_vertex_temp.txt";
        if(exists_test(temp_filename)){
            outfile.open(temp_filename, std::ios_base::app); // append instead of overwrite     
        }
        else {
            outfile.open(temp_filename.c_str(), ios_base::trunc | ios_base::out | ios_base::in);
            outfile << "# This is a temporary file created by TEG version 2.0 \n";
            outfile << "# It contains the sampled production vertex (x,y,z) for given neutrino energy \n\n";
            outfile << "# x (m)    y (m)    z (m) \n";
        }
        outfile << std::fixed << std::setprecision(6) << production_vertex[0] << "    ";
        outfile << std::fixed << std::setprecision(6) << production_vertex[1] << "    ";
        outfile << std::fixed << std::setprecision(6) << production_vertex[2] << "\n";
        outfile.close();
    }
    return production_vertex;
}

//*************************************************
// Generate unweighted events
//*************************************************
   
void GenerateEvents(){
      
      double random;
      int i = 0;
      
      errorcounter=0;
      eventcounter=0;
      reweightcounter=0; 
      
      float progress = 0.0;
      int barWidth = 50;
      
      std::cout << "\n\n";
      LoadEvsZhist(flux_file);

      while(i<Nevents){
		
      //progress bar
      if (eventcounter % 500000 == 0){
      progress = (i + 0.1)/Nevents;
      std::cout << "generating events:  [";
      int pos = barWidth * progress;
      for (int j = 0; j < barWidth; ++j) {
        if (j < pos) std::cout << "=";
        else if (j == pos) std::cout << ">";
        else std::cout << " ";
      }
      std::cout << "]  " << int(progress * 100.0) << " %\r";
      std::cout.flush();
      }

      eventcounter=eventcounter+1;
      bool get_vertex = true;
      GenerateEvent(get_vertex);
      DetermineWeight();
      
      if(weight>2*maxweight){
	 reweightcounter=reweightcounter+1;
//	 std::cout << "\n\n";
//	 std::cout << "Warning: weight of an event is larger than the maximum weight \n";
//	 std::cout << "weight / maxweight = " << weight/2/maxweight << "\n\n";
      }
      
      random = realdistribution(generator);
      
      if(random < weight/maxweight/2.0){
	
        ::data[28*i] = event[0][0];
        ::data[28*i+1] = event[0][1];
        ::data[28*i+2] = event[0][2];
        ::data[28*i+3] = event[0][3];
        ::data[28*i+4] = event[1][0];
        ::data[28*i+5] = event[1][1];
        ::data[28*i+6] = event[1][2];
        ::data[28*i+7] = event[1][3];
        ::data[28*i+8] = event[2][0];
        ::data[28*i+9] = event[2][1];
        ::data[28*i+10] = event[2][2];
        ::data[28*i+11] = event[2][3];
        ::data[28*i+12] = event[3][0];
        ::data[28*i+13] = event[3][1];
        ::data[28*i+14] = event[3][2];
        ::data[28*i+15] = event[3][3];
        ::data[28*i+16] = event[4][0];
        ::data[28*i+17] = event[4][1];
        ::data[28*i+18] = event[4][2];
        ::data[28*i+19] = event[4][3];
        ::data[28*i+20] = event[5][0];
        ::data[28*i+21] = event[5][1];
        ::data[28*i+22] = event[5][2];
        ::data[28*i+23] = event[5][3];
        ::data[28*i+24] = event[6][0];
        ::data[28*i+25] = event[6][1];
        ::data[28*i+26] = event[6][2];
        ::data[28*i+27] = event[6][3];
	
	i++;}}
    
    std::cout << "\n\n";
    std::cout << "To generate the " << Nevents << " events, " << eventcounter << " weighted events were computed \n\n";
    
    if (errorcounter > 0){
    std::cout << "Out of the " << eventcounter << " weighted events, " << errorcounter << " events with unphysical kinematics were ignored \n\n";}
    
    if (reweightcounter > 0){
    std::cout << "Out of the " << eventcounter << " weighted events, " << reweightcounter << " events had weights larger than the maximum weight determined earlier \n\n";}
	
    return;

}