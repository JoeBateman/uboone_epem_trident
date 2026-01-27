float KEfromMom(float mom[3],float mass){
  float energy = sqrt(mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2] + mass*mass);
  float KE = energy - mass;
  return KE;
  
}

void make_2d_flux_hists(){
  int nfiles = 500; // Number of flux files to iterate over
  string particleNames[4] = {"#nu_{#mu}","#bar{#nu}_{#mu}","#nu_{e}","#bar{#nu}_{e}"};
  string particleSymbols[4] = {"numu","numubar","nue","nuebar"};
  int pdgID[4] = {14,-14,12,-12};

  int nHists = sizeof(particleNames)/sizeof(particleNames[0]); // Number of histograms to create
  TH1D* histsE[nHists]; // Energy histograms
  TH1D* histsz[nHists]; // z histograms
  TH2D* hists_Ez[nHists]; // E vs z histograms

  double cumulativePOT = 0.0;

  for(int iHist=0;iHist<nHists;iHist++){
      string histName = "hE"+particleSymbols[iHist]+"_cv";
      TH1D* temp = new TH1D(histName.c_str(),histName.c_str(),200,0,10); // 200 bins from 0 to 10 GeV (50 MeV each)
      temp->GetXaxis()->SetTitle("Energy (GeV)");
      temp->GetYaxis()->SetTitle("#nu / bin / POT / cm^{2}"); // Replace X with actual POT later
      histsE[iHist] = temp;

      string histName_z = "hz_"+particleSymbols[iHist]+"_cv";
      TH1D* temp_z = new TH1D(histName_z.c_str(),histName_z.c_str(),50,-100,5000); // 50 bins from -100 to 5000 cm
      temp_z->GetXaxis()->SetTitle("z (cm)");
      temp_z->GetYaxis()->SetTitle("#nu / bin / POT / cm^{2}"); // Replace X with actual POT later
      histsz[iHist] = temp_z;

      string histName_Ez = "hE_vs_z_"+particleSymbols[iHist];
      TH2D* temp_Ez = new TH2D(histName_Ez.c_str(),histName_Ez.c_str(),200,0,10,50,-100,5000);
      temp_Ez->GetXaxis()->SetTitle("Energy (GeV)");
      temp_Ez->GetYaxis()->SetTitle("z (cm)");
      hists_Ez[iHist] = temp_Ez;
    }

  TH1D * pot = new TH1D("POT","POT",1,0,1); // Dummy histogram to store POT info
  pot->SetBinContent(1, cumulativePOT);

  //Loop through files
  string basePath = "/pnfs/uboone/persistent/uboonebeam/bnb_gsimple/bnb_gsimple_fluxes_01.09.2019_463/converted_beammc_wincorr_";
  for(int iFile=0;iFile<nfiles;iFile++){
    string fileIndex = to_string(iFile);
    int n_zero = 4;
    fileIndex = string(n_zero - std::min(n_zero, (int)fileIndex.length()), '0') + fileIndex;
    string inFile = basePath + fileIndex + ".root";
    std::cout << "Processing file: " << inFile << std::endl;
    
    TFile* f = TFile::Open(inFile.c_str());
    
    // Add after TFile::Open
    if (!f || f->IsZombie()) {
      std::cout << "Error: Could not open file " << inFile << std::endl;
      continue;
    }

    TTree * t = (TTree*) f->Get("flux");
    TTree * meta = (TTree*) f->Get("meta");

    double filePOT;
    meta->SetBranchAddress("meta.protons",&filePOT);
    meta->GetEntry(0);
    cumulativePOT += filePOT;

    if (!t) {
      std::cout << "Error: Could not find 'flux' tree in " << inFile << std::endl;
      f->Close();
      continue;
    }

    Long64_t n_entries = t->GetEntries();
    int pdg;
    double E_nu;
    double wgt;
    double vtx_z;

    t->SetBranchAddress("numi.vz",&vtx_z);
    t->SetBranchAddress("entry.pdg",&pdg);
    t->SetBranchAddress("entry.E",&E_nu);
    t->SetBranchAddress("entry.wgt",&wgt);

    std::cout << "Number of entries in tree: " << n_entries << std::endl;
    int n_particles = 0;
    for(int iEntry = 0; iEntry<n_entries;iEntry++){
      t->GetEntry(iEntry);
        //Find the matching particle ID histogram
      for(int iFind=0;iFind<nHists;iFind++){
        if(pdg==pdgID[iFind]){
          histsE[iFind]->Fill(E_nu,wgt);
          histsz[iFind]->Fill(vtx_z,wgt);
          hists_Ez[iFind]->Fill(E_nu,vtx_z,wgt);
          n_particles++;
        }
      }
    }
    std::cout << "Grabbed a total of " << n_particles << " from this file." << std::endl;
    f->Close();
  }
  
  std::cout << "Cumulative POT: " << cumulativePOT << std::endl;

  


  TFile* outFile = new TFile("/exp/uboone/app/users/jbateman/workdir/DarkNews/Trident/data/flux/bnb/MCC9_FluxHist_volTPCActive_w2D_hists_alt.root","RECREATE");

  // Normalize histograms to per POT per cm2

  // 256.35*233.
  double activeVolFace = 256.35*233.; // in cm^2
  double normFactor = 1 / (cumulativePOT * 256.35*233.) ;
  for(int iOut=0;iOut<nHists;iOut++){
    histsE[iOut]->Scale(normFactor);
    histsz[iOut]->Scale(normFactor);
    hists_Ez[iOut]->Scale(normFactor);
    histsE[iOut]->Write();
    histsz[iOut]->Write();
    hists_Ez[iOut]->Write();
  }
  pot->SetBinContent(1, cumulativePOT);
  pot->Write();
  outFile->Close();  // Add this line
}


// void plot_Parents(){
  

  
//   //Set up variabls for TTree. Largest array in files is 30 long, so make sure they are longer
//   int nmult;
//   int id[40];
//   float mom[40][3];
//   float wgt[40];
  
//   //Loop through files
//   for(int iFile=0;iFile<1000;iFile++){

//     //Create filename string, open file and get tree
//     string inFile = "/exp/sbnd/data/users/jpaton/G4BNB_ProductionNtuple/NuBeam_production_Production_BooNE_50m_I174000A_"+to_string(iFile)+".dk2nu.root";
//     TFile* f = TFile::Open(inFile.c_str());
    
//     // Add after TFile::Open
//     if (!f || f->IsZombie()) {
//       cout << "Error: Could not open file " << inFile << endl;
//       continue;
//     }


//     TTree * t = (TTree*) f->Get("production");
    
//     // Add after getting the tree
//     if (!t) {
//       cout << "Error: Could not find 'production' tree in " << inFile << endl;
//       f->Close();
//       continue;
//     }

//     //Set branch addresses and get number of entries
//     t->SetBranchAddress("nmult",&nmult);
//     t->SetBranchAddress("id",&id);
//     t->SetBranchAddress("mom",&mom);
//     t->SetBranchAddress("wgt",&wgt);
//     Long64_t n_entries = t->GetEntries();

//     //Loop through entries
//     for(int iEntry = 0; iEntry<n_entries;iEntry++){
//       t->GetEntry(iEntry);
      
//       //Loop through the length of the array (nmult)
//       for(int i=0;i<nmult;i++){     
	
// 	//Find the matching particle ID histogram
// 	for(int iFind=0;iFind<nHists;iFind++){
// 	  if(pdg[i]==pdgID[iFind]){
	    
// 	    //Calculate KE and fill relevant histogram
// 	    float KE = KEfromMom(mom[i],particleMasses[iFind]);
// 	    hists[iFind]->Fill(KE,wgt[i]);
	    
// 	    //Fill 2D histogram with momentum and
// 	    float MomMag = sqrt(mom[i][0]*mom[i][0] + mom[i][1]*mom[i][1] + mom[i][2]*mom[i][2]);
// 	    float cosTheta = mom[i][2]/MomMag;
// 	    histsMom[iFind]->Fill(MomMag,cosTheta,wgt[i]);
	    
// 	  }
// 	}
//       }
      
//       //Fill the ID with 0s to avoid double counting (arrays do not fully reset for each entry)
//       for(int iNull=0;iNull<40;iNull++){
// 	id[iNull]=0;
//       }
//     }
//   }
//   TFile* outFile = new TFile("pi0_eta_flux.root","RECREATE");
//   for(int iOut=0;iOut<nHists;iOut++){
//     hists[iOut]->Write();
//     histsMom[iOut]->Write();
//   }
//   outFile->Close();  // Add this line
// }