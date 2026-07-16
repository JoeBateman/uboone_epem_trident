float KEfromMom(float mom[3],float mass){
  float energy = sqrt(mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2] + mass*mass);
  float KE = energy - mass;
  return KE;
  
}

bool IntersectsDetector(double vx, double vy, double vz,
                         double px, double py, double pz)
{
  // Detector box (detector coordinates, cm)
  static const double box_min[3] = {   0.0, -116.0,    0.0};
  static const double box_max[3] = { 256.0,  116.0, 1036.0};

  // Beam -> detector offset (cm)
  static const double offset[3] = {128.175, -0.93, -46336.3525};

  // Beam -> detector rotation (identity here, but written generally)
  static const double rot[3][3] = {
    {1.0, 0.0, 0.0},
    {0.0, 1.0, 0.0},
    {0.0, 0.0, 1.0}
  };

  // Rotate + translate vertex; rotate-only the momentum direction
  double pos[3], dir[3];
  for (int i = 0; i < 3; ++i) {
    pos[i] = rot[i][0]*vx + rot[i][1]*vy + rot[i][2]*vz + offset[i];
    dir[i] = rot[i][0]*px + rot[i][1]*py + rot[i][2]*pz;
  }

  double tmin = -std::numeric_limits<double>::infinity();
  double tmax =  std::numeric_limits<double>::infinity();

  for (int i = 0; i < 3; ++i) {
    if (dir[i] == 0.0) {
      // Ray parallel to this pair of faces: must already lie within the slab
      if (pos[i] < box_min[i] || pos[i] > box_max[i]) return false;
      continue;
    }
    double t1 = (box_min[i] - pos[i]) / dir[i];
    double t2 = (box_max[i] - pos[i]) / dir[i];
    if (t1 > t2) std::swap(t1, t2);

    tmin = std::max(tmin, t1);
    tmax = std::min(tmax, t2);
    if (tmax < tmin) return false;
  }

  // Require the intersection to be ahead of the vertex (forward ray, not full line)
  tmin = std::max(tmin, 0.0);
  return tmax >= tmin && tmax >= 0.0;
}

//Loop through files
string basePath = "/pnfs/uboone/persistent/uboonebeam/bnb_gsimple/bnb_gsimple_fluxes_01.09.2019_463/converted_beammc_wincorr_";
string outFileName = "/exp/uboone/app/users/jbateman/workdir/DarkNews/Trident/data/flux/bnb/MCC9_FluxHist_volTPCActive_w2D_hists_test.root";

// string basePath = "/pnfs/uboone/persistent/uboonebeam/numi_gsimple/numi_gsimple_fluxes_12.17.2015_470/gsimple_microboone-numi_mn000z200i_rp11_bs1.1_pnut_lowth_f112c0f093bbird_";
// string outFileName = "/exp/uboone/app/users/jbateman/workdir/DarkNews/Trident/data/flux/numi/MCC9_FluxHist_volTPCActive_w2D_hists.root";

void make_2d_flux_hists_BNB(){

  int nfiles = 0;
  // std::filesystem::path dir = "/pnfs/uboone/persistent/uboonebeam/bnb_gsimple/bnb_gsimple_fluxes_01.09.2019_463/";
  // std::string prefix = "converted_beammc_wincorr_";
  
  std::filesystem::path dir = basePath.substr(0, basePath.find_last_of("/"));
  string prefix = basePath.substr(basePath.find_last_of("/") + 1);
  std::cout << "Searching in directory: " << dir << std::endl;
  std::cout << "Looking for files with prefix: " << prefix << std::endl;


  for (const auto& entry : std::filesystem::directory_iterator(dir)) {
      if (!entry.is_regular_file()) continue;

      std::string name = entry.path().filename().string();
      if (name.rfind(prefix, 0) == 0 && entry.path().extension() == ".root") {
          ++nfiles;
      }
  }


  std::cout << "Found " << nfiles << " files matching the pattern." << std::endl;  

  string particleNames[4] = {"#nu_{#mu}","#bar{#nu}_{#mu}","#nu_{e}","#bar{#nu}_{e}"};
  string particleSymbols[4] = {"numu","numubar","nue","nuebar"};
  int pdgID[4] = {14,-14,12,-12};

  int nHists = sizeof(particleNames)/sizeof(particleNames[0]); // Number of histograms to create
  TH1D* histsE[nHists]; // Energy histograms
  TH1D* histsz[nHists]; // z histograms
  TH2D* hists_Ez[nHists]; // E vs z histograms

  double cumulativePOT = 0.0;

  for(int iHist=0;iHist<nHists;iHist++){
      string histName = "hE_"+particleSymbols[iHist]+"_cv";
      TH1D* temp = new TH1D(histName.c_str(),histName.c_str(),400,0,10); // 400 bins from 0 to 10 GeV (25 MeV each)
      temp->GetXaxis()->SetTitle("Energy (GeV)");
      temp->GetYaxis()->SetTitle("#nu / bin / POT / cm^{2}"); // Replace X with actual POT later
      histsE[iHist] = temp;

      string histName_z = "hz_"+particleSymbols[iHist]+"_cv";
      // TH1D* temp_z = new TH1D(histName_z.c_str(),histName_z.c_str(),50,-100,5000); // 50 bins from -100 to 5000 cm (BNB)
      TH1D* temp_z = new TH1D(histName_z.c_str(),histName_z.c_str(),400,-1000,75000); // 400 bins from -1000 to 75000 cm (NuMI)


      temp_z->GetXaxis()->SetTitle("z (cm)");
      temp_z->GetYaxis()->SetTitle("#nu / bin / POT / cm^{2}"); // Replace X with actual POT later
      histsz[iHist] = temp_z;

      string histName_Ez = "hE_vs_z_"+particleSymbols[iHist]+"_cv";
      // TH2D* temp_Ez = new TH2D(histName_Ez.c_str(),histName_Ez.c_str(),200,0,10,50,-100,5000); // 200 bins from 0 to 10 GeV and 50 bins from -100 to 5000 cm (BNB)
      TH2D* temp_Ez = new TH2D(histName_Ez.c_str(),histName_Ez.c_str(),400,0,10,400,-1000,75000); // 400 bins from 0 to 10 GeV and 400 bins from -1000 to 75000 cm (NuMI)
      temp_Ez->GetXaxis()->SetTitle("Energy (GeV)");
      temp_Ez->GetYaxis()->SetTitle("z (cm)");
      hists_Ez[iHist] = temp_Ez;
    }

  TH1D * pot = new TH1D("POT","POT",1,0,1); // Dummy histogram to store POT info
  pot->SetBinContent(1, cumulativePOT);

  
  std::cout << "Starting to loop over " << nfiles << " files..." << std::endl;
  int N_loaded_hists = 0;
  for(int iFile=0;iFile<nfiles;iFile++){
    string fileIndex = to_string(iFile);
    int n_zero = 4;
    fileIndex = string(n_zero - std::min(n_zero, (int)fileIndex.length()), '0') + fileIndex;
    
    string inFile = basePath + fileIndex + ".root";
    // string inFile = basePath + fileIndex +"_0" + fileIndex + ".root"; // Adjusted for numi file naming convention

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
    double vtx_x;
    double vtx_y;
    double vtx_z;
    double px;
    double py;
    double pz;

    t->SetBranchAddress("numi.vx",&vtx_x);
    t->SetBranchAddress("numi.vy",&vtx_y);
    t->SetBranchAddress("numi.vz",&vtx_z);
    t->SetBranchAddress("entry.pdg",&pdg);
    t->SetBranchAddress("entry.px",&px);
    t->SetBranchAddress("entry.py",&py);
    t->SetBranchAddress("entry.pz",&pz);
    t->SetBranchAddress("entry.E",&E_nu);
    t->SetBranchAddress("entry.wgt",&wgt);

    std::cout << "Number of entries in tree: " << n_entries << std::endl;
    int n_particles = 0;
    for(int iEntry = 0; iEntry<n_entries;iEntry++){
      t->GetEntry(iEntry);



      // Make sure the neutrino intersects with the TPC
      if(!IntersectsDetector(vtx_x,vtx_y,vtx_z,px,py,pz)) continue;


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
    N_loaded_hists++;
  }
  
  std::cout << "Cumulative POT: " << cumulativePOT << std::endl;
  // float alt_POT = 5e8 * N_loaded_hists;
  // std::cout << "Alt POT: " << alt_POT << std::endl;
  
  TFile* outFile = new TFile(outFileName.c_str(),"RECREATE");
  // Normalize histograms to per POT per cm2

  // 256.35*233.
  double activeVolFace = 256.35*233.; // in cm^2
  // double POT_factor = 4997*5e8; // Total POT in the flux files
  double normFactor = 1 / (cumulativePOT * activeVolFace); // Normalize to nu/POT/bin/cm2

  // std::cout<< "Using alt POT for normalization: " << alt_POT << std::endl;
  // double normFactor = 1 / (alt_POT * 256.35*233.) ;
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
