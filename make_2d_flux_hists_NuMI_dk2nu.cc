#include <iostream>
#include <string>
#include <cmath>
#include <limits>
#include <algorithm>
#include <filesystem>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"

#include "dk2nu/tree/dk2nu.h"
#include "dk2nu/tree/dkmeta.h"

using std::string;

// -----------------------------------------------------------------------------
// Configuration
// -----------------------------------------------------------------------------

string basePath =
  "/pnfs/uboone/persistent/uboonebeam/numi_dk2nu_g4_10_4_zero_threshold/FHC/"
  "g4numi_minervame_me000z200i_";

string outFileName =
  "/exp/uboone/app/users/jbateman/workdir/DarkNews/Trident/data/flux/numi/"
  "NuMI_g4_10_4_zero_threshold_FHC_FluxHist_w2D_hists.root";


// string basePath =
//   "/pnfs/uboone/persistent/uboonebeam/numi_dk2nu_g4_10_4_zero_threshold/RHC/"
//   "g4numi_minervame_me000z-200i_";

// string outFileName =
//   "/exp/uboone/app/users/jbateman/workdir/DarkNews/Trident/data/flux/numi/"
//   "NuMI_g4_10_4_zero_threshold_RHC_FluxHist_w2D_hists.root";


// MicroBooNE location index in dk2nu metadata/nuray vector.
static constexpr int kMicroBooNENuRayIdx = 6;

// Set this negative to use all matching files.
static constexpr int kMaxFilesToUse = -1;

// Set the number of bins per histogram dimension. The energy range is 0-10 GeV, and the z range is -1000 to 75000 cm.
static constexpr int kNumBinsE = 400; // 100 -> per 100 MeV, 400 -> per 25 MeV
static constexpr int kNumBinsZ = 400; // 100 -> per 760 cm, 400 -> per 190 cm

// -----------------------------------------------------------------------------
// Unit conversion for nuray[i].wgt.
//
// Per dk2nu/tree/calcLocationWeights.{h,cxx} (LBNEAnalysis::CalcLocationWeights,
// which builds nuray[i] for registered locations like MicroBooNE), the weight
// stored in nuray[i].wgt comes from bsim::calcEnuWgt(). Two independent
// downstream consumers of that same function document its output as:
//
//     nuray[i].wgt = flux / (pi * m^2)   [neutrinos / POT / (pi * m^2)]
//
//   - SBN DK2NuInterface.cxx:
//       "weight for the window area and divide by pi
//        (since wgt returned by calcEnuWgt function is flux/(pi*m^2))"
//   - DUNE ppfx FillIMapHists.cpp:
//       weight = nwtnear*nimpwt; hists->_h_nuflux->Fill(enu, weight/pival);
//
// So to get neutrinos/POT/cm^2 you must:
//   1) divide by pi  (both references above do this explicitly)
//   2) convert m^2 -> cm^2: since 1 m^2 = 1e4 cm^2, a per-m^2 density is
//      LARGER than the equivalent per-cm^2 density by 1e4, so multiply by 1e-4.
//
// This was NOT independently re-derived from the raw calcEnuWgt() arithmetic
// (GitHub would not serve the raw source file to this fetcher) -- it is taken
// from two separate production consumers of that function's output, which is
// strong but not airtight evidence. If you have a known total-flux number to
// check against, it's worth a sanity check.
static constexpr double kM2ToCm2 = 1.0e-4;

// -----------------------------------------------------------------------------
// Optional legacy helper: not used by the nuray method.
// -----------------------------------------------------------------------------

float KEfromMom(float mom[3], float mass)
{
  float energy = std::sqrt(
    mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2] + mass*mass
  );
  return energy - mass;
}


// -----------------------------------------------------------------------------
// Optional legacy geometry cut: kept for comparison/debugging only.
// NOT used for nuray-based flux construction -- nuray[kMicroBooNENuRayIdx]
// is already evaluated at the registered MicroBooNE location, so no separate
// geometric cut is needed.
//
// Note: nuray[i].px/py/pz (and decay.vx/vy/vz, decay.pdpx/pdpy/pdpz) are all
// in BEAM coordinates, per LBNEAnalysis::CalcLocationWeights() -- they are
// NOT already expressed in detector coordinates. If you ever re-enable this
// cut, the momentum direction must be rotated the same way the vertex
// position is (rotation only, no translation, since it's a direction vector).
// -----------------------------------------------------------------------------

bool IntersectsDetector(double vx, double vy, double vz,
                        double px, double py, double pz)
{
  // Detector box in detector coordinates, cm.
  static const double box_min[3] = {   0.0, -116.0,    0.0};
  static const double box_max[3] = { 256.0,  116.0, 1036.0};

  // Beam -> detector offset, cm.
  static const double offset[3] = {-31388, -3316, -60100};

  // Beam -> detector rotation.
  static const double rot[3][3] = {
    { 0.9210385380402568,       0.0227135048039241207,  0.38880857519374290},
    { 0.0000462540012621546684, 0.99829162468141475,   -0.0584279894529063024},
    {-0.38947144863934974,      0.0538324139386641073,  0.91946400794392302}
  };

  double pos[3];
  double dir[3];

  for (int i = 0; i < 3; ++i) {
    pos[i] = rot[i][0]*vx + rot[i][1]*vy + rot[i][2]*vz + offset[i];
    // Momentum is a direction vector: rotate only, do not translate.
    dir[i] = rot[i][0]*px + rot[i][1]*py + rot[i][2]*pz;
  }

  double tmin = -std::numeric_limits<double>::infinity();
  double tmax =  std::numeric_limits<double>::infinity();

  for (int i = 0; i < 3; ++i) {
    if (dir[i] == 0.0) {
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

  tmin = std::max(tmin, 0.0);
  return tmax >= tmin && tmax >= 0.0;
}


// -----------------------------------------------------------------------------
// Helper: collect matching ROOT flux files.
// -----------------------------------------------------------------------------

std::vector<string> CollectInputFiles(const string& base_path)
{
  std::filesystem::path dir =
    base_path.substr(0, base_path.find_last_of("/"));

  string prefix =
    base_path.substr(base_path.find_last_of("/") + 1);

  std::vector<string> files;

  std::cout << "Searching in directory: " << dir << std::endl;
  std::cout << "Looking for files with prefix: " << prefix << std::endl;

  for (const auto& entry : std::filesystem::directory_iterator(dir)) {
    if (!entry.is_regular_file()) continue;

    const std::filesystem::path path = entry.path();
    const string name = path.filename().string();

    if (name.rfind(prefix, 0) == 0 && path.extension() == ".root") {
      files.push_back(path.string());
    }
  }

  std::sort(files.begin(), files.end());

  std::cout << "Found " << files.size()
            << " files matching the pattern." << std::endl;

  if (kMaxFilesToUse > 0 && static_cast<int>(files.size()) > kMaxFilesToUse) {
    files.resize(kMaxFilesToUse);
    std::cout << "Only using first " << files.size()
              << " files for now." << std::endl;
  }

  return files;
}


// -----------------------------------------------------------------------------
// Main macro.
// -----------------------------------------------------------------------------

void make_2d_flux_hists_NuMI_dk2nu()
{
  // Print the dkmeta location table once, from the first available file,
  // so kMicroBooNENuRayIdx can be visually confirmed rather than assumed.
  {
    std::vector<string> checkFiles = CollectInputFiles(basePath);
    if (!checkFiles.empty()) {
      TFile* checkF = TFile::Open(checkFiles.front().c_str(), "READ");
      if (checkF && !checkF->IsZombie()) {
        TTree* checkMeta = dynamic_cast<TTree*>(checkF->Get("dkmetaTree"));
        if (checkMeta && checkMeta->GetEntries() > 0) {
          bsim::DkMeta* checkDkMeta = nullptr;
          checkMeta->SetBranchAddress("dkmeta", &checkDkMeta);
          checkMeta->GetEntry(0);
          std::cout << "\ndkmeta location table (first file):" << std::endl;
          for (size_t i = 0; i < checkDkMeta->location.size(); ++i) {
            std::cout << "  [" << i << "] " << checkDkMeta->location[i].name
                      << std::endl;
          }
          std::cout << "Using index " << kMicroBooNENuRayIdx
                    << " as MicroBooNE -- confirm this matches the table above.\n"
                    << std::endl;
        }
      }
      if (checkF) { checkF->Close(); delete checkF; }
    }
  }

  const std::vector<string> inputFiles = CollectInputFiles(basePath);

  if (inputFiles.empty()) {
    std::cerr << "ERROR: no input files found." << std::endl;
    return;
  }

  string particleNames[4]   = {"#nu_{#mu}", "#bar{#nu}_{#mu}", "#nu_{e}", "#bar{#nu}_{e}"};
  string particleSymbols[4] = {"numu",      "numubar",         "nue",     "nuebar"};
  int    pdgID[4]           = {14,          -14,                12,        -12};

  const int nHists = sizeof(particleNames) / sizeof(particleNames[0]);

  TH1D* histsE[nHists];
  TH1D* histsz[nHists];
  TH2D* hists_Ez[nHists];

  for (int iHist = 0; iHist < nHists; ++iHist) {
    string histName = "hE_" + particleSymbols[iHist] + "_cv";
    histsE[iHist] = new TH1D(
      histName.c_str(),
      histName.c_str(),
      kNumBinsE, 0.0, 10.0
    );
    histsE[iHist]->GetXaxis()->SetTitle("Energy at MicroBooNE location (GeV)");
    histsE[iHist]->GetYaxis()->SetTitle("#nu / bin / POT / cm^{2}");

    string histName_z = "hz_" + particleSymbols[iHist] + "_cv";
    histsz[iHist] = new TH1D(
      histName_z.c_str(),
      histName_z.c_str(),
      kNumBinsZ, -1000.0, 75000.0
    );
    histsz[iHist]->GetXaxis()->SetTitle("Decay z position, beam coordinates (cm)");
    histsz[iHist]->GetYaxis()->SetTitle("#nu / bin / POT / cm^{2}");

    string histName_Ez = "hE_vs_z_" + particleSymbols[iHist] + "_cv";
    hists_Ez[iHist] = new TH2D(
      histName_Ez.c_str(),
      histName_Ez.c_str(),
      kNumBinsE, 0.0, 10.0,
      kNumBinsZ, -1000.0, 75000.0
    );
    hists_Ez[iHist]->GetXaxis()->SetTitle("Energy at MicroBooNE location (GeV)");
    hists_Ez[iHist]->GetYaxis()->SetTitle("Decay z position, beam coordinates (cm)");
    hists_Ez[iHist]->GetZaxis()->SetTitle("#nu / bin / POT / cm^{2}");
  }

  TH1D* pot = new TH1D("POT", "POT", 1, 0.0, 1.0);

  double cumulativePOT = 0.0;
  int nLoadedFiles = 0;

  std::cout << "Starting to loop over " << inputFiles.size()
            << " files..." << std::endl;

  for (const string& inFile : inputFiles) {
    std::cout << "\nProcessing file: " << inFile << std::endl;

    TFile* f = TFile::Open(inFile.c_str(), "READ");
    if (!f || f->IsZombie()) {
      std::cout << "  Warning: could not open file, skipping." << std::endl;
      if (f) {
        f->Close();
        delete f;
      }
      continue;
    }

    TTree* fluxTree = dynamic_cast<TTree*>(f->Get("dk2nuTree"));
    TTree* metaTree = dynamic_cast<TTree*>(f->Get("dkmetaTree"));

    if (!fluxTree || !metaTree) {
      std::cout << "  Warning: missing dk2nuTree/dkmetaTree, skipping."
                << std::endl;
      f->Close();
      delete f;
      continue;
    }

    bsim::Dk2Nu*  dk2nu  = nullptr;
    bsim::DkMeta* dkmeta = nullptr;

    fluxTree->SetBranchAddress("dk2nu",  &dk2nu);
    metaTree->SetBranchAddress("dkmeta", &dkmeta);

    // -------------------------------------------------------------------------
    // Sum POT from metadata.
    // -------------------------------------------------------------------------

    double filePOT = 0.0;
    const Long64_t nMeta = metaTree->GetEntries();

    for (Long64_t iMeta = 0; iMeta < nMeta; ++iMeta) {
      metaTree->GetEntry(iMeta);
      filePOT += dkmeta->pots;
    }

    cumulativePOT += filePOT;

    std::cout << "  File POT: " << filePOT << std::endl;

    // -------------------------------------------------------------------------
    // Event loop.
    //
    // With nuray, do NOT require the generated decay momentum to intersect the
    // TPC. nuray[kMicroBooNENuRayIdx] is the ray already evaluated at the
    // registered MicroBooNE location -- both its momentum direction and its
    // weight are computed specifically for that point (beam-coordinate frame,
    // per LBNEAnalysis::CalcLocationWeights()).
    // -------------------------------------------------------------------------

    const Long64_t nEntries = fluxTree->GetEntries();

    Long64_t nFilledThisFile = 0;
    Long64_t nBadNuRay = 0;

    for (Long64_t iEntry = 0; iEntry < nEntries; ++iEntry) {
      fluxTree->GetEntry(iEntry);

      const int pdg = dk2nu->decay.ntype;
      const double vtx_z = dk2nu->decay.vz;

      if (kMicroBooNENuRayIdx < 0 ||
          kMicroBooNENuRayIdx >= static_cast<int>(dk2nu->nuray.size())) {
        ++nBadNuRay;
        continue;
      }

      const bsim::NuRay& ray = dk2nu->nuray[kMicroBooNENuRayIdx];

      const double E_nu = ray.E;
      const double nwtnear = ray.wgt;
      const double nimpwt = dk2nu->decay.nimpwt;

      if (!std::isfinite(E_nu) ||
          !std::isfinite(nwtnear) ||
          !std::isfinite(nimpwt)) {
        ++nBadNuRay;
        continue;
      }

      if (E_nu <= 0.0) {
        ++nBadNuRay;
        continue;
      }

      // Combined dk2nu weight for this detector location.
      //
      // nwtnear: location/geometric ray weight for nuray[kMicroBooNENuRayIdx],
      //          documented (via calcEnuWgt's downstream consumers) as
      //          neutrinos / POT / (pi * m^2).
      // nimpwt : decay-level importance-sampling correction.
      //
      // Divide by pi and convert m^2 -> cm^2 to land on neutrinos/POT/cm^2.
      // Final POT normalization is applied after all files are processed.
      const double weight =
        (nwtnear * nimpwt / TMath::Pi()) * kM2ToCm2;

      if (!std::isfinite(weight) || weight == 0.0) {
        continue;
      }

      for (int iFind = 0; iFind < nHists; ++iFind) {
        if (pdg == pdgID[iFind]) {
          histsE[iFind]->Fill(E_nu, weight);
          histsz[iFind]->Fill(vtx_z, weight);
          hists_Ez[iFind]->Fill(E_nu, vtx_z, weight);
          ++nFilledThisFile;
          break;
        }
      }
    }

    std::cout << "  Filled " << nFilledThisFile
              << " weighted neutrino entries from this file." << std::endl;

    if (nBadNuRay > 0) {
      std::cout << "  Skipped " << nBadNuRay
                << " entries with missing/non-finite/invalid nuray information."
                << std::endl;
    }

    f->Close();
    delete f;

    ++nLoadedFiles;
  }

  std::cout << "\nLoaded files: " << nLoadedFiles << std::endl;
  std::cout << "Cumulative POT: " << cumulativePOT << std::endl;

  if (cumulativePOT <= 0.0) {
    std::cerr << "ERROR: cumulative POT is <= 0. Histograms will not be "
              << "meaningfully normalized." << std::endl;
  }

  const double normFactor =
    (cumulativePOT > 0.0) ? 1.0 / cumulativePOT : 0.0;

  std::cout << "Applying final normalization factor: "
            << normFactor << std::endl;

  TFile* outFile = new TFile(outFileName.c_str(), "RECREATE");

  for (int iOut = 0; iOut < nHists; ++iOut) {
    histsE[iOut]->Scale(normFactor);
    histsz[iOut]->Scale(normFactor);
    hists_Ez[iOut]->Scale(normFactor);

    histsE[iOut]->Write();
    histsz[iOut]->Write();
    hists_Ez[iOut]->Write();
  }

  pot->SetBinContent(1, cumulativePOT);
  pot->Write();

  outFile->Close();
  delete outFile;

  std::cout << "Wrote output file: " << outFileName << std::endl;
}