#include <cstdlib>
#include <iostream>
#include <vector>
#include <chrono>

#include <TApplication.h>
#include <TCanvas.h>
#include <TString.h>
#include <TH1F.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>

#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Random.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/ViewFEMesh.hh"
#include "Garfield/ViewDrift.hh"
#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/TrackHeed.hh"
#include "Garfield/ViewSignal.hh"
#include "Garfield/RandomEngineRoot.hh"


using namespace Garfield;

// Helper function to save a canvas to a folder
void SaveCanvasToFolder(TCanvas* canvas, const TString& folder, const TString& name, const TString& ext = "png") {
    TString filename = folder + "/" + name + "." + ext;
    canvas->SaveAs(filename);
    std::cout << "Saved " << filename << "\n";
}

// Planes height [cm]
const double ytop = 0.51;
const double ymesh = 0.01;
const double ybottom = 0.0;

// Dimensions of the sensor [cm]
const double xmin = -0.15; double xmax = 0.15;
const double zmin = 0.; double zmax = 0.3;  

int main(int argc, char * argv[]) {
  // Check command line arguments
  if(argc<4){
      fprintf(stderr,"Use %s DeltaV, rPenning, event\n",argv[0]);
      exit(1);
  }
  
  // Read command line arguments and create branch variables
  int DeltaV = atoi(argv[1]);       // Voltage difference in V 
  double rPenning = atof(argv[2]);  // Penning transfer efficiency
  int event = atoi(argv[3]);        // Event number
  //std::string particle = argv[4]; // Particle type (muon or photon)
  std::string particle = "muon";    // Fixed to muon for this example
  int nPrimaryElectrons = 0;        // Number of primary electrons
  int nFinalElectrons = 0;          // Number of final electrons after the avalanche

  /* Check if the particle type is valid
  if (particle != "muon" && particle != "photon") {
    std::cerr << "Invalid particle type. Use 'muon' or 'photon'.\n";
    return 1;
  }*/

  // Initialize random number generator
  Garfield::RandomEngineRoot rndEngine;
  rndEngine.SetSeed(event);

  // Set the voltage values
  int Vbottom = 0;        // Voltage on the bottom plane in V 
  int Vmesh = -DeltaV;    // Voltage on the mesh in V
  int Vtop = Vmesh - 250; // Voltage on the top plane in V

  // Create output file name
  TString deltaV = Form("%d", DeltaV);          // Voltage difference
  TString penning = Form("%.2f", rPenning);     // Penning transfer efficiency
  TString ev = Form("%d", event);               // Event number
  TString Particle = "muon";                    // Particle type (fixed to muon for this example)
  TString fileName = "plots/prog_"+Particle+"_"+deltaV+"V_"+penning+"_"+ev+".root";
  TFile file(fileName, "RECREATE");
  std::cout << "Particle: " << "muon" << " DeltaV: " << deltaV << " V, Penning: " << penning << ", Event: " << event << "\n"
            << "Output file: " << fileName << "\n";
  
  // TTree creation
  TTree tree("tree", "tree");
  tree.Branch("event", &event);
  tree.Branch("particle", &particle);
  tree.Branch("DeltaV", &DeltaV);
  tree.Branch("rPenning", &rPenning);
  tree.Branch("nPrimaryElectrons", &nPrimaryElectrons, "nPrimaryElectrons/I");
  tree.Branch("nFinalElectrons", &nFinalElectrons, "nFinalElectrons/I");

  // Setup------------------------------------------------------------------
  // Initialize ROOT application
  TApplication app("app", &argc, argv);
  gROOT->SetBatch(kTRUE);

  // Gas
  MediumMagboltz gas("ar", 93., "co2", 5., "c4h10", 2.);
  //MediumMagboltz gas("ar", 80., "co2", 20.);
  gas.SetTemperature(293.15);
  gas.SetPressure(760.);
  gas.Initialise(true);
  gas.LoadIonMobility("IonMobility_Ar+_Ar.txt");
  //gas.LoadIonMobility("IonMobility_C8Hn+_iC4H10.txt");

  // Set the Penning transfer efficiency.
  //constexpr double rPenning = 0.51;
  constexpr double lambdaPenning = 0.;
  gas.EnablePenningTransfer(rPenning, lambdaPenning, "ar");

  // Component 
  // Drift region (top to mesh)
  ComponentAnalyticField driftRegion;
  driftRegion.SetMedium(&gas);
  driftRegion.AddPlaneY(ytop, Vtop, "top");      // Top plane
  driftRegion.AddPlaneY(ymesh, Vmesh, "mesh");   // Mesh plane

  // Amplification region (mesh to bottom)
  ComponentAnalyticField ampRegion;
  ampRegion.SetMedium(&gas);
  ampRegion.AddPlaneY(ymesh, Vmesh, "mesh");        // Mesh plane
  ampRegion.AddPlaneY(ybottom, Vbottom, "bottom");  // Bottom plane

  std::cout << "Vtop: " << Vtop << " V, Vmesh: " << Vmesh << " V, Vbottom: " << Vbottom << " V\n" 
            << "ytop: " << ytop << " cm, ymesh: " << ymesh << " cm, ybottom: " << ybottom << " cm\n";

  // Strips in the amplification region
  ampRegion.AddStripOnPlaneY('z', ybottom, -0.15, -0.06, "strip1", 0.015);
  ampRegion.AddStripOnPlaneY('z', ybottom, -0.045, 0.045, "strip2", 0.015);
  ampRegion.AddStripOnPlaneY('z', ybottom, 0.06, 0.15, "strip3", 0.015);
  const double x1 = -(0.15 + 0.06)/2; // Center of strip1
  const double x2 = 0.0;              // Center of strip2
  const double x3 = -x1;              // Center of strip3

  // Sensor
  Sensor sensor;
  sensor.SetArea(xmin, ybottom, zmin, xmax, ytop, zmax);
  sensor.AddComponent(&driftRegion);
  sensor.AddComponent(&ampRegion);

  // Add electrodes
  sensor.AddElectrode(&ampRegion, "strip1");
  sensor.AddElectrode(&ampRegion, "strip2");
  sensor.AddElectrode(&ampRegion, "strip3");
  
  // Track muon/photon
  TrackHeed track;
  track.SetParticle("muon") ;
  track.SetEnergy(17.e9); // 17 GeV for muon
  /*
  track.SetParticle(particle);
  if(particle == "muon") {
    track.SetEnergy(17.e9); // 17 GeV
  }else if(particle == "photon") {
    const double p = 167. * RndmUniform();
    const double egamma = p < 100. ? 5898.8 : p < 150. ? 5887.6 : 6490.4; // Energy of the photon
    track.SetEnergy(egamma);
  }*/
  track.SetSensor(&sensor);

  // View for drift lines
  ViewDrift driftView;

  // Avalanche for electrons
  AvalancheMicroscopic avalanche;
  avalanche.SetSensor(&sensor);
  avalanche.EnableMagneticField(true);
  avalanche.EnableSignalCalculation(true);
  //avalanche.EnablePlotting(&driftView);

  // Drift for ions
  AvalancheMC drift;
  drift.SetSensor(&sensor);
  drift.SetDistanceSteps(2e-5); // original 0.00001 cm
  drift.EnableSignalCalculation(true);
  //drift.EnablePlotting(&driftView);

  // Set time window and number of bins
  sensor.ClearSignal();
  const double tStart = 0.; const double tStep = 0.05; const double tEnd = 400; // ns original step 0.01 ns
  const unsigned int nBins = (tEnd-tStart)/tStep;
  sensor.SetTimeWindow(tStart, tStep, nBins);

  // Charge view
  ViewSignal *chargeView1 = new ViewSignal(&sensor);
  TCanvas *cCharge1 = new TCanvas("cCharge1", "", 600, 600);
  chargeView1->SetCanvas(cCharge1);

  ViewSignal *chargeView2 = new ViewSignal(&sensor);
  TCanvas *cCharge2 = new TCanvas("cCharge2", "", 600, 600);
  chargeView2->SetCanvas(cCharge2);

  ViewSignal *chargeView3 = new ViewSignal(&sensor);
  TCanvas *cCharge3 = new TCanvas("cCharge3", "", 600, 600);
  chargeView3->SetCanvas(cCharge3);
  
  // Other Canvas
  TCanvas* c = nullptr;
  TCanvas* czoom = nullptr;
  TCanvas* cSignal1 = nullptr;
  TCanvas* cSignal2 = nullptr;
  TCanvas* cSignal3 = nullptr;
  TCanvas* cSignals = nullptr;  

  // Drift------------------------------------------------------------------
  // Set random starting point
  std::vector<double> centers = {x1, x2, x3}; // Vector of centers for the strips
  int rnd = static_cast<int>(RndmUniform() * centers.size());   // Generate a random index between 0 and vec.size()-1
  
  //const double edge = 0.1*(xmax-xmin);
  //const double x0 = xmin+edge+RndmUniform()*(xmax-xmin-2*edge);
  const double x0 = centers[rnd]; // Random x0 from the centers vector
  const double y0 = ytop;
  const double z0 = zmin+0.5*(zmax-zmin);
  const double t0 = 0.0;

  // Start a new track
  //track.EnableDebugging();
  //track.EnablePlotting(&driftView);
  track.DisableDeltaElectronTransport();
  track.NewTrack(x0, y0, z0, t0, 0.0, -1.0, 0.0);
  
  unsigned int nTotal = 0, nBF = 0;         // Total number of ions and back-flowing ions
  double xc, yc, zc, tc, ec, extra; int ne; // Cluster properties
  double xe, ye, ze, te, ee, dx, dy, dz;    // Electron properties

  // Loop initial time
  auto loopStart = std::chrono::high_resolution_clock::now();


  // Loop over the clusters
  for (const auto& cluster : track.GetClusters()) {
    xc = cluster.x; yc = cluster.y; zc = cluster.z;
    ec = cluster.energy; ne = cluster.electrons.size();
    std::cout << "Cluster at (" << xc << ", " << yc << ", " << zc << ") with " << ne << " electrons, energy " << ec << " eV.\n";

    // Loop over the electrons in the cluster
    for (const auto& electron : cluster.electrons) {
      /*
      avalanche = AvalancheMicroscopic();
      avalanche.SetSensor(&sensor);
      avalanche.EnableMagneticField(true);
      avalanche.EnableSignalCalculation(true);
      */

      // Set the electron properties
      xe = electron.x; ye = electron.y; ze = electron.z;
      te = electron.t; ee = electron.e;
      dx = electron.dx; dy = electron.dy; dz = electron.dz;
      std::cout << "    Electron at (" << xe << ", " << ye << ", " << ze << ") with " << ee << " eV.\n";
      avalanche.AvalancheElectron(xe, ye, ze, te, ee, dx, dy, dz);

      const auto& electrons = avalanche.GetElectrons();

      // 1st Avalanche size
      int nel1 = 0, ni1 = 0;
      avalanche.GetAvalancheSize(nel1, ni1);
      std::cout << "    1st Avalanche size: " << nel1 << " electrons, " << ni1
                << " ions.\n";
      nPrimaryElectrons += nel1; // Count primary electrons

      // Loop over the electrons in the avalanche
      for (const auto& electron : electrons) { 
        // Endpoint of the electrons
        const auto& endpoint = electron.path.back();
        std::cout << "        Electron endpoint at mesh: (" << endpoint.x << ", " << endpoint.y << ", " << endpoint.z << ")\n";
        // Check if the electron is in the amplification region
        if (endpoint.y <= 0.01) {
          // Electron after the mesh
          const double xafter = endpoint.x;   const double yafter = endpoint.y; const double zafter = endpoint.z;  
          const double tafter = endpoint.t;   const double eafter = endpoint.energy;
          const double dxafter = endpoint.kx; const double dyafter = endpoint.ky; const double dzafter = endpoint.kz;
        
          avalanche.AvalancheElectron(xafter, yafter, zafter, tafter, eafter, dxafter, dyafter, dzafter);

          // 2nd Avalanche size
          int nel2 = 0, ni2 = 0;
          avalanche.GetAvalancheSize(nel2, ni2);
          std::cout << "        2nd Avalanche size: " << nel2 << " electrons, " << ni2 << " ions.\n";
          nFinalElectrons += nel2; // Count final electrons after the avalanche

          // Drift ions in the amplification region
          for (const auto& avalElectron : avalanche.GetElectrons()) {
            const auto& p0 = avalElectron.path[0];
            drift.DriftIon(p0.x, p0.y, p0.z, p0.t);
            ++nTotal;
            const auto& ionEndpoint = drift.GetIons().front().path.back();
            if (ionEndpoint.y > ymesh) ++nBF;
        }
        } else {
          std::cout << "Electron ended before the mesh.\n";
        }
      }
    }   
  }

  // Print results
  std::cout << "Total number of ions: " << nTotal << std::endl;
  std::cout << "Number of back-flowing ions: " << nBF << std::endl;
  std::cout << "Fraction of back-flowing ions: " << double(nBF) / double(nTotal) << "\n";
  std::cout << "Number of primary electrons: " << nPrimaryElectrons << "\n";
  std::cout << "Number of final electrons after the avalanche: " << nFinalElectrons << "\n";
  std::cout << "Gain: " << double(nFinalElectrons) / double(nPrimaryElectrons) << "\n";  

  // Loop final time
  auto loopEnd = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> loopDuration = loopEnd - loopStart;
  
  // Plot initial time 
  auto plotStart = std::chrono::high_resolution_clock::now();

  // Drift Plotting------------------------------------------------------------------
  constexpr bool PlotDrift = false;
  if (PlotDrift) {
    // Zoomed area
    double xzoom_min = std::min(x0-0.05,xmin); double xzoom_max = std::max(x0+0.05, xmax);
  
    // Large area
    c = new TCanvas("cDrift", "Drift", 800, 600);
    driftView.SetArea(xmin, ybottom, zmin, xmax, ytop, zmax);
    driftView.SetPlaneXY();
    driftView.SetCanvas(c);
    constexpr bool twod = true;
    driftView.Plot(twod);

    // Zoomed area
    czoom = new TCanvas("cDriftZoom", "Drift (zoomed)", 800, 600);
    driftView.SetArea(xzoom_min, ybottom, zmin, xzoom_max, ymesh, zmax);
    driftView.SetPlaneXY();
    driftView.SetCanvas(czoom);
    driftView.Plot(twod);
    
    // Update the canvases
    c->Update();
    czoom->Update();

    // Save the canvases
    SaveCanvasToFolder(c, "plots", "drift");
    SaveCanvasToFolder(czoom, "plots", "drift_zoom");
  }


  // Signal Plotting------------------------------------------------------------------
  constexpr bool PlotSignal = true;
  // Original plot with PlotSignal
  if (PlotSignal){
      // Canvases for signals
      cSignal1 = new TCanvas("cSignal1", "Signal1", 800, 600);
      cSignal2 = new TCanvas("cSignal2", "Signal2", 800, 600);
      cSignal3 = new TCanvas("cSignal3", "Signal3", 800, 600);

      cSignal1->SetTitle("Signal on Strip 1;Time [ns];Signal [fC]");
      cSignal2->SetTitle("Signal on Strip 2;Time [ns];Signal [fC]");
      cSignal3->SetTitle("Signal on Strip 3;Time [ns];Signal [fC]");

      sensor.PlotSignal("strip1", cSignal1);
      sensor.PlotSignal("strip2", cSignal2);
      sensor.PlotSignal("strip3", cSignal3);

      // Update and save the canvases
      cSignal1->Update(); SaveCanvasToFolder(cSignal1, "plots", "signal1");
      cSignal2->Update(); SaveCanvasToFolder(cSignal2, "plots", "signal2");
      cSignal3->Update(); SaveCanvasToFolder(cSignal3, "plots", "signal3");
  }

  // Alternative signal plotting in 1 canvas with GetSignal
  constexpr bool PlotSignal2 = true;
  if (PlotSignal2) {

    // Create vectors for time and signals
    std::vector<double> times(nBins);
    std::vector<double> signal1(nBins);
    std::vector<double> signal2(nBins);
    std::vector<double> signal3(nBins);

    // Fill time vector and get signals
    for (unsigned int i = 0; i < nBins; ++i) {
      times[i] = tStart + i * tStep;
      signal1[i] = sensor.GetSignal("strip1", i);
      signal2[i] = sensor.GetSignal("strip2", i);
      signal3[i] = sensor.GetSignal("strip3", i);
    }

    // Create a canvas to visualize the signals
    cSignals = new TCanvas("cSignals", "Signals", 1200, 900);
    cSignals->Divide(1, 3);

    // Create histograms
    TH1F* hSignal1 = new TH1F("hSignal1", "Signal on Strip 1;Time [ns];Signal [fC]", nBins, tStart, tEnd);
    TH1F* hSignal2 = new TH1F("hSignal2", "Signal on Strip 2;Time [ns];Signal [fC]", nBins, tStart, tEnd);
    TH1F* hSignal3 = new TH1F("hSignal3", "Signal on Strip 3;Time [ns];Signal [fC]", nBins, tStart, tEnd);

    // Fill histograms
    for (unsigned int i = 0; i < nBins; ++i) {
      hSignal1->SetBinContent(i + 1, signal1[i]);
      hSignal2->SetBinContent(i + 1, signal2[i]);
      hSignal3->SetBinContent(i + 1, signal3[i]);
    }

    // Draw histograms
    cSignals->cd(1);
    hSignal1->Draw();
    gPad->SetGridx();
    gPad->SetGridy();

    cSignals->cd(2);
    hSignal2->Draw();
    gPad->SetGridx();
    gPad->SetGridy();

    cSignals->cd(3);
    hSignal3->Draw();
    gPad->SetGridx();
    gPad->SetGridy();

    cSignals->Update();
    // Save the canvas as a PNG file
    SaveCanvasToFolder(cSignals, "plots", "signals");
  }
  
  // Charge plotting with Integrate
  constexpr bool PlotCharge = false;
  if (PlotCharge) {
    sensor.IntegrateSignal("strip1");
    chargeView1->PlotSignal("strip1");
    cCharge1->Update();
    gSystem->ProcessEvents();
    //sensor.ExportSignal("strip1", "charge1");

    sensor.IntegrateSignal("strip2");
    chargeView2->PlotSignal("strip2");
    cCharge2->Update();
    gSystem->ProcessEvents();
    //sensor.ExportSignal("strip2", "charge2");

    sensor.IntegrateSignal("strip3");
    chargeView3->PlotSignal("strip3");
    cCharge3->Update();
    gSystem->ProcessEvents();
    //sensor.ExportSignal("strip3", "charge3");

    SaveCanvasToFolder(cCharge1, "plots", "charge1");
    SaveCanvasToFolder(cCharge2, "plots", "charge2");
    SaveCanvasToFolder(cCharge3, "plots", "charge3");
  }
  // End of plotting---------------------------------------------------------
  // Plot final time
  auto plotEnd = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> plotDuration = plotEnd - plotStart;

  // Print execution times
  std::cout << "Execution times:\n";
  std::cout << "Principal loop execution time: " << loopDuration.count() << " seconds\n";
  std::cout << "Plotting execution time: " << plotDuration.count() << " seconds\n";
  
  // Print total execution time
  auto totalEnd = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> totalDuration = totalEnd - loopStart;
  std::cout << "Total execution time: " << totalDuration.count() << " seconds\n";

  // Tree filling and file writing
  file.cd();
  tree.Fill();
  tree.Write();
  //c->Write(); czoom->Write(); 
  cSignals->Write();
  cSignal1->Write(); cSignal2->Write(); cSignal3->Write();
  //cCharge1->Write(); cCharge2->Write(); cCharge3->Write();

  file.Close();

  // Run the ROOT event loop
  //app.Run(kTRUE);
  app.Terminate();
  return 0;
}


  /*
  const int nClusters = track.GetClusters().size();
  std::cout << "Number of clusters: " << nClusters << std::endl;

  for (int i = 0; i < nClusters; ++i) {
    double xc, yc, zc, tc, ec, extra;
    int ne, ni, np;
    track.GetCluster(xc, yc, zc, tc, ne, ni, np, ec, extra);

    for (int j = 0; j < ne; ++j) {
      double xe, ye, ze, te, ee, dx, dy, dz;
      track.GetElectron(j, xe, ye, ze, te, ee, dx, dy, dz);

      // Avalanche for each electron
      avalanche.AvalancheElectron(xe, ye, ze, te, ee, dx, dy, dz);

      int neAval, niAval;
      avalanche.GetAvalancheSize(neAval, niAval);

      const auto& electrons = avalanche.GetElectrons();
      for (const auto& e : electrons) {
        const auto& p0 = e.path[0];  // last position of the electron
        drift.DriftIon(p0.x, p0.y, p0.z, p0.t);
        ++nTotal;

        const auto& ions = drift.GetIons();
        if (!ions.empty()) {
          const auto& endpoint = ions.front().path.back();
          if (endpoint.z > 0.005) ++nBF;
        }
      }
    }
  }
  */


  /*
  DriftLineRKF line(&sensor);

  for (const auto& cluster : track.GetClusters()) { 
    // Loop over the electrons in the cluster. 
    for (const auto& electron : cluster.electrons) { 
      line.DriftElectron(electron.x, electron.y, electron.z, electron.t); 
    }
  }

  if (ni != 0) {
          for (const auto& electron : avalanche.GetElectrons()) {
          const auto& p0 = electron.path[0];
          drift.DriftIon(p0.x, p0.y, p0.z, p0.t);
          ++nTotal;
          const auto& endpoint = drift.GetIons().front().path.back();
          if (endpoint.z > 0.002) ++nBF;
        }
      }
    
  */

  /*
for (const auto& electron : cluster.electrons) {
  xe = electron.x; ye = electron.y; ze = electron.z;
  te = electron.t; ee = electron.e;
  dx = electron.dx; dy = electron.dy; dz = electron.dz;
  std::cout << "Electron at (" << xe << ", " << ye << ", " << ze << ") with "
  << ee << " eV.\n";
  avalanche.AvalancheElectron(xe, ye, ze, te, ee, dx, dy, dz);

  const auto& electrons = avalanche.GetElectrons();
  const auto& endpoint = electrons.back().path.back();
  std::cout << "Electron endpoint: (" << endpoint.x << ", " << endpoint.y << ", " << endpoint.z << ")\n";

  int nel = 0, ni = 0;
  avalanche.GetAvalancheSize(nel, ni);
  std::cout << "Avalanche size: " << nel << " electrons, " << ni
            << " ions.\n";

  // Drift ions
  for (const auto& electron : avalanche.GetElectrons()) {
    const auto& p0 = electron.path[0];
    drift.DriftIon(p0.x, p0.y, p0.z, p0.t);
    ++nTotal;
    const auto& endpoint = drift.GetIons().front().path.back();
    if (endpoint.z > 0.002) ++nBF;
  }
}
  */