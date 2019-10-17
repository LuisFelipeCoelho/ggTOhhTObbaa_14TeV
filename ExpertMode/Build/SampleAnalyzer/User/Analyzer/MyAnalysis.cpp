#include "SampleAnalyzer/User/Analyzer/MyAnalysis.h"
#include "SampleAnalyzer/Commons/Vector/MABoost.h"
//#include "tools/delphes/classes/DelphesModule.h"
#include <iostream>
using namespace MA5;
using namespace std;

// -----------------------------------------------------------------------------
// Initialize
// function called one time at the beginning of the analysis
// -----------------------------------------------------------------------------
bool MyAnalysis::Initialize(const MA5::Configuration& cfg, const std::map<std::string,std::string>& parameters)
{
  cout << "BEGIN Initialization" << endl;
  // initialize variables, histos
  PHYSICS->mcConfig().Reset();
  // Initializing the histogram
  //myHisto1 = new TH1F("my Histo1", "#theta (diagrama da caixa)", 80, 0, 3.25);
  //myHisto2 = new TH1F("my Histo2", "#theta (depois do boost, diagrama da caixa)", 80, 0, 3.25);
  //myHisto1 = new TH1F("my Histo1", "#theta (diagrama do triangulo)", 80, 0, 3.25);
  //myHisto2 = new TH1F("my Histo2", "#theta (depois do boost, diagrama do triangulo)", 80, 0, 3.25);
  myHisto1 = new TH1F("my Histo1", "#theta (caixa e triangulo)", 80, 0, 3.25);
  myHisto2 = new TH1F("my Histo2", "#theta (depois do boost, caixa e triangulo)", 80, 0, 3.25);



  cout << "END   Initialization" << endl;
  return true;
}

// -----------------------------------------------------------------------------
// Finalize
// function called one time at the end of the analysis
// -----------------------------------------------------------------------------
void MyAnalysis::Finalize(const SampleFormat& summary, const std::vector<SampleFormat>& files)
{
  cout << "BEGIN Finalization" << endl;
  // saving histos
  // Color of the canvas background
  gStyle->SetCanvasColor(0);
  // Turning off the border lines of the canvas
  gStyle->SetCanvasBorderMode(0);
  // Configuring statistics printing
  gStyle->SetOptStat(111110);
  // Creating the output root file
  TCanvas* myCanvas1 = new TCanvas("myCanvas1", "");
  // Setting background color
  myHisto1->SetFillColor(kBlue-2);
  // Setting axis title
  myHisto1->GetXaxis()->SetTitle("#theta");
  // Drawing histogram
  myHisto1->Draw();
  // Saving plot
  //myCanvas1->SaveAs(("theta_caixa.pdf"));
  //myCanvas1->SaveAs(("theta_caixa.epf"));
  //myCanvas1->SaveAs(("theta_triangulo.pdf"));
  //myCanvas1->SaveAs(("theta_triangulo.epf"));
  myCanvas1->SaveAs(("theta.pdf"));
  myCanvas1->SaveAs(("theta.epf"));

  TCanvas* myCanvas2 = new TCanvas("myCanvas2", "");
  myHisto2->SetFillColor(kBlue-2);
  myHisto2->GetXaxis()->SetTitle("#theta");
  myHisto2->Draw();
  //myCanvas2->SaveAs(("theta_caixa_restframe.png"));
  //myCanvas2->SaveAs(("theta_caixa_restframe.epf"));
  //myCanvas2->SaveAs(("theta_triangulo_restframe.png"));
  //myCanvas2->SaveAs(("theta_triangulo_restframe.epf"));
  myCanvas2->SaveAs(("theta_restframe.png"));
  myCanvas2->SaveAs(("theta_restframe.epf"));

  cout << "END   Finalization" << endl;
}

// -----------------------------------------------------------------------------
// Execute
// function called each time one event is read
// -----------------------------------------------------------------------------
bool MyAnalysis::Execute(SampleFormat& sample, const EventFormat& event)
{
  // Event weight
  double myEventWeight;
  if(Configuration().IsNoEventWeight()) myEventWeight=1.;
  else if(event.mc()->weight()!=0.) myEventWeight=event.mc()->weight();
  else
  {
    //WARNING << "Found one event with a zero weight. Skipping..." << endmsg;
    return false;
  }
  Manager()->InitializeForNewEvent(myEventWeight);
  // Empty event
  if (event.mc()==0) {return true;}

  // Initialization
  MALorentzVector pBottom, pBottomBar, pPhoton1, pPhoton2, pHiggspp, pHiggsbb;
  MALorentzVector V_cross_bottom, V_cross_photon;
  const MCParticleFormat* higgspp = 0;
  const MCParticleFormat* higgsbb = 0;
  const MCParticleFormat* Bottom1 = 0;
  const MCParticleFormat* Bottom2 = 0;
  const MCParticleFormat* Photon1 = 0;
  const MCParticleFormat* Photon2 = 0;
  double nPhoton = 0;
  double nBottom = 0;

  // Identification of the particles of interest
  for (unsigned i=0; i<event.mc()->particles().size(); i++)
  {
    const MCParticleFormat* part = &(event.mc()->particles()[i]);

    // The bottom is a final-state particle
    if (!PHYSICS->Id->IsFinalState(part)) continue;

    // Getting the mother
    if (abs(part->pdgid())!=5 && abs(part->pdgid()!=22)) continue;
    const MCParticleFormat* mother = part->mothers()[0];
    if (mother==0) continue;
    if (abs(mother->pdgid())!=25) continue;

    // Bottom selection based on the PDG-id code
    if (part->pdgid()==5 && nBottom==0) //  && mother->pdgid()==25
    {
      pBottom = part->momentum();
      Bottom1 = part;
      nBottom += 1;
      higgsbb = mother;
      pHiggsbb = mother->momentum();
      continue;
    }
    // anti-Bottom selection based on the PDG-id code
    if (part->pdgid()==-5 && nBottom==1)
    {
      pBottomBar = part->momentum();
      Bottom2 = part;
      nBottom += 1;
      // higgsbb = mother;
      continue;
    }
    // Photon1 selection based on the PDG-id code
    if (part->pdgid()==22 && nPhoton==0)
    {
      pPhoton1 = part->momentum();
      Photon1 = part;
      nPhoton = 1;
      higgspp = mother;
      pHiggspp = mother->momentum();
      continue;
    }
    // Photon2 selection based on the PDG-id code
    if (part->pdgid()==22 && nPhoton==1)
    {
      pPhoton2 = part->momentum();
      Photon2 = part;
      nPhoton = 2;
      // higgspp = mother;
      continue;
    }
    // In case ther is more than 2 photons
    if (nPhoton>=2)
    {
      WARNING << "More than 2 photons:" << endmsg;
      nPhoton += 1;
      MALorentzVector pPhotons3 = part->momentum();
      cout << nPhoton << "\t" << pPhotons3.Pt()  << "\t" << pPhoton1.Pt()  << "\t" << pPhoton2.Pt() << endl;
    }
  }

  // -----------------Cross Product-------------------
  double Vbx = (pBottom.Py()*pBottomBar.Pz())-(pBottom.Pz()*pBottomBar.Py());
  double Vby = (pBottom.Pz()*pBottomBar.Px())-(pBottom.Px()*pBottomBar.Pz());
  double Vbz = (pBottom.Px()*pBottomBar.Py())-(pBottom.Py()*pBottomBar.Px());

  double Vpx = (pPhoton1.Py()*pPhoton2.Pz())-(pPhoton1.Pz()*pPhoton2.Py());
  double Vpy = (pPhoton1.Pz()*pPhoton2.Px())-(pPhoton1.Px()*pPhoton2.Pz());
  double Vpz = (pPhoton1.Px()*pPhoton2.Py())-(pPhoton1.Py()*pPhoton2.Px());

  V_cross_bottom.SetXYZT(Vbx,Vby,Vbz,0);
  V_cross_photon.SetXYZT(Vpx,Vpy,Vpz,0);

  // Angle between V_cross_bottom and V_cross_photon
  myHisto1->Fill(V_cross_bottom.Angle(V_cross_photon));

  /*
  cout << "----------------------"<< endl;
  cout << "n:      " << nBottom << endl;
  cout << "p:      " << pBottom.Px() << "\t" << pBottomBar.Px() << "\t" << pPhoton1.Px() << "\t" << pPhoton2.Px() << endl;
  cout << "Vp:     " << Vpx << "\t" << Vpy << "\t" << Vpz << "\t" << endl;
  cout << "Vb:     " << Vbx << "\t" << Vby << "\t" << Vbz << "\t" << endl;
  */

  // ------------Virtual particle and Lorentz Boost-----------------------
  MALorentzVector new_pPhoton1, new_pPhoton2, new_pBottom1, new_pBottom2;
  MCParticleFormat new_bottom1 = *Bottom1;
  MCParticleFormat new_bottom2 = *Bottom2;
  MCParticleFormat new_photon1 = *Photon1;
  MCParticleFormat new_photon2 = *Photon2;

  MCParticleFormat Virtual, Cross_photon, Cross_bottom;
  MALorentzVector pVirtual = pHiggsbb + pHiggspp;
  Virtual.setMomentum(pVirtual);

  // Lorentz Boost
  new_photon1.ToRestFrame(Virtual);
  new_photon2.ToRestFrame(Virtual);
  new_bottom1.ToRestFrame(Virtual);
  new_bottom2.ToRestFrame(Virtual);

  new_pPhoton1 = new_photon1.momentum();
  new_pPhoton2 = new_photon2.momentum();
  new_pBottom1 = new_bottom1.momentum();
  new_pBottom2 = new_bottom2.momentum();

  MALorentzVector new_V_cross_photon = Cross_photon.momentum();
  MALorentzVector new_V_cross_bottom = Cross_bottom.momentum();

  double new_Vbx = (new_pBottom1.Py()*new_pBottom2.Pz())-(new_pBottom1.Pz()*new_pBottom2.Py());
  double new_Vby = (new_pBottom1.Pz()*new_pBottom2.Px())-(new_pBottom1.Px()*new_pBottom2.Pz());
  double new_Vbz = (new_pBottom1.Px()*new_pBottom2.Py())-(new_pBottom1.Py()*new_pBottom2.Px());

  double new_Vpx = (new_pPhoton1.Py()*new_pPhoton2.Pz())-(new_pPhoton1.Pz()*new_pPhoton2.Py());
  double new_Vpy = (new_pPhoton1.Pz()*new_pPhoton2.Px())-(new_pPhoton1.Px()*new_pPhoton2.Pz());
  double new_Vpz = (new_pPhoton1.Px()*new_pPhoton2.Py())-(new_pPhoton1.Py()*new_pPhoton2.Px());

  new_V_cross_bottom.SetXYZT(new_Vbx, new_Vby, new_Vbz,0);
  new_V_cross_photon.SetXYZT(new_Vpx, new_Vpy, new_Vpz,0);

  // Angle between Cross_bottom and Cross_photon in the rest frame of the Virtual particle
  myHisto2->Fill(new_V_cross_bottom.Angle(new_V_cross_photon));

/*
  // ------------Virtual particle and Lorentz Boost-----------------------
  pVirtual = pHiggsbb + pHiggspp;
  MABoost boost1 = MABoost(pVirtual.Px(), pVirtual.Py(), pVirtual.Pz());
  boost1.boost(V_cross_bottom);
  boost1.boost(V_cross_photon);

  // Angle between Cross_bottom and Cross_photon in the rest frame of the Virtual particle
  myHisto2->Fill(V_cross_bottom.Angle(V_cross_photon));
*/

  return true;
}
