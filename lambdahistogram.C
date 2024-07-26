#include "Pythia8/Pythia.h"
//#include "Pythia8Plugins/ColourReconnectionHooks.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH3D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THnSparse.h"
#include "TList.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TNtuple.h"
#include "TProfile.h"
#include "TRandom3.h"
#include "TString.h"
#include "TTree.h"
#include "TVector3.h"
//#include "fastjet/PseudoJet.hh"
#include <algorithm>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <vector>
//#include "fastjet/ClusterSequence.hh"
//#include "fastjet/ClusterSequenceArea.hh"
#include <cstdio> // needed for io
#include <ctime>
#include <iostream> // needed for io
#include <time.h>   /* time */
#include <valarray>
//#include <yaml.h>
// include <stdio.h>
// include <glib.h>
//#include <yaml-cpp/yaml.h>

using namespace Pythia8;

void lambdahistogram() {
    Pythia pythia;
    // Generate event
    pythia.readString("Beams:eCM = 13000.");
    pythia.readString("SoftQCD:inelastic = on");
    pythia.init();

    int nevents = 10000;
    int lambda_PDG = 3122;

    int nbins = 100;
    double multV0A_min = 0.0, multV0A_max = 100.0;
    double pt_min = 1.0, pt_max = 3.0;
    double cos_theta_min = -1.0, cos_theta_max = 1.0;
    double plan_event_min = -100.0, plan_event_max = 100.0;

    // Note the change from THnD to THnSparseD
    int ndim = 5;
    int nbins_array[5] = { nbins, nbins, nbins, nbins, nbins };
    double xmin[5] = { multV0A_min, pt_min, cos_theta_min, plan_event_min, plan_event_min };
    double xmax[5] = { multV0A_max, pt_max, cos_theta_max, plan_event_max, plan_event_max };
    THnSparseD* histogram = new THnSparseD("histogram", "Lambda histogram;multiplicity;transverse momentum;cosine theta;plan_event_x;plan_event_y", ndim, nbins_array, xmin, xmax);
    TH1D* projection = new TH1D("projection", "projection lambda", 100, 0, 10);

    for (int event = 0; event < nevents; ++event) {
        if (!pythia.next()) continue;

        double multV0A = 0;
        double multV0C = 0;
        double multiplicity = 0;
        double plan_event_x = 0;
        double plan_event_y = 0;
        double plan_norm = 0;

        for (int i = 0; i < pythia.event.size(); ++i) {
            Float_t eta = pythia.event[i].eta();
            Float_t pt = pythia.event[i].pT();
            Int_t pid = pythia.event[i].id();

            if (!pythia.event[i].isFinal() || !pythia.event[i].isCharged()) continue;
            if (eta > -3.7 && eta < -1.7) {
                multV0A += 1;
                plan_event_x += TMath::Cos(2.0 * pythia.event[i].phi());
                plan_event_y += TMath::Sin(2.0 * pythia.event[i].phi());
                plan_norm += 1.0;
            }
            if (eta > 2.8 && eta < 5.1) {
                multV0C += 1;
                plan_event_x += pt * TMath::Cos(2.0 * pythia.event[i].phi());
                plan_event_y += pt * TMath::Sin(2.0 * pythia.event[i].phi());
                plan_norm += 1.0;
            }
            if (TMath::Abs(eta) < 0.5) multiplicity += 1;
        }

        for (int i = 0; i < pythia.event.size(); ++i) {
            Float_t eta = pythia.event[i].eta();
            Float_t pt = pythia.event[i].pT();
            Int_t pid = pythia.event[i].id();

            if (TMath::Abs(pid) == lambda_PDG) {
                Vec4 lambdaMomentum = pythia.event[i].p(); // getting the lambda momentum

                // std::cout << lambdaMomentum << std::endl;
                int daughter1 = pythia.event[i].daughter1();
                int daughter2 = pythia.event[i].daughter2();

        
 		if (!((TMath::Abs(pythia.event[daughter1].id()) == 2212 && (TMath::Abs(pythia.event[daughter2].id()) == 211 || TMath::Abs(pythia.event[daughter2].id()) == 111)) || 
                      ((TMath::Abs(pythia.event[daughter1].id()) == 211 || TMath::Abs(pythia.event[daughter1].id()) == 111) && TMath::Abs(pythia.event[daughter2].id()) == 2212))) {
                    continue;
                }
                
                int proton = (TMath::Abs(pythia.event[daughter1].id()) == 2212) ? daughter1 : daughter2;
                if (proton < 0) continue;
                Vec4 protonMomentum = pythia.event[proton].p();

                double betaX = lambdaMomentum.px() / lambdaMomentum.e(); // Creating the Lorentz boost beta/E
                double betaY = lambdaMomentum.py() / lambdaMomentum.e();
                double betaZ = lambdaMomentum.pz() / lambdaMomentum.e();
                Vec4 boostedProtonMomentum = protonMomentum;
                boostedProtonMomentum.bst(-betaX, -betaY, -betaZ); // boosting the Proton Momentum

                double cosine_theta = TMath::Cos(boostedProtonMomentum.theta());

                double multV0M = multV0A + multV0C; // total event multiplicity

                double values[5] = { multV0M, pt, cosine_theta, plan_event_x, plan_event_y };
                histogram->Fill(values);
                projection->Fill(pt);
            }
        }
    }

    // Save the histogram to a ROOT file
    std::unique_ptr<TFile> file(TFile::Open("Mult_pt_cos_Lambda.root", "RECREATE"));
    histogram->Write();
    file->Close();

    std::unique_ptr<TFile> arquivo(TFile::Open("1dpt.root", "RECREATE"));
    projection->Write();
    arquivo->Close();

    delete histogram;
}

int main() {
    lambdahistogram();
    return 0;
}
