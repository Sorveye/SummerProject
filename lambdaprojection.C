#include <TFile.h>
#include <TDirectoryFile.h>
#include <THnSparse.h>
#include <TAxis.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <memory>

void lambdaprojection() {
    
    // Defining centrality and pt ranges
    double ptmin[] = {1.0, 1.5, 2.0, 2.5};
    double ptmax[] = {1.5, 2.0, 2.5, 3.0};
    
    double cosmin[] = {-1., -0.8, -0.6, -0.4, -0.2, 0., 0.2, 0.4, 0.6, 0.8};
    double cosmax[] = {-0.8, -0.6, -0.4, -0.2, 0., 0.2, 0.4, 0.6, 0.8, 1.};
    
    int nset = 4;
    int cosnset = 10;
    
    // Opening the file
    TFile* fin = new TFile("Mult_pt_cos_Lambda.root", "READ");
    
    // Extracting the histogram
    THnSparseD* h_lambda = dynamic_cast<THnSparseD*>(fin->Get("histogram"));
    
    TAxis* axis_mult = h_lambda->GetAxis(0);
    TAxis* axis_pt = h_lambda->GetAxis(1);
    TAxis* axis_cossine = h_lambda->GetAxis(2);
    TAxis* axis_planx = h_lambda->GetAxis(3);
    TAxis* axis_plany = h_lambda->GetAxis(4);
    
    double sigma[nset];
    double mean[nset];
    
    double asigma[nset];
    double amean[nset];
    
    TH1D* projection[nset];
    TH1D* aprojection[nset];
    
    TF1* fitFunc = new TF1("fitFunc", "[0] * exp(-0.5*((x-[1])/[2])^2) + [3] + [4]*x + [5]*x*x", 1.0, 2.0);
    
    fitFunc->SetParameters(1000, 1.11, 0.006, 1.0, 1.0, 1.0);
    fitFunc->SetParLimits(1, 1.113, 1.118);
    fitFunc->SetParLimits(2, 0.001, 0.01);
    
    for (int j = 0; j < cosnset; j++) {
        int binmin = h_lambda->GetAxis(2)->FindBin(cosmin[j]);
        int binmax = h_lambda->GetAxis(2)->FindBin(cosmax[j]);
        h_lambda->GetAxis(2)->SetRange(binmin, binmax); //Selecting the cossine range
        
        int nBins[5] = {axis_mult->GetNbins(), axis_pt->GetNbins(), axis_cossine->GetNbins()};
        double xMin[5] = {axis_mult->GetXmin(), axis_pt->GetXmin(), axis_cossine->GetXmin()};
        double xMax[5] = {axis_mult->GetXmax(), axis_pt->GetXmax(), axis_cossine->GetXmax()};
        
        THnSparseD* histogram = new THnSparseD(Form("histogram_cos_%.1f_%.1f", cosmin[j], cosmax[j]), "Lambda Projection", 5, nBins, xMin, xMax);//Saving the normalization parameter, the pt range and the cossine range
        
        double normalization[nset];
       
        for (int i = 0; i < nset; i++) {
            int bin_min = axis_pt->FindBin(ptmin[i]);
            int bin_max = axis_pt->FindBin(ptmax[i]);
            h_lambda->GetAxis(1)->SetRange(bin_min, bin_max); //Selecting the pt range
		
            projection[i] = h_lambda->Projection(0); //Projecting over the multiplicity (should it be the invariant mass?)
            projection[i]->Sumw2();
            projection[i]->Fit(fitFunc);
		
            mean[i] = projection[i]->GetMean();
            sigma[i] = projection[i]->GetStdDev();
		
            normalization[i] = fitFunc->GetParameter(0);
		
            histogram->Fill(normalization[i], (ptmin[i] + ptmax[i])/2, (cosmin[j] + cosmax[j])/2);
        }
        std::string filename = Form("projection_cos_%.1f_%.1f.root", cosmin[j], cosmax[j]);
        std::unique_ptr<TFile> file(TFile::Open(filename.c_str(), "RECREATE"));
        histogram->Write();
        file->Close();
    }
}

