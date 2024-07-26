#include <TFile.h>
#include <THnSparse.h>
#include <TAxis.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TF1.h>
#include <vector>
#include <string>
#include <iostream>

void lambdafitting() {
    // Defining centrality and pt ranges
    double ptmin[] = {1.0, 1.5, 2.0, 2.5};
    double ptmax[] = {1.5, 2.0, 2.5, 3.0};
    
    double cosmin[] = {-1., -0.8, -0.6, -0.4, -0.2, 0., 0.2, 0.4, 0.6, 0.8};
    double cosmax[] = {-0.8, -0.6, -0.4, -0.2, 0., 0.2, 0.4, 0.6, 0.8, 1.};
    
    int nset = 4;
    int cosnset = 10;

    std::vector<double> cos_values(cosnset);
    std::vector<double> normalization_values(cosnset);

    TF1* fitFunc = new TF1("fitFunc", "[0]*(1 + [1]*x)", -1.0, 1.0);
    fitFunc->SetParameters(1.0, 0.1);
    
    for (int j = 0; j < cosnset; j++) {
        std::string filename = Form("projection_cos_%.1f_%.1f.root", cosmin[j], cosmax[j]);
        TFile* fin = TFile::Open(filename.c_str(), "READ");
      
        THnSparseD* h_lambda = dynamic_cast<THnSparseD*>(fin->Get(Form("histogram_cos_%.1f_%.1f", cosmin[j], cosmax[j])));
        if (!h_lambda) {
            std::cerr << "Histogram not found in file: " << filename << std::endl;
            continue;
        }

        TAxis* axis_pt = h_lambda->GetAxis(1);

        double total_normalization = 0.0;
	
        for (int i = 0; i < nset; i++) {
            int bin_min = axis_pt->FindBin(ptmin[i]);
            int bin_max = axis_pt->FindBin(ptmax[i]);
            h_lambda->GetAxis(1)->SetRange(bin_min, bin_max);

            TH1D* proj = h_lambda->Projection(0);  // Projection on the mass axis
            proj->Fit(fitFunc, "Q");  // Fit the projection with the defined function
            total_normalization += fitFunc->GetParameter(0);

            delete proj;
        }
        
        normalization_values[j] = total_normalization;
        cos_values[j] = 0.5 * (cosmin[j] + cosmax[j]);

        fin->Close();
    }

    TGraph* graph = new TGraph(cosnset, cos_values.data(), normalization_values.data());
    graph->SetTitle("Normalization vs Cosine Theta;Cosine Theta;Normalization");

    TCanvas* canvas = new TCanvas("canvas", "Normalization vs Cosine Theta", 800, 600);
    graph->Draw("AP");
    graph->Fit(fitFunc, "R");

    canvas->SaveAs("normalization_vs_cos_theta.png");

    delete graph;
    delete canvas;
    delete fitFunc;
}

int main() {
    lambdafitting();
    return 0;
}

