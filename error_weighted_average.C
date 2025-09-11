#include<TMath.h>
#include<TFile.h>
#include<TTree.h>
#include<TCanvas.h>
#include<TH1D.h>
#include<TCut.h>
#include<TPad.h>
#include<TStyle.h>
#include<TROOT.h>
#include<TF1.h>
#include<Math/Minimizer.h>

#include "ROOT/RDataFrame.hxx"

#include<iostream>
#include<fstream>

void error_weighted_average(TString particle, TString t_bin) {
    TCanvas* wa_canv = new TCanvas("wa_canv", "weighted average canvas", 500, 500);

    // TString particle = "pion";
    // TString t_bin = "-t0.060750-0.110250";
    
    TFile* center_f = new TFile("output/Asym_"+particle+"_center_"+t_bin+".root");
    TFile* right_f = new TFile("output/Asym_"+particle+"_right_"+t_bin+".root");
    TFile* left_f = new TFile("output/Asym_"+particle+"_left_"+t_bin+".root");

    TCanvas* c_canv = (TCanvas*)center_f->Get("c4");
    TCanvas* r_canv = (TCanvas*)right_f->Get("c4");
    TCanvas* l_canv = (TCanvas*)left_f->Get("c4");

    TH1F* center_h = (TH1F*)c_canv->FindObject("asym_hist");
    TH1F* right_h = (TH1F*)r_canv->FindObject("asym_hist");
    TH1F* left_h = (TH1F*)l_canv->FindObject("asym_hist");

    Int_t nbins = center_h->GetNbinsX();
    TH1F* asym = new TH1F("asym",t_bin,nbins,-3.14,3.14);
    
    for (Int_t i = 1; i <= nbins; ++i) {
        Double_t nc = center_h->GetBinContent(i),
                 nr = right_h->GetBinContent(i), 
                 nl = left_h->GetBinContent(i);
        
        Double_t sc = center_h->GetBinError(i),
                 sr = right_h->GetBinError(i),
                 sl = left_h->GetBinError(i);

        Double_t wc, wr, wl;
        
        if (sc==0){wc=0;} else {wc = 1/(sc*sc);}
        if (sr==0){wr=0;} else {wr = 1/(sr*sr);}
        if (sl==0){wl=0;} else {wl = 1/(sl*sl);}

        Double_t avg = (nl*wl +nc*wc +nr*wr)/(wl+wc+wr);
        Double_t avg_err = sqrt(1/(wl+wc+wr));

        asym->SetBinContent(i, avg);
        asym->SetBinError(i, avg_err);
    }

    wa_canv->cd();
    //plot nicely
    asym->SetStats(0);
	asym->SetMarkerStyle(21);
	asym->SetMarkerSize(1.5);
	asym->SetLineWidth(3);
	//asym->SetTitleSize(0.1);
	asym->GetXaxis()->SetTitle("#phi");
	asym->GetYaxis()->SetTitle("BSA");
	asym->SetTitleSize(0.05,"y");
	asym->SetLabelSize(0.04,"y");
	asym->SetLabelSize(0.04,"x");
	asym->SetTitleSize(0.05,"x");
	asym->SetTitleSize(0.07,"t");
	asym->GetYaxis()->SetRangeUser(-0.5,0.5);
	asym->GetYaxis()->SetTitleOffset(0.8);
	asym->GetXaxis()->SetTitleOffset(0.7);
	asym->Draw("E0");
	gStyle->SetTitleH(0.09);

    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2", "Migrad"); // Attempted to change minimizer to get a better fit

	//fit with sin curve and 3 parameter full functional form curve
	// TF1 *threepar = new TF1("threepar","[0]*sin(x)/(1+[1]*cos(x)+[2]*cos(2*x))");
	TF1 *sinfit = new TF1("sinfit","[0]*sin(x)");
	sinfit->SetLineStyle(kDashed);
	sinfit->SetLineWidth(3);
	// threepar->SetLineWidth(3);
    // threepar->SetParLimits(1, -0.9999, 0.9999);
    // threepar->SetParLimits(2, -0.95, 0.95);

	// asym->Fit("threepar", "M");

    // std::cout << "\nChi2: " << threepar->GetChisquare() << std::endl;
    // std::cout << "NDF: " << threepar->GetNDF() << std::endl;
    // std::cout << "Reduced Chi2: " << (threepar->GetChisquare()/threepar->GetNDF()) << "\n" << std::endl;

	asym->Fit("sinfit");

    std::cout << "\n Chi2: " << sinfit->GetChisquare() << std::endl;
    std::cout << "NDF: " << sinfit->GetNDF() << std::endl;
    std::cout << "Reduced Chi2: " << (sinfit->GetChisquare()/sinfit->GetNDF()) << std::endl;

    asym->SetTitle(t_bin + " rChi2:" + std::to_string(sinfit->GetChisquare()/sinfit->GetNDF()).c_str());

    wa_canv->SaveAs("output/WA_"+particle+"_"+t_bin+".root");
}