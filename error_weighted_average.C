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
#include<TGraphAsymmErrors.h>

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

    auto center_h = (TGraphAsymmErrors*)c_canv->FindObject("asym_hist");
    auto right_h = (TGraphAsymmErrors*)r_canv->FindObject("asym_hist");
    auto left_h = (TGraphAsymmErrors*)l_canv->FindObject("asym_hist");

    Int_t nbins = center_h->GetN();
    auto asym = new TGraphAsymmErrors();
    
    for (Int_t i = 1; i <= nbins; ++i) {
        Double_t nc = center_h->GetPointY(i),
                 nr = right_h->GetPointY(i), 
                 nl = left_h->GetPointY(i);

        Double_t sc = center_h->GetErrorY(i),
                 sr = right_h->GetErrorY(i),
                 sl = left_h->GetErrorY(i);

        Double_t wc, wr, wl;
        
        if (sc==0){wc=0;} else {wc = 1/(sc*sc);}
        if (sr==0){wr=0;} else {wr = 1/(sr*sr);}
        if (sl==0){wl=0;} else {wl = 1/(sl*sl);}

        Double_t avg = (nl*wl +nc*wc +nr*wr)/(wl+wc+wr);
        Double_t avg_err = sqrt(1/(wl+wc+wr));

        asym->SetPoint(i, center_h->GetPointX(i), avg);
        asym->SetPointError(i, 0, 0, avg_err, avg_err);
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
	asym->GetYaxis()->SetRangeUser(-0.5,0.5);
	asym->GetYaxis()->SetTitleOffset(0.8);
	asym->GetXaxis()->SetTitleOffset(0.7);
	asym->Draw("AEP");

    gStyle->SetOptFit(1111);
    gStyle->SetTitleFontSize(0.05);
    gStyle->SetTitleX(0.5);
    gStyle->SetTitleAlign(23);
    gStyle->SetTitleW(0.6);
    gStyle->SetTitleY(0.95);
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