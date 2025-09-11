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
#include<TGraph.h>
#include<TGraphErrors.h>
#include<TSystem.h>

#include "ROOT/RDataFrame.hxx"

#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<cstdio>

void csv_to_tree(std::string particle, std::string shms_pos, TFile *root_file) {
    std::string csv_name = "output/kinav_"+particle+"_"+shms_pos+".csv";
    std::string root_name = "output/kinav_"+particle+"_"+shms_pos+".root";

    
    std::ifstream kinav_csv(csv_name);
    std::string line;
    
    if (!kinav_csv.is_open()) {
        std::cerr << "Error: Could not open file " << csv_name << std::endl;
        return;
    }

    std::getline(kinav_csv, line); // Skip header lines

    Double_t Q2, W, xB, epsilon, T;
    Double_t Q2_err, W_err, xB_err, epsilon_err, T_err;
    Double_t t_min, t_max;

    TTree *kinav_tree = new TTree(shms_pos.c_str(), "Tree for KinAv CSV");

    kinav_tree->Branch("Q2", &Q2, "Q2/D");
    kinav_tree->Branch("Q2_err", &Q2_err, "Q2_err/D");
    kinav_tree->Branch("W", &W, "W/D");
    kinav_tree->Branch("W_err", &W_err, "W_err/D");
    kinav_tree->Branch("xB", &xB, "xB/D");
    kinav_tree->Branch("xB_err", &xB_err, "xB_err/D");
    kinav_tree->Branch("epsilon", &epsilon, "epsilon/D");
    kinav_tree->Branch("epsilon_err", &epsilon_err, "epsilon_err/D");
    kinav_tree->Branch("t", &T, "t/D");
    kinav_tree->Branch("t_err", &T_err, "t_err/D");
    kinav_tree->Branch("t_min", &t_min, "t_min/D");
    kinav_tree->Branch("t_max", &t_max, "t_max/D");

    while (std::getline(kinav_csv, line)) {
        std::stringstream ss(line);
        std::string val;

        std::getline(ss, val, ','); Q2 = std::stod(val);
        std::getline(ss, val, ','); Q2_err = std::stod(val);
        std::getline(ss, val, ','); W = std::stod(val);
        std::getline(ss, val, ','); W_err = std::stod(val);
        std::getline(ss, val, ','); xB = std::stod(val);
        std::getline(ss, val, ','); xB_err = std::stod(val);
        std::getline(ss, val, ','); epsilon = std::stod(val);
        std::getline(ss, val, ','); epsilon_err = std::stod(val);
        std::getline(ss, val, ','); T = std::stod(val);
        std::getline(ss, val, ','); T_err = std::stod(val);
        std::getline(ss, val, ','); t_min = std::stod(val);
        std::getline(ss, val, ','); t_max = std::stod(val);
        
        kinav_tree->Fill();
    }
    kinav_csv.close();

    root_file->cd();
    kinav_tree->Write(shms_pos.c_str(), TObject::kOverwrite);
    std::cout << shms_pos + " tree written." << std::endl;
}

void root_to_plot(const std::string particle, TFile *root_file) {
    std::string root_name_src = "output/kinav_"+particle+".root";
    std::string root_name_out = "output/kinav_"+particle+"_plots.root";
    
    TFile *root_file_out = TFile::Open(root_name_out.c_str(), "RECREATE");

    ROOT::RDataFrame kinav_tree_df("WATree", root_name_src);

    // Iterate over Q2, W, xB, epsilon to get plots vs T
    std::vector<std::string> variables = {"Q2", "W", "xB", "epsilon", "t"};

    for (const auto &var : variables) {    
        std::string x_var = "t", y_var = var;


        auto x_vals = kinav_tree_df.Take<Double_t>(x_var); if (x_var == "t") for (Double_t &v: *x_vals) {v = -v;} // Negate t values
        auto x_errs = kinav_tree_df.Take<Double_t>(x_var+"_err"); if (x_var == "t") for (Double_t &v: *x_errs) {v = -v;} // Negate t errors

        auto y_vals = kinav_tree_df.Take<Double_t>(y_var);
        auto y_errs = kinav_tree_df.Take<Double_t>(y_var+"_err");

        TGraphErrors *kinav_graph = new TGraphErrors(x_vals->size(), x_vals->data(), y_vals->data(), x_errs->data(), y_errs->data());
        kinav_graph->SetTitle((y_var + " vs " + x_var).c_str());
        kinav_graph->GetXaxis()->SetTitle(x_var.c_str());
        kinav_graph->GetYaxis()->SetTitle(y_var.c_str());
        kinav_graph->SetMarkerStyle(20);
        kinav_graph->SetMarkerColor(kBlue);
        kinav_graph->SetLineColor(kBlue);
        kinav_graph->SetLineWidth(2);

        TCanvas *kinav_canvas = new TCanvas(("kinav_" + particle + "_" + var).c_str(), ("kinav_" + particle + "_" + var).c_str(), 800, 600);
        kinav_canvas->SetGrid();
        kinav_graph->Draw("APL");

        root_file_out->cd();
        kinav_canvas->Write();
    }

    root_file_out->Close();
}

void kinav_reader(std::string particle, Int_t t_bins) {
    ROOT::EnableImplicitMT(8);
    // const std::string particle = "kaonL";
    // const Int_t t_bins = 3; // Number of t bins, change as needed
    std::string root_file_name = "output/kinav_"+particle+".root";

    TFile *root_file = TFile::Open((root_file_name).c_str(), "RECREATE");

    csv_to_tree(particle, "center", root_file);
    csv_to_tree(particle, "right", root_file);
    csv_to_tree(particle, "left", root_file);

    root_file->Close();
    
    // Can possibly be written as a function but not necessary
    root_file = TFile::Open((root_file_name).c_str(), "UPDATE");

    ROOT::RDataFrame center_df("center", root_file_name);
    ROOT::RDataFrame right_df("right", root_file_name);
    ROOT::RDataFrame left_df("left", root_file_name);

    std::vector<std::string> variables = {"Q2", "W", "xB", "epsilon", "t"};

    Double_t Q2, W, xB, epsilon, T;
    Double_t Q2_err, W_err, xB_err, epsilon_err, T_err;

    TTree *kinav_tree = new TTree("WATree", "Tree for KinAv CSV");

    kinav_tree->Branch("Q2", &Q2, "Q2/D");
    kinav_tree->Branch("Q2_err", &Q2_err, "Q2_err/D");
    kinav_tree->Branch("W", &W, "W/D");
    kinav_tree->Branch("W_err", &W_err, "W_err/D");
    kinav_tree->Branch("xB", &xB, "xB/D");
    kinav_tree->Branch("xB_err", &xB_err, "xB_err/D");
    kinav_tree->Branch("epsilon", &epsilon, "epsilon/D");
    kinav_tree->Branch("epsilon_err", &epsilon_err, "epsilon_err/D");
    kinav_tree->Branch("t", &T, "t/D");
    kinav_tree->Branch("t_err", &T_err, "t_err/D");

    std::vector<std::vector<Double_t>> t_var(t_bins, std::vector<Double_t>(10, -1));

    for (const auto &var : variables) {
        auto nc = center_df.Take<Double_t>(var);
        auto nr = right_df.Take<Double_t>(var);
        auto nl = left_df.Take<Double_t>(var);
        
        auto sc = center_df.Take<Double_t>(var + "_err");
        auto sr = right_df.Take<Double_t>(var + "_err");
        auto sl = left_df.Take<Double_t>(var + "_err");

        for (Int_t i = 0; i < nc->size(); ++i) {
            if (std::isnan((*nc)[i]) || std::isnan((*nr)[i]) || std::isnan((*nl)[i]) ||
                std::isnan((*sc)[i]) || std::isnan((*sr)[i]) || std::isnan((*sl)[i])) {
                std::cerr << "Warning: NaN value encountered in " << var << " at index " << i << ". Skipping." << std::endl;
                continue;
            }

            Double_t wc, wr, wl;
            
            if ((*sc)[i]==0){wc=0;} else {wc = 1/((*sc)[i]*(*sc)[i]);}
            if ((*sr)[i]==0){wr=0;} else {wr = 1/((*sr)[i]*(*sr)[i]);}
            if ((*sl)[i]==0){wl=0;} else {wl = 1/((*sl)[i]*(*sl)[i]);}

            Double_t avg = ((*nl)[i]*wl + (*nc)[i]*wc + (*nr)[i]*wr)/(wl+wc+wr);
            Double_t avg_err = sqrt(1/(wl+wc+wr));

            auto it = std::find(variables.begin(), variables.end(), var);
            auto index = std::distance(variables.begin(), it);

            t_var[i][index*2] = avg;
            t_var[i][index*2+1] = avg_err;
        }
    }

    for (Int_t i = 0; i < t_var.size(); ++i) {
        for (Int_t j = 0; j < t_var[0].size(); ++j) {
            if (t_var[i][j] == -1) {
                std::cerr << "Warning: No valid data for index " << i << ", variable " << j << ". Skipping." << std::endl;
                continue;
            }

            switch (j) {
                case 0: Q2 = t_var[i][j]; Q2_err = t_var[i][j+1]; break;
                case 2: W = t_var[i][j]; W_err = t_var[i][j+1]; break;
                case 4: xB = t_var[i][j]; xB_err = t_var[i][j+1]; break;
                case 6: epsilon = t_var[i][j]; epsilon_err = t_var[i][j+1]; break;
                case 8: T = t_var[i][j]; T_err = t_var[i][j+1]; break;
            }
        }
        kinav_tree->Fill();
    }

    root_file->cd();
    kinav_tree->Write("WATree", TObject::kOverwrite);
    // End of possible weighted average function.
    root_file->Close();

    root_to_plot(particle, root_file);

    std::remove(root_file_name.c_str()); // Remove the original root file
}