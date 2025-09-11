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
#include "ROOT/RDataFrame.hxx"
#include "asym_head.h"

#include<iostream>
#include<fstream>
#include<sys/stat.h>

TH1D *k_pi_sub(TH1D *k_hist, TH1D *pi_hist, TH1D *k_dum, TH1D *pi_dum, TH1D *k_rnd, TH1D *pi_rnd, std::string shms_pos) {
    ConfigManager& config = ConfigManager::getInstance();

    Double_t prompt = 2;
    Double_t rand = 12;
    Double_t rnd_scale = prompt/rand;

    Double_t dum_scale = 0;
    if (shms_pos == "center") {dum_scale = 2803415/(4.8579*298159);}
    if (shms_pos == "right") {dum_scale = 1876555/(4.8579*341334);}
    if (shms_pos == "left") {dum_scale = 2611276/(4.8579*305498);}

    TH1D *k_sub = (TH1D*)k_hist->Clone();
    TH1D *pi_sub = (TH1D*)pi_hist->Clone();

    k_dum->Scale(dum_scale);
    pi_dum->Scale(dum_scale);

    k_rnd->Scale(rnd_scale);
    rnd_scale = rnd_scale * 1.5;
    pi_rnd->Scale(rnd_scale);

    k_sub->Add(k_dum, -1);
    k_sub->Add(k_rnd, -1);

    pi_sub->Add(pi_dum, -1);
    pi_sub->Add(pi_rnd, -1);

    pi_sub->Scale(config.getPI_sub_scale());
    k_sub->Add(pi_sub, -1);

    // backscale dum and rand
    Double_t og_dum_sc = 1/dum_scale;
    Double_t og_PIrnd_sc = 1/(rnd_scale * 1.5);
    Double_t og_Krnd_sc = 1/rnd_scale;

    k_dum->Scale(og_dum_sc);
    pi_dum->Scale(og_dum_sc);
    k_rnd->Scale(og_Krnd_sc);
    pi_rnd->Scale(og_PIrnd_sc);

    return (TH1D*)k_sub;
}

TH1D *pi_clean(TH1D* pi_hist, TH1D* pi_dum, TH1D* pi_rnd, std::string shms_pos) {
    Double_t prompt = 2;
    Double_t rand = 12;
    Double_t rnd_scale = prompt/rand;

    Double_t dum_scale = 0;
    if (shms_pos == "center") {dum_scale = 2803415/(4.8579*298159);}
    if (shms_pos == "right") {dum_scale = 1876555/(4.8579*341334);}
    if (shms_pos == "left") {dum_scale = 2611276/(4.8579*305498);}

    TH1D *pi_sub = (TH1D*)pi_hist->Clone();

    pi_dum->Scale(dum_scale);
    rnd_scale = rnd_scale * 1.5;
    pi_rnd->Scale(rnd_scale);

    pi_sub->Add(pi_dum, -1);
    pi_sub->Add(pi_rnd, -1);

    return (TH1D*)pi_sub;
}

void kin_to_csv(TH1D *Q2, TH1D *W, TH1D *xB, TH1D *epsilon, TH1D *T, Double_t t_min, Double_t t_max,  std::string filename) {
    //const Double_t m_p = 0.938; // Mass of proton

    Double_t Q2_p = Q2->GetMean();
    Double_t Q2_p_err = Q2->GetMeanError();
    Double_t W_p =  W->GetMean();
    Double_t W_p_err = W->GetMeanError();
    // Double_t xB_p = Q2_p / (Q2_p + pow(W_p, 2) - pow(0.938, 2));
    // Double_t xB_p_err = (1/pow((Q2_p + pow(W_p, 2) - pow(m_p, 2)), 4)) * (pow((pow(m_p, 2) + pow(W_p, 2)), 2) * pow(Q2_p_err, 2) + 4 * pow(Q2_p, 2) * pow(W_p, 2) * pow(W_p_err, 2));
    Double_t xB_p = xB->GetMean();
    Double_t xB_p_err = xB->GetMeanError();
    Double_t epsilon_p = epsilon->GetMean();
    Double_t epsilon_p_err = epsilon->GetMeanError();
    Double_t t_p = T->GetMean();
    Double_t t_p_err = T->GetMeanError();
    
    std::fstream csv_file(filename, std::ios::out | std::ios::app);

    bool header_written = false;
    // Check if header has been written
    struct stat file_stat;
    if (stat(filename.c_str(), &file_stat) != 0 || file_stat.st_size == 0) {
        header_written = true;
    } else {
        header_written = false;
    }

    if (header_written) {
        std::cout<< "Header written" << std::endl;
        csv_file << "Q2, Q2 err, W, W err, xB, xB err, epsilon, epsilon err, t, t err, t min, t max\n";
    } 


    csv_file << Q2_p << ", " << Q2_p_err << ", "
            << W_p << ", " << W_p_err << ", "
            << xB_p << ", " << xB_p_err << ", "
            << epsilon_p << ", " << epsilon_p_err << ", "
            << t_p << ", " << t_p_err << ", "
            << t_min << ", " << t_max << "\n";

    csv_file.close();
}

int main(int argc, char* argv[]) {
    ROOT::EnableImplicitMT(8);

    // const std::string particle = "kaonL";    // Select particle
    const std::string root_naming = "Nacer"; // Nacer's/Alicia's root file naming scheme
    // const std::string shms_pos = "left";

    ConfigManager& config = ConfigManager::getInstance();
    config.loadFromJSON("config/cut_config.json");
    
    std::string particle = argv[1]; // Particle type from command line argument
    std::string shms_pos = argv[2]; // SHMS position from command line argument
    Double_t t_min = std::stod(argv[3]); // Minimum t value from command line argument
    Double_t t_max = std::stod(argv[4]); // Maximum t value from command line argument
    
    
    // const Int_t t_bins = 3;
    // Double_t t_edges[t_bins+1] = {0.00675, 0.02025, 0.03075, 0.04275, 0.06075, 0.11025}; // 5 pion bins
    // Double_t t_edges[t_bins+1] = {0.06225, 0.09975, 0.12375, 0.19125}; // 3 kaon bins
    // Double_t t_edges[t_bins+1] = {0.070,0.086,0.094,0.100,0.106,0.114,0.126,0.152}; // nacer bins

    auto start = std::chrono::high_resolution_clock::now();

    std::string csv_file_name = "output/kinav_"+particle+"_"+shms_pos+".csv";


    // Get a tree from root file
    ROOT::RDataFrame posdata_df("Pos_Data", "input/Q0p5W2p40"+shms_pos+"_highe_kaon.root");
    ROOT::RDataFrame negdata_df("Neg_Data", "input/Q0p5W2p40"+shms_pos+"_highe_kaon.root");
    ROOT::RDataFrame posdum_df("Pos_Dummy", "input/Q0p5W2p40"+shms_pos+"_highe_kaon.root");
    ROOT::RDataFrame negdum_df("Neg_Dummy", "input/Q0p5W2p40"+shms_pos+"_highe_kaon.root");

    // Naming convention selection
    std::string H_cal_etottracknorm, H_cer_npeSum, P_hgcer_npeSum, P_aero_npeSum, P_RF_Dist, P_gtr_xp, CTime_eKCoinTime_ROC1, MMK, ph_q;
    std::string Q2, W, epsilon, MandelT;

    if (root_naming == "Nacer" || root_naming == "nacer") {
        H_cal_etottracknorm = "H.cal.etottracknorm";
        H_cer_npeSum = "H.cer.npeSum";
        P_hgcer_npeSum = "P.hgcer.npeSum";
        P_aero_npeSum = "P.aero.npeSum";
        P_RF_Dist = "RFTime.SHMS_RFtimeDist";
        P_gtr_xp = "P.gtr.th";
        CTime_eKCoinTime_ROC1 = "CTime.eKCoinTime_ROC1";
        MMK = "P.kin.secondary.MMK";
        ph_q = "P.kin.secondary.ph_xq";
        MandelT = "P.kin.secondary.MandelT";
        Q2 = "H.kin.primary.Q2";
        W = "H.kin.primary.W";
        epsilon = "H.kin.primary.epsilon";
    }
    if (root_naming == "Alicia" || root_naming == "alicia") {
        H_cal_etottracknorm = "H_cal_etottracknorm";
        H_cer_npeSum = "H_cer_npeSum";
        P_hgcer_npeSum = "P_hgcer_npeSum";
        P_aero_npeSum = "P_aero_npeSum";
        P_RF_Dist = "P_RF_Dist";
        P_gtr_xp = "P_gtr_xp";
        CTime_eKCoinTime_ROC1 = "CTime_eKCoinTime_ROC1";
        MMK = "MMK";
        ph_q = "ph_q";
        MandelT = "MandelT";
        Q2 = "Q2";
        W = "W";
        epsilon = "epsilon";
    }
    
    // RF corection
    std::string RF_corected = P_RF_Dist+"+(3.05*"+P_gtr_xp + ")"; //P_RF_Dist+(3.05*P_gtr_xp)
    // RF cut limits
    std::string RF_high = std::to_string(config.getRF_high());
    std::string RF_low = std::to_string(config.getRF_low());
    // RF cuts
    std::string RF_Kaon_PID = "!("+RF_corected+ ">"+RF_low+" && "+RF_corected+"<"+RF_high+")";
    std::string RF_Pion_PID = RF_corected+">"+RF_low+" && "+RF_corected+"<"+RF_high;


    //HMS Cut limits
    std::string H_cal_low = std::to_string(config.getH_cal_low());
    std::string H_cer_low = std::to_string(config.getH_cer_low());

    //HMS cuts
    std::string HMS = H_cal_etottracknorm+">"+H_cal_low+" && "+H_cer_npeSum+">"+H_cer_low;

    // PID with hgcer and aero
    std::string Kaon_PID = P_hgcer_npeSum+"<1.4 && "+P_aero_npeSum+">2.3";
    std::string Pion_PID = P_hgcer_npeSum+">=0 && "+P_aero_npeSum+">2.3";
    
    //Coin time cut limits
    std::string CT_Kaon_low =  std::to_string(config .getCT_Kaon_low());
    std::string CT_Kaon_high = std::to_string(config.getCT_Kaon_high());
    std::string CT_Pion_low =  std::to_string(config .getCT_Pion_low());
    std::string CT_Pion_high = std::to_string(config.getCT_Pion_high());
    std::string CT_Rand_low_right =  std::to_string(config .getCT_Rand_low_right());
    std::string CT_Rand_low_left =   std::to_string(config  .getCT_Rand_low_left());
    std::string CT_Rand_high_right = std::to_string(config.getCT_Rand_high_right());
    std::string CT_Rand_high_left =  std::to_string(config .getCT_Rand_high_left());

    //Coin Time cuts
    std::string CTime_Kaon = CTime_eKCoinTime_ROC1+">"+CT_Kaon_low+" && "+CTime_eKCoinTime_ROC1+"<"+CT_Kaon_high;
    std::string CTime_Pion = CTime_eKCoinTime_ROC1+"<"+CT_Pion_high+" && "+CTime_eKCoinTime_ROC1+">"+CT_Pion_low;
    std::string CTime_Rand = "!("+CTime_eKCoinTime_ROC1+">"+CT_Rand_low_left+" && "+CTime_eKCoinTime_ROC1+"<"+CT_Rand_low_right+") && "
                            +CTime_eKCoinTime_ROC1+">"+CT_Rand_high_left+" && "+CTime_eKCoinTime_ROC1+"<"+CT_Rand_high_right;

    std::string gen_Kaon_cut = CTime_Kaon + " && " + HMS + " && " + RF_Kaon_PID;
    std::string gen_Pion_cut = CTime_Pion + " && " + HMS + " && " + RF_Pion_PID;

    // MMK Cut limits
    std::string kcut_low = "1.1";
    std::string kcut_high = "1.15";

    // MMK Cuts
    std::string Kaon_MMK = MMK + ">"+kcut_low+" && "+MMK+"<"+kcut_high;


    // for (Int_t i = 1; i<=t_bins; ++i) {
    std::string t_cut = MandelT + ">" + std::to_string(-t_max) + " && " + MandelT + "<" + std::to_string(-t_min);

    auto check1 = std::chrono::high_resolution_clock::now();
    std::cout << "Opening Histograms..." << std::endl;
    
    // Histogram models
    ROOT::RDF::TH1DModel model_Q2("Q2 Average", Q2.c_str(), 200, 0, 1);
    ROOT::RDF::TH1DModel model_W("W Average", W.c_str(), 200, 2, 2.5);
    ROOT::RDF::TH1DModel model_xB("xB Average", "xB",200, 0, 1);
    ROOT::RDF::TH1DModel model_epsilon("epsilon Average", epsilon.c_str(), 200, 0, 1);
    ROOT::RDF::TH1DModel model_T("t Average", MandelT.c_str(), 200, -0.3, 0.1);

    //General combined cut
    std::string Kaon_cut = gen_Kaon_cut + " && " + t_cut;
    std::string Pion_cut = gen_Pion_cut + " && " + t_cut;
    std::string Kaon_r_cut = HMS  + " && " + RF_Kaon_PID + " && " + CTime_Rand + " && " + t_cut;
    std::string Pion_r_cut = HMS  + " && " + RF_Pion_PID + " && " + CTime_Rand + " && " + t_cut;

    // "Lazy" generated histogram. Through RDF these are not type specified and need to be converted to root hist
    // Takes the longes, faster drive=faster execution
    // Positive
    auto h_k_p_Q2 = posdata_df.Filter(Kaon_cut).Histo1D(model_Q2, Q2);
    auto h_k_p_W = posdata_df.Filter(Kaon_cut).Histo1D(model_W, W);
    auto h_k_p_xB = posdata_df.Define("xB", Q2 + "/(" + Q2 + " + (" + W + "*"+ W +") - (0.938*0.938))").Filter(Kaon_cut).Histo1D(model_xB, "xB");
    auto h_k_p_epsilon = posdata_df.Filter(Kaon_cut).Histo1D(model_epsilon, epsilon);
    auto h_k_p_T = posdata_df.Filter(Kaon_cut).Histo1D(model_T, MandelT);
    //Random
    auto h_kr_p_Q2 = posdata_df.Filter(Kaon_r_cut).Histo1D(model_Q2, Q2);
    auto h_kr_p_W = posdata_df.Filter(Kaon_r_cut).Histo1D(model_W, W);
    auto h_kr_p_xB = posdata_df.Define("xB", Q2 + "/(" + Q2 + " + (" + W + "*"+ W +") - (0.938*0.938))").Filter(Kaon_r_cut).Histo1D(model_xB, "xB");
    auto h_kr_p_epsilon = posdata_df.Filter(Kaon_r_cut).Histo1D(model_epsilon, epsilon);
    auto h_kr_p_T = posdata_df.Filter(Kaon_r_cut).Histo1D(model_T, MandelT);
    //Dummy
    auto h_kd_p_Q2 = posdum_df.Filter(Kaon_cut).Histo1D(model_Q2, Q2);
    auto h_kd_p_W = posdum_df.Filter(Kaon_cut).Histo1D(model_W, W);
    auto h_kd_p_xB = posdum_df.Define("xB", Q2 + "/(" + Q2 + " + (" + W + "*"+ W +") - (0.938*0.938))").Filter(Kaon_cut).Histo1D(model_xB, "xB");
    auto h_kd_p_epsilon = posdum_df.Filter(Kaon_cut).Histo1D(model_epsilon, epsilon);
    auto h_kd_p_T = posdum_df.Filter(Kaon_cut).Histo1D(model_T, MandelT);


    //Negative
    auto h_k_n_Q2 = negdata_df.Filter(Kaon_cut).Histo1D(model_Q2, Q2);
    auto h_k_n_W = negdata_df.Filter(Kaon_cut).Histo1D(model_W, W);
    auto h_k_n_xB = negdata_df.Define("xB", Q2 + "/(" + Q2 + " + (" + W + "*"+ W +") - (0.938*0.938))").Filter(Kaon_cut).Histo1D(model_xB, "xB");
    auto h_k_n_epsilon = negdata_df.Filter(Kaon_cut).Histo1D(model_epsilon, epsilon);
    auto h_k_n_T = negdata_df.Filter(Kaon_cut).Histo1D(model_T, MandelT);
    //Random
    auto h_kr_n_Q2 = negdata_df.Filter(Kaon_r_cut).Histo1D(model_Q2, Q2);
    auto h_kr_n_W = negdata_df.Filter(Kaon_r_cut).Histo1D(model_W, W);
    auto h_kr_n_xB = negdata_df.Define("xB", Q2 + "/(" + Q2 + " + (" + W + "*"+ W +") - (0.938*0.938))").Filter(Kaon_r_cut).Histo1D(model_xB, "xB");
    auto h_kr_n_epsilon = negdata_df.Filter(Kaon_r_cut).Histo1D(model_epsilon, epsilon);
    auto h_kr_n_T = negdata_df.Filter(Kaon_r_cut).Histo1D(model_T, MandelT);
    //Dummy
    auto h_kd_n_Q2 = negdum_df.Filter(Kaon_cut).Histo1D(model_Q2, Q2);
    auto h_kd_n_W = negdum_df.Filter(Kaon_cut).Histo1D(model_W, W);
    auto h_kd_n_xB = negdum_df.Define("xB", Q2 + "/(" + Q2 + " + (" + W + "*"+ W +") - (0.938*0.938))").Filter(Kaon_cut).Histo1D(model_xB, "xB");
    auto h_kd_n_epsilon = negdum_df.Filter(Kaon_cut).Histo1D(model_epsilon, epsilon);
    auto h_kd_n_T = negdum_df.Filter(Kaon_cut).Histo1D(model_T, MandelT);

    // Pions
    // Positive
    auto h_p_p_Q2 = posdata_df.Filter(Pion_cut).Histo1D(model_Q2, Q2);
    auto h_p_p_W = posdata_df.Filter(Pion_cut).Histo1D(model_W, W);
    auto h_p_p_xB = posdata_df.Define("xB", Q2 + "/(" + Q2 + " + (" + W + "*"+ W +") - (0.938*0.938))").Filter(Pion_cut).Histo1D(model_xB, "xB");
    auto h_p_p_epsilon = posdata_df.Filter(Pion_cut).Histo1D(model_epsilon, epsilon);
    auto h_p_p_T = posdata_df.Filter(Pion_cut).Histo1D(model_T, MandelT);
    //Random
    auto h_pr_p_Q2 = posdata_df.Filter(Pion_r_cut).Histo1D(model_Q2, Q2);
    auto h_pr_p_W = posdata_df.Filter(Pion_r_cut).Histo1D(model_W, W);
    auto h_pr_p_xB = posdata_df.Define("xB", Q2 + "/(" + Q2 + " + (" + W + "*"+ W +") - (0.938*0.938))").Filter(Pion_r_cut).Histo1D(model_xB, "xB");
    auto h_pr_p_epsilon = posdata_df.Filter(Pion_r_cut).Histo1D(model_epsilon, epsilon);
    auto h_pr_p_T = posdata_df.Filter(Pion_r_cut).Histo1D(model_T, MandelT);
    //Dummy
    auto h_pd_p_Q2 = posdum_df.Filter(Pion_cut).Histo1D(model_Q2, Q2);
    auto h_pd_p_W = posdum_df.Filter(Pion_cut).Histo1D(model_W, W);
    auto h_pd_p_xB = posdum_df.Define("xB", Q2 + "/(" + Q2 + " + (" + W + "*"+ W +") - (0.938*0.938))").Filter(Pion_cut).Histo1D(model_xB, "xB");
    auto h_pd_p_epsilon = posdum_df.Filter(Pion_cut).Histo1D(model_epsilon, epsilon);
    auto h_pd_p_T = posdum_df.Filter(Pion_cut).Histo1D(model_T, MandelT);

    //Negative
    auto h_p_n_Q2 = negdata_df.Filter(Pion_cut).Histo1D(model_Q2, Q2);
    auto h_p_n_W = negdata_df.Filter(Pion_cut).Histo1D(model_W, W);
    auto h_p_n_xB = negdata_df.Define("xB", Q2 + "/(" + Q2 + " + (" + W + "*"+ W +") - (0.938*0.938))").Filter(Pion_cut).Histo1D(model_xB, "xB");
    auto h_p_n_epsilon = negdata_df.Filter(Pion_cut).Histo1D(model_epsilon, epsilon);
    auto h_p_n_T = negdata_df.Filter(Pion_cut).Histo1D(model_T, MandelT);
    //Random
    auto h_pr_n_Q2 = negdata_df.Filter(Pion_r_cut).Histo1D(model_Q2, Q2);
    auto h_pr_n_W = negdata_df.Filter(Pion_r_cut).Histo1D(model_W, W);
    auto h_pr_n_xB = negdata_df.Define("xB", Q2 + "/(" + Q2 + " + (" + W + "*"+ W +") - (0.938*0.938))").Filter(Pion_r_cut).Histo1D(model_xB, "xB");
    auto h_pr_n_epsilon = negdata_df.Filter(Pion_r_cut).Histo1D(model_epsilon, epsilon);
    auto h_pr_n_T = negdata_df.Filter(Pion_r_cut).Histo1D(model_T, MandelT);
    //Dummy
    auto h_pd_n_Q2 = negdum_df.Filter(Pion_cut).Histo1D(model_Q2, Q2);
    auto h_pd_n_W = negdum_df.Filter(Pion_cut).Histo1D(model_W, W);
    auto h_pd_n_xB = negdum_df.Define("xB", Q2 + "/(" + Q2 + " + (" + W + "*"+ W +") - (0.938*0.938))").Filter(Pion_cut).Histo1D(model_xB, "xB");
    auto h_pd_n_epsilon = negdum_df.Filter(Pion_cut).Histo1D(model_epsilon, epsilon);
    auto h_pd_n_T = negdum_df.Filter(Pion_cut).Histo1D(model_T, MandelT);


    std::cout << "Loading histograms..." << std::endl;
    //Defining actual histograms.
    //Kaon
    //Positive
    TH1D *hist_k_Q2_p = new TH1D(*h_k_p_Q2);
    TH1D *hist_k_W_p = new TH1D(*h_k_p_W);
    TH1D *hist_k_xB_p = new TH1D(*h_k_p_xB);
    TH1D *hist_k_epsilon_p = new TH1D(*h_k_p_epsilon);
    TH1D *hist_k_T_p = new TH1D(*h_k_p_T);
    //Random
    TH1D *hist_kr_Q2_p = new TH1D(*h_kr_p_Q2);
    TH1D *hist_kr_W_p = new TH1D(*h_kr_p_W);
    TH1D *hist_kr_xB_p = new TH1D(*h_kr_p_xB);
    TH1D *hist_kr_epsilon_p = new TH1D(*h_kr_p_epsilon);
    TH1D *hist_kr_T_p = new TH1D(*h_kr_p_T);
    //Dummy
    TH1D *hist_kd_Q2_p = new TH1D(*h_kd_p_Q2);
    TH1D *hist_kd_W_p = new TH1D(*h_kd_p_W);
    TH1D *hist_kd_xB_p = new TH1D(*h_kd_p_xB);
    TH1D *hist_kd_epsilon_p = new TH1D(*h_kd_p_epsilon);
    TH1D *hist_kd_T_p = new TH1D(*h_kd_p_T);

    //Negative
    TH1D *hist_k_Q2_n = new TH1D(*h_k_n_Q2);
    TH1D *hist_k_W_n = new TH1D(*h_k_n_W);
    TH1D *hist_k_xB_n = new TH1D(*h_k_n_xB);
    TH1D *hist_k_epsilon_n = new TH1D(*h_k_n_epsilon);
    TH1D *hist_k_T_n = new TH1D(*h_k_n_T);
    //Random
    TH1D *hist_kr_Q2_n = new TH1D(*h_kr_n_Q2);
    TH1D *hist_kr_W_n = new TH1D(*h_kr_n_W);
    TH1D *hist_kr_xB_n = new TH1D(*h_kr_n_xB);
    TH1D *hist_kr_epsilon_n = new TH1D(*h_kr_n_epsilon);
    TH1D *hist_kr_T_n = new TH1D(*h_kr_n_T);
    //Dummy 
    TH1D *hist_kd_Q2_n = new TH1D(*h_kd_n_Q2);
    TH1D *hist_kd_W_n = new TH1D(*h_kd_n_W);
    TH1D *hist_kd_xB_n = new TH1D(*h_kd_n_xB);
    TH1D *hist_kd_epsilon_n = new TH1D(*h_kd_n_epsilon);
    TH1D *hist_kd_T_n = new TH1D(*h_kd_n_T);

    //Pion
    //Positive
    TH1D *hist_p_Q2_p = new TH1D(*h_p_p_Q2);
    TH1D *hist_p_W_p = new TH1D(*h_p_p_W);
    TH1D *hist_p_xB_p = new TH1D(*h_p_p_xB);
    TH1D *hist_p_epsilon_p = new TH1D(*h_p_p_epsilon);
    TH1D *hist_p_T_p = new TH1D(*h_p_p_T);
    //Random
    TH1D *hist_pr_Q2_p = new TH1D(*h_pr_p_Q2);
    TH1D *hist_pr_W_p = new TH1D(*h_pr_p_W);
    TH1D *hist_pr_xB_p = new TH1D(*h_pr_p_xB);
    TH1D *hist_pr_epsilon_p = new TH1D(*h_pr_p_epsilon);
    TH1D *hist_pr_T_p = new TH1D(*h_pr_p_T);
    //Dummy
    TH1D *hist_pd_Q2_p = new TH1D(*h_pd_p_Q2);
    TH1D *hist_pd_W_p = new TH1D(*h_pd_p_W);
    TH1D *hist_pd_xB_p = new TH1D(*h_pd_p_xB);
    TH1D *hist_pd_epsilon_p = new TH1D(*h_pd_p_epsilon);
    TH1D *hist_pd_T_p = new TH1D(*h_pd_p_T);

    //Negative
    TH1D *hist_p_Q2_n = new TH1D(*h_p_n_Q2);
    TH1D *hist_p_W_n = new TH1D(*h_p_n_W);
    TH1D *hist_p_xB_n = new TH1D(*h_p_n_xB);
    TH1D *hist_p_epsilon_n = new TH1D(*h_p_n_epsilon);
    TH1D *hist_p_T_n = new TH1D(*h_p_n_T);
    //Random
    TH1D *hist_pr_Q2_n = new TH1D(*h_pr_n_Q2);
    TH1D *hist_pr_W_n = new TH1D(*h_pr_n_W);
    TH1D *hist_pr_xB_n = new TH1D(*h_pr_n_xB);
    TH1D *hist_pr_epsilon_n = new TH1D(*h_pr_n_epsilon);
    TH1D *hist_pr_T_n = new TH1D(*h_pr_n_T);
    //Dummy
    TH1D *hist_pd_Q2_n = new TH1D(*h_pd_n_Q2);
    TH1D *hist_pd_W_n = new TH1D(*h_pd_n_W);
    TH1D *hist_pd_xB_n = new TH1D(*h_pd_n_xB);
    TH1D *hist_pd_epsilon_n = new TH1D(*h_pd_n_epsilon);
    TH1D *hist_pd_T_n = new TH1D(*h_pd_n_T);

    // Combine positive and negative
    //Kaon
    hist_k_Q2_p->Add(hist_k_Q2_n);
    hist_k_W_p->Add(hist_k_W_n);
    hist_k_xB_p->Add(hist_k_xB_n);
    hist_k_epsilon_p->Add(hist_k_epsilon_n);
    hist_k_T_p->Add(hist_k_T_n);
    // Random
    hist_kr_Q2_p->Add(hist_kr_Q2_n);
    hist_kr_W_p->Add(hist_kr_W_n);
    hist_kr_xB_p->Add(hist_kr_xB_n);
    hist_kr_epsilon_p->Add(hist_kr_epsilon_n);
    hist_kr_T_p->Add(hist_kr_T_n);
    // Dummy
    hist_kd_Q2_p->Add(hist_kd_Q2_n);
    hist_kd_W_p->Add(hist_kd_W_n);
    hist_kd_xB_p->Add(hist_kd_xB_n);
    hist_kd_epsilon_p->Add(hist_kd_epsilon_n);
    hist_kd_T_p->Add(hist_kd_T_n);

    //Pion
    hist_p_Q2_p->Add(hist_p_Q2_n);
    hist_p_W_p->Add(hist_p_W_n);
    hist_p_xB_p->Add(hist_p_xB_n);
    hist_p_epsilon_p->Add(hist_p_epsilon_n);
    hist_p_T_p->Add(hist_p_T_n);
    // Random
    hist_pr_Q2_p->Add(hist_pr_Q2_n);
    hist_pr_W_p->Add(hist_pr_W_n);
    hist_pr_xB_p->Add(hist_pr_xB_n);
    hist_pr_epsilon_p->Add(hist_pr_epsilon_n);
    hist_pr_T_p->Add(hist_pr_T_n);
    // Dummy
    hist_pd_Q2_p->Add(hist_pd_Q2_n);
    hist_pd_W_p->Add(hist_pd_W_n);
    hist_pd_xB_p->Add(hist_pd_xB_n);
    hist_pd_epsilon_p->Add(hist_pd_epsilon_n);
    hist_pd_T_p->Add(hist_pd_T_n);

    // Positive histograms are now the combination of positive and negative helicities
    std::cout << "Cleaning kaon sample..." << std::endl;

    TH1D *Q2_h = nullptr, *W_h = nullptr, *xB_h = nullptr, *epsilon_h = nullptr, *T_h = nullptr;

    if (particle == "kaonL") {
        Q2_h = k_pi_sub(hist_k_Q2_p, hist_p_Q2_p, hist_kd_Q2_p, hist_pd_Q2_p, hist_kr_Q2_p, hist_pr_Q2_p, shms_pos);
        W_h = k_pi_sub(hist_k_W_p, hist_p_W_p, hist_kd_W_p, hist_pd_W_p, hist_kr_W_p, hist_pr_W_p, shms_pos);
        xB_h = k_pi_sub(hist_k_xB_p, hist_p_xB_p, hist_kd_xB_p, hist_pd_xB_p, hist_kr_xB_p, hist_pr_xB_p, shms_pos);
        epsilon_h = k_pi_sub(hist_k_epsilon_p, hist_p_epsilon_p, hist_kd_epsilon_p, hist_pd_epsilon_p, hist_kr_epsilon_p, hist_pr_epsilon_p, shms_pos);
        T_h = k_pi_sub(hist_k_T_p, hist_p_T_p, hist_kd_T_p, hist_pd_T_p, hist_kr_T_p, hist_pr_T_p, shms_pos);
    }
    if (particle == "pion") {
        Q2_h = pi_clean(hist_p_Q2_p, hist_pd_Q2_p, hist_pr_Q2_p, shms_pos);
        W_h = pi_clean(hist_p_W_p, hist_pd_W_p, hist_pr_W_p, shms_pos);
        xB_h = pi_clean(hist_p_xB_p, hist_pd_xB_p, hist_pr_xB_p, shms_pos);
        epsilon_h = pi_clean(hist_p_epsilon_p, hist_pd_epsilon_p, hist_pr_epsilon_p, shms_pos);
        T_h = pi_clean(hist_p_T_p, hist_pd_T_p, hist_pr_T_p, shms_pos);
    }
    std::cout << "Calculating variables... " << std::endl;

    kin_to_csv(Q2_h, W_h, xB_h, epsilon_h, T_h, t_min, t_max, csv_file_name);

    std::cout << "t bin" << " complete\n" << "---------------------------------------------------------" << std::endl;
    
    auto check2 = std::chrono::high_resolution_clock::now();
    std::cout<< "Bin Time: " << std::chrono::duration<double>(check2 - check1).count() << "s" << std::endl;


    delete Q2_h;
    delete W_h;
    delete epsilon_h;
    delete T_h;

    delete hist_k_Q2_p;
    delete hist_k_W_p;
    delete hist_k_epsilon_p;
    delete hist_k_T_p;

    delete hist_kr_Q2_p;
    delete hist_kr_W_p;
    delete hist_kr_epsilon_p;
    delete hist_kr_T_p;

    delete hist_kd_Q2_p;
    delete hist_kd_W_p;
    delete hist_kd_epsilon_p;
    delete hist_kd_T_p;

    delete hist_k_Q2_n;
    delete hist_k_W_n;
    delete hist_k_epsilon_n;
    delete hist_k_T_n;

    delete hist_kr_Q2_n;
    delete hist_kr_W_n;
    delete hist_kr_epsilon_n;
    delete hist_kr_T_n;

    delete hist_kd_Q2_n;
    delete hist_kd_W_n;
    delete hist_kd_epsilon_n;
    delete hist_kd_T_n;

    delete hist_p_Q2_p;
    delete hist_p_W_p;
    delete hist_p_epsilon_p;
    delete hist_p_T_p;

    delete hist_pr_Q2_p;
    delete hist_pr_W_p;
    delete hist_pr_epsilon_p;
    delete hist_pr_T_p;

    delete hist_pd_Q2_p;
    delete hist_pd_W_p;
    delete hist_pd_epsilon_p;
    delete hist_pd_T_p;

    delete hist_p_Q2_n;
    delete hist_p_W_n;
    delete hist_p_epsilon_n;
    delete hist_p_T_n;

    delete hist_pr_Q2_n;
    delete hist_pr_W_n;
    delete hist_pr_epsilon_n;
    delete hist_pr_T_n;

    delete hist_pd_Q2_n;
    delete hist_pd_W_n;
    delete hist_pd_epsilon_n;
    delete hist_pd_T_n;
    // }
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Elapsed time: " << std::chrono::duration<double>(end - start).count() << " s\n";
}