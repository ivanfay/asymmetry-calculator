#include "asym_head.h"

std::string H_cal_etottracknorm;
std::string H_cer_npeSum;
std::string P_hgcer_npeSum;
std::string P_aero_npeSum;
std::string P_RF_Dist;
std::string P_gtr_xp;
std::string CTime_eKCoinTime_ROC1;
std::string MMK;
std::string ph_q;
std::string MandelT;

int main(int argc, char* argv[]) {
    ROOT::EnableImplicitMT(6);
    
    // Timing benchmark
    auto start = std::chrono::high_resolution_clock::now();

    std::string infile = argv[1];
    std::string particle = argv[2];
    std::string shms_pos = argv[3];
    Double_t t_high = -1*std::stod(argv[4]);
    Double_t t_low = -1*std::stod(argv[5]);


    FormatInput(infile, particle, shms_pos);

    if (infile == "") {
        std::cout << "Make sure to check your particle and SHMS phi coverage." << std::endl;
        return 0;
    }

    ROOT::RDataFrame posdata_df("Pos_Data", infile);
    ROOT::RDataFrame negdata_df("Neg_Data", infile);
    ROOT::RDataFrame posdum_df("Pos_Dummy", infile);
    ROOT::RDataFrame negdum_df("Neg_Dummy", infile);

    // CheckNaming("Nacer");
    CheckNaming("Alicia");
    // RF corection
    std::string RF_corected = P_RF_Dist+"+(3.05*"+P_gtr_xp + ")"; //P_RF_Dist+(3.05*P_gtr_xp)
    // RF cut limits
    std::string RF_high = "1.15";
    std::string RF_low = "0.15";
    // RF cuts
    std::string RF_Kaon_PID = "!("+RF_corected+ ">"+RF_low+" && "+RF_corected+"<"+RF_high+")";
    std::string RF_Pion_PID = RF_corected+">"+RF_low+" && "+RF_corected+"<"+RF_high;


    //HMS Cut limits
    std::string H_cal_low = "0.5";
    std::string H_cer_low = "2";

    //HMS cuts
    std::string HMS = H_cal_etottracknorm+">"+H_cal_low+" && "+H_cer_npeSum+">"+H_cer_low;

    // PID with hgcer and aero
    std::string Kaon_PID = P_hgcer_npeSum+"<1.4 && "+P_aero_npeSum+">2.3";
    std::string Pion_PID = P_hgcer_npeSum+">=0 && "+P_aero_npeSum+">2.3";
    
    //Coin time cut limits
    // std::string CT_Kaon_low = "-0.5";
    // std::string CT_Kaon_high = "1.5";
    // std::string CT_Pion_low = "-2.5";
    // std::string CT_Pion_high = "-0.5";
    // std::string CT_Rand_low_right = "9.5";
    // std::string CT_Rand_low_left = "-8.5";
    // std::string CT_Rand_high_right = "15.5";
    // std::string CT_Rand_high_left = "-14.5";


    std::string CT_Kaon_low = "-3";
    std::string CT_Kaon_high = "3";
    std::string CT_Pion_low = "-3";
    std::string CT_Pion_high = "3";
    std::string CT_Rand_low_right = "9";
    std::string CT_Rand_low_left = "-9";
    std::string CT_Rand_high_right = "15";
    std::string CT_Rand_high_left = "-15";

    //Coin Time cuts
    std::string CTime_Kaon = CTime_eKCoinTime_ROC1+">"+CT_Kaon_low+" && "+CTime_eKCoinTime_ROC1+"<"+CT_Kaon_high;
    std::string CTime_Pion = CTime_eKCoinTime_ROC1+"<"+CT_Pion_high+" && "+CTime_eKCoinTime_ROC1+">"+CT_Pion_low;
    std::string CTime_Rand = "!("+CTime_eKCoinTime_ROC1+">"+CT_Rand_low_left+" && "+CTime_eKCoinTime_ROC1+"<"+CT_Rand_low_right+") && " 
    +CTime_eKCoinTime_ROC1+">"+CT_Rand_high_left+" && "+CTime_eKCoinTime_ROC1+"<"+CT_Rand_high_right;

    // T Cuts
    std::string t_cut = MandelT+">"+std::to_string(t_low)+" && "+MandelT+"<"+std::to_string(t_high);

    // std::cout << t_cut << std::endl;
    // return 0;

    // std::string gen_Kaon_cut = HMS + " && " + CTime_Kaon + " && " + RF_Kaon_PID;
    // std::string gen_Kaon_rand_cut = HMS + " && " + CTime_Rand + " && " + RF_Kaon_PID;
    // std::string gen_Pion_cut = HMS + " && " + CTime_Pion + " && " + RF_Pion_PID;
    // std::string gen_Pion_rand_cut = HMS + " && " + CTime_Rand + " && " + RF_Pion_PID; 
    std::string gen_Kaon_cut = HMS + " && " + CTime_Kaon + " && " + Kaon_PID;
    std::string gen_Kaon_rand_cut = HMS + " && " + CTime_Rand + " && " + Kaon_PID;
    std::string gen_Pion_cut = HMS + " && " + CTime_Pion + " && " + Pion_PID;
    std::string gen_Pion_rand_cut = HMS + " && " + CTime_Rand + " && " + Pion_PID; 


    std::string csv_name = "output/Yields_"+particle+"_"+shms_pos+"_"+"-t"+std::to_string(t_high).substr(1)+"-"+std::to_string(t_low).substr(1) + ".csv";
    // Assymetry histogram -pi to +pi
    TH1D *asym_hist = new TH1D("asym_hist","Asym",phi_bins,-TMath::Pi(),TMath::Pi());

    // phi bin partitioning
    Double_t single_cut_val = (2 * TMath::Pi()) / phi_bins;
    Double_t cut_sum = -TMath::Pi();

    for (Int_t i = 1; i <= phi_bins; ++i) {
        auto check1 = std::chrono::high_resolution_clock::now();
        // phi uper bin limit
        cut_sum += single_cut_val;

        // Phi bining (dynamic) cuts for data
        std::string phi_cut_upper = ph_q+"<" + std::to_string(cut_sum);
        std::string phi_cut_lower = ph_q+">" + std::to_string(cut_sum - single_cut_val);
        std::string phi_cut = phi_cut_upper + " && " + phi_cut_lower;

        std::cout << "loading histograms..." << std::endl;

        // Names for iteration of histograms
        std::string nameKaon_hist_p = "p_h_kaon" + std::to_string(i); // Makes unique name for interations
        std::string namePion_hist_p = "p_h_pion" + std::to_string(i);
        std::string nameRnd_p = "p_rnd" + std::to_string(i);
        std::string nameRnd_pion_p = "p_rnd_pion" + std::to_string(i);
        std::string nameKaon_hist_n = "n_h_kaon" + std::to_string(i);
        std::string namePion_hist_n = "n_h_pion" + std::to_string(i);
        std::string nameRnd_n = "n_rnd" + std::to_string(i);
        std::string nameRnd_pion_n = "n_rnd_pion" + std::to_string(i);
        std::string nameDum_p = "p_dum" + std::to_string(i);
        std::string nameDum_pion_p = "p_dum_pion" + std::to_string(i);
        std::string nameDum_n = "n_dum" + std::to_string(i);
        std::string nameDum_pion_n = "n_dum" + std::to_string(i);
        

        Int_t bin_num = 200;
        Double_t mmk_low = 0.7;
        Double_t mmk_high = 1.4;
        // Histogram models for intialization of histograms in RDF. Maybe there is a better way, but this works.
        ROOT::RDF::TH1DModel model_kaon_hist_p(nameKaon_hist_p.c_str(), MMK.c_str(),    bin_num, mmk_low, mmk_high);// create histogram model for RDF, Its picky like that
        ROOT::RDF::TH1DModel model_pion_hist_p(namePion_hist_p.c_str(), MMK.c_str(),    bin_num, mmk_low, mmk_high);
        ROOT::RDF::TH1DModel model_Rnd_p(nameRnd_p.c_str(), MMK.c_str(),                bin_num, mmk_low, mmk_high);
        ROOT::RDF::TH1DModel model_Rnd_pion_p(nameRnd_pion_p.c_str(), MMK.c_str(),      bin_num, mmk_low, mmk_high);
        ROOT::RDF::TH1DModel model_kaon_hist_n(nameKaon_hist_n.c_str(), MMK.c_str(),    bin_num, mmk_low, mmk_high);
        ROOT::RDF::TH1DModel model_pion_hist_n(namePion_hist_n.c_str(), MMK.c_str(),    bin_num, mmk_low, mmk_high);
        ROOT::RDF::TH1DModel model_Rnd_n(nameRnd_n.c_str(), MMK.c_str(),                bin_num, mmk_low, mmk_high);
        ROOT::RDF::TH1DModel model_Rnd_pion_n(nameRnd_pion_n.c_str(), MMK.c_str(),      bin_num, mmk_low, mmk_high);
        ROOT::RDF::TH1DModel model_Dum_p(nameDum_p.c_str(), MMK.c_str(),                bin_num, mmk_low, mmk_high);
        ROOT::RDF::TH1DModel model_Dum_pion_p(nameDum_pion_p.c_str(), MMK.c_str(),      bin_num, mmk_low, mmk_high);
        ROOT::RDF::TH1DModel model_Dum_n(nameDum_n.c_str(), MMK.c_str(),                bin_num, mmk_low, mmk_high);
        ROOT::RDF::TH1DModel model_Dum_pion_n(nameDum_pion_n.c_str(), MMK.c_str(),      bin_num, mmk_low, mmk_high);

        auto kaon_h_p =     posdata_df.Filter(t_cut + " && " + phi_cut + " && " + gen_Kaon_cut).Histo1D(model_kaon_hist_p, MMK);
        auto pion_h_p =     posdata_df.Filter(t_cut + " && " + phi_cut + " && " + gen_Pion_cut).Histo1D(model_pion_hist_p, MMK);
        auto rnd_p =        posdata_df.Filter(t_cut + " && " + phi_cut + " && " + gen_Kaon_rand_cut).Histo1D(model_Rnd_p, MMK);
        auto rnd_pion_p =   posdata_df.Filter(t_cut + " && " + phi_cut + " && " + gen_Pion_rand_cut).Histo1D(model_Rnd_pion_p, MMK);
        auto kaon_h_n =     negdata_df.Filter(t_cut + " && " + phi_cut + " && " + gen_Kaon_cut).Histo1D(model_kaon_hist_n, MMK);
        auto pion_h_n =     negdata_df.Filter(t_cut + " && " + phi_cut + " && " + gen_Pion_cut).Histo1D(model_pion_hist_n, MMK);
        auto rnd_n =        negdata_df.Filter(t_cut + " && " + phi_cut + " && " + gen_Kaon_rand_cut).Histo1D(model_Rnd_n, MMK);
        auto rnd_pion_n =   negdata_df.Filter(t_cut + " && " + phi_cut + " && " + gen_Pion_rand_cut).Histo1D(model_Rnd_pion_n, MMK);
        auto dum_p =        posdum_df.Filter(t_cut + " && " + phi_cut + " && " + gen_Kaon_cut).Histo1D(model_Dum_p, MMK);
        auto dum_pion_p =   posdum_df.Filter(t_cut + " && " + phi_cut + " && " + gen_Pion_cut).Histo1D(model_Dum_pion_p, MMK);
        auto dum_n =        negdum_df.Filter(t_cut + " && " + phi_cut + " && " + gen_Kaon_cut).Histo1D(model_Dum_n, MMK);
        auto dum_pion_n =   negdum_df.Filter(t_cut + " && " + phi_cut + " && " + gen_Pion_cut).Histo1D(model_Dum_pion_n, MMK);

        std::cout << "filling histograms..." << std::endl;
        
        TH1D *kaon_hist_p = nullptr;
        TH1D *pion_hist_p = nullptr;
        TH1D *rnd_hist_p = nullptr;
        TH1D *rnd_hist_pion_p = nullptr;
        TH1D *kaon_hist_n = nullptr;
        TH1D *pion_hist_n = nullptr;
        TH1D *rnd_hist_n = nullptr;
        TH1D *rnd_hist_pion_n = nullptr;
        TH1D *dum_hist_p = nullptr;
        TH1D *dum_hist_pion_p = nullptr;
        TH1D *dum_hist_n = nullptr;
        TH1D *dum_hist_pion_n = nullptr;

        // Filing ROOT histograms for calculations. Auto initialized don't have a defined type

        if (particle == "kaonL" || particle == "kaonS") {
            kaon_hist_p = new TH1D(*kaon_h_p);
            rnd_hist_p = new TH1D(*rnd_p);

            kaon_hist_n = new TH1D(*kaon_h_n);
            rnd_hist_n = new TH1D(*rnd_n);

            dum_hist_p = new TH1D(*dum_p);
            dum_hist_n = new TH1D(*dum_n);
        }

        pion_hist_p = new TH1D(*pion_h_p);
        rnd_hist_pion_p = new TH1D(*rnd_pion_p);
        pion_hist_n = new TH1D(*pion_h_n);
        rnd_hist_pion_n = new TH1D(*rnd_pion_n);
        dum_hist_pion_p = new TH1D(*dum_pion_p);
        dum_hist_pion_n = new TH1D(*dum_pion_n);

        std::cout << "applying corrections and calculating asymmetry..." << std::endl;


        TH1D *kaon_hist_sub_p = nullptr;
        TH1D *kaon_hist_sub_n = nullptr;

        // Create a copy for multiple subtractions. Preserve the original kaonhistogram
        if (particle == "kaonL" || particle == "kaonS") {
            kaon_hist_sub_p = (TH1D*)kaon_hist_p->Clone();
            kaon_hist_sub_n = (TH1D*)kaon_hist_n->Clone();
        }

        // Rnd scaling
        // Double_t prompt = 2;
        Double_t prompt = 6;
        Double_t rand = 12;
        Double_t s = prompt/rand;

        // Dummy charge scaling scaling
        Double_t s_dum = 0;
        if (shms_pos == "center") {s_dum = 2803415/(4.8579*298159);}
        if (shms_pos == "right") {s_dum = 1876555/(4.8579*341334);}
        if (shms_pos == "left") {s_dum = 2611276/(4.8579*305498);}

        // Subtract the random from kaons
        if (particle == "kaonL" || particle == "kaonS") {
            kaon_hist_sub_p->Add(ptrScale(rnd_hist_p, s), -1);
            kaon_hist_sub_n->Add(ptrScale(rnd_hist_n, s), -1);

            // Subtract the dummy data from kaons
            kaon_hist_sub_p->Add(ptrScale(dum_hist_p, s_dum), -1);
            kaon_hist_sub_n->Add(ptrScale(dum_hist_n, s_dum), -1);
        }


        // Correct the pion for rnd
        pion_hist_p->Add(ptrScale(rnd_hist_pion_p, s*1.5), -1);
        pion_hist_n->Add(ptrScale(rnd_hist_pion_n, s*1.5), -1);
        // Correct the pion for dummy
        pion_hist_p->Add(ptrScale(dum_hist_pion_p, s_dum), -1);
        pion_hist_n->Add(ptrScale(dum_hist_pion_n, s_dum), -1);


        // Scale the pion MMK to Kaon MMK
        if (particle == "kaonL" || particle == "KaonS") {
            kaon_hist_sub_p->Add(IntegralScale(kaon_hist_p, pion_hist_p, 0.92, 0.905), -1);
            kaon_hist_sub_n->Add(IntegralScale(kaon_hist_n, pion_hist_n, 0.92, 0.905), -1);
        }
    

        if (particle == "kaonL" || particle == "kaonS") {
            SaveAsymetry(particle, kaon_hist_sub_p, kaon_hist_sub_n, asym_hist, i, csv_name);
        } else if (particle == "pion") {
            SaveAsymetry(particle, pion_hist_p, pion_hist_n, asym_hist, i, csv_name);
        } else if (particle == "dummyK") {
            SaveAsymetry(particle, dum_hist_p, dum_hist_n, asym_hist, i, csv_name);
        } else if (particle == "dummyP") {
            SaveAsymetry(particle, dum_hist_pion_p, dum_hist_pion_n, asym_hist, i, csv_name);
        }

        // Clean memory and speed up execution.
        delete kaon_hist_p;
        delete pion_hist_p;
        delete rnd_hist_p;
        delete rnd_hist_pion_p;
        delete kaon_hist_n;
        delete pion_hist_n;
        delete rnd_hist_n;
        delete rnd_hist_pion_n;
        delete dum_hist_p;
        delete dum_hist_pion_p;
        delete dum_hist_n;
        delete dum_hist_pion_n;
        delete kaon_hist_sub_p;
        delete kaon_hist_sub_n;


        auto check2 = std::chrono::high_resolution_clock::now();
        std::cout << "Phi bin #" << i << " time: " << std::chrono::duration<double>(check2 - check1).count() << "s" << std::endl;
    }

    // Draw the asymmetry histogram on a separate canvas
    TCanvas *asym_canv = new TCanvas("c4","asymmetry",1600,950);
    asym_canv->cd();
    asym_hist->DrawClone("E");
    // Fit a sin function for later value extraction
    TF1 *sinfit = new TF1("fit1","[0]*sin(x)");
    asym_hist->Fit(sinfit);

    SaveCanvas(asym_canv, particle, shms_pos, "-t"+std::to_string(t_high).substr(1)+"-"+std::to_string(t_low).substr(1)); // substr to remove '-'

    // Timing benchmark end
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "\nTotal elapsed time: " 
    << std::chrono::duration<double>(end - start).count()
    << " s\n\n";

    delete asym_hist;
    return 0;
}