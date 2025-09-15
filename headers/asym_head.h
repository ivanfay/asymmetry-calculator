#pragma once
#ifndef ASYM_HEAD_H
#define ASYM_HEAD_H

#include<TROOT.h>
#include "ROOT/RDataFrame.hxx"

#include<TStyle.h>
#include<TCanvas.h>
#include<TPad.h>

#include<TFile.h>

#include<TMath.h>

#include<TTree.h>
#include<TCut.h>

#include<TH1D.h>
#include<TF1.h>

#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<nlohmann/json.hpp>


class ConfigManager {
    public:
        ConfigManager(const ConfigManager&) = delete;
        ConfigManager& operator=(const ConfigManager&) = delete;

        static ConfigManager& getInstance();

        // variable accecsors
        const Double_t& getRF_high() const;
        const Double_t& getRF_low() const;
        const Double_t& getH_cal_low() const;
        const Double_t& getH_cer_low() const;
        const Double_t& getCT_Kaon_low() const;
        const Double_t& getCT_Kaon_high() const;
        const Double_t& getCT_Pion_low() const;
        const Double_t& getCT_Pion_high() const;
        const Double_t& getCT_Rand_low_right() const;
        const Double_t& getCT_Rand_high_right() const;
        const Double_t& getCT_Rand_low_left() const;
        const Double_t& getCT_Rand_high_left() const;
        const Double_t& getPI_sub_scale() const;
        const Double_t& getMMK_low() const;
        const Double_t& getMMK_high() const;

        // JSON save/load
        void saveToJSON(const std::string& filename) const;
        void loadFromJSON(const std::string& filename);
    
    private:
        ConfigManager();

        // Member variables
        Double_t RF_high;
        Double_t RF_low;
        Double_t H_cal_low;
        Double_t H_cer_low;
        Double_t CT_Kaon_low;
        Double_t CT_Kaon_high;
        Double_t CT_Pion_low;
        Double_t CT_Pion_high;
        Double_t CT_Rand_low_right;
        Double_t CT_Rand_high_right;
        Double_t CT_Rand_low_left;
        Double_t CT_Rand_high_left;
        Double_t PI_sub_scale;
        Double_t MMK_low;
        Double_t MMK_high;

};


void CngPad(TCanvas* canv, int n);

TH1D* IntegralScale(TH1D* hist1, TH1D* hist2, Double_t high, Double_t low);
TH1D* ptrScale(TH1D* hist, Double_t scale);

Double_t assymetry(std::vector<TH1D*> histP_vec, std::vector<TH1D*> histN_vec, Double_t low, Double_t high, Double_t pol, std::string ptcl);
std::array<Double_t, 2> assymetry_error(std::vector<TH1D*> histP_vec, std::vector<TH1D*> histN_vec, Double_t low, Double_t high, Double_t pol, std::string ptcl);

void yield_to_CSV(std::vector<TH1D*> histP_vec, std::vector<TH1D*> histN_vec, Double_t low, Double_t high, Int_t phi_bin, std::string ptcl, const std::string filename);

void FormatInput(std::string &input, std::string &particle, std::string &shms_pos);

void CheckNaming(std::string who_named);

void SaveAsymmetry(std::string ptcl, std::vector<TH1D*> histP_vec, std::vector<TH1D*> histN_vec, TGraphAsymmErrors* storage, Double_t phi_bin_mid, Int_t phi_bin, std::string filename);

void SaveCanvas(TCanvas* canv, std::string particle, std::string shms_pos, std::string t_bin);

const int phi_bins = 15;

extern std::string H_cal_etottracknorm;
extern std::string H_cer_npeSum;
extern std::string P_hgcer_npeSum;
extern std::string P_aero_npeSum;
extern std::string P_RF_Dist;
extern std::string P_gtr_xp;
extern std::string CTime_eKCoinTime_ROC1;
extern std::string MMK;
extern std::string ph_q;
extern std::string MandelT;


#endif