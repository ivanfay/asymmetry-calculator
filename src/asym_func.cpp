#include "asym_head.h"

// Function to change pads on a canvas with margins to prevent title overlap
void CngPad(TCanvas* canv, Int_t n) {
    canv->cd(n);
    canv->cd(n)->SetRightMargin(0.2);
	canv->cd(n)->SetLeftMargin(0.1);
}

// Scale 2 histograms by integrating over a specific region
// Histograms nust have same binning and range
// good for subtraction
TH1D* IntegralScale(TH1D* hist1, TH1D* hist2, Double_t high, Double_t low){
    ConfigManager& config = ConfigManager::getInstance();
    // Arbitrary variables
    Double_t x,y,z;

    // Find low and high bin numbers. Based on histogram 2
    Int_t mmkLow = hist2->FindBin(low);
    Int_t mmkHigh = hist2->FindBin(high);
    
    // Count the regions on both histogram 
    y = hist2->Integral(mmkLow, mmkHigh);
    x = hist1->Integral(mmkLow, mmkHigh);
    // Conversion factor. Ratio of hist1 to hist2
    z = x/y;
    
    z = config.getPI_sub_scale(); // This is main scaling, everything above is optional for integral. 
    // Safety cloning
    TH1D* hist2C = (TH1D*)hist2->Clone();

    // scale and return histogram 2 fit to histogram 1
    hist2C->Scale(z);
    return hist2C;
}

// One line, safe, scaler
TH1D* ptrScale(TH1D* hist, Double_t scale){
    TH1D* histC = (TH1D*)hist->Clone();
    histC->Scale(scale);
    return histC;
}

// Calculate assymetry between positive and negative histograms, make sure to set polarity
// The asymmetry is calculated in a specific region of interest
// Histogram should have same binning and range
Double_t assymetry(std::vector<TH1D*> histP_vec, std::vector<TH1D*> histN_vec, Double_t low, Double_t high, Double_t pol, std::string ptcl){
    // Find bins for the range selected
    Int_t bin_low = histP_vec[0]->FindBin(low);
    Int_t bin_high = histP_vec[0]->FindBin(high);

    TH1D* gross_histP = (TH1D*)histP_vec[0]->Clone(); 
    TH1D* gross_histN = (TH1D*)histN_vec[0]->Clone();

    if (ptcl == "kaonL" || ptcl == "kaonS") {
        gross_histP->Add(histP_vec[1], -1);
        gross_histN->Add(histN_vec[1], -1);

        gross_histP->Add(histP_vec[2], -1);
        gross_histN->Add(histN_vec[2], -1);
    }

    // Find total number of counts in range
    Double_t Y_p = gross_histP->Integral(bin_low, bin_high);
    Double_t Y_n = gross_histN->Integral(bin_low, bin_high);
    
    // Calculate and return assymetry in region
    Double_t asym;
    
    if (Y_n + Y_p == 0) {
        asym = 0;
    } else {
        asym = (1/pol) * ((Y_p - Y_n) / (abs(Y_p + Y_n)));
    }
    return asym;
}

// Same function as assymetry, but for the error of assymetry
std::array<Double_t, 2> assymetry_error(std::vector<TH1D*> histP_vec, std::vector<TH1D*> histN_vec, Double_t low, Double_t high, Double_t pol, std::string ptcl){
    // Find bins for the range selected
    Int_t bin_low = histP_vec[0]->FindBin(low);
    Int_t bin_high = histP_vec[0]->FindBin(high);

    std::array<Double_t, 2> err;
    Double_t pol_err[2] = {pol*0.03, pol*0.01}; // i = 0 is -3% polarization error, i = 1 is +1% polarization error
    Double_t asym;

    if (ptcl == "kaonL" || ptcl == "kaonS") {
        Double_t rawY_p, pionY_p, rndY_p, rawY_n, pionY_n, rndY_n;

        rawY_p = histP_vec[0]->Integral(bin_low, bin_high);
        pionY_p = histP_vec[1]->Integral(bin_low, bin_high);
        rndY_p = histP_vec[2]->Integral(bin_low, bin_high);

        rawY_n = histN_vec[0]->Integral(bin_low, bin_high);
        pionY_n = histN_vec[1]->Integral(bin_low, bin_high);
        rndY_n = histN_vec[2]->Integral(bin_low, bin_high);

        Double_t Y_p = rawY_p - pionY_p - rndY_p;
        Double_t Y_n = rawY_n - pionY_n - rndY_n;
        Double_t sum_p = rawY_p + pionY_p + rndY_p;
        Double_t sum_n = rawY_n + pionY_n + rndY_n;
        // error equation parts

        if (Y_n + Y_p == 0) {
            asym = 0;
        } else {
            asym = (1/pol) * ((Y_p - Y_n) / (abs(Y_p + Y_n)));
        }

        Double_t term1 = Y_n*Y_n*sum_p + Y_p*Y_p*sum_n;
        Double_t term1_coef = 4/(pol*pol*pow(Y_p + Y_n, 4));
        Double_t term2_H_pol_err = pow(asym, 2)*pow(pol_err[1], 2);
        Double_t term2_L_pol_err = pow(asym, 2)*pow(pol_err[0], 2);
        
        if ((pol*pow(Y_p + Y_n, 2)) == 0 || ((Y_n*Y_n)*sum_p + (Y_p*Y_p)*sum_n) < 0) {
            err[0] = 10;
            err[1] = 10;
        } else {
            err[0] = sqrt(term1_coef*term1 + term2_L_pol_err);
            err[1] = sqrt(term1_coef*term1 + term2_H_pol_err);
        }

        
    } else {
        Double_t Y_p = histP_vec[0]->Integral(bin_low, bin_high);
        Double_t Y_n = histN_vec[0]->Integral(bin_low, bin_high);

        if (Y_p*Y_n == 0 || Y_n + Y_p == 0) {
            err[0] = 10;
            err[1] = 10;
        } else {
            err[0] = (2/pol) * ((sqrt(abs(Y_p*Y_n)))/sqrt(pow(abs(Y_p+Y_n),3)));
            err[1] = (2/pol) * ((sqrt(abs(Y_p*Y_n)))/sqrt(pow(abs(Y_p+Y_n),3)));
        }
    }

    return err;
}

// Function to send yields (Y_p, Y_n) per phi bin to a CSV
// also adds assymetry 
// Histograms passed should have same binning and range
void yield_to_CSV(std::vector<TH1D*> histP_vec, std::vector<TH1D*> histN_vec, Double_t low, Double_t high, Int_t phi_bin, std::string ptcl, const std::string filename){
    // Open or create a CSV file. Append to the existing file
    // Best to create a new file before using this function
    std::ofstream out(filename, std::ios::app);

    // Check if header has been written
    static bool header_written = false;

    Int_t bin_low = histP_vec[0]->FindBin(low);
    Int_t bin_high = histP_vec[0]->FindBin(high);

    Double_t pol = 0.38;
    Double_t pol_err[2] = {pol*0.03, pol*0.01}; // i = 0 is -3% polarization error, i = 1 is +1% polarization error
    Double_t asym;
    Double_t err[2];

    if (ptcl == "kaonL" || ptcl == "kaonS") {
        Double_t rawY_p, pionY_p, rndY_p, rawY_n, pionY_n, rndY_n;

        rawY_p = histP_vec[0]->Integral(bin_low, bin_high);
        pionY_p = histP_vec[1]->Integral(bin_low, bin_high);
        rndY_p = histP_vec[2]->Integral(bin_low, bin_high);

        rawY_n = histN_vec[0]->Integral(bin_low, bin_high);
        pionY_n = histN_vec[1]->Integral(bin_low, bin_high);
        rndY_n = histN_vec[2]->Integral(bin_low, bin_high);

        Double_t Y_p = rawY_p - pionY_p - rndY_p;
        Double_t Y_n = rawY_n - pionY_n - rndY_n;
        Double_t sum_p = rawY_p + pionY_p + rndY_p;
        Double_t sum_n = rawY_n + pionY_n + rndY_n;

        // error equation parts

        if (Y_n + Y_p == 0) {
            asym = 0;
        } else {
            asym = (1/pol) * ((Y_p - Y_n) / (abs(Y_p + Y_n)));
        }

        Double_t term1 = Y_n*Y_n*sum_p + Y_p*Y_p*sum_n;
        Double_t term1_coef = 4/(pol*pol*pow(Y_p + Y_n, 4));
        Double_t term2_H_pol_err = pow(asym, 2)*pow(pol_err[1], 2);
        Double_t term2_L_pol_err = pow(asym, 2)*pow(pol_err[0], 2);

        if ((pol*pow(Y_p + Y_n, 2)) == 0 || Y_p + Y_n == 0 || ((Y_n*Y_n)*sum_p + (Y_p*Y_p)*sum_n) < 0) {
            err[0] = 10;
            err[1] = 10;
        } else {
            err[0] = sqrt(term1_coef*term1 + term2_L_pol_err);
            err[1] = sqrt(term1_coef*term1 + term2_H_pol_err);
        }

        if (!header_written){
            out << "Phi_bin, raw Yield p, raw Yield n, pion Yield p, pion Yield n, rnd Yield p, rnd Yield n, net Yield p, net Yield n, Asym, Asym_errL, Asym_errH  \n";
            header_written = true;
        }

        // Write 1 phi bin entry to the CSV file
        out << phi_bin << ", " << rawY_p << ", " << rawY_n << ", " << pionY_p << ", " << pionY_n << ", " << rndY_p << ", " << rndY_n << ", " <<  Y_p << ", " << Y_n << ", " << asym << ", " << err[0] << ", " << err[1] << "\n";
        out.close();

    } else {
        // Asymetry and error calculation
        Double_t Y_p = histP_vec[0]->Integral(bin_low, bin_high);
        Double_t Y_n = histN_vec[0]->Integral(bin_low, bin_high);


        if (Y_p*Y_n == 0 || Y_n + Y_p == 0) {
            asym = 0;
            err[0] = 0;
            err[1] = 0;
        } else {
            asym = 1/pol * ((Y_p - Y_n) / (abs(Y_p + Y_n)));
            err[0] = (2/pol)*((sqrt(abs(Y_p*Y_n)))/sqrt(pow(abs(Y_p+Y_n),3)));
            err[1] = (2/pol)*((sqrt(abs(Y_p*Y_n)))/sqrt(pow(abs(Y_p+Y_n),3)));
        }

        if (!header_written){
            out << "Phi_bin, Yield_p, Yield_n, Asym, Asym_errL, Asym_errH \n";
            header_written = true;
        }

        // Write 1 phi bin entry to the CSV file
        out << phi_bin << ", " << Y_p << ", " << Y_n << ", " << asym << ", " << err[0] << ", " << err[1] << "\n";
        out.close();
    }
}


void FormatInput(std::string &infile, std::string &particle, std::string &shms_pos) {
    infile = "input/" + infile + ".root";

    if (particle == "Kaon"
        || particle == "kaon"
        || particle == "KaonL" 
        || particle == "kaonL" 
        || particle == "Kaonl" 
        || particle == "kaonl" 
        || particle == "Kaon Lambda"
        || particle == "kaon_lambda"
        || particle == "kaon lambda") {particle = "kaonL";}
    else if (particle == "KaonS"
        || particle == "kaonS" 
        || particle == "Kaons" 
        || particle == "kaons" 
        || particle == "Kaon Sigma"
        || particle == "kaon_sigma"
        || particle == "kaon sigma") {particle = "kaonS";}
    else if (particle == "Pion"
        || particle == "pion") {particle = "pion";}
    else if (particle == "DummyK"
        || particle == "dummyK"
        || particle == "Dummyk"
        || particle == "dummyk") {particle = "dummyK";}
    else if (particle == "DummyP"
        || particle == "dummyP"
        || particle == "Dummyp"
        || particle == "dummyp") {particle = "dummyP";}
    else {
        std::cout << "Particle " + particle + " is not valid. Try: kaonL, kaonS, pion" << std::endl;
        infile = "";
    }

    if (shms_pos == "Center"
        || shms_pos == "center"
        || shms_pos == "C"
        || shms_pos == "c") {shms_pos = "center";}
    else if (shms_pos == "Right"
        || shms_pos == "right"
        || shms_pos == "R"
        || shms_pos == "r") {shms_pos = "right";}
    else if (shms_pos == "Left"
        || shms_pos == "left"
        || shms_pos == "L"
        || shms_pos == "l") {shms_pos = "left";}
    else {
        std::cout << "SHMS "+shms_pos+" postion is not valid. Try: center, right, left" << std::endl;
        infile = "";
    }

}

// TODO: Make this into a config file that dynamicly assigns user's root file names
void CheckNaming(std::string who_named) {
    std::string root_naming = who_named;
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
    }
    if (root_naming == "Alicia" || root_naming == "alicia") {
        H_cal_etottracknorm = "H_cal_etottracknorm";
        H_cer_npeSum = "H_cer_npeSum";
        P_hgcer_npeSum = "P_hgcer_npeSum";
        P_aero_npeSum = "P_aero_npeSum";
        P_RF_Dist = "P_RF_Dist";
        P_gtr_xp = "P_gtr_xp";
        CTime_eKCoinTime_ROC1 = "CTime_eKCoinTime_ROC1";
        // CTime_eKCoinTime_ROC1 = "CTime_ePiCoinTime_ROC1";
        MMK = "MMK";
        // MMK = "MMpi";
        ph_q = "ph_q";
        MandelT = "MandelT";
    }
}

void SaveAsymmetry(std::string ptcl, std::vector<TH1D*> histP_vec, std::vector<TH1D*> histN_vec, TGraphAsymmErrors* storage, Double_t phi_bin_mid, Int_t phi_bin, std::string filename) { 
    ConfigManager& config = ConfigManager::getInstance();

    Double_t low, high, polarity = 0.38;
    if (ptcl == "kaonL") {
        // shifted +0.002 because data was moved by the shift amount to fit the simC
        low = config.getMMK_low();
        high = config.getMMK_high();
        // low = 1.12;
        // high = 1.168;

        std::array<Double_t, 2> asym_err = assymetry_error(histP_vec, histN_vec, low, high, polarity, ptcl);

        // Integrate kaon output and combine positive and negative to calculate the assimetry. Save in histogram
        storage->SetPoint(phi_bin, phi_bin_mid, assymetry(histP_vec, histN_vec, low, high, polarity, ptcl));
        storage->SetPointError(phi_bin, 0, 0, asym_err[0], asym_err[1]);

        yield_to_CSV(histP_vec, histN_vec, low, high, phi_bin, ptcl,filename);

    } else if (ptcl == "kaonS") {
        low = config.getMMK_low();
        high = config.getMMK_high();

        std::array<Double_t, 2> asym_err = assymetry_error(histP_vec, histN_vec, low, high, polarity, ptcl);

        // Integrate kaon output and combine positive and negative to calculate the assimetry. Save in histogram
        storage->AddPoint(phi_bin_mid, assymetry(histP_vec, histN_vec, low, high, polarity, ptcl));
        storage->SetPointError(phi_bin_mid, 0, 0, asym_err[0], asym_err[1]);

        yield_to_CSV(histP_vec, histN_vec, low, high, phi_bin, ptcl, filename);

    } else if (ptcl == "pion") {
        low = config.getMMK_low();
        high = config.getMMK_high();

        std::array<Double_t, 2> asym_err = assymetry_error(histP_vec, histN_vec, low, high, polarity, ptcl);

        storage->AddPoint(phi_bin, assymetry(histP_vec, histN_vec, low, high, polarity, ptcl));
        storage->SetPointError(phi_bin, 0, 0, asym_err[0], asym_err[1]);

        yield_to_CSV(histP_vec, histN_vec, low, high, phi_bin, ptcl, filename);
    } else if (ptcl == "dummyK" || ptcl == "dummyP") {
        low = config.getMMK_low();
        high = config.getMMK_high();

        std::array<Double_t, 2> asym_err = assymetry_error(histP_vec, histN_vec, low, high, polarity, ptcl);

        storage->AddPoint(phi_bin, assymetry(histP_vec, histN_vec, low, high, polarity, ptcl));
        storage->SetPointError(phi_bin, 0, 0, asym_err[0], asym_err[1]);

        yield_to_CSV(histP_vec, histN_vec, low, high, phi_bin, ptcl, filename);

    }
    else {
        std::cout << "Particle not supported. Try: kaonL, kaonS, pion" << std::endl;
    }

}

void SaveCanvas(TCanvas* canv, std::string particle, std::string shms_pos, std::string t_bin) {
    std::string filename = "output/Asym_"+particle+"_"+shms_pos+"_"+t_bin+".root";
    canv->SaveAs(filename.c_str());

    std::cout << "Compilation of "+particle+" "+shms_pos+" data for "+t_bin+" is complete." << std::endl;
}