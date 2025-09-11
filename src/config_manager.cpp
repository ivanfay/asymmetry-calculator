#include "asym_head.h"

ConfigManager& ConfigManager::getInstance() {
    static ConfigManager instance;
    return instance;
}

ConfigManager::ConfigManager():
    RF_high(1.15),
    RF_low(0.2),
    H_cal_low(0.5),
    H_cer_low(2),
    CT_Kaon_low(-0.5),
    CT_Kaon_high(1.5),
    CT_Pion_low(-2.5),
    CT_Pion_high(-0.5),
    CT_Rand_low_right(9.5),
    CT_Rand_high_right(15.5),
    CT_Rand_low_left(-8.5),
    CT_Rand_high_left(-14.5),
    PI_sub_scale(0.0169), 
    MMK_low(1.080),
    MMK_high(1.1803){}

const Double_t& ConfigManager::getRF_high()           const{return RF_high;}
const Double_t& ConfigManager::getRF_low()            const{return RF_low;}
const Double_t& ConfigManager::getH_cal_low()         const{return H_cal_low;}
const Double_t& ConfigManager::getH_cer_low()         const{return H_cer_low;}
const Double_t& ConfigManager::getCT_Kaon_low()       const{return CT_Kaon_low;}
const Double_t& ConfigManager::getCT_Kaon_high()      const{return CT_Kaon_high;}
const Double_t& ConfigManager::getCT_Pion_low()       const{return CT_Pion_low;}
const Double_t& ConfigManager::getCT_Pion_high()      const{return CT_Pion_high;}
const Double_t& ConfigManager::getCT_Rand_low_right() const{return CT_Rand_low_right;}
const Double_t& ConfigManager::getCT_Rand_high_right()const{return CT_Rand_high_right;}
const Double_t& ConfigManager::getCT_Rand_low_left()  const{return CT_Rand_low_left;}
const Double_t& ConfigManager::getCT_Rand_high_left() const{return CT_Rand_high_left;}
const Double_t& ConfigManager::getPI_sub_scale()      const{return PI_sub_scale;}
const Double_t& ConfigManager::getMMK_low()          const{return MMK_low;}
const Double_t& ConfigManager::getMMK_high()         const{return MMK_high;}

// Save to JSON file
void ConfigManager::saveToJSON(const std::string& filename) const {
    nlohmann::json j;
    j["RF_high"] = RF_high;
    j["RF_low"] = RF_low;
    j["H_cal_low"] = H_cal_low;
    j["H_cer_low"] = H_cer_low;
    j["CT_Kaon_low"] = CT_Kaon_low;
    j["CT_Kaon_high"] = CT_Kaon_high;
    j["CT_Pion_low"] = CT_Pion_low;
    j["CT_Pion_high"] = CT_Pion_high;
    j["CT_Rand_low_right"] = CT_Rand_low_right;
    j["CT_Rand_high_right"] = CT_Rand_high_right;
    j["CT_Rand_low_left"] = CT_Rand_low_left;
    j["CT_Rand_high_left"] = CT_Rand_high_left;
    j["PI_sub_scale"] = PI_sub_scale;
    j["MMK_low"] = MMK_low;
    j["MMK_high"] = MMK_high;

    std::ofstream file(filename);
    if (file.is_open()) {
        file << j.dump(4); // Pretty print with 4 spaces
        file.close();
    } else {
        std::cerr << "Error opening file for writing: " << filename << std::endl;
    }
}

// Load from JSON file
void ConfigManager::loadFromJSON(const std::string& filename) {
    nlohmann::json j;
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file for reading: " << filename << std::endl;
        return;
    }
    try {
        file >> j;

        if (j.contains("RF_high")) {RF_high = j["RF_high"].get<Double_t>();}
        if (j.contains("RF_low")) {RF_low = j["RF_low"].get<Double_t>();}
        if (j.contains("H_cal_low")) {H_cal_low = j["H_cal_low"].get<Double_t>();}
        if (j.contains("H_cer_low")) {H_cer_low = j["H_cer_low"].get<Double_t>();}
        if (j.contains("CT_Kaon_low")) {CT_Kaon_low = j["CT_Kaon_low"].get<Double_t>();}
        if (j.contains("CT_Kaon_high")) {CT_Kaon_high = j["CT_Kaon_high"].get<Double_t>();}
        if (j.contains("CT_Pion_low")) {CT_Pion_low = j["CT_Pion_low"].get<Double_t>();}
        if (j.contains("CT_Pion_high")) {CT_Pion_high = j["CT_Pion_high"].get<Double_t>();}
        if (j.contains("CT_Rand_low_right")) {CT_Rand_low_right = j["CT_Rand_low_right"].get<Double_t>();}
        if (j.contains("CT_Rand_high_right")) {CT_Rand_high_right = j["CT_Rand_high_right"].get<Double_t>();}
        if (j.contains("CT_Rand_low_left")) {CT_Rand_low_left = j["CT_Rand_low_left"].get<Double_t>();}
        if (j.contains("CT_Rand_high_left")) {CT_Rand_high_left = j["CT_Rand_high_left"].get<Double_t>();}
        if (j.contains("PI_sub_scale")) {PI_sub_scale = j["PI_sub_scale"].get<Double_t>();}
        if (j.contains("MMK_low")) {MMK_low = j["MMK_low"].get<Double_t>();}
        if (j.contains("MMK_high")) {MMK_high = j["MMK_high"].get<Double_t>();}

    } catch (nlohmann::json::parse_error& e) {
        std::cerr << "JSON parse error: " << e.what() << std::endl;
        return; 
    } catch (nlohmann::json::type_error& e) {
        std::cerr << "JSON type error: " << e.what() << std::endl;
        return;
    }
}