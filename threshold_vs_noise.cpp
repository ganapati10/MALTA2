#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <algorithm>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>

#include "TROOT.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TSystem.h"

#include "AtlasStyle.h"

struct DataRow {
    double ITHR;
    double IDB;
    double ICASN;
    double TH_MEAN;
    double NOISE_MEAN;
    double TH_SD;
    double NOISE_SD;
};

bool isValid(const std::string& s) {
    return !(s.empty() || s == "m" || s == " ");
}

std::vector<DataRow> readCSV(std::string filename)
{
    std::ifstream file(filename);
    std::vector<DataRow> data;

    if (!file.is_open()) {
        std::cerr << "Cannot open file " << filename << std::endl;
        return data;
    }

    std::string line;
    std::getline(file, line); // skip header

    while (std::getline(file, line)) {

        std::stringstream ss(line);
        std::string item;
        std::vector<std::string> tokens;

        while (std::getline(ss, item, ',')) tokens.push_back(item);

        if (tokens.size() < 7) continue;
        if (!isValid(tokens[3]) || !isValid(tokens[4]) || !isValid(tokens[5]) || !isValid(tokens[6])) continue;

        DataRow r;

        r.ITHR       = std::stod(tokens[0]);
        r.IDB        = std::stod(tokens[1]);
        r.ICASN      = std::stod(tokens[2]);
        r.TH_MEAN    = std::stod(tokens[3]);
        r.TH_SD      = std::stod(tokens[4]);
        r.NOISE_MEAN = std::stod(tokens[5]);
        r.NOISE_SD   = std::stod(tokens[6]);

        data.push_back(r);
    }

    file.close();
    return data;
}

// Helper to get the scanned parameter value
double getScanValue(const DataRow& r, const std::string& scanVar) {
    if (scanVar == "ITHR")  return r.ITHR;
    if (scanVar == "IDB")   return r.IDB;
    if (scanVar == "ICASN") return r.ICASN;
    return 0.0;
}

void drawAtlasPlot(const std::vector<DataRow>& irr,
                   const std::vector<DataRow>& unir,
                   const std::string& fixedLabel,
                   const std::string& filename,
                   const std::string& scanVar)
{
    if (irr.empty() && unir.empty()) return;

    TCanvas* c = new TCanvas(filename.c_str(), filename.c_str(), 800, 700);
    c->SetTicks();

    double xmin = 1e9, xmax = -1e9;
    double ymin = 1e9, ymax = -1e9;

    auto updateRange = [&](const std::vector<DataRow>& rows) {
        for (auto& r : rows) {
            xmin = std::min(xmin, r.TH_MEAN - r.TH_SD);
            xmax = std::max(xmax, r.TH_MEAN + r.TH_SD);
            ymin = std::min(ymin, r.NOISE_MEAN - r.NOISE_SD);
            ymax = std::max(ymax, r.NOISE_MEAN + r.NOISE_SD);
        }
    };

    updateRange(irr);
    updateRange(unir);

    if (xmin >= xmax || ymin >= ymax) {
        std::cerr << "Invalid plotting range for " << filename << std::endl;
        return;
    }

    double dx = 0.05 * (xmax - xmin);
    double dy = 0.05 * (ymax - ymin);

    //TH1F* frame = c->DrawFrame(xmin - dx, ymin - dy, xmax + dx, ymax * 1.35);
    double ymin_fixed = 0.0;
    TH1F* frame = c->DrawFrame(xmin - dx, ymin_fixed, xmax + dx, ymax * 2.0);
    frame->SetTitle("");
    frame->GetXaxis()->SetTitle("Threshold (e^{-})");
    frame->GetYaxis()->SetTitle("Noise (e^{-})");

    std::vector<int> colors = {kBlue+1, kRed+1, kGreen+2, kOrange+7, kViolet+1, kCyan+2, kMagenta+1};

    std::vector<TGraphErrors*> irrGraphs;
    std::vector<TGraphErrors*> unirGraphs;

    // -------- Irradiated --------
    for (size_t i = 0; i < irr.size(); i++) {

        TGraphErrors* g = new TGraphErrors(1);

        g->SetPoint(0, irr[i].TH_MEAN, irr[i].NOISE_MEAN);
        g->SetPointError(0, irr[i].TH_SD, irr[i].NOISE_SD);

        int color = colors[i % colors.size()];

        g->SetMarkerStyle(20); // filled circle
        g->SetMarkerSize(1.5);
        g->SetMarkerColor(color);
        g->SetLineColor(color);

        g->Draw("P SAME");
        irrGraphs.push_back(g);
    }

    // -------- Unirradiated --------
    for (size_t i = 0; i < unir.size(); i++) {

        TGraphErrors* g = new TGraphErrors(1);

        g->SetPoint(0, unir[i].TH_MEAN, unir[i].NOISE_MEAN);
        g->SetPointError(0, unir[i].TH_SD, unir[i].NOISE_SD);

        int color = colors[i % colors.size()];

        g->SetMarkerStyle(22); // triangle
        g->SetMarkerSize(1.5);
        g->SetMarkerColor(color);
        g->SetLineColor(color);

        g->Draw("P SAME");
        unirGraphs.push_back(g);
    }

    // -------- Legend 1: Irradiated / Unirradiated (vertical) --------
TLegend* leg1 = new TLegend(0.18, 0.70, 0.38, 0.88);
leg1->SetBorderSize(0);
leg1->SetFillStyle(0);
leg1->SetTextFont(42);
leg1->SetTextSize(0.032);

TGraph* irrMarker = new TGraph();
irrMarker->SetMarkerStyle(20);
irrMarker->SetMarkerSize(1.5);

TGraph* unirMarker = new TGraph();
unirMarker->SetMarkerStyle(22);
unirMarker->SetMarkerSize(1.5);

leg1->AddEntry(irrMarker, "Irradiated", "p");
leg1->AddEntry(unirMarker, "Unirradiated", "p");

leg1->Draw();


// -------- Legend 2: Scan variable (beside it) --------
TLegend* leg2 = new TLegend(0.42, 0.55, 0.72, 0.88);
leg2->SetBorderSize(0);
leg2->SetFillStyle(0);
leg2->SetTextFont(42);
leg2->SetTextSize(0.032);

// Add scan values (colors)
for (size_t i = 0; i < irr.size(); i++) {
    double val = getScanValue(irr[i], scanVar);
    leg2->AddEntry(irrGraphs[i],
                   Form("%s = %.0f", scanVar.c_str(), val),
                   "p");
}

leg2->Draw();

    // -------- Label --------
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.03);
    latex.DrawLatex(0.18, 0.92, fixedLabel.c_str());

    c->SaveAs(filename.c_str());
    delete c;
}

void threshold_vs_noise()
{
    gROOT->SetBatch(kFALSE);

    gROOT->LoadMacro("AtlasStyle.C");
    // SetAtlasStyle();

    gSystem->mkdir("plots_thvsnoise", kTRUE);

    auto irr  = readCSV("iitm_MALTA2testing-threshold_noise_irrediated.csv");
    auto unir = readCSV("iitm_MALTA2testing-threshold_noise_unirrediated.csv");

    // =========================================================
    // 1) ITHR scan  -> fix IDB, ICASN
    // =========================================================
    std::map<std::pair<double,double>, std::vector<DataRow>> irrMapITHR, unirMapITHR;

    for (auto& r : irr)  irrMapITHR[{r.IDB, r.ICASN}].push_back(r);
    for (auto& r : unir) unirMapITHR[{r.IDB, r.ICASN}].push_back(r);

    for (auto& entry : irrMapITHR) {

        auto key = entry.first;

        auto rowsIrr  = entry.second;
        auto rowsUnir = unirMapITHR[key];

        std::sort(rowsIrr.begin(), rowsIrr.end(),
                  [](auto& a, auto& b) { return a.ITHR < b.ITHR; });

        std::sort(rowsUnir.begin(), rowsUnir.end(),
                  [](auto& a, auto& b) { return a.ITHR < b.ITHR; });

        drawAtlasPlot(rowsIrr,
                      rowsUnir,
                      Form("IDB = %.0f, ICASN = %.0f", key.first, key.second),
                      Form("plots_thvsnoise/Scan_ITHR_IDB%.0f_ICASN%.0f.pdf",
                           key.first, key.second),
                      "ITHR");
    }

    // =========================================================
    // 2) IDB scan  -> fix ITHR, ICASN
    // =========================================================
    std::map<std::pair<double,double>, std::vector<DataRow>> irrMapIDB, unirMapIDB;

    for (auto& r : irr)  irrMapIDB[{r.ITHR, r.ICASN}].push_back(r);
    for (auto& r : unir) unirMapIDB[{r.ITHR, r.ICASN}].push_back(r);

    for (auto& entry : irrMapIDB) {

        auto key = entry.first;

        auto rowsIrr  = entry.second;
        auto rowsUnir = unirMapIDB[key];

        std::sort(rowsIrr.begin(), rowsIrr.end(),
                  [](auto& a, auto& b) { return a.IDB < b.IDB; });

        std::sort(rowsUnir.begin(), rowsUnir.end(),
                  [](auto& a, auto& b) { return a.IDB < b.IDB; });

        drawAtlasPlot(rowsIrr,
                      rowsUnir,
                      Form("ITHR = %.0f, ICASN = %.0f", key.first, key.second),
                      Form("plots_thvsnoise/Scan_IDB_ITHR%.0f_ICASN%.0f.pdf",
                           key.first, key.second),
                      "IDB");
    }

    // =========================================================
    // 3) ICASN scan  -> fix ITHR, IDB
    // =========================================================
    std::map<std::pair<double,double>, std::vector<DataRow>> irrMapICASN, unirMapICASN;

    for (auto& r : irr)  irrMapICASN[{r.ITHR, r.IDB}].push_back(r);
    for (auto& r : unir) unirMapICASN[{r.ITHR, r.IDB}].push_back(r);

    for (auto& entry : irrMapICASN) {

        auto key = entry.first;

        auto rowsIrr  = entry.second;
        auto rowsUnir = unirMapICASN[key];

        std::sort(rowsIrr.begin(), rowsIrr.end(),
                  [](auto& a, auto& b) { return a.ICASN < b.ICASN; });

        std::sort(rowsUnir.begin(), rowsUnir.end(),
                  [](auto& a, auto& b) { return a.ICASN < b.ICASN; });

        drawAtlasPlot(rowsIrr,
                      rowsUnir,
                      Form("ITHR = %.0f, IDB = %.0f", key.first, key.second),
                      Form("plots_thvsnoise/Scan_ICASN_ITHR%.0f_IDB%.0f.pdf",
                           key.first, key.second),
                      "ICASN");
    }

    std::cout << "Plots created successfully\n";
}
