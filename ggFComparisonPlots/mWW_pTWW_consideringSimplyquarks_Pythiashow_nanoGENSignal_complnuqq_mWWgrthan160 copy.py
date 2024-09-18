import ROOT
import argparse
import time

# Make global style change
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetTextFont(42)
ROOT.gStyle.SetTextSize(0.05)

# Parse the files using argparse 
parser = argparse.ArgumentParser(description="Provide root files")
parser.add_argument("--rootfile1", type=str, help="Path to the first root file")
parser.add_argument("--rootfile2", type=str, help="Path to the second root file")
args = parser.parse_args()
rootfile1 = args.rootfile1
rootfile2 = args.rootfile2

# Start the timer
start = time.time()

# Open the first root file and get the TTree
file1 = ROOT.TFile(rootfile1, "READ")
tree1 = file1.Get("Events")

# Open the second root file and get the TTree
file2 = ROOT.TFile(rootfile2, "READ")
tree2 = file2.Get("Events")

# Creating 1D histograms to store the mWW distributions for both samples
hist_mWW_1 = ROOT.TH1F("hist_mWW_1", "Distribution of mWW - Sample 1", 100, 0, 1000)
hist_mWW_2 = ROOT.TH1F("hist_mWW_2", "Distribution of mWW - Sample 2", 100, 0, 1000)

# Loop over events in the first TTree
num_events_tree1 = tree1.GetEntries()
for i in range(num_events_tree1):
    tree1.GetEntry(i)
    
    # Collect quarks
    quark_coll = []
    nu_coll = []
    lep_coll = []
    lep, nu, jet1, jet2, qfat, met = [ROOT.TLorentzVector() for _ in range(6)]
    
    for part in range(tree1.nGenPart):
        if abs(tree1.GenPart_pdgId[tree1.GenPart_genPartIdxMother[part]]) == 24:
            if abs(tree1.GenPart_pdgId[part]) in [11, 13, 15]:
                lep_coll.append(part)
            elif abs(tree1.GenPart_pdgId[part]) <= 4:
                quark_coll.append(part)
            elif abs(tree1.GenPart_pdgId[part]) in [12, 14, 16]:
                nu_coll.append(part)
    
    # Require exactly 1 lepton and 1 neutrino
    if len(nu_coll) != 1 or len(lep_coll) != 1 or len(quark_coll) != 2:
        continue

    # Calculate and fill histograms
    # Build the TLorentzVectors
    lep.SetPtEtaPhiM(tree1.GenPart_pt[lep_coll[0]], tree1.GenPart_eta[lep_coll[0]], tree1.GenPart_phi[lep_coll[0]], tree1.GenPart_mass[lep_coll[0]])
    nu.SetPtEtaPhiM(tree1.GenPart_pt[nu_coll[0]], tree1.GenPart_eta[nu_coll[0]], tree1.GenPart_phi[nu_coll[0]], tree1.GenPart_mass[nu_coll[0]])
    jet1.SetPtEtaPhiM(tree1.GenPart_pt[quark_coll[0]], tree1.GenPart_eta[quark_coll[0]], tree1.GenPart_phi[quark_coll[0]], tree1.GenPart_mass[quark_coll[0]])
    jet2.SetPtEtaPhiM(tree1.GenPart_pt[quark_coll[1]], tree1.GenPart_eta[quark_coll[1]], tree1.GenPart_phi[quark_coll[1]], tree1.GenPart_mass[quark_coll[1]])
    met.SetPtEtaPhiM(tree1.GenMET_pt, 0.0, tree1.GenMET_phi, 0.0)

    # Calculate quantities
    lnu = lep + nu
    jj = jet1 + jet2
    lnujj = lnu + jj 
    
    # Fill histograms
    hist_mWW_1.Fill(lnujj.M())
    

# Loop over events in the second TTree
num_events_tree2 = tree2.GetEntries()
for i in range(num_events_tree2):
    tree2.GetEntry(i)
    
    # Collect quarks
    quark_coll_2 = []
    nu_coll_2 = []
    lep_coll_2 = []
    lep, nu, jet1, jet2, qfat, met = [ROOT.TLorentzVector() for _ in range(6)]

    for part in range(tree2.nGenPart):
        if abs(tree2.GenPart_pdgId[tree2.GenPart_genPartIdxMother[part]]) == 24 :
            if abs(tree2.GenPart_pdgId[part]) in [11, 13, 15]:
                lep_coll_2.append(part)
            elif abs(tree2.GenPart_pdgId[part]) <= 4:
                quark_coll_2.append(part)
            elif abs(tree2.GenPart_pdgId[part]) in [12, 14, 16]:
                nu_coll_2.append(part)
    
    # Require exactly 1 lepton and 1 neutrino
    if len(nu_coll_2) != 1 or len(lep_coll_2) != 1 or len(quark_coll_2) != 2:
        continue

    # Calculate and fill histograms
    # Build the TLorentzVectors
    lep.SetPtEtaPhiM(tree2.GenPart_pt[lep_coll_2[0]], tree2.GenPart_eta[lep_coll_2[0]], tree2.GenPart_phi[lep_coll_2[0]], tree2.GenPart_mass[lep_coll_2[0]])
    nu.SetPtEtaPhiM(tree2.GenPart_pt[nu_coll_2[0]], tree2.GenPart_eta[nu_coll_2[0]], tree2.GenPart_phi[nu_coll_2[0]], tree2.GenPart_mass[nu_coll_2[0]])
    jet1.SetPtEtaPhiM(tree2.GenPart_pt[quark_coll_2[0]], tree2.GenPart_eta[quark_coll_2[0]], tree2.GenPart_phi[quark_coll_2[0]], tree2.GenPart_mass[quark_coll_2[0]])
    jet2.SetPtEtaPhiM(tree2.GenPart_pt[quark_coll_2[1]], tree2.GenPart_eta[quark_coll_2[1]], tree2.GenPart_phi[quark_coll_2[1]], tree2.GenPart_mass[quark_coll_2[1]])
    met.SetPtEtaPhiM(tree2.GenMET_pt, 0.0, tree2.GenMET_phi, 0.0)

    # Calculate quantities
    lnu = lep + nu
    jj = jet1 + jet2
    lnujj = lnu + jj 
    
    # Fill histograms
    hist_mWW_2.Fill(lnujj.M())

# Normalize histograms
hist_mWW_1.Scale(1.0 / hist_mWW_1.Integral())
hist_mWW_2.Scale(1.0 / hist_mWW_2.Integral())

# Set histogram styles
hist_mWW_1.SetLineColor(ROOT.kBlue)
hist_mWW_2.SetLineColor(ROOT.kRed)

# Draw histograms on a canvas
canvas = ROOT.TCanvas("canvas_mWW_comparison", "mWW Comparison", 1000, 800)
canvas.SetLogy(1)  # Set logarithmic scale for y-axis
hist_mWW_1.Draw()
hist_mWW_1.GetXaxis().SetTitle("mWW [GeV]")
hist_mWW_1.GetYaxis().SetTitle("Normalized Events")
hist_mWW_2.Draw("same")

# Add legend
#legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
#legend.AddEntry(hist_mWW_1, "Madgraph (Central)", "l")
#legend.AddEntry(hist_mWW_2, "Madgraph+Pythia (WW decay) ", "l")
#legend.Draw()

# Add legend without the outer box and with symbols
legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
legend.SetBorderSize(0)  # Remove outer box
legend.AddEntry(hist_mWW_1, "Madgraph(Central)", "l")
legend.AddEntry(hist_mWW_2, "Madgraph+Pythia(WW decay) ", "l")
#legend.SetFillStyle(0) # Remove legend fill
legend.Draw()


# Save canvas as PNG
canvas.SaveAs("mWW_comparison_steve_ggWW.pdf")

# Close the ROOT files
file1.Close()
file2.Close()

# Print processing time
print("Processing time:", time.time() - start)
