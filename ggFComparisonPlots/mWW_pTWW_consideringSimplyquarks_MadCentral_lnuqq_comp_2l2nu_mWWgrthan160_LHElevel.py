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
hist_mWW_1 = ROOT.TH1F("hist_mWW_1", "Distribution of mWW - Sample 1 (lnuqq)", 100, 0, 1000)
hist_mWW_2 = ROOT.TH1F("hist_mWW_2", "Distribution of mWW - Sample 2 (2l2nu)", 100, 0, 1000)

hist_mWW_1 = ROOT.TH1F("hist_mWW_1", "Distribution of mWW - Sample 1 (lnuqq)", 100, 0, 1000)
hist_mWW_2 = ROOT.TH1F("hist_mWW_2", "Distribution of mWW - Sample 2 (2l2nu)", 100, 0, 1000)

hist_lep_1 = ROOT.TH1F("hist_lnu_1", "Distribution of lep Pt of lnuqq vs 2l2nu ", 50,0,500 )
hist_lep_2 = ROOT.TH1F("hist_lnu_2", "Distribution of lep Pt of lnuqq vs 2l2nu ", 50,0,500 )


# Loop over events in the first TTree (lnuqq events)
num_events_tree1 = tree1.GetEntries()
for i in range(num_events_tree1):
    tree1.GetEntry(i)
    
    # Collect quarks, leptons, and neutrinos from LHE particles
    quark_coll = []
    nu_coll = []
    lep_coll = []
    lep, nu, jet1, jet2 = [ROOT.TLorentzVector() for _ in range(4)]
    
    for part in range(tree1.nLHEPart):
        if tree1.LHEPart_status[part] == 1:  # outgoing particles
            if abs(tree1.LHEPart_pdgId[part]) in [11, 13, 15]:  # Leptons
                lep_coll.append(part)
            elif abs(tree1.LHEPart_pdgId[part]) in [12, 14, 16]:  # Neutrinos
                nu_coll.append(part)
            elif abs(tree1.LHEPart_pdgId[part]) <= 4:  # Quarks
                quark_coll.append(part)
    
    # Require exactly 1 lepton, 1 neutrino, and 2 quarks
    if len(lep_coll) != 1 or len(nu_coll) != 1 or len(quark_coll) != 2:
        continue

    # Build the TLorentzVectors
    lep.SetPtEtaPhiM(tree1.LHEPart_pt[lep_coll[0]], tree1.LHEPart_eta[lep_coll[0]], tree1.LHEPart_phi[lep_coll[0]], tree1.LHEPart_mass[lep_coll[0]])
    nu.SetPtEtaPhiM(tree1.LHEPart_pt[nu_coll[0]], tree1.LHEPart_eta[nu_coll[0]], tree1.LHEPart_phi[nu_coll[0]], tree1.LHEPart_mass[nu_coll[0]])
    jet1.SetPtEtaPhiM(tree1.LHEPart_pt[quark_coll[0]], tree1.LHEPart_eta[quark_coll[0]], tree1.LHEPart_phi[quark_coll[0]], tree1.LHEPart_mass[quark_coll[0]])
    jet2.SetPtEtaPhiM(tree1.LHEPart_pt[quark_coll[1]], tree1.LHEPart_eta[quark_coll[1]], tree1.LHEPart_phi[quark_coll[1]], tree1.LHEPart_mass[quark_coll[1]])

    # Calculate mWW (lnuqq system)
    lnu = lep + nu
    jj = jet1 + jet2
    lnujj = lnu + jj
    
    # Fill histogram with mWW (lnuqq system)
    #hist_mWW_1.Fill(lnujj.M())
    hist_mWW_1.Fill(lnujj.Pt())
    hist_lep_1.Fill(lep.Pt())
    print(lnujj.M())
    print("lnu.Pt()", "lnu.M()", lnu.Pt(), lnu.M())
    print("jj.Pt()", "jj.M()", jj.Pt(), jj.M())
    print("q")
    print(lnujj.Pt(),lnujj.Px(), lnujj.Py())

# Loop over events in the second TTree (2l2nu events)
num_events_tree2 = tree2.GetEntries()
for i in range(num_events_tree2):
    tree2.GetEntry(i)
    
    # Collect leptons and neutrinos from LHE particles (2 leptons and 2 neutrinos)
    lep_coll_2 = []
    nu_coll_2 = []
    lep1, nu1, lep2, nu2 = [ROOT.TLorentzVector() for _ in range(4)]

    for part in range(tree2.nLHEPart):
        if tree2.LHEPart_status[part] == 1:  # outgoing particles
            if abs(tree2.LHEPart_pdgId[part]) in [11, 13, 15]:  # Leptons
                lep_coll_2.append(part)
                print(part)
            elif abs(tree2.LHEPart_pdgId[part]) in [12, 14, 16]:  # Neutrinos
                nu_coll_2.append(part)
    
    print("did we reach at step 2")
    # Require exactly 2 leptons and 2 neutrinos
    if len(lep_coll_2) != 2 or len(nu_coll_2) != 2:
        continue

    # Build the TLorentzVectors for the two leptons and two neutrinos
    lep1.SetPtEtaPhiM(tree2.LHEPart_pt[lep_coll_2[0]], tree2.LHEPart_eta[lep_coll_2[0]], tree2.LHEPart_phi[lep_coll_2[0]], tree2.LHEPart_mass[lep_coll_2[0]])
    lep2.SetPtEtaPhiM(tree2.LHEPart_pt[lep_coll_2[1]], tree2.LHEPart_eta[lep_coll_2[1]], tree2.LHEPart_phi[lep_coll_2[1]], tree2.LHEPart_mass[lep_coll_2[1]])
    nu1.SetPtEtaPhiM(tree2.LHEPart_pt[nu_coll_2[0]], tree2.LHEPart_eta[nu_coll_2[0]], tree2.LHEPart_phi[nu_coll_2[0]], tree2.LHEPart_mass[nu_coll_2[0]])
    nu2.SetPtEtaPhiM(tree2.LHEPart_pt[nu_coll_2[1]], tree2.LHEPart_eta[nu_coll_2[1]], tree2.LHEPart_phi[nu_coll_2[1]], tree2.LHEPart_mass[nu_coll_2[1]])
    
    print(tree2.LHEPart_pt[lep_coll_2[0]])
    print(tree2.LHEPart_pt[lep_coll_2[1]])
    print(tree2.LHEPart_pt[nu_coll_2[0]])
    print(tree2.LHEPart_pt[nu_coll_2[1]])

    print(tree2.LHEPart_phi[lep_coll_2[0]])
    print(tree2.LHEPart_phi[lep_coll_2[1]])
    print(tree2.LHEPart_pt[nu_coll_2[0]])
    print(tree2.LHEPart_pt[nu_coll_2[1]])

    print(tree2.LHEPart_phi[lep_coll_2[0]])
    print(tree2.LHEPart_phi[lep_coll_2[1]])
    print(tree2.LHEPart_phi[nu_coll_2[0]])
    print(tree2.LHEPart_phi[nu_coll_2[1]])

    print(tree2.LHEPart_mass[lep_coll_2[0]])
    print(tree2.LHEPart_mass[lep_coll_2[1]])
    print(tree2.LHEPart_mass[nu_coll_2[0]])
    print(tree2.LHEPart_mass[nu_coll_2[1]])


    # Calculate mWW for the 2l2nu system (ll + nn)
    lnu1 = lep1 + nu1
    lnu2 = lep2 + nu2
    llnunu = lnu1 + lnu2
    
    print("ll.pt", lnu1.Pt())
    print("nn.Pt()",lnu2.Pt())
    print("lnu1.Pt()", "lnu1.M())", lnu1.Pt(), lnu1.M())
    print("lnu2.Pt()", "lnu2.M())", lnu2.Pt(), lnu2.M())
    # Fill histogram with mWW (2l2nu system)
    #hist_mWW_2.Fill(llnunu.Pt())
    
    hist_mWW_2.Fill(llnunu.Pt())
    hist_lep_2.Fill(lep1.Pt())

    print("llnunu.Pt()",llnunu.Pt())



# Normalize histograms
hist_mWW_1.Scale(1.0 / hist_mWW_1.Integral())
hist_mWW_2.Scale(1.0 / hist_mWW_2.Integral())

hist_lep_1.Scale(1.0 / hist_mWW_1.Integral())
hist_lep_2.Scale(1.0 / hist_mWW_2.Integral())

# Set histogram styles
hist_mWW_1.SetLineColor(ROOT.kBlue)
hist_mWW_2.SetLineColor(ROOT.kRed)


hist_lep_1.SetLineColor(ROOT.kBlue)
hist_lep_2.SetLineColor(ROOT.kRed)


# Draw histograms on a canvas
canvas_mWW = ROOT.TCanvas("canvas_mWW_comparison", "mWW Comparison", 1000, 800)
canvas_mWW.SetLogy(1)  # Set logarithmic scale for y-axis
hist_mWW_1.Draw()
hist_mWW_1.GetXaxis().SetTitle("mWW [GeV]")
hist_mWW_1.GetYaxis().SetTitle("Normalized Events")
hist_mWW_2.Draw("same")

# Add legend without the outer box and with symbols
legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
legend.SetBorderSize(0)  # Remove outer box
legend.AddEntry(hist_mWW_1, "lnuqq Sample", "l")
legend.AddEntry(hist_mWW_2, "2l2nu Sample", "l")
legend.Draw()

canvas_mWW.SaveAs("mWW_comparison_lnuqq_2l2nu_lhe.pdf")

# Draw histograms on a canvas
canvas_lep_Pt = ROOT.TCanvas("canvas_mWW_comparison", "mWW Comparison", 1000, 800)
canvas_lep_Pt.SetLogy(1)  # Set logarithmic scale for y-axis
hist_lep_1.Draw()
hist_lep_1.GetXaxis().SetTitle("lep Pt[GeV]")
hist_lep_1.GetYaxis().SetTitle("Normalized Events")
hist_lep_2.Draw("same")

# Add legend without the outer box and with symbols
legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
legend.SetBorderSize(0)  # Remove outer box
legend.AddEntry(hist_lep_1, "lnuqq Sample", "l")
legend.AddEntry(hist_lep_2, "2l2nu Sample", "l")
legend.Draw()



# Save canvas as PDF
canvas_lep_Pt.SaveAs("lepPt_comparison_lnuqq_2l2nu_lhe.pdf")

# Close the ROOT files
file1.Close()
file2.Close()

# Print processing time
print("Processing time:", time.time() - start)
