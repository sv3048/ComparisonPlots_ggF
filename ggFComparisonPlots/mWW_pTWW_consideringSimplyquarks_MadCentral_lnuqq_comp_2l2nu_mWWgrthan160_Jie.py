#python3 mWW_pTWW_consideringSimplyquarks_MadCentral_lnuqq_comp_2l2nu_mWWgrthan160.py --rootfile1 --rootfile2 


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

# Creating 1D histograms to store the mWW distributions for both samples
hist_pTWW_1 = ROOT.TH1F("hist_pTWW_1", "Distribution of pTWW ", 100, 0, 1000)
hist_pTWW_2 = ROOT.TH1F("hist_pTWW_2", "Distribution of pTWW ", 100, 0, 1000)


hist_mWW_1 = ROOT.TH1F("hist_mWW_1", "Distribution of mWW", 100, 0, 1000)
hist_mWW_2 = ROOT.TH1F("hist_mWW_2", "Distribution of mWW", 100, 0, 1000)



# Creating 1D histograms to store the mWW distributions for both samples
hist_PTWW_1_new = ROOT.TH1F("hist_pTWW_1_new", "Distribution of pTWW_new - Sample 1", 100, 0, 1000)
hist_pTWW_2_new = ROOT.TH1F("hist_pTWW_2_new", "Distribution of pTWW_new - Sample 2", 100, 0, 1000)


hist_mWW_1_new = ROOT.TH1F("hist_mWW_1_new", "Distribution of mWW_new ", 100, 0, 1000)
hist_mWW_2_new = ROOT.TH1F("hist_mWW_2_new", "Distribution of mWW_new ", 100, 0, 1000)

print ("did we reach at line 52")
# Loop over events in the first TTree
num_events_tree1 = tree1.GetEntries()
for i in range(num_events_tree1):
    tree1.GetEntry(i)
    
    # Collect quarks
    quark_coll = []
    nu_coll = []
    lep_coll = []
    W_coll = []
   
    #lep, nu, jet1, jet2, qfat, met = [ROOT.TLorentzVector() for _ in range(6)]
    lep, nu, jet1, jet2, qfat, met, W1, W2 = [ROOT.TLorentzVector() for _ in range(8)]

    
    for part in range(tree1.nGenPart):
        #is_first_copy = (tree1.GenPart_statusFlags[part] >>12 & 1)
        is_first_copy = (tree1.GenPart_statusFlags[part] >>12 & 1)

        if is_first_copy and abs(tree1.GenPart_pdgId[tree1.GenPart_genPartIdxMother[part]]) == 24:
            W_coll.append(part)
            #print (len(W_coll))
            if abs(tree1.GenPart_pdgId[part]) in [11, 13, 15]:
                lep_coll.append(part)
            elif abs(tree1.GenPart_pdgId[part]) <= 4:
                quark_coll.append(part)
                #print(len(quark_coll))
            elif abs(tree1.GenPart_pdgId[part]) in [12, 14, 16]:
                nu_coll.append(part)
    
    # Require exactly 1 lepton and 1 neutrino
    if len(nu_coll) != 1 or len(lep_coll) != 1 or len(quark_coll) != 2 or len(W_coll) != 2 :

    #if len(nu_coll) != 1 or len(lep_coll) != 1 or len(quark_coll) != 2 or len(W_coll) != 2 :
        print(f"Neutrinos: {len(nu_coll)}, Leptons: {len(lep_coll)}, Quarks: {len(quark_coll)}, Ws: {len(W_coll)}")
        continue
    print ("did we reach here ")

    print ("did we reach here at line 85")
    # Calculate and fill histograms
    # Build the TLorentzVectors
    lep.SetPtEtaPhiM(tree1.GenPart_pt[lep_coll[0]], tree1.GenPart_eta[lep_coll[0]], tree1.GenPart_phi[lep_coll[0]], tree1.GenPart_mass[lep_coll[0]])
    nu.SetPtEtaPhiM(tree1.GenPart_pt[nu_coll[0]], tree1.GenPart_eta[nu_coll[0]], tree1.GenPart_phi[nu_coll[0]], tree1.GenPart_mass[nu_coll[0]])
    jet1.SetPtEtaPhiM(tree1.GenPart_pt[quark_coll[0]], tree1.GenPart_eta[quark_coll[0]], tree1.GenPart_phi[quark_coll[0]], tree1.GenPart_mass[quark_coll[0]])
    jet2.SetPtEtaPhiM(tree1.GenPart_pt[quark_coll[1]], tree1.GenPart_eta[quark_coll[1]], tree1.GenPart_phi[quark_coll[1]], tree1.GenPart_mass[quark_coll[1]])
    met.SetPtEtaPhiM(tree1.GenMET_pt, 0.0, tree1.GenMET_phi, 0.0)
    W1.SetPtEtaPhiM(tree1.GenPart_pt[W_coll[0]], tree1.GenPart_eta[W_coll[0]], tree1.GenPart_phi[W_coll[0]], tree1.GenPart_mass[W_coll[0]])
    W2.SetPtEtaPhiM(tree1.GenPart_pt[W_coll[1]], tree1.GenPart_eta[W_coll[1]], tree1.GenPart_phi[W_coll[1]], tree1.GenPart_mass[W_coll[1]])

    print (W1.M())
    # Calculate quantities
    lnu = lep + nu
    jj = jet1 + jet2
    lnujj = lnu + jj 
    WW = W1+W2

    if lnujj.M() > 160:
        hist_pTWW_1.Fill(lnujj.Pt())
        hist_mWW_1.Fill(lnujj.M())


    print ("hist_pTWW_1", lnujj.Pt())
    hist_mWW_1_new.Fill(WW.M())   
    hist_pTWW_1_new.Fill(WW.Pt())
    

# Loop over events in the second TTree
# Loop over events in the second TTree
num_events_tree2 = tree2.GetEntries()
for i in range(num_events_tree2):
    tree2.GetEntry(i)
    
    # Collect quarks
    quark_coll_2 = []
    nu_coll_2 = []
    lep_coll_2 = []
    W_coll_2 = []
    #lep1, nu1, lep2,nu2, qfat, met = [ROOT.TLorentzVector() for _ in range(6)]
    lep1, nu1, lep2,nu2, qfat, met, W1, W2 = [ROOT.TLorentzVector() for _ in range(8)]

    for part in range(tree2.nGenPart):
        is_first_copy = (tree2.GenPart_statusFlags[part]>>12 & 1)
        
        if is_first_copy and abs(tree2.GenPart_pdgId[tree2.GenPart_genPartIdxMother[part]]) == 24 :
            W_coll_2.append(part)
            if abs(tree2.GenPart_pdgId[part]) in [11, 13, 15]:
                lep_coll_2.append(part)
            #elif abs(tree2.GenPart_pdgId[part]) <= 4:
            #    quark_coll_2.append(part)
            elif abs(tree2.GenPart_pdgId[part]) in [12, 14, 16]:
                nu_coll_2.append(part)
    
    # Require exactly 1 lepton and 1 neutrino
    #if len(nu_coll_2) != 1 or len(lep_coll_2) != 1 or len(quark_coll_2) != 2:
    if(len(lep_coll_2)!=2 or len(nu_coll_2)!=2)  :

    #if(len(lep_coll_2)!=2 or len(nu_coll_2)!=2)  or len(W_coll_2) != 2:
        continue

    # Calculate and fill histograms
    # Build the TLorentzVectors
    lep1.SetPtEtaPhiM(tree2.GenPart_pt[lep_coll_2[0]], tree2.GenPart_eta[lep_coll_2[0]], tree2.GenPart_phi[lep_coll_2[0]], tree2.GenPart_mass[lep_coll_2[0]])
    nu1.SetPtEtaPhiM(tree2.GenPart_pt[nu_coll_2[0]], tree2.GenPart_eta[nu_coll_2[0]], tree2.GenPart_phi[nu_coll_2[0]], tree2.GenPart_mass[nu_coll_2[0]])
    #jet1.SetPtEtaPhiM(tree2.GenPart_pt[quark_coll_2[0]], tree2.GenPart_eta[quark_coll_2[0]], tree2.GenPart_phi[quark_coll_2[0]], tree2.GenPart_mass[quark_coll_2[0]])
    #jet2.SetPtEtaPhiM(tree2.GenPart_pt[quark_coll_2[1]], tree2.GenPart_eta[quark_coll_2[1]], tree2.GenPart_phi[quark_coll_2[1]], tree2.GenPart_mass[quark_coll_2[1]])
    lep2.SetPtEtaPhiM(tree2.GenPart_pt[lep_coll_2[1]], tree2.GenPart_eta[lep_coll_2[1]], tree2.GenPart_phi[lep_coll_2[1]], tree2.GenPart_mass[lep_coll_2[1]])
    nu2.SetPtEtaPhiM(tree2.GenPart_pt[nu_coll_2[1]], tree2.GenPart_eta[nu_coll_2[1]], tree2.GenPart_phi[nu_coll_2[1]], tree2.GenPart_mass[nu_coll_2[1]])
    met.SetPtEtaPhiM(tree2.GenMET_pt, 0.0, tree2.GenMET_phi, 0.0)
    W1.SetPtEtaPhiM(tree2.GenPart_pt[W_coll_2[0]], tree2.GenPart_eta[W_coll_2[0]], tree2.GenPart_phi[W_coll_2[0]], tree2.GenPart_mass[W_coll_2[0]])
    W2.SetPtEtaPhiM(tree2.GenPart_pt[W_coll_2[1]], tree2.GenPart_eta[W_coll_2[1]], tree2.GenPart_phi[W_coll_2[1]], tree2.GenPart_mass[W_coll_2[1]])

            
    ll = lep1 + lep2
    nn = nu1 + nu2
    ll_met = lep1 + lep2 + met
    llnunu = ll + nn 
    WW= W1+W2


    
    # # Calculate quantities
    # print("lep_coll_2[0] GenPart_pt:", tree2.GenPart_pt[lep_coll_2[0]])
    # print("lep_coll_2[1] GenPart_pt:", tree2.GenPart_pt[lep_coll_2[1]])
    # print("nu_coll_2[0] GenPart_pt:", tree2.GenPart_pt[nu_coll_2[0]])
    # print("nu_coll_2[1] GenPart_pt:", tree2.GenPart_pt[nu_coll_2[1]])

    # print("lep_coll_2[0] GenPart_phi:", tree2.GenPart_phi[lep_coll_2[0]])
    # print("lep_coll_2[1] GenPart_phi:", tree2.GenPart_phi[lep_coll_2[1]])
    # print("nu_coll_2[0] GenPart_pt:", tree2.GenPart_pt[nu_coll_2[0]])
    # print("nu_coll_2[1] GenPart_pt:", tree2.GenPart_pt[nu_coll_2[1]])

    # print("lep_coll_2[0] GenPart_phi:", tree2.GenPart_phi[lep_coll_2[0]])
    # print("lep_coll_2[1] GenPart_phi:", tree2.GenPart_phi[lep_coll_2[1]])
    # print("nu_coll_2[0] GenPart_phi:", tree2.GenPart_phi[nu_coll_2[0]])
    # print("nu_coll_2[1] GenPart_phi:", tree2.GenPart_phi[nu_coll_2[1]])

    # print("lep_coll_2[0] GenPart_mass:", tree2.GenPart_mass[lep_coll_2[0]])
    # print("lep_coll_2[1] GenPart_mass:", tree2.GenPart_mass[lep_coll_2[1]])
    # print("nu_coll_2[0] GenPart_mass:", tree2.GenPart_mass[nu_coll_2[0]])
    # print("nu_coll_2[1] GenPart_mass:", tree2.GenPart_mass[nu_coll_2[1]])


    
    if llnunu.M() > 160:
        hist_pTWW_2.Fill(llnunu.Pt())
        hist_mWW_2.Fill(llnunu.M())
    
    hist_pTWW_2_new.Fill(WW.Pt())
    hist_mWW_2_new.Fill(WW.M())



# Normalize histograms
hist_pTWW_1.Scale(1.0 / hist_pTWW_1.Integral())
hist_pTWW_2.Scale(1.0 / hist_pTWW_2.Integral())

hist_mWW_1.Scale(1.0 / hist_mWW_1.Integral())
hist_mWW_2.Scale(1.0 / hist_mWW_2.Integral())

# Normalize histograms
hist_pTWW_1_new.Scale(1.0 / hist_pTWW_1_new.Integral())
hist_pTWW_2_new.Scale(1.0 / hist_pTWW_2_new.Integral())

hist_mWW_1_new.Scale(1.0 / hist_mWW_1_new.Integral())
hist_mWW_2_new.Scale(1.0 / hist_mWW_2_new.Integral())

# Set histogram styles
hist_pTWW_1.SetLineColor(ROOT.kBlue)
hist_pTWW_2.SetLineColor(ROOT.kRed)

hist_mWW_1.SetLineColor(ROOT.kBlue)
hist_mWW_2.SetLineColor(ROOT.kRed)

# Set histogram styles
hist_pTWW_1_new.SetLineColor(ROOT.kBlue)
hist_pTWW_2_new.SetLineColor(ROOT.kRed)

hist_mWW_1_new.SetLineColor(ROOT.kBlue)
hist_mWW_2_new.SetLineColor(ROOT.kRed)


# Draw histograms on a canvas
canvas_pTWW = ROOT.TCanvas("canvas_PTWW", "pTWW Comparison", 1000, 800)
canvas_pTWW.SetLogy(1)  # Set logarithmic scale for y-axis
hist_pTWW_1.Draw()
hist_pTWW_1.GetXaxis().SetTitle("pTWW [GeV]")
hist_pTWW_1.GetYaxis().SetTitle("Normalized Events")
hist_pTWW_2.Draw("same")




# Add legend without the outer box and with symbols
legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
legend.SetBorderSize(0)  # Remove outer box
legend.AddEntry(hist_pTWW_1, "Madgraph(Central)", "l")
#legend.AddEntry(hist_mWW_2, "Madgraph+Pythia(WW decay) ", "l")
legend.AddEntry(hist_pTWW_2, "MCFM 2l2nu Sample", "l")
#legend.SetFillStyle(0) # Remove legend fill
legend.Draw()



#canvas.SaveAs("mWW_comparison_lnuqq_madCental_2l2nu_mcfm_trywithFIRSTCOPY.pdf")
canvas_pTWW.SaveAs("PTWW_comparison_lnuqq_madCental_2l2nu_mcfm_trywithFIRSTCOPY_mWWGRT160_NEW.pdf")


# Draw histograms on a canvas
canvas_mWW = ROOT.TCanvas("canvas_mWW_comparison", "mWW Comparison", 1000, 800)
canvas_mWW.SetLogy(1)  # Set logarithmic scale for y-axis
hist_mWW_1.Draw()
hist_mWW_1.GetXaxis().SetTitle("mWW [GeV]")
hist_mWW_1.GetYaxis().SetTitle("Normalized Events")
hist_mWW_2.Draw("same")


# Add legend
legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
legend.AddEntry(hist_mWW_1, "Madgraph (Central)", "l")
legend.AddEntry(hist_mWW_2, " MCFM 2l2nu ", "l")
legend.Draw()

canvas_mWW.SaveAs("mWW_comparison_lnuqq_madCental_2l2nu_mcfm_trywithFIRSTCOPY_mWWGRT160_NEW.pdf")

############### Canvas 3 ###################

# Draw histograms on a canvas
canvas_pTWW_new = ROOT.TCanvas("canvas_mWW_comparison", "mWW Comparison", 1000, 800)
canvas_pTWW_new.SetLogy(1)  # Set logarithmic scale for y-axis
hist_pTWW_.Draw()
hist_pTWW_1_new.GetXaxis().SetTitle("mWW [GeV]")
hist_pTWW_1_new.GetYaxis().SetTitle("Normalized Events")
hist_pTWW_2_new.Draw("same")


# Add legend
legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
legend.AddEntry(hist_pTWW_1_new, "Madgraph (Central)", "l")
legend.AddEntry(hist_pTWW_2_new, " MCFM 2l2nu ", "l")
legend.Draw()

canvas_pTWW_new.SaveAs("mWW_comparison_lnuqq_madCentral_2l2nu_mcfm_tryBeforeRadiationOfW_mWWGRT160_NEW.pdf")

############### Canvas 4 ###################

#canvas.SaveAs("mWW_comparison_lnuqq_madCental_2l2nu_mcfm_trywithFIRSTCOPY.pdf")


# Draw histograms on a canvas
canvas_mWW_new = ROOT.TCanvas("canvas_mWW_comparison", "mWW Comparison", 1000, 800)
canvas_mWW_new.SetLogy(1)  # Set logarithmic scale for y-axis
hist_mWW_1_new.Draw()
hist_mWW_1_new.GetXaxis().SetTitle("mWW [GeV]")
hist_mWW_1_new.GetYaxis().SetTitle("Normalized Events")
hist_mWW_2_new.Draw("same")


# Add legend
legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
legend.AddEntry(hist_mWW_1_new, "Madgraph (Central)", "l")
legend.AddEntry(hist_mWW_2_new, " MCFM 2l2nu ", "l")
legend.Draw()

canvas_mWW_new.SaveAs("mWW_comparison_lnuqq_madCentral_2l2nu_mcfm_tryBeforeRadiationOfW_mWWGRT160_NEW.pdf")



# Close the ROOT files
file1.Close()
file2.Close()

# Print processing time
print("Processing time:", time.time() - start)
