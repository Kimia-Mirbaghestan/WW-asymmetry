import ROOT
import ctypes
import copy
import math
from Samples import samp

class MyAnalysis(object):
   
    def __init__(self, sample):

        """ The Init() function is called when an object MyAnalysis is initialised
        The tree corresponding to the specific sample is picked up
        and histograms are booked.
        """
        self._tree = ROOT.TTree()
        self.A = 0
        self.AA = 0
        self.un_AA = 0
        self.charge_neg = 0
        self.charge_pos = 0
        self.Triggered_pos = 0  
        self.Triggered_neg = 0 
        self.bac_neg = 0
        self.bac_pos = 0
        
        if(sample not in samp.keys() and sample != "data"):
            print ("Error") #RuntimeError("Sample %s not valid. please, choose among these: %s" % (sample, str(samp.keys())) )
            exit()
        self.histograms = {}
        self.sample = sample
        self._file = ROOT.TFile("files/"+sample+".root")
        self._file.cd()
        tree = self._file.Get("events")
        self._tree = tree
        self.acceptance_pos = 0
        self.acceptance_neg = 0
        self.integral_aftercut_pos = 0
        self.integral_aftercut_neg = 0
        self.integral_beforecut_pos = 0
        self.integral_beforecut_neg = 0
        self.triggered_pass_neg = 0
        self.triggered_pass_pos = 0
        self.acceptance_uncertainty_neg = 0
        self.acceptance_uncertainty_pos = 0
        self.Trigger_Efficiency_neg = 1.0
        self.Trigger_Efficiency_pos = 1.0
        self.Trigger_Efficiency_uncertainty_neg = 0
        self.Trigger_Efficiency_uncertainty_pos = 0
        self.nEvents = self._tree.GetEntries()
        print ("-------------------------------------------------------")
        print ("Number of entries for " + self.sample + ": " + str(self.nEvents))
        ### Book histograms
        self.bookHistos()

    def getTree(self):
        return self._tree

    def getHistos(self):
        return self.histograms

    def bookHistos(self):
        h_nJet = ROOT.TH1F("NJet","#of jets", 6, -0.5, 6.5)
        h_nJet.SetXTitle("%# of jets")
        self.histograms["NJet"] = h_nJet

        h_nJetFinal = ROOT.TH1F("NJetFinal","#of jets", 6, -0.5, 6.5) #Plotted
        h_nJetFinal.SetXTitle("%# of jets")
        self.histograms["NJetFinal"] = h_nJetFinal

        h_MuonIso = ROOT.TH1F("Muon_Iso","Muon Isolation", 25, 0., 3.)
        h_MuonIso.SetXTitle("Muon Isolation")
        self.histograms["Muon_Iso"] = h_MuonIso

        h_NIsoMu = ROOT.TH1F("NIsoMu","Number of isolated muons", 5, 0.5, 5.5) #Plotted
        h_NIsoMu.SetXTitle("Number of isolated muons")
        self.histograms["NIsoMu"] = h_NIsoMu

        h_MuonPt = ROOT.TH1F("Muon_Pt","Muon P_T", 50, 0., 200.) #Plotted
        h_MuonPt.SetXTitle("Muon P_T")
        self.histograms["Muon_Pt"] = h_MuonPt

        h_METpt = ROOT.TH1F("MET_Pt","MET P_T", 25, 0., 300.) #Plotted
        h_METpt.SetXTitle("MET P_T")
        self.histograms["MET_Pt"] = h_METpt
        
        h_METp = ROOT.TH1F("MET_p","MET_p", 25, 0., 500.) #Plotted
        h_METp.SetXTitle("MET_p")
        self.histograms["MET_p"] = h_METp
        
        h_METn = ROOT.TH1F("MET_n","MET_n", 25, 0., 500.) #Plotted
        h_METn.SetXTitle("MET_n")
        self.histograms["MET_n"] = h_METn

        h_JetPt = ROOT.TH1F("Jet_Pt","Jet P_T", 50, 0., 200.) #Plotted
        h_JetPt.SetXTitle("Jet P_T")
        self.histograms["Jet_Pt"] = h_JetPt

        h_JetBtag = ROOT.TH1F("Jet_Btag","Jet B tag", 10, 1., 6.)
        h_JetBtag.SetXTitle("Jet B tag")
        self.histograms["Jet_btag"] = h_JetBtag

        h_NBtag = ROOT.TH1F("NBtag","Jet B tag", 4, 0.5, 4.5) #Plotted
        h_NBtag.SetXTitle("Number of B tagged jets")
        self.histograms["NBtag"] = h_NBtag
        
        
        h_nJet_nocut = ROOT.TH1F("NJet_NoCut","# of Jet", 12, -0.5, 12.5) #To Check Integral without Cut    
        h_nJet_nocut.SetXTitle("Number of Jets without any Cut")
        self.histograms["NJet_NoCut"] =  h_nJet_nocut
        
        h_nJet_triggered_pos = ROOT.TH1F("Triggered_pos","# of Jet", 100, -0.5, 199.5) #To Check Integral without Cut
        h_nJet_triggered_pos.SetXTitle("Number of Jets with Trigger")
        self.histograms["Triggered_pos"] =  h_nJet_triggered_pos
        
        h_nJet_triggered_neg = ROOT.TH1F("Triggered_neg","# of Jet", 100, -0.5, 199.5) #To Check Integral without Cut
        h_nJet_triggered_neg.SetXTitle("Number of Jets with Trigger")
        self.histograms["Triggered_neg"] =  h_nJet_triggered_neg
        
        self.histograms["Mass_W_Positive"] = ROOT.TH1F("Mass_W_Positive", "Mass of W (positive)", 100,  0., 150.)
        self.histograms["Mass_W_Positive"].SetXTitle("Mass of W (positive)")
        
        self.histograms["Mass_W_Negative"] = ROOT.TH1F("Mass_W_Negative", "Mass of W (negative)", 100,  0., 150.)
        self.histograms["Mass_W_Negative"].SetXTitle("Mass of W (negative)")
        
        self.histograms["Mass_W2"] = ROOT.TH1F("Mass_W2", "Mass of W (positive)", 100,  0., 150.)
        self.histograms["Mass_W2"].SetXTitle("Mass of W (2)")
        
        h_Muon_charge_positive = ROOT.TH1F("Muon_charge_positive","# of Jet", 12, -0.5, 12.5) #To Check Integral without Cut
        h_Muon_charge_positive.SetXTitle("charge of the W positive")
        self.histograms["Muon_charge_positive"] =  h_Muon_charge_positive
        
        h_Muon_charge_negative = ROOT.TH1F("Muon_charge_negative","# of Jet", 12, -0.5, 12.5) #To Check Integral without Cut
        h_Muon_charge_negative.SetXTitle("charge of the W negative")
        self.histograms["Muon_charge_negative"] =  h_Muon_charge_negative
        
        h_Muon_charge = ROOT.TH1F("Muon_charge","# of Jet", 12, -0.5, 12.5) #To Check Integral without Cut
        h_Muon_charge.SetXTitle("charge")
        self.histograms["Muon_charge"] =  h_Muon_charge
        
        

    def saveHistos(self):
        outfilename = "Histogram_Root/"+self.sample + "_histos.root"
        outfile = ROOT.TFile(outfilename, "RECREATE")
        outfile.cd()
        for h in self.histograms.values():
            h.Write()
        outfile.Close()

    ### processEvent function implements the actions to perform on each event
    ### This is the place where to implement the analysis strategy: study of most sensitive variables
    ### and signal-like event selection
    
    def processEvent2(self, entry):
        tree = self.getTree()
        tree.GetEntry(entry)
        w = tree.EventWeight
        metPx2=tree.MET_px
        metPy2=tree.MET_py
        metPt2=math.sqrt(metPx2**2 + metPy2**2)
        for m in range(tree.NMuon):
            if tree.Muon_Charge[m]:
              MET2 = ROOT.TLorentzVector(tree.MET_px, tree.MET_py, 0., metPt2)
              muon_T2 = ROOT.TLorentzVector(tree.Muon_Px[m], tree.Muon_Py[m], 0., math.sqrt(tree.Muon_E[m]**2-tree.Muon_Pz[m]**2)) 
              Mass_W2 = (MET2 + muon_T2).M()
              self.histograms["Mass_W2"].Fill(Mass_W2, w)
        
        

    def processEvent(self, entry):
        tree = self.getTree()
        tree.GetEntry(entry)
        w = tree.EventWeight

        ### Muon selection - Select events with at least 1 isolated muon
        ### with pt>25 GeV to match trigger requirements
        muonPtCut = 30    #25.
        muonRelIsoCut = 0.05
        nIsoMu = 0
        muon_eta_cut = 0.4
        muon_eta_cut_min = 1.5  # -0.4
        muon_eta_cut_max = 1.8
        metPx=tree.MET_px
        metPy=tree.MET_py
        metPt=math.sqrt(metPx**2 + metPy**2)
        # Jet and b-jet counts
        nJets = 0
        nBjets = 0
                    
        # Object selection: selecting the muon
        for m in range(tree.NMuon):
            muon = ROOT.TLorentzVector(tree.Muon_Px[m],tree.Muon_Py[m],tree.Muon_Pz[m],tree.Muon_E[m])
            self.histograms["Muon_Iso"].Fill(tree.Muon_Iso[m], w)
            if muon.Pt() > muonPtCut and (tree.Muon_Iso[m]/muon.Pt()) < muonRelIsoCut:
                if (muon_eta_cut_min < abs(muon.Eta()) < muon_eta_cut_max):
                #if (abs(muon.Eta()) < muon_eta_cut):
                   nIsoMu += 1

        
        # Event Selection (Apply cuts)
        if nIsoMu >= 1  and tree.triggerIsoMu24 and metPt > 30:
            for m in range(tree.NMuon):
                muon = ROOT.TLorentzVector(tree.Muon_Px[m],tree.Muon_Py[m],tree.Muon_Pz[m],tree.Muon_E[m])
                if muon.Pt() > muonPtCut and (tree.Muon_Iso[m]/muon.Pt()) < muonRelIsoCut and muon_eta_cut_min < abs(muon.Eta()) < muon_eta_cut_max:
                    MET = ROOT.TLorentzVector(tree.MET_px, tree.MET_py, 0., metPt)
                    muon_T = ROOT.TLorentzVector(tree.Muon_Px[m], tree.Muon_Py[m], 0., math.sqrt(tree.Muon_E[m]**2-tree.Muon_Pz[m]**2))
                    Mass_W = (MET + muon_T).M()
                    if tree.Muon_Charge[m] > 0:
                        self.histograms["Mass_W_Positive"].Fill(Mass_W, w)
                        self.histograms["MET_p"].Fill(metPt, w)
                        self.charge_pos +=1
                        if tree.triggerIsoMu24 :
                            self.Triggered_pos += 1
                            self.histograms["Triggered_pos"].Fill(Mass_W,w)
                    if tree.Muon_Charge[m] < 0:
                        self.histograms["Mass_W_Negative"].Fill(Mass_W, w)
                        self.histograms["MET_n"].Fill(metPt, w)
                        self.charge_neg +=1
                        if tree.triggerIsoMu24 :
                            self.Triggered_neg += 1
                            self.histograms["Triggered_neg"].Fill(Mass_W,w)

        
    # computing acceptance and triggered efficiency  
    def processEvents(self):
        #W_negative
        minBin_beforecut_neg = self.histograms["Mass_W2"].GetXaxis().GetFirst()
        maxBin_beforecut_neg = self.histograms["Mass_W2"].GetXaxis().GetLast()
        nDataErr_beforecut_neg = ctypes.c_double(0.0)
        minBin_aftercut_neg = self.histograms["MET_n"].GetXaxis().GetFirst()
        maxBin_aftercut_neg = self.histograms["MET_n"].GetXaxis().GetLast()
        nDataErr_aftercut_neg = ctypes.c_double(0.0)
        minBin_triggered_neg = self.histograms["Triggered_neg"].GetXaxis().GetFirst()
        maxBin_triggered_neg = self.histograms["Triggered_neg"].GetXaxis().GetLast()
        nDataErr_triggered_neg = ctypes.c_double(0.0)
        
        #W_positive
        minBin_beforecut_pos = self.histograms["NJet_NoCut"].GetXaxis().GetFirst()
        maxBin_beforecut_pos = self.histograms["NJet_NoCut"].GetXaxis().GetLast()
        nDataErr_beforecut_pos = ctypes.c_double(0.0)
        minBin_aftercut_pos = self.histograms["MET_p"].GetXaxis().GetFirst()
        maxBin_aftercut_pos = self.histograms["MET_p"].GetXaxis().GetLast()
        nDataErr_aftercut_pos = ctypes.c_double(0.0)
        minBin_triggered_pos = self.histograms["Triggered_pos"].GetXaxis().GetFirst()
        maxBin_triggered_pos = self.histograms["Triggered_pos"].GetXaxis().GetLast()
        nDataErr_triggered_pos = ctypes.c_double(0.0)
        
        nevts = self.nEvents
        self.charge_pos = 0
        self.charge_neg = 0
        for i in range(nevts):
            self.processEvent(i)
            if self.sample == "wjets":
                self.processEvent2(i)
                
        print ("Events Passed for negative: " + str(self.charge_neg))
        print ("Events Passed for positive: " + str(self.charge_pos))
                
        #W_negative
        print ("For negative W:  ")
        self.integral_beforecut_neg = self.histograms["Mass_W2"].IntegralAndError(minBin_beforecut_neg,maxBin_beforecut_neg,nDataErr_beforecut_neg)
        self.integral_aftercut_neg = self.histograms["MET_n"].IntegralAndError(minBin_aftercut_neg,maxBin_aftercut_neg,nDataErr_aftercut_neg)
        #self.integral_aftercut_neg = self.histograms["Mass_W_Negative"].IntegralAndError(minBin_aftercut_neg,maxBin_aftercut_neg,nDataErr_aftercut_neg)
        print ("Integral without Cut negative: "+str(self.integral_beforecut_neg)+" Integral uncertainty without Cut : "+str(nDataErr_beforecut_neg.value))
        print ("Integral with Cut negative: "+str(self.integral_aftercut_neg)+" Integral uncertainty with Cut : "+str(nDataErr_aftercut_neg.value))
        #W_positive
        print ("For positive W:  ")
        self.integral_beforecut_pos = self.histograms["Mass_W2"].IntegralAndError(minBin_beforecut_pos,maxBin_beforecut_pos,nDataErr_beforecut_pos)
        self.integral_aftercut_pos = self.histograms["MET_p"].IntegralAndError(minBin_aftercut_pos,maxBin_aftercut_pos,nDataErr_aftercut_pos)
        #self.integral_aftercut_pos = self.histograms["Mass_W_Positive"].IntegralAndError(minBin_aftercut_pos,maxBin_aftercut_pos,nDataErr_aftercut_pos)
        print ("Integral without Cut : "+str(self.integral_beforecut_pos)+" Integral uncertainty without Cut : "+str(nDataErr_beforecut_pos.value))
        print ("Integral with Cut : "+str(self.integral_aftercut_pos)+" Integral uncertainty with Cut : "+str(nDataErr_aftercut_pos.value))
  
        
        
        #if self.sample == "wjets":
        if False:
            #W_negative
            print ("For negative W:  ")
            self.integral_beforecut_neg = self.histograms["Mass_W2"].IntegralAndError(minBin_beforecut_neg,maxBin_beforecut_neg,nDataErr_beforecut_neg)
            self.integral_aftercut_neg = self.histograms["Mass_W_Negative"].IntegralAndError(minBin_aftercut_neg,maxBin_aftercut_neg,nDataErr_aftercut_neg)
            self.acceptance_neg = (self.integral_aftercut_neg/self.integral_beforecut_neg) 
            self.acceptance_uncertainty_neg = (nDataErr_aftercut_neg.value/self.integral_beforecut_neg)**2 +((nDataErr_beforecut_neg.value*self.integral_aftercut_neg)/(self.integral_beforecut_neg**2))**2
            print ("Integral without Cut : "+str(self.integral_beforecut_neg)+" Integral uncertainty without Cut : "+str(nDataErr_beforecut_neg.value))
            print ("Integral with Cut : "+str(self.integral_aftercut_neg)+" Integral uncertainty with Cut : "+str(nDataErr_aftercut_neg.value))
            print ("Acceptance : "+str(self.acceptance_neg)+" acceptance uncertanty : " +str(self.acceptance_uncertainty_neg))
            
            #W_positive
            print ("For positive W:  ")
            self.integral_beforecut_pos = self.histograms["Mass_W2"].IntegralAndError(minBin_beforecut_pos,maxBin_beforecut_pos,nDataErr_beforecut_pos)
            self.integral_aftercut_pos = self.histograms["Mass_W_Positive"].IntegralAndError(minBin_aftercut_pos,maxBin_aftercut_pos,nDataErr_aftercut_pos)
            self.acceptance_pos = (self.integral_aftercut_pos/self.integral_beforecut_neg) 
            self.acceptance_uncertainty_pos = (nDataErr_aftercut_pos.value/self.integral_beforecut_pos)**2 +((nDataErr_beforecut_pos.value*self.integral_aftercut_pos)/(self.integral_beforecut_pos**2))**2
            print ("Integral without Cut : "+str(self.integral_beforecut_pos)+" Integral uncertainty without Cut : "+str(nDataErr_beforecut_pos.value))
            print ("Integral with Cut : "+str(self.integral_aftercut_pos)+" Integral uncertainty with Cut : "+str(nDataErr_aftercut_pos.value))
            print ("Acceptance : "+str(self.acceptance_pos)+" acceptance uncertanty : " +str(self.acceptance_uncertainty_pos))
            
        #if self.sample == "wjets":
        if False:
            #W_negative
            print ("For negative W:  ")
            self.triggered_pass_neg = self.histograms["Triggered_neg"].IntegralAndError(minBin_triggered_neg,maxBin_triggered_neg,nDataErr_triggered_neg)
            print ("Events Passed with Trigger : " + str(self.triggered_pass_neg))
            self.Trigger_Efficiency_neg = (self.triggered_pass_neg/self.integral_aftercut_neg) 
            self.Trigger_Efficiency_uncertainty_neg = (nDataErr_triggered_neg.value/self.integral_aftercut_neg)**2 +((nDataErr_triggered_neg.value*self.triggered_pass_neg)/(self.integral_aftercut_neg**2))**2
            print ("Trigger Efficiency : " + str(self.Trigger_Efficiency_neg)+" Trigger Efficiency uncertainty : " + str(self.Trigger_Efficiency_uncertainty_neg))
            
            #W_positive
            print ("For positive W:  ")
            self.triggered_pass_pos = self.histograms["Triggered_pos"].IntegralAndError(minBin_triggered_pos,maxBin_triggered_pos,nDataErr_triggered_pos)
            print ("Events Passed with Trigger : " + str(self.triggered_pass_pos))
            self.Trigger_Efficiency_pos = (self.triggered_pass_pos/self.integral_aftercut_pos) 
            self.Trigger_Efficiency_uncertainty_pos = (nDataErr_triggered_pos.value/self.integral_aftercut_pos)**2 +((nDataErr_triggered_pos.value*self.triggered_pass_pos)/(self.integral_aftercut_pos**2))**2
            print ("Trigger Efficiency : " + str(self.Trigger_Efficiency_pos)+" Trigger Efficiency uncertainty : " + str(self.Trigger_Efficiency_uncertainty_pos))
            self.Asymmetry()
            
        if self.sample == "data":
            self.A = (self.charge_pos - self.charge_neg)/(self.charge_pos + self.charge_neg)
            self.AA = (self.integral_aftercut_pos - self.integral_aftercut_neg)/(self.integral_aftercut_pos + self.integral_aftercut_neg)
            self.un_AA = math.sqrt(((nDataErr_aftercut_pos.value)/(self.integral_aftercut_pos + self.integral_aftercut_neg)+((self.integral_aftercut_pos - self.integral_aftercut_neg)*nDataErr_aftercut_pos.value)/(self.integral_aftercut_pos + self.integral_aftercut_neg)**2)**2+((-nDataErr_aftercut_neg.value)/(self.integral_aftercut_pos + self.integral_aftercut_neg)+((self.integral_aftercut_pos - self.integral_aftercut_neg)*nDataErr_aftercut_neg.value)/(self.integral_aftercut_pos + self.integral_aftercut_neg)**2)**2)
            print("I think asymmetry1 : ", str(self.A))
            print("I think asymmetry2 : ", str(self.AA))
            print("I think uncertainty2 : ", str(self.un_AA))
            
            
        print ("Cut Efficiency for negative: " + str((self.charge_neg/nevts)*100))
        print ("Cut Efficiency for positive: " + str((self.charge_pos/nevts)*100))
        self.saveHistos()
               
    def getNTriggers(self):
        return self.n_triggers

    
    
    def Asymmetry(self):
        # Assuming the number of observed events is stored in the histogram "NJetFinal"
        N_observed_pos = self.charge_neg     #14034   #397          
        N_observed_neg = self.charge_neg     #11331   #664         

        # Integrated luminosity (in pb^-1)
        luminosity = 50  # Example value, replace with the actual integrated luminosity
        print ("Acceptance Neg "+str(self.acceptance_neg))
        print ("Acceptance Pos "+str(self.acceptance_pos))
        print ("Trigger Neg "+str(self.Trigger_Efficiency_neg))
        print ("Trigger Pos "+str(self.Trigger_Efficiency_pos))

        # Compute cross-section
        sigma_observed_pos = N_observed_pos / (luminosity * self.acceptance_pos * self.Trigger_Efficiency_pos)
        sigma_observed_neg = N_observed_neg / (luminosity * self.acceptance_neg * self.Trigger_Efficiency_neg)

        un_N_observed_pos = math.sqrt(N_observed_pos)
        un_N_observed_neg = math.sqrt(N_observed_neg)
        un_luminosity = 0.1
        un_Trigger_Efficiency_pos = self.Trigger_Efficiency_uncertainty_pos
        un_Trigger_Efficiency_neg = self.Trigger_Efficiency_uncertainty_neg
        un_acceptance_pos = self.acceptance_uncertainty_pos
        un_acceptance_neg = self.acceptance_uncertainty_neg

        total_uncertainty_pos = math.sqrt((un_N_observed_pos / luminosity * self.acceptance_pos * self.Trigger_Efficiency_pos) ** 2 + \
                                ((un_luminosity * N_observed_pos) / ((luminosity ** 2) * self.acceptance_pos *
                                                               self.Trigger_Efficiency_pos)) ** 2 + \
                                ((un_Trigger_Efficiency_pos * N_observed_pos) / (
                                        (self.acceptance_pos ** 2) * self.Trigger_Efficiency_pos *
                                    luminosity * luminosity)) ** 2 + \
                                ((un_acceptance_pos * N_observed_pos) / (
                                        (self.Trigger_Efficiency_pos) * self.acceptance_pos * luminosity)) ** 2)

        total_uncertainty_neg = math.sqrt((un_N_observed_neg / luminosity * self.acceptance_neg * self.Trigger_Efficiency_neg) ** 2 + \
                                ((un_luminosity * N_observed_neg) / ((luminosity ** 2) * self.acceptance_neg *
                                                               self.Trigger_Efficiency_neg)) ** 2 + \
                                ((un_Trigger_Efficiency_neg * N_observed_neg) / (
                                        (self.acceptance_neg ** 2) * self.Trigger_Efficiency_neg *
                                        luminosity * luminosity)) ** 2 + \
                                ((un_acceptance_neg * N_observed_neg) / (
                                        (self.Trigger_Efficiency_neg) * self.acceptance_neg * luminosity)) ** 2)
        Asymmetry = (sigma_observed_pos - sigma_observed_neg)/(sigma_observed_pos + sigma_observed_pos)
        
        # Calculate the partial derivatives
        partial_sigma_pos = 2 * sigma_observed_neg / (sigma_observed_pos + sigma_observed_neg)**2
        partial_sigma_neg = -2 * sigma_observed_pos / (sigma_observed_pos + sigma_observed_neg)**2
        # Calculate the uncertainty in the asymmetry using error propagation
        delta_asymmetry = math.sqrt((partial_sigma_pos * total_uncertainty_pos)**2 + (partial_sigma_neg * total_uncertainty_neg)**2)
        print("Observed cross-section for positive W:", sigma_observed_pos, "pb")  
        print("Total uncertainty for positive W:", total_uncertainty_pos, "pb")

        print("Observed cross-section for negative W:", sigma_observed_neg, "pb")
        print("Total uncertainty for negative W:", total_uncertainty_neg, "pb")
        
        print("Asymmetry: " + str(Asymmetry) + " Asymmetry uncertainty: " + str(delta_asymmetry))



