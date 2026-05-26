import ROOT
import math
import os
import glob

plot_dir = "/home/aidangardnerokearny/mucoll/analysisCode/efficiencyProblem/efficiencyPlots"
os.makedirs(plot_dir, exist_ok=True)

efficiencyHisto = ROOT.TH1D("efficiency", "", 100, 0, 100)


def getTLV(obj):
    obj_p = obj.getMomentum()
    obj_e = obj.getEnergy()
    obj_tlv = ROOT.TLorentzVector()
    obj_tlv.SetPxPyPzE(obj_p[0], obj_p[1], obj_p[2], obj_e)
    return obj_tlv


reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
reader.open("/home/aidangardnerokearny/mucoll/reco/electronGun_pT_0_50/electronGun_pT_0_50_reco_0.slcio")
for event in reader:
    mcps = event.getCollection("MCParticle")
    pfos = event.getCollection("PandoraPFOs")
    trks = event.getCollection("SiTracks_Refitted")
    clusters = event.getCollection("PandoraClusters")

    mcp_electrons = []
    for mcp in mcps:
        if mcp.getGeneratorStatus() != 1 or abs(mcp.getPDG()) != 11: continue
        tlv = getTLV(mcp)
        print(tlv.E())
 
