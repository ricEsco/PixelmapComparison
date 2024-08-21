import ROOT
import math

chip = "12"
histdirectory = "Detector/Board_0/OpticalGroup_0/Hybrid_0/Chip_"
histnameprefix = "/D_B(0)_O(0)_H(0)_"

forwardroot = ROOT.TFile("/home/kalib/Analysis/inputroot/Run000005_SCurve.root", "READ")
reverseroot = ROOT.TFile("/home/kalib/Analysis/inputroot/Run000004_SCurve.root", "READ")

forward_thr_canvas = forwardroot.Get(histdirectory+chip+histnameprefix+"Threshold2D_Chip("+chip+")")
h_forward_thr = None
primitives = forward_thr_canvas.GetListOfPrimitives()
for primitive in primitives:
    if isinstance(primitive, ROOT.TH2):
        h_forward_thr = primitive
        break

reverse_thr_canvas = reverseroot.Get(histdirectory+chip+histnameprefix+"Threshold2D_Chip("+chip+")")
h_reverse_thr = None
primitives = reverse_thr_canvas.GetListOfPrimitives()
for primitive in primitives:
    if isinstance(primitive, ROOT.TH2):
        h_reverse_thr = primitive
        break

forward_ns_canvas = forwardroot.Get(histdirectory+chip+histnameprefix+"Noise2D_Chip("+chip+")")
h_forward_ns = None
primitives = forward_ns_canvas.GetListOfPrimitives()
for primitive in primitives:
    if isinstance(primitive, ROOT.TH2):
        h_forward_ns = primitive
        break

reverse_ns_canvas = reverseroot.Get(histdirectory+chip+histnameprefix+"Noise2D_Chip("+chip+")")
h_reverse_ns = None
primitives = reverse_ns_canvas.GetListOfPrimitives()
for primitive in primitives:
    if isinstance(primitive, ROOT.TH2):
        h_reverse_ns = primitive
        break

delta_thr = h_forward_thr.Clone("delta_thr")
delta_ns = h_forward_ns.Clone("delta_ns")

delta_thr.Add(h_reverse_thr, -1)
delta_ns.Add(h_reverse_ns, -1)

# Creating missing map
missing_map = delta_thr.Clone("missing_map")
for i in range(1, missing_map.GetNbinsX() + 1):
    for j in range(1, missing_map.GetNbinsY() + 1):
        missing_map.SetBinContent(i, j, 0)

for i in range(1, missing_map.GetNbinsX() + 1):
    for j in range(1, missing_map.GetNbinsY() + 1):
        if (math.sqrt(((delta_thr.GetBinContent(i, j))**2) + ((delta_ns.GetBinContent(i, j))**2)) <= 5.0):
            missing_map.SetBinContent(i, j, 1)

c = ROOT.TCanvas("c", "Canvas", missing_map.GetNbinsX(), missing_map.GetNbinsY())
c.Clear()
missing_map.Draw("missing_map")
c.SaveAs("results/thr5missing_map.png")

output_file = ROOT.TFile("results/frbias/histograms.root", "RECREATE")
delta_thr.Write()
delta_ns.Write()
missing_map.Write()
output_file.Close()
forwardroot.Close()
reverseroot.Close()
