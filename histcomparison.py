import ROOT
from array import array

def compare_xtalk_xray(xray_file_path, xtalk_file_path, frbias_file_path, output_file_path):
    # Load the ROOT files
    xray_file = ROOT.TFile(xray_file_path, "READ")
    xtalk_file = ROOT.TFile(xtalk_file_path, "READ")
    frbias_file = ROOT.TFile(frbias_file_path, "READ")
    
    h_xray = xray_file.Get("MissingMap")
    h_xtalk = xtalk_file.Get("h_confirmed2D")
    h_frbias = frbias_file.Get("missing_map")

    # Clone the X-ray histogram to create a result histogram
    h_xray_exclusive = h_xray.Clone("h_xray_exclusive")
    h_xtalk_exclusive = h_xray.Clone("h_xtalk_exclusive")
    h_frbias_exclusive = h_xray.Clone("h_frbias_exclusive")
    h_xray_xtalk = h_xray.Clone("h_xray_xtalk")
    h_xray_frbias = h_xray.Clone("h_xray_frbias")
    h_xtalk_frbias = h_xray.Clone("h_xtalk_frbias")
    h_xray_xtalk_frbias = h_xray.Clone("h_xray_xtalk_frbias")

    histograms = [h_xray_exclusive,
                  h_xtalk_exclusive,
                  h_frbias_exclusive,
                  h_xray_xtalk,
                  h_xray_frbias,
                  h_xtalk_frbias,
                  h_xray_xtalk_frbias
                  ]

    # Filling histograms with 0's
    for histogram in histograms:
        for i in range(1, histogram.GetNbinsX()+1):
            for j in range(1, histogram.GetNbinsY()+1):
                histogram.SetBinContent(i, j, 0)

    for i in range(1, h_xray_exclusive.GetNbinsX() + 1):
        for j in range(1, h_xray_exclusive.GetNbinsY() +1):
            xraybin = h_xray.GetBinContent(i, j)
            xtalkbin = h_xtalk.GetBinContent(i, j)
            frbiasbin = h_frbias.GetBinContent(i, j)
            if ((xraybin) and (not xtalkbin) and (not frbiasbin)):
                h_xray_exclusive.SetBinContent(i, j, 1)
            elif ((not xraybin) and (xtalkbin) and (not frbiasbin)):
                h_xtalk_exclusive.SetBinContent(i, j, 1)
            elif ((not xraybin) and (not xtalkbin) and (frbiasbin)):
                h_frbias_exclusive.SetBinContent(i, j, 1)
            elif ((xraybin) and (xtalkbin) and (not frbiasbin)):
                h_xray_xtalk.SetBinContent(i, j, 1)
            elif ((xraybin) and (not xtalkbin) and (frbiasbin)):
                h_xray_frbias.SetBinContent(i, j, 1)
            elif ((not xraybin) and (xtalkbin) and (frbiasbin)):
                h_xtalk_frbias.SetBinContent(i, j, 1)
            elif ((xraybin) and (xtalkbin) and (frbiasbin)):
                h_xray_xtalk_frbias.SetBinContent(i, j, 1)

    # Set custom color palettes
    def set_palette(histogram, color_0, color_1):
        ROOT.gStyle.SetNumberContours(2)
        ROOT.gStyle.SetOptStat(0)
        colors = array('i', [color_0, color_1])
        ROOT.gStyle.SetPalette(len(colors), colors)
        histogram.SetContour(2, array('d', [0.0, 1.0]))
        histogram.SetMinimum(0)
        histogram.SetMaximum(1)

    # Define colors
    color_sets = {
        "h_xray_exclusive": (ROOT.TColor.GetColor("#FFFFFF"), ROOT.TColor.GetColor("#0000FF")),  # White for 0, Blue for 1
        "h_xtalk_exclusive": (ROOT.TColor.GetColor("#FFFFFF"), ROOT.TColor.GetColor("#FF0000")),  # Red for 1
        "h_frbias_exclusive": (ROOT.TColor.GetColor("#FFFFFF"), ROOT.TColor.GetColor("#00FF00")),   # Green for 1
        "h_xray_xtalk": (ROOT.TColor.GetColor("#FFFFFF"), ROOT.TColor.GetColor("#FF00FF")),  # Magenta for 1
        "h_xray_frbias": (ROOT.TColor.GetColor("#FFFFFF"), ROOT.TColor.GetColor("#00FFFF")),  # Cyan for 1
        "h_xtalk_frbias": (ROOT.TColor.GetColor("#FFFFFF"), ROOT.TColor.GetColor("#FFFF00")),  # Yellow for 1
        "h_xray_xtalk_frbias": (ROOT.TColor.GetColor("#FFFFFF"), ROOT.TColor.GetColor("#000000")),  # Black for 1
    }

    # Create a canvas
    c = ROOT.TCanvas("c", "Canvas", h_xray_exclusive.GetNbinsX(), h_xray_exclusive.GetNbinsY())
    
    # Draw and save histograms

    histnames = ["h_xray_exclusive", "h_xtalk_exclusive", "h_frbias_exclusive", "h_xray_xtalk", "h_xray_frbias", "h_xtalk_frbias", "h_xray_xtalk_frbias"]
    
    for hist, name in zip(histograms, histnames):
        # Set the palette for the current histogram
        set_palette(hist, *color_sets[f"{name}"])
        c.Clear()
        hist.Draw("COLZ")
        c.SaveAs(f"{png_file_path_base}_{name}.png")

    # Create an output file to save the result histogram
    output_file = ROOT.TFile(output_file_path, "RECREATE")
    h_xray_exclusive.Write()
    h_xray_exclusive.Write()
    h_xtalk_exclusive.Write()
    h_xray_xtalk.Write()
    h_xray_frbias.Write()
    h_xtalk_frbias.Write()
    h_xray_xtalk_frbias.Write()
    output_file.Close()
    
    # Close the input files
    xray_file.Close()
    xtalk_file.Close()
    frbias_file.Close()

    print("Comparison complete. Result saved to: ", output_file_path)

# Example usage
xray_file_path = "/home/kalib/Analysis/outputroot/xray/xrayroot12.root"
xtalk_file_path = "/home/kalib/Analysis/outputroot/xtalk/h_missing2dC12.root"
frbias_file_path = "/home/kalib/Analysis/outputroot/frbias/histograms.root"
output_file_path = "/home/kalib/Analysis/outputroot/RH0026C12Comparison.root"
png_file_path_base = "/home/kalib/Analysis/results/histogram"

compare_xtalk_xray(xray_file_path, xtalk_file_path, frbias_file_path, output_file_path)
