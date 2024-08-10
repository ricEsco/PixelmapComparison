import ROOT
from array import array

def compare_xtalk_xray(xray_file_path, xtalk_file_path, output_file_path):
    # Load the ROOT files
    xray_file = ROOT.TFile(xray_file_path, "READ")
    xtalk_file = ROOT.TFile(xtalk_file_path, "READ")

    
    h_xray = xray_file.Get("MissingMap")
    h_xtalk = xtalk_file.Get("h_confirmed2D")

    # Clone the X-ray histogram to create a result histogram
    h_both = h_xray.Clone("h_both")
    h_xray_exclusive = h_xray.Clone("h_xray_exclusive")
    h_xtalk_exclusive = h_xray.Clone("h_xtalk_exclusive")
    
    for i in range(1, h_both.GetNbinsX() + 1):
        for j in range(1, h_both.GetNbinsY() +1):
            xraybin = h_xray.GetBinContent(i, j)
            xtalkbin = h_xtalk.GetBinContent(i, j)
            h_both.SetBinContent(i, j, int(xraybin or xtalkbin))

    for i in range(1, h_both.GetNbinsX() + 1):
        for j in range(1, h_both.GetNbinsY() +1):
            xraybin = h_xray.GetBinContent(i, j)
            xtalkbin = h_xtalk.GetBinContent(i, j)
            h_xray_exclusive.SetBinContent(i, j, int(xraybin or not(xtalkbin)))

    for i in range(1, h_both.GetNbinsX() + 1):
        for j in range(1, h_both.GetNbinsY() +1):
            xraybin = h_xray.GetBinContent(i, j)
            xtalkbin = h_xtalk.GetBinContent(i, j)
            h_xtalk_exclusive.SetBinContent(i, j, int(xtalkbin or not(xraybin)))

    for i in range(1, h_both.GetNbinsX() + 1):
        for j in range(1, h_both.GetNbinsY() +1):
            if not h_both.GetBinContent(i, j):
                h_both.SetBinContent(i, j, 0.001)

    for i in range(1, h_xray_exclusive.GetNbinsX() + 1):
        for j in range(1, h_xray_exclusive.GetNbinsY() +1):
            if not h_xray_exclusive.GetBinContent(i, j):
                h_xray_exclusive.SetBinContent(i, j, 0.001)

    for i in range(1, h_xtalk_exclusive.GetNbinsX() + 1):
        for j in range(1, h_xtalk_exclusive.GetNbinsY() +1):
            if not h_xtalk_exclusive.GetBinContent(i, j):
                h_xtalk_exclusive.SetBinContent(i, j, 0.001)


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
        "h_both": (ROOT.TColor.GetColor("#FF00FF"), ROOT.TColor.GetColor("#FFFFFF")),  # Magenta for 0, White for 1
        "h_xray_exclusive": (ROOT.TColor.GetColor("#0000FF"), ROOT.TColor.GetColor("#FFFFFF")),  # Blue for 0, White for 1
        "h_xtalk_exclusive": (ROOT.TColor.GetColor("#FF0000"), ROOT.TColor.GetColor("#FFFFFF"))   # Red for 0, White for 1
    }


    # Create a canvas
    c = ROOT.TCanvas("c", "Canvas", h_both.GetNbinsX(), h_both.GetNbinsY())
    
    # Draw and save histograms
    histograms = [h_both, h_xray_exclusive, h_xtalk_exclusive]
    histogram_names = ["h_both", "h_xray_exclusive", "h_xtalk_exclusive"]
    
    for hist, name in zip(histograms, histogram_names):
        # Set the palette for the current histogram
        set_palette(hist, *color_sets[name])
        c.Clear()
        hist.Draw("COLZ")
        c.SaveAs(f"{png_file_path_base}_{name}.png")

    # Create an output file to save the result histogram
    output_file = ROOT.TFile(output_file_path, "RECREATE")
    h_both.Write()
    h_xray_exclusive.Write()
    h_xtalk_exclusive.Write()
    output_file.Close()
    
    # Close the input files
    xray_file.Close()
    xtalk_file.Close()

    print("Comparison complete. Result saved to: ", output_file_path)

# Example usage
xray_file_path = "/home/kalib/Analysis/outputroot/xrayroot12.root"
xtalk_file_path = "/home/kalib/Analysis/outputroot/h_missing2dC12.root"
output_file_path = "/home/kalib/Analysis/outputroot/RH0026C12Comparison.root"
png_file_path_base = "/home/kalib/Analysis/results/histogram"

compare_xtalk_xray(xray_file_path, xtalk_file_path, output_file_path)
