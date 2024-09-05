import ROOT
import numpy as np

# Create a ROOT file to store the histograms
output_file = ROOT.TFile("toy_histograms.root", "RECREATE")

# Define histogram parameters
nbins_x = 432
nbins_y = 336

# Create three TH2F histograms
# The arguments for the TH2F constructor are: name, title, number of bins in x, xlow, xhigh, number of bins in y, ylow, yhigh
hist1 = ROOT.TH2F("hist1", "Histogram 1", nbins_x, 0, nbins_x, nbins_y, 0, nbins_y)
hist2 = ROOT.TH2F("hist2", "Histogram 2", nbins_x, 0, nbins_x, nbins_y, 0, nbins_y)
hist3 = ROOT.TH2F("hist3", "Histogram 3", nbins_x, 0, nbins_x, nbins_y, 0, nbins_y)

# Fill the histograms
for i in range(1, nbins_x + 1):
    for j in range(1, nbins_y + 1):
        # Fill hist1 with 1s in the first half and 0s in the second half
        value1 = 1 if i <= nbins_x // 2 else 0
        hist1.SetBinContent(i, j, value1)

        # Fill hist2 with 0s in the first half and 1s in the second half
        value2 = 1 if j > nbins_y // 2 else 0
        hist2.SetBinContent(i, j, value2)

        # Fill hist3 with a circle in the middle
        value3 = 1 if (i - nbins_x // 2) ** 2 + (j - nbins_y // 2) ** 2 <= 150 ** 2 else 0
        hist3.SetBinContent(i, j, value3)

# Write histograms to the file and close it
hist1.Write()
hist2.Write()
hist3.Write()
output_file.Close()

print("Toy histograms have been created and saved to toy_histograms.root")
