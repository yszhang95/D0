Out-dated
# General Description
 This is a repo to optimize the mva selection with the existing TTree files.
 Several steps are implemented.
 1. Read MVA trees to get the histograms prepared, specifically, 
 TH3 object, defined as DCA(z) vs Mass(x) and Mva(y). 
 All the ouput would be stored in the file `prefix_hists.root`
 2. Read the output file from the first step as the input.
 It would extract the DCA distribution of signal by fitting the invariant mass in bins of DCA.
 Also it would project 
 And the histograms would be written into the output file, `prefix_dca_hists.root`.
 3. Read the histograms created by the step 2. Then fit the fraction of non-prompt (and prompt) D0 out.
Taking the fraction into account, combined with the invariant mass fitting result, 
give the significance curve.

# Usage
