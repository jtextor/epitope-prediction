# EpitopePrediction

Predicts binding of 9- to 12-mer peptides to MHC class I molecules using the [stabilized matrix method](http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-10-394) (B. Peters and colleagues). Can predict binding for several alleles from humans (HLA A/B), mice (H-2), chimpanzees (Patr A/B), and rhesus macaques (Mamu A/B).

# License

This package is licensed under a CC BY-NC-SA 4.0 license. This means it is free for use *for non-commercial purposes* (such as academic research). Companies interested in using this package should contact [Bjoern Peters](http://www.liai.org/pages/faculty-peters). 

# Example


	# Install the package
	devtools::install_github("jtextor/epitope-prediction")

	# Load the package
	library( EpitopePrediction )

	# This is the CORE protein from the Hepatitis C virus reference sequence available
	# at https://hcv.lanl.gov/content/sequence/LOCATE/locate.html
	hcv.core <- paste("MSTNPKPQRKTKRNTNRRPQDVKFPGGGQIVGGVYLLPRRGPRLGVRATRKTSERSQPRGRR",
		"QPIPKARRPEGRTWAQPGYPWPLYGNEGCGWAGWLLSPRGSRPSWGPTDPRRRSRNLGKVIDTLTCGFADLMGYIP",
		"LVGAPLGGAARALAHGVRVLEDGVNYATGNLPGCSFSIFLLALLSCLTVPASA",sep="")

	# Predict 9-mer binders to human HLA-A02:01
	binders( hcv.core, "HLA-A-02:01", 9 )
