

# A treemap R script produced by the REVIGO server at http://revigo.irb.hr/
# If you found REVIGO useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from REVIGO. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","freqInDbPercent","abslog10pvalue","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0007049","cell cycle",1.410,9.8539,0.894,0.000,"cell cycle"),
c("GO:0051301","cell division",1.357,4.0410,0.894,0.118,"cell cycle"),
c("GO:0042590","antigen processing and presentation of exogenous peptide antigen via MHC class I",0.001,3.3799,0.931,0.000,"antigen processing and presentation of exogenous peptide antigen via MHC class I"),
c("GO:0044419","interspecies interaction between organisms",1.784,3.3497,0.974,0.000,"interspecies interaction between organisms"),
c("GO:0051641","cellular localization",1.638,4.0438,0.928,0.000,"cellular localization"),
c("GO:1902580","single-organism cellular localization",0.296,4.3372,0.869,0.285,"cellular localization"),
c("GO:0033036","macromolecule localization",1.969,3.1385,0.928,0.352,"cellular localization"),
c("GO:0071840","cellular component organization or biogenesis",5.425,3.0223,0.975,0.000,"cellular component organization or biogenesis"),
c("GO:0048285","organelle fission",0.071,5.4498,0.833,0.025,"organelle fission"),
c("GO:0006996","organelle organization",0.932,4.1568,0.830,0.506,"organelle fission"),
c("GO:0051128","regulation of cellular component organization",1.020,3.7100,0.666,0.651,"organelle fission"),
c("GO:0007032","endosome organization",0.002,3.5884,0.856,0.518,"organelle fission"),
c("GO:0045898","regulation of RNA polymerase II transcriptional preinitiation complex assembly",0.000,3.3507,0.698,0.307,"organelle fission"),
c("GO:0044265","cellular macromolecule catabolic process",1.109,4.4737,0.728,0.032,"cellular macromolecule catabolism"),
c("GO:0007264","small GTPase mediated signal transduction",0.234,3.9136,0.733,0.651,"cellular macromolecule catabolism"),
c("GO:0043065","positive regulation of apoptotic process",0.021,3.9101,0.646,0.622,"cellular macromolecule catabolism"),
c("GO:0045862","positive regulation of proteolysis",0.012,4.0645,0.542,0.133,"cellular macromolecule catabolism"),
c("GO:0051340","regulation of ligase activity",0.002,3.8097,0.751,0.208,"cellular macromolecule catabolism"),
c("GO:0031331","positive regulation of cellular catabolic process",0.013,3.2692,0.591,0.639,"cellular macromolecule catabolism"),
c("GO:0000209","protein polyubiquitination",0.010,3.0467,0.701,0.243,"cellular macromolecule catabolism"),
c("GO:0035556","intracellular signal transduction",2.718,4.0600,0.693,0.273,"cellular macromolecule catabolism"),
c("GO:0009057","macromolecule catabolic process",1.643,3.8386,0.823,0.500,"cellular macromolecule catabolism"),
c("GO:0031145","anaphase-promoting complex-dependent proteasomal ubiquitin-dependent protein catabolic process",0.003,4.4237,0.711,0.547,"cellular macromolecule catabolism"),
c("GO:1903047","mitotic cell cycle process",0.077,8.0550,0.558,0.088,"mitotic cell cycle process"),
c("GO:0000278","mitotic cell cycle",0.087,8.0600,0.657,0.591,"mitotic cell cycle process"),
c("GO:0022402","cell cycle process",0.527,7.1798,0.623,0.684,"mitotic cell cycle process"),
c("GO:0045786","negative regulation of cell cycle",0.022,3.5317,0.528,0.539,"mitotic cell cycle process"),
c("GO:0000075","cell cycle checkpoint",0.020,3.6144,0.657,0.665,"mitotic cell cycle process"),
c("GO:0044770","cell cycle phase transition",0.030,6.9136,0.632,0.683,"mitotic cell cycle process"),
c("GO:0031577","spindle checkpoint",0.006,4.0429,0.652,0.616,"mitotic cell cycle process"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$abslog10pvalue <- as.numeric( as.character(stuff$abslog10pvalue) );
stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="revigo_treemap.pdf", width=16, height=9 ) # width and height are in inches

# check the tmPlot command documentation for all possible parameters - there are a lot more
tmPlot(
	stuff,
	index = c("representative","description"),
	vSize = "abslog10pvalue",
	type = "categorical",
	vColor = "representative",
	title = "REVIGO Gene Ontology treemap",
	inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
	lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
	bg.labels = "#CCCCCCAA",     # define background color of group labels
												       # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
	position.legend = "none"
)

dev.off()
