
# IMPORT REQUIRED PACKAGES #
# After loading into the correct directory

install.packages("ggseqlogo")
devtools::install_github("omarwagih/ggseqlogo")

require(ggplot2)
require(ggseqlogo)

library(ggplot2)
library(ggseqlogo)


# BUILD THE SEQUENCE LOGO #

B = read.delim("B_L100_K101_NN.txt") #file for each subtype and presence/absence of treatment for all 16 NNIBP residues with the correspondent nucleotide sequence of that residue(s)

# Set the parameters:
# annotate is going to create a rectangle filled with whatever colour you want to label the residue you're working with
# geom_logo is going to plots the sequence logo as a layer on ggplot, 
  # as, in this case, we're working with nucleotide sequences, we want the sequence type (seq_type) to be "dna", but amino acid ("aa"), rna ("rna") and other ("other") sequences types can be used
  # method is the height method you want to use, which can be bits or probability
  # stack_width sets the width of letter stack and must be between 0 and 1
# scale_x_continuous positions scales for continuous data - in this case, the x axis
  # breaks sets a numeric vector of positions, meaning how many nucleotides (in this case) you want to show in the sequence logo
  # labels, in this instance, are representing the nucleotide position number in the entire RT (reverse transcriptase) genome sequence

ggplot() + 
  annotate('rect', xmin = 2.5, xmax = 5.5, ymin = -0.05, ymax = 1.05, alpha = .1, col='black', fill='yellow') +
  geom_logo(B, seq_type="dna",method="prob", stack_width = 0.90) + 
  scale_x_continuous(breaks=1:11, labels=296:306) +
  theme_logo()

# (this process is then repeated for all the files, meaning all the set conditions and residues)