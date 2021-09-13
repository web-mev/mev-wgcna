suppressMessages(suppressWarnings(library(optparse)))
suppressMessages(suppressWarnings(library(WGCNA)))

# Enables parallelizable calculations with WGCNA
enableWGCNAThreads()

option_list <- list(
    make_option(
        c('-f','--input_file'),
        help='Path to the count matrix input.'
    ),
    make_option(
        c('-m','--max_var_genes'),
        default=5000, # Maybe we want this by default to be bigger
        help="Max number of top most variable genes to include."
    )
)

# Load the pre-normalized expression
# Can be RNAseq, microarray, doesn't matter. Just normalize in WebMev, then
# use the normalized count matrix as input.
exprs <- read.table(
    opt$input_file,
    sep="\t",
    header=T,
    row.names=1,
    stringsAsFactors=F
)
exprs.flt <- exprs[
    order(
        apply(exprs, 1, mad), decreasing = T
    )[1:opt$max_var_genes], 
]

# Note: do we want to drop mitochondrial genes first?

# First, filter the expression matrix for top variable genes
# Then transpose the matrix
wgcna_matrix <- t(exprs.flt)

# Create a correlation network
# Use biweight midcorrelation. Median based, so less sensitive to outliers.
# Absolute value to reduce everything to positive 0-1 interpretation
s <- abs(bicor(wgcna_matrix))

# Create powers fomr 1 - 25
powers = c(c(1:30))
# Calculate the soft thresholds
sft = pickSoftThreshold(wgcna_matrix, powerVector = powers, verbose = 5)
# Identify the first max delta delta
# Estimate maximum acceleration from the i+1 diff
x <- sft$fitIndices[,1]
y <- -sign(sft$fitIndices[,3])*sft$fitIndices[,2]
d.first <- y[-1] - y[-length(y)]
d.second <- d.first[-1] - d.first[-length(d.first)]
d.strength <- d.second - d.first[-1]

# Recalculate the corr matrix with exponentiation by chosen beta
beta = which.max(d.strength) + 2 # Add two to compensate for shift
a = s^beta
w = 1-a

# Munge the exprs.flt into numeric
vec <- seq(1, length(colnames(exprs.flt)))
exprs.flt[, vec] <- apply(exprs.flt[, vec, drop=F], 2, as.numeric)
# Then transpose and convert to matrix
exprs.flt <- as.matrix(t(exprs.flt))

# Parameters
# maxBlockSize:
#   2000 for 2G?
#   5000 for 4G
#   20000 for 16G
#   30000 for 32G
# power: calculated beta value from above
# corType: bidweight midcorrelation ("pearson" is also an option)
# pearsonFallback: if MAD is zero, fallback to pearson for that individual gene
#   "all" if you want to fallback to Pearson for all genes if a single 0 MAD
# minModuleSize: opt?
# TOMType & networkType : WGCNA recommends signed (or signed hybrid)
#   (despite the default being unsigned)
#   opt?
# mergeCutHeight: set to 0.15, opt? (0.15 is a common cutoff)
#   Manages the dendrogram cut height for merging
#   No real perfect single value
# nThreads: 0 = dynamically determine threading
# verbose: 0 = silent (higher for more verbosity)
# There are many more parameters, but this should be sufficient
bwnet = blockwiseModules(
    exprs.flt, 
    maxBlockSize = 2000,
    corType = "bicor"
    pearsonFallback = "individual",
    power = beta, 
    networkType = "signed",
    TOMType = "signed",
    minModuleSize = 30,
    reassignThreshold = 1e-6, 
    mergeCutHeight = 0.15,
    nThreads = 0,
    numericLabels = TRUE,
    saveTOMs = FALSE,
    verbose = 0
)

# Colors denotes the label / module the gene was assigned.
# Export as JSON
bwnet$colors