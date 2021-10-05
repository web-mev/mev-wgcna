suppressMessages(suppressWarnings(library(optparse)))
suppressMessages(suppressWarnings(library(WGCNA)))
suppressMessages(suppressWarnings(library("rjson", character.only=T, warn.conflicts = F, quietly = T)))

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
    ),
    make_option(
        c('-b','--beta'),
        type='numeric',
        help="Beta value for setting network connectivity. If not set, will use an auto-selection."
    )
)

opt <- parse_args(OptionParser(option_list=option_list))

if('beta' %in% names(opt)){
    # if beta is given as something that can't be parsed with as.numeric,
    # then it is assigned NA. NA will ultimately trigger the heuristic auto-detection
    # option below
    user_beta = as.numeric(opt$beta)
} else {
    user_beta = NA
}

# change the working directory to co-locate with the counts file:
working_dir <- dirname(opt$input_file)
setwd(working_dir)

# Load the expression data
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

# Create powers for 1 - 30. WGCNA blocks choosing of powers over 30 anyway.
powers = c(1:30)
# Calculate the soft thresholds
sft = pickSoftThreshold(wgcna_matrix, powerVector = powers, verbose = 5)

# Identify the first max delta delta
# Estimate maximum acceleration from the i+1 diff
x <- sft$fitIndices[,1]
y <- -sign(sft$fitIndices[,3])*sft$fitIndices[,2]

if(is.na(user_beta)){

    # Instead of trying to calculate numerical derivatives, the code below
    # outlines a 'voting' heuristic. In essense, we expect a plot that has a
    # logistic-like shape and we want to pick the elbow where we get diminishing
    # returns. To find this elbow, we calculate linear regressions using different
    # numbers of points. Most of those lines coming from the lower and upper half 
    # will intersect *roughly* near the elbow. We bin all those intersection points
    # and pick the densest. This gives us an x-coordinate/beta which can be used
    # to set the scale of the network.

    # Calculate the linear fits. lh=left-hand, rh=right-hand
    lh_fits = c()
    rh_fits = c()
    L = length(y)
    for(i in 2:L){
        xvals_lh = x[i:L]
        yvals_lh = y[i:L]
        xvals_rh = x[L-i:L]
        yvals_rh = y[L-i:L]
        tryCatch({
            lh_f = lm(yvals_lh ~ xvals_lh)
            rh_f = lm(yvals_rh ~ xvals_rh)
            lh_fits = c(lh_fits, lh_f)
            rh_fits = c(rh_fits, rh_f)
        }, error=function(err){
        }
        )
    }

    # Now that we have a bunch of fitted regression lines, find all the pairwise
    # intersections of those lines and store them
    intersections_x = c()
    intersections_y = c()
    for(i in 1:length(lh_fits)){
        for(j in i:length(rh_fits)){
            lc = lh_fits[i]
            rc = rh_fits[j]
            BL_0 = lc$coefficients[[1]]
            BL_1 = lc$coefficients[[2]]
            BR_0 = rc$coefficients[[1]]
            BR_1 = rc$coefficients[[2]]
            x0 = (BR_0 - BL_0)/(BL_1 - BR_1)
            y0 = BL_0 + BL_1*x0
            intersections_x = c(intersections_x, x0)
            intersections_y = c(intersections_y, y0)
        }
    }

    # make a dataframe so we can line-up those intersections up and remove NAs
    intersection_df = na.omit(data.frame(x=intersections_x, y=intersections_y))
    nbins = 20
    xbins = seq(min(x), max(x), length.out=nbins)
    ybins = seq(min(y), max(y), length.out=nbins)
    cuts_x = cut(intersection_df$x, xbins)
    cuts_y = cut(intersection_df$y, ybins)

    # create an abundance table based on the binned data
    tt = table(cuts_x, cuts_y)

    # Finding the max index flattens the matrix, so we have to use 
    # integer division +  modulus to get the actual row/col index
    max_idx = which.max(tt)
    row_idx = max_idx %% (dim(tt)[1])
    col_idx = 1 + max_idx %/% (dim(tt)[1])

    # double check that our indexing operation actually found the max location
    if(tt[row_idx, col_idx] != max(tt)){
        message('An error occurred during the auto-selection of the thresholding parameter. Try setting a manual default.')
        quit(status=1)
    }
    beta = ceiling(xbins[row_idx])
} else {
    beta = user_beta
}

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

# Munge the wgcna matrix into numeric (in case it's an int)
vec <- seq(1, dim(wgcna_matrix)[2])
wgcna_matrix[, vec] <- apply(wgcna_matrix[, vec, drop=F], 2, as.numeric)


bwnet = blockwiseModules(
    wgcna_matrix, 
    maxBlockSize = 5000,
    corType = "bicor",
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
modules_df = data.frame(group=bwnet$colors, gene=names(bwnet$colors))
summary_df = aggregate(x=modules_df$gene, by=list(modules_df$group), FUN=paste, collapse=',')
q = apply(summary_df, 1, function(r){
            list(
                    module=r['Group.1'][[1]],
                    genes=r['x'][[1]]
                )
        }
    )
results_json_str <- toJSON(q)
output_filename = 'wgcna_modules.json'
results_json_file <- paste(working_dir, output_filename, sep='/')
write(results_json_str, results_json_file)

# save the data for the network topology so users can evaluate
# the choice of beta
thresholding <- list(x=x,y=y,beta=beta)
output_filename = 'network_threshold_metadata.json'
thresholds_json_file <- paste(working_dir, output_filename, sep='/')
write(toJSON(thresholding), thresholds_json_file)
# for WebMEV compatability, need to create an outputs.json file.
json_str = paste0(
       '{"module_results":"', 
       results_json_file, 
       '",',
       '"network_connectivity_thresholds":"',
        thresholds_json_file,
       '"}'
)
output_json <- paste(working_dir, 'outputs.json', sep='/')
write(json_str, output_json)
