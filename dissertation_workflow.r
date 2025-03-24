# --------------------------------------------------------------------------------TABLE OF CONTENTS------------------------------------------------------------------------------------ #
# 1 Data Collection (Digitisation) ---------------------------------------------------------------- Line 30
# 2 Data Processing ------------------------------------------------------------------------------- Line 54
#   2.1 Correct File Name & Examine Data Structure ------------------------------------------------ Line 60
#   2.2 shapRead + shapFix ------------------------------------------------------------------------ Line 148
#   2.3 GPA --------------------------------------------------------------------------------------- Line 203
# 3 Data Analyses --------------------------------------------------------------------------------- Line 347
#   3.1 PCA --------------------------------------------------------------------------------------- Line 353
#       3.1.1 Taxonomic Order Mapping (6 Orders) -------------------------------------------------- Line 356
#       3.1.2 Temperal Interval Mapping (6 Timeseries) -------------------------------------------- Line 704
#       3.1.3 Palaeogeographic Mapping (4 Palaeocontinents) --------------------------------------- Line 900
#       3.1.4 Thin-Plate Spline Deformation Grid -------------------------------------------------- Line 1140
#       3.1.5 Total PC Variances ------------------------------------------------------------------ Line 1270
#       3.1.6 SoV + Rate Plot (Based on Three Classification Schemes) ----------------------------- Line 1363
#       3.1.7 SoV + Rate Plot (Only Combined Palaeocontinents + Time) ----------------------------- Line 2017
#       3.1.8 Statistical Analyses Of SoV + Rate -------------------------------------------------- Line 2520
#   3.2 Phylogenetic Tree ------------------------------------------------------------------------- Line 3480
#       3.2.1 Generate Morphological Tree (Using PCA Scores) -------------------------------------- Line 3483
#       3.2.2 Phylogenetic Tree With Evolutionary Rate -------------------------------------------- Line 3565
# 4 Variability Analyses -------------------------------------------------------------------------- Line 3785
# 5 Validation Analyses --------------------------------------------------------------------------- Line 4190
#   5.1 Validation Analyses Of Substitution Locus (Landmark 11-12) -------------------------------- Line 4195
#   5.2 Validation Analyses Of Substitution Facial Suture (Connected by LM10 + LM12) -------------- Line 4444
#   5.3 Validation Analyses Of Substitution Scheme For Suborder Olenellina In PCA ----------------- Line 4680



# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 1 Data Collection (Digitisation)
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Clean environment and set up workspace
rm(list=ls(all=TRUE))
graphics.off()

# Install and load required packages
install.packages("geomorph")
install.packages("StereoMorph")
library(geomorph)
library(StereoMorph)

# Set working directory
setwd("~/Desktop/dissertation_project_trilobite")

# Digitisation of Trilobites' 2D Figures
digitizeImage(image.file='trilobite_figure_dataset', shapes.file='xml_output_dataset',landmarks.ref='landmarks.txt', curves.ref='curves.txt')



# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 2 Data Processing
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------------------
# 2.1 Correct File Names & Examine Data Structures
# -------------------------------------------------------------------------------------------------
# ----------------------------------------
# Process file names
# ----------------------------------------
# Get file list
files <- list.files("xml_output_dataset", pattern = "\\.txt$", full.names = TRUE)
cat("Found", length(files), "text files\n")

# Process file names (Remove extra .txt extensions)
file_names <- basename(files)
clean_names <- gsub("\\.txt\\.txt$", ".txt", file_names)
clean_names <- gsub("\\.txt$", "", clean_names)

# ----------------------------------------
# Examine data structure of all specimens
# ----------------------------------------
# Import TriloMorph functions (Serra et al. 2023) (WARNING!!! You must download Serra's R script (TriloMorph_Functions.r) and ensure correct file path)
source("~/Downloads/TriloMorph_Functions.r")

# Read all specimens
all_specimens <- shapRead(clean_names, subdir = "xml_output_dataset")

# Create storage structures
dimensions <- numeric()
fixed_landmarks <- numeric()
curve_points <- list()

# Examine each specimen
for(specimen_name in names(all_specimens)) {
  specimen <- all_specimens[[specimen_name]]
  
  # Store basic information
  dimensions <- c(dimensions, specimen$k)
  fixed_landmarks <- c(fixed_landmarks, specimen$plm)
  
  # Store curve points
  curve_points[[specimen_name]] <- specimen$cv.pts
}

# ----------------------------------------
# Examination results
# ----------------------------------------
# Create matrix to store all curve points
curve_matrix <- matrix(NA, 
                      nrow = length(curve_points), 
                      ncol = max(sapply(curve_points, length)),
                      dimnames = list(names(curve_points), 
                                    c("glabella", "suture", "anterior", "posterior")))

# Fill the matrix
for(i in names(curve_points)) {
  curve_matrix[i, 1:length(curve_points[[i]])] <- curve_points[[i]]
}

# ----------------------------------------
# Output examination results
# ----------------------------------------
cat("\nData Dimension Analysis:\n")
cat("Dimensions (k):", unique(dimensions), "\n")
cat("Number of fixed landmarks:", unique(fixed_landmarks), "\n")

cat("\nPoint count statistics for each curve:\n")
for(curve in 1:ncol(curve_matrix)) {
  cat("\nCurve", colnames(curve_matrix)[curve], ":\n")
  cat("Minimum:", min(curve_matrix[,curve], na.rm=TRUE), "\n")
  cat("Maximum:", max(curve_matrix[,curve], na.rm=TRUE), "\n")
  cat("Mean:", round(mean(curve_matrix[,curve], na.rm=TRUE), 2), "\n")
}

# ----------------------------------------
# Generate nlms template
# ----------------------------------------
nlms <- c(unique(dimensions)[1],                    # Dimensions
          unique(fixed_landmarks)[1],               # Number of fixed landmarks
          min(curve_matrix[,1], na.rm=TRUE),        # Minimum points for glabella curve
          min(curve_matrix[,2], na.rm=TRUE),        # Minimum points for facial suture curve
          min(curve_matrix[,3], na.rm=TRUE),        # Minimum points for anterior margin curve
          min(curve_matrix[,4], na.rm=TRUE))        # Minimum points for posterior margin curve

cat("\nRecommended nlms template:\n")
cat("nlms <- c(", paste(nlms, collapse=", "), ")\n")

# Save examination results
write.csv(curve_matrix, "landmark_examination_results.csv")


# -------------------------------------------------------------------------------------------------
# 2.2 Standardise Landmarks Dataset By Utilising shapRead & shapFix (Serra et al. 2023) (It can be considered as the first step in GPA)
# -------------------------------------------------------------------------------------------------
# ----------------------------------------
# Set up the correct landmark template
# ----------------------------------------
nlms <- c(2,    # Dimensions
          16,   # Fixed landmarks
          95,   # Minimum points for glabella curve
          199,  # Minimum points for facial suture curve
          173,  # Minimum points for anterior margin curve
          134)  # Minimum points for posterior margin curve

# ----------------------------------------
# shapRead
# ----------------------------------------
# Read the landmark data for all specimens
all_specimens <- shapRead(clean_names, subdir = "xml_output_dataset")

# Print summary information about all loaded specimens
cat("\nLoaded Specimens Summary:\n")
cat("Number of specimens:", length(all_specimens), "\n")

# Check for any loading issues
if(length(attr(all_specimens, "specimens with missing datafile")) > 0) {
    cat("\nWarning: Missing datafiles for:\n")
    print(attr(all_specimens, "specimens with missing datafile"))
}
if(length(attr(all_specimens, "specimens with failed loading")) > 0) {
    cat("\nWarning: Failed loading for:\n")
    print(attr(all_specimens, "specimens with failed loading"))
}

# ----------------------------------------
# shapFix
# ----------------------------------------
# Process and standardise the landmark data
cat("\nProcessing landmarks with shapFix...\n")
processed_specimens <- shapFix(all_specimens, nlms, lm.scale = TRUE)

# Print summary information about the processed data
cat("\nProcessed Data Summary:\n")
cat("Array dimensions:", paste(dim(processed_specimens), collapse=" x "), "\n")

# Check for removed specimens
if(!is.null(attr(processed_specimens, "removed"))) {
    cat("\nSpecimens removed during processing:\n")
    print(attr(processed_specimens, "removed"))
}

# Save the processed data
save(processed_specimens, nlms, file = "processed_landmarks.RData")
cat("\nProcessing complete. Data saved to 'processed_landmarks.RData'\n")


# -------------------------------------------------------------------------------------------------
# 2.3 GPA
# -------------------------------------------------------------------------------------------------
# Load saved processed landmarks (if needed)
if(!exists("processed_specimens")) {
    load("processed_landmarks.RData")
}

# -----------------------------------------------------------------------------
# Perform GPA analysis
# -----------------------------------------------------------------------------
# Perform Generalised Procrustes Analysis
gpan <- geomorph::gpagen(processed_specimens, 
                        Proj = TRUE,      # Project onto tangent space
                        PrinAxes = FALSE) # Don't align by principal axes

# Calculate mean shape and Procrustes distances
mean_shape <- mshape(gpan$coords)
Pd <- c() # Pd Stands for Procrustes Distances
for(i in 1:dim(gpan$coords)[3]) {
    Pd[i] <- sqrt(sum((gpan$coords[,,i] - mean_shape)^2))
}

# Calculate statistics for Procrustes distances
mean_pd <- mean(Pd)
sd_pd <- sd(Pd)
upper_threshold <- mean_pd + 2*sd_pd

# -----------------------------------------------------------------------------
# Create visualisation (Two plots side by side)
# -----------------------------------------------------------------------------
# Create new plotting window
dev.new(width = 12, height = 6)
par(mfrow = c(1,2))

# ----------------------------------------
# Plot 1: Superimposition Plot
# ----------------------------------------
# Set up empty plot
plot(NA, NA, 
     xlim = range(gpan$coords[,1,]), 
     ylim = range(gpan$coords[,2,]),
     xlab = "X coordinate (Procrustes aligned)",
     ylab = "Y coordinate (Procrustes aligned)",
     main = "Procrustes Superimposition of Trilobite Cephalon")

# Plot individual specimens
for(i in 1:dim(gpan$coords)[3]) {
    points(gpan$coords[,,i], 
           pch = 20,       
           col = "gray50", 
           cex = 0.3)      
}

# Plot consensus shape
points(mean_shape, 
       pch = 16,        
       cex = 1.2,       
       col = "black")   

# Add sample size
mtext(paste0("n = ", dim(gpan$coords)[3]), 
      side = 3, 
      adj = 1, 
      font = 3)

# Add legend at top right
legend("bottomleft", 
       legend = c("Consensus shape", "Individual specimens"),
       col = c("black", "gray50"), 
       pch = c(16, 20),
       pt.cex = c(1.2, 0.3),
       cex = 0.8,
       bty = "n",
       inset = 0.02)

# ----------------------------------------
# Plot 2: Procrustes Distances Histogram
# ----------------------------------------
# Create histogram
hist(Pd, 
     main = "Procrustes Distances from Consensus",
     xlab = "Procrustes Distance",
     ylab = "Frequency",
     ylim = c(0, 20),     # Set y-axis limit from 0 to 25
     xlim = c(0.1, 0.8),  # Adjust x-axis to fit all data
     breaks = sqrt(dim(gpan$coords)[3]),
     col = "lightblue",
     border = "white")

# Add vertical lines to maximum height
segments(x0 = mean_pd, y0 = 0, x1 = mean_pd, y1 = 19, col = "red", lwd = 2)
segments(x0 = upper_threshold, y0 = 0, x1 = upper_threshold, y1 = 19, 
         col = "red", lty = 2, lwd = 2)

# Add text annotations
# Mean label on the left side of its line
text(mean_pd, 18.7, paste("Mean =", round(mean_pd, 3)), 
     pos = 2, # Left side
     col = "red", 
     cex = 0.8)

# Outlier threshold label on the right side of its line
text(upper_threshold, 18.7, paste("Outlier threshold =", round(upper_threshold, 3)), 
     pos = 4,  # Right side
     col = "red", 
     cex = 0.8)

# SD between mean and threshold lines (SD Stands for Standard Deviation)
text((mean_pd + upper_threshold)/2, 18.8,  # Position between the two lines
     paste("SD =", round(sd_pd, 3)),
     col = "black", 
     cex = 0.8)

# Return to default plotting parameters
par(mfrow = c(1,1))

# Print Analysis Results
cat("\nGPA Analysis Summary:\n")
cat("----------------------------------------\n")
cat("Number of specimens:", dim(gpan$coords)[3], "\n")
cat("Procrustes Sum of Squares:", gpan$procD, "\n")
cat("Tangent Sum of Squares:", gpan$tangentD, "\n")

cat("\nProcrustes Distances from Consensus:\n")
cat("Mean:", round(mean_pd, 4), "\n")
cat("SD:", round(sd_pd, 4), "\n")
cat("Range:", round(range(Pd), 4), "\n")

# Identify potential outliers
potential_outliers <- which(Pd > upper_threshold)
if(length(potential_outliers) > 0) {
    cat("\nPotential outliers (>2SD from mean):\n")
    cat("Specimen IDs:\n")
    print(dimnames(gpan$coords)[[3]][potential_outliers])
    cat("\nTheir Procrustes distances:\n")
    print(round(Pd[potential_outliers], 4))
}

# Save results
save(gpan, Pd, mean_shape, file = "gpa_analysis_results.RData")


# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 3 Data Analyses
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------------------
# 3.1 PCA
# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------
# 3.1.1 Taxonomic Order Mapping (6 Orders)
# -------------------------------------------------------------------------------------------------
# Set working directory
setwd("~/Desktop/dissertation_project_trilobite")
# Load required packages
library(geomorph)
# Load the GPA results
load("gpa_analysis_results.RData")

# ----------------------------------------
# Create mapping dictionary
# ----------------------------------------
taxa_to_order <- list(
   # EODISCIDA
    "Hebediscina_longispinus" = "EODISCIDA",
    "Neocobboldia_chinlinica" = "EODISCIDA",
    "Pagetia_ocellata" = "EODISCIDA",
    "Tsunyidiscus_longquanensis" = "EODISCIDA",
    
    # REDLICHIIDA
    "Agraulos_ceticephalus" = "REDLICHIIDA",
    "Balcoracania_dailyi" = "REDLICHIIDA",
    "Bigotina_bivallata" = "REDLICHIIDA",
    "Choubertella_spinosa" = "REDLICHIIDA",
    "Daguinaspis_ambroggii" = "REDLICHIIDA",
    "Dolerolenus_zoppii" = "REDLICHIIDA",
    "Eccaparadoxides_pusillus" = "REDLICHIIDA",
    "Ellipsocephalus_gripi" = "REDLICHIIDA",
    "Elliptocephala_asaphoides" = "REDLICHIIDA",
    "Emuella_polymera" = "REDLICHIIDA",
    "Eoredlichia_intermediata" = "REDLICHIIDA",
    "Estaingia_bilobata" = "REDLICHIIDA",
    "Fallotaspis_bondoni" = "REDLICHIIDA",
    "Hamatolenus_vincenti" = "REDLICHIIDA",
    "Holmia_kjerulfi" = "REDLICHIIDA",
    "Judomia_tera" = "REDLICHIIDA",
    "Lemdadella_antarcticae" = "REDLICHIIDA",
    "Megapharanaspis_nedini" = "REDLICHIIDA",
    "Metadoxides_armatus" = "REDLICHIIDA",
    "Montezumaspis_parallela" = "REDLICHIIDA",
    "Nephrolenellus_geniculatus" = "REDLICHIIDA",
    "Nevadia_weeksi" = "REDLICHIIDA",
    "Olenellus_gilberti" = "REDLICHIIDA",
    "Parabadiella_huoi" = "REDLICHIIDA",
    "Peachella_iddingsi" = "REDLICHIIDA",
    "Protolenus_densigranulatus" = "REDLICHIIDA",
    "Redlichia_takooensis" = "REDLICHIIDA",
    "Strenuaeva_inflata" = "REDLICHIIDA",
    "Wanneria_walcottana" = "REDLICHIIDA",
    "Xystridura_templetonensis" = "REDLICHIIDA",
    "Yunnanocephalus_yunnanensis" = "REDLICHIIDA",
    "Zhangshania_typica" = "REDLICHIIDA",
    
    # CORYNEXOCHIDA
    "Corynexochus_plumula" = "CORYNEXOCHIDA",
    "Dinesus_ida" = "CORYNEXOCHIDA",
    "Fuchouia_fecunda" = "CORYNEXOCHIDA",
    "Guangxiaspis_guangxiensis" = "CORYNEXOCHIDA",
    "Jakutus_primigenius" = "CORYNEXOCHIDA",
    "Olenoides_serratus" = "CORYNEXOCHIDA",
    "Oryctocephalus_indicus" = "CORYNEXOCHIDA",
    "Shergoldia_laevigata" = "CORYNEXOCHIDA",
    "Wanshania_wanshanensis" = "CORYNEXOCHIDA",
    
    # AULACOPLEURIDA
    "Coosella_kieri" = "AULACOPLEURIDA",
    "Elrathia_kingii" = "AULACOPLEURIDA",
    "Modocia_kohli" = "AULACOPLEURIDA",
    "Sao_hirsuta" = "AULACOPLEURIDA",
    "Tricrepicephalus_texanus" = "AULACOPLEURIDA",
    
    # OLENIDA
    "Asaphiscus_wheeleri" = "OLENIDA",
    "Burnetiella_leechi" = "OLENIDA",
    "Cedaria_minor" = "OLENIDA",
    "Haniwa_quadrata" = "OLENIDA",
    "Iveria_iverensis" = "OLENIDA",
    "Labiostria_westropi" = "OLENIDA",
    "Loganellus_logani" = "OLENIDA",
    "Olenus_wahlenbergi" = "OLENIDA",
    "Orygmaspis_billingsi" = "OLENIDA",
    "Pterocephalia_norfordi" = "OLENIDA",
    
    # UNCERTAIN
    "Bathynotus_kueichouensis" = "UNCERTAIN",
    "Bolaspidella_housensis" = "UNCERTAIN",
    "Brachyaspidion_microps" = "UNCERTAIN",
    "Catillicephala_shawi" = "UNCERTAIN",
    "Cheilocephalus_brevilobus" = "UNCERTAIN",
    "Dikelocephalus_minnesotensis" = "UNCERTAIN",
    "Elvinia_roemeri" = "UNCERTAIN",
    "Eurekia_rintintini" = "UNCERTAIN",
    "Hungaia_puelchana" = "UNCERTAIN",
    "Jiulongshania_longispina" = "UNCERTAIN",
    "Lisania_paibiensis" = "UNCERTAIN",
    "Myopsolenites_altus" = "UNCERTAIN",
    "Neodrepanura_premesnili" = "UNCERTAIN",
    "Norwoodia_boninoi" = "UNCERTAIN",
    "Palaeadotes_hunanensis" = "UNCERTAIN",
    "Papyriaspis_lanceola" = "UNCERTAIN",
    "Penarosa_netenta" = "UNCERTAIN",
    "Plethopeltis_armatus" = "UNCERTAIN",
    "Ptychoparia_striata" = "UNCERTAIN",
    "Tasmacephalus_platypus" = "UNCERTAIN",
    "Xingrenaspis_xingrenensis" = "UNCERTAIN"
)

# ----------------------------------------
# Create specimen to order mapping
# ----------------------------------------
# Get specimen names
specimens <- dimnames(gpan$coords)[[3]]

# Create mapping data frame
specimen_orders <- data.frame(
    specimen = specimens,
    order = NA_character_,
    stringsAsFactors = FALSE
)

# Map specimens to orders
for(i in 1:nrow(specimen_orders)) {
    # Extract species name from specimen ID
    spec_name <- specimens[i]
    # Find matching taxon in our dictionary
    for(taxon in names(taxa_to_order)) {
        if(grepl(taxon, spec_name, ignore.case = TRUE)) {
            specimen_orders$order[i] <- taxa_to_order[[taxon]]
            break
        }
    }
}

# ----------------------------------------
# Define colors for orders
# ----------------------------------------
orders <- c(
    "EODISCIDA", "REDLICHIIDA", "CORYNEXOCHIDA", 
    "AULACOPLEURIDA", "OLENIDA", "UNCERTAIN"
)

order_colors <- c(
    "EODISCIDA" = "#E41A1C",      # Red
    "REDLICHIIDA" = "#4DAF4A",    # Green
    "CORYNEXOCHIDA" = "#377EB8",  # Blue
    "AULACOPLEURIDA" = "#FF7F00", # Orange
    "OLENIDA" = "#984EA3",        # Purple
    "UNCERTAIN" = "#999999"       # Gray
)

# Create transparent versions for hulls
hull_colors <- sapply(order_colors, function(x) {
    rgb(t(col2rgb(x))/255, alpha = 0.2)
})
# Function to draw convex hulls
draw_hull <- function(x, y, order) {
    order_points <- specimen_orders$order == order
    if(sum(order_points) > 2) {
        hpts <- chull(x[order_points], y[order_points])
        hpts <- c(hpts, hpts[1])
        polygon(x[order_points][hpts], y[order_points][hpts], 
               col = hull_colors[order], border = NA)
    }
}
# ----------------------------------------
# Check fault & save
# ----------------------------------------
# Print mapping results
cat("\nNumber of specimens per order:\n")
print(table(specimen_orders$order))

# Print unmapped specimens
unmapped <- specimen_orders[is.na(specimen_orders$order), ]
if(nrow(unmapped) > 0) {
    cat("\nWarning: The following specimens couldn't be mapped to orders:\n")
    print(unmapped)
}

# Save the mapping
save(specimen_orders, orders, order_colors, hull_colors, file = "specimen_order_mapping.RData")

# -----------------------------------------------------------------------------
# Performe PCA analysis & visualise PCA morphospace
# -----------------------------------------------------------------------------
# Load required packages
library(geomorph)
# Load GPA results and specimen-order mapping
load("gpa_analysis_results.RData")
load("specimen_order_mapping.RData")


# Perform PCA analysis
pcan <- geomorph::gm.prcomp(gpan$coords)
variances <- pcan$sdev^2                 # Calculate the standard deviations of the PC
prop_var <- variances/sum(variances)     # Calculate the proportion of variance explained by each PC
cum_var <- cumsum(prop_var)              # Calculate the cumulative proportion of variance explained

# Save PCA results
# Create a list to store all PCA-related results
pca_results <- list(
    # Store the PCA object itself
    pca = pcan,
    
    # Store the scores (coordinates in PC space)
    scores = pcan$x,
    
    # Store the loadings (eigenvectors)
    loadings = pcan$rotation,
    
    # Store variance information
    variance = list(
        eigenvalues = variances,
        proportion = prop_var,
        cumulative = cum_var
    ),
    
    # Store the center (mean shape)
    center = pcan$center,
    
    # Store the scale
    scale = pcan$scale
)

# Add metadata
pca_results$metadata <- list(
    date = Sys.Date(),
    n_specimens = dim(gpan$coords)[3],
    n_landmarks = dim(gpan$coords)[1],
    n_dimensions = dim(gpan$coords)[2]
)

# Save the results
save(pca_results, file = "pca_analysis_results.RData")

# Print summary to confirm
cat("\nPCA Results saved successfully.\n")
cat("Summary of saved data:\n")
cat("Number of specimens:", pca_results$metadata$n_specimens, "\n")
cat("Number of landmarks:", pca_results$metadata$n_landmarks, "\n")
cat("Number of PCs:", length(pca_results$variance$eigenvalues), "\n")
cat("\nVariance explained by first 4 PCs:\n")
for(i in 1:4) {
    cat(sprintf("PC%d: %.1f%% (Cumulative: %.1f%%)\n",
                i, 
                pca_results$variance$proportion[i]*100,
                pca_results$variance$cumulative[i]*100))
}

# -----------------------------------------------------------------------------
# Visualise PCA morphospace
# -----------------------------------------------------------------------------
# Set up plotting window
dev.new(width = 12, height = 6)
par(mfrow = c(1,2))

# ----------------------------------------
# Plot 1: PC1 vs PC2
# ----------------------------------------
# Set up empty plot
plot(NA, NA,
     xlim = range(pcan$x[,1]),
     ylim = range(pcan$x[,2]),
     xlab = paste0("PC1 (", round(prop_var[1]*100, 1), "%)"),
     ylab = paste0("PC2 (", round(prop_var[2]*100, 1), "%)"),
     main = "PC1 vs PC2",
     asp = 1)

# Add zero lines
abline(h = 0, v = 0, lty = 2, col = "gray")

# Draw hulls for each order
for(order in orders) {
    draw_hull(pcan$x[,1], pcan$x[,2], order)
}

# Add points
points(pcan$x[,1], pcan$x[,2],
       pch = 19,  # Solid circle
       col = order_colors[specimen_orders$order],
       cex = 1.2)

# Add legend
legend("topright", 
       legend = orders,
       col = order_colors[orders],
       pch = 19,
       cex = 0.8,
       bty = "n")

# Add sample size
mtext(paste0("n = ", nrow(pcan$x)), side = 3, adj = 1, font = 3)

# ----------------------------------------
# Plot 2: PC3 vs PC4
# ----------------------------------------
# Set up empty plot
plot(NA, NA,
     xlim = range(pcan$x[,3]),
     ylim = range(pcan$x[,4]),
     xlab = paste0("PC3 (", round(prop_var[3]*100, 1), "%)"),
     ylab = paste0("PC4 (", round(prop_var[4]*100, 1), "%)"),
     main = "PC3 vs PC4",
     asp = 1)

# Add zero lines
abline(h = 0, v = 0, lty = 2, col = "gray")

# Draw hulls for each order
for(order in orders) {
    draw_hull(pcan$x[,3], pcan$x[,4], order)
}

# Add points
points(pcan$x[,3], pcan$x[,4],
       pch = 19,  # Solid circle
       col = order_colors[specimen_orders$order],
       cex = 1.2)

# Add legend
legend("topright", 
       legend = orders,
       col = order_colors[orders],
       pch = 19,
       cex = 0.8,
       bty = "n")

# Add sample size
mtext(paste0("n = ", nrow(pcan$x)), side = 3, adj = 1, font = 3)

# Reset plotting parameters
par(mfrow = c(1,1))


# Print analysis results
cat("\nPCA Analysis Summary:\n")
cat("----------------------------------------\n")
cat("Number of specimens:", nrow(pcan$x), "\n")
cat("\nVariance explained by first 4 PCs:\n")
for(i in 1:4) {
    cat(sprintf("PC%d: %.1f%% (Cumulative: %.1f%%)\n",
                i, prop_var[i]*100, sum(prop_var[1:i]*100)))
}

# Save the plots
dev.copy2pdf(file = "pca_morphospace.pdf")


# -------------------------------------------------------------------------------------------------
# 3.1.2 Temperal Interval Mapping (6 Timeseries)
# -------------------------------------------------------------------------------------------------
# Set working directory
setwd("~/Desktop/dissertation_project_trilobite")
# Load required packages
library(geomorph)
# Load the GPA results
load("gpa_analysis_results.RData")
load("pca_analysis_results.RData")

# Create time series mapping
taxa_to_time <- list(
    # Stage 3
    "Bigotina_bivallata" = "Stage_3",
    "Choubertella_spinosa" = "Stage_3",
    "Daguinaspis_ambroggii" = "Stage_3",
    "Dolerolenus_zoppii" = "Stage_3",
    "Eoredlichia_intermediata" = "Stage_3",
    "Fallotaspis_bondoni" = "Stage_3",
    "Hebediscina_longispinus" = "Stage_3",
    "Jakutus_primigenius" = "Stage_3",
    "Judomia_tera" = "Stage_3",
    "Lemdadella_antarcticae" = "Stage_3",
    "Metadoxides_armatus" = "Stage_3",
    "Montezumaspis_parallela" = "Stage_3",
    "Nevadia_weeksi" = "Stage_3",
    "Parabadiella_huoi" = "Stage_3",
    "Tsunyidiscus_longquanensis" = "Stage_3",
    "Yunnanocephalus_yunnanensis" = "Stage_3",
    "Zhangshania_typica" = "Stage_3",
    
    # Stage 4
    "Balcoracania_dailyi" = "Stage_4",
    "Bathynotus_kueichouensis" = "Stage_4",
    "Ellipsocephalus_gripi" = "Stage_4",
    "Elliptocephala_asaphoides" = "Stage_4",
    "Emuella_polymera" = "Stage_4",
    "Estaingia_bilobata" = "Stage_4",
    "Hamatolenus_vincenti" = "Stage_4",
    "Holmia_kjerulfi" = "Stage_4",
    "Megapharanaspis_nedini" = "Stage_4",
    "Myopsolenites_altus" = "Stage_4",
    "Neocobboldia_chinlinica" = "Stage_4",
    "Nephrolenellus_geniculatus" = "Stage_4",
    "Olenellus_gilberti" = "Stage_4",
    "Peachella_iddingsi" = "Stage_4",
    "Protolenus_densigranulatus" = "Stage_4",
    "Redlichia_takooensis" = "Stage_4",
    "Strenuaeva_inflata" = "Stage_4",
    "Wanneria_walcottana" = "Stage_4",
    
    # Wuliuan
    "Dinesus_ida" = "Wuliuan",
    "Olenoides_serratus" = "Wuliuan",
    "Oryctocephalus_indicus" = "Wuliuan",
    "Pagetia_ocellata" = "Wuliuan",
    "Xingrenaspis_xingrenensis" = "Wuliuan",
    "Xystridura_templetonensis" = "Wuliuan",
    
    # Drumian
    "Agraulos_ceticephalus" = "Drumian",
    "Asaphiscus_wheeleri" = "Drumian",
    "Bolaspidella_housensis" = "Drumian",
    "Brachyaspidion_microps" = "Drumian",
    "Eccaparadoxides_pusillus" = "Drumian",
    "Elrathia_kingii" = "Drumian",
    "Fuchouia_fecunda" = "Drumian",
    "Modocia_kohli" = "Drumian",
    "Papyriaspis_lanceola" = "Drumian",
    "Penarosa_netenta" = "Drumian",
    "Ptychoparia_striata" = "Drumian",
    "Sao_hirsuta" = "Drumian",
    
    # Guzhangian
    "Catillicephala_shawi" = "Guzhangian",
    "Cedaria_minor" = "Guzhangian",
    "Coosella_kieri" = "Guzhangian",
    "Jiulongshania_longispina" = "Guzhangian",
    "Lisania_paibiensis" = "Guzhangian",
    "Neodrepanura_premesnili" = "Guzhangian",
    "Norwoodia_boninoi" = "Guzhangian",
    "Palaeadotes_hunanensis" = "Guzhangian",
    "Tasmacephalus_platypus" = "Guzhangian",
    "Tricrepicephalus_texanus" = "Guzhangian",
    "Wanshania_wanshanensis" = "Guzhangian",
    
    # Furongian
    "Corynexochus_plumula" = "Furongian",
    "Cheilocephalus_brevilobus" = "Furongian",
    "Olenus_wahlenbergi" = "Furongian",
    "Burnetiella_leechi" = "Furongian",
    "Elvinia_roemeri" = "Furongian",
    "Eurekia_rintintini" = "Furongian",
    "Guangxiaspis_guangxiensis" = "Furongian",
    "Haniwa_quadrata" = "Furongian",
    "Iveria_iverensis" = "Furongian",
    "Labiostria_westropi" = "Furongian",
    "Loganellus_logani" = "Furongian",
    "Orygmaspis_billingsi" = "Furongian",
    "Pterocephalia_norfordi" = "Furongian",
    "Shergoldia_laevigata" = "Furongian",
    "Dikelocephalus_minnesotensis" = "Furongian",
    "Hungaia_puelchana" = "Furongian",
    "Plethopeltis_armatus" = "Furongian"
)

# Create color scheme for time periods
time_colors <- c(
    "Stage_3" = "#1B9E77",      # Dark teal
    "Stage_4" = "#D95F02",      # Orange
    "Wuliuan" = "#7570B3",      # Purple
    "Drumian" = "#E7298A",      # Pink
    "Guzhangian" = "#66A61E",   # Green
    "Furongian" = "#666666"     # Grey
)

# Create mapping for specimens
specimens <- dimnames(gpan$coords)[[3]]
specimen_times <- data.frame(
    specimen = specimens,
    time_period = NA_character_,
    stringsAsFactors = FALSE
)

# Map specimens to time periods
for(i in 1:nrow(specimen_times)) {
    spec_name <- specimens[i]
    for(taxon in names(taxa_to_time)) {
        if(grepl(taxon, spec_name, ignore.case = TRUE)) {
            specimen_times$time_period[i] <- taxa_to_time[[taxon]]
            break
        }
    }
}

# Hull drawing function
draw_hull <- function(x, y, period) {
    period_points <- specimen_times$time_period == period
    period_points[is.na(period_points)] <- FALSE
    if(sum(period_points, na.rm = TRUE) > 2) {
        hpts <- chull(x[period_points], y[period_points])
        hpts <- c(hpts, hpts[1])
        lines(x[period_points][hpts], y[period_points][hpts], 
              col = time_colors[period],
              lwd = 1.5)
    }
}

# Save the mapping
save(specimen_times, time_colors, file = "specimen_time_mapping.RData")


# Perform PCA
pcan <- geomorph::gm.prcomp(gpan$coords)
variances <- pcan$sdev^2
prop_var <- variances/sum(variances)
cum_var <- cumsum(prop_var)

# Plot PCA results (PC1 vs PC2)
dev.new(width = 8, height = 6)
plot(NA, NA,
     xlim = range(pcan$x[,1]),
     ylim = range(pcan$x[,2]),
     xlab = paste0("PC1 (", round(prop_var[1]*100, 1), "%)"),
     ylab = paste0("PC2 (", round(prop_var[2]*100, 1), "%)"),
     main = "PC1 vs PC2",
     asp = 1)

abline(h = 0, v = 0, lty = 2, col = "gray")

for(period in names(time_colors)) {
    draw_hull(pcan$x[,1], pcan$x[,2], period)
}

points(pcan$x[,1], pcan$x[,2],
       pch = 19,
       col = time_colors[specimen_times$time_period],
       cex = 1.2)

legend("topright", 
       legend = names(time_colors),
       col = time_colors,
       pch = 19,
       cex = 0.8,
       bty = "n")

mtext(paste0("n = ", nrow(pcan$x)), side = 3, adj = 1, font = 3)

# Save results
save(pcan, prop_var, cum_var, file = "timeseries_pca_results.RData")

# Save plot
dev.copy2pdf(file = "timeseries_pca_morphospace.pdf")


# -------------------------------------------------------------------------------------------------
# 3.1.3 Palaeogeographic Mapping (4 Palaeocontinents)
# -------------------------------------------------------------------------------------------------
# Set working directory
setwd("~/Desktop/dissertation_project_trilobite")

# Load required packages
library(geomorph)

# Load the GPA and PCA results 
load("pca_analysis_results.RData")

# Create palaeocontinent mapping
taxa_to_palaeocontinent <- list(
    # BALTICA
    "Ellipsocephalus_gripi" = "BALTICA",
    "Holmia_kjerulfi" = "BALTICA",
    "Strenuaeva_inflata" = "BALTICA", 
    "Olenus_wahlenbergi" = "BALTICA",
    "Judomia_tera" = "BALTICA",
    "Jakutus_primigenius" = "BALTICA",
    
    # CHINA
    "Hebediscina_longispinus" = "CHINA",
    "Neocobboldia_chinlinica" = "CHINA",
    "Tsunyidiscus_longquanensis" = "CHINA",
    "Eoredlichia_intermediata" = "CHINA",
    "Parabadiella_huoi" = "CHINA",
    "Redlichia_takooensis" = "CHINA",
    "Yunnanocephalus_yunnanensis" = "CHINA",
    "Zhangshania_typica" = "CHINA",
    "Guangxiaspis_guangxiensis" = "CHINA",
    "Shergoldia_laevigata" = "CHINA",
    "Wanshania_wanshanensis" = "CHINA",
    "Haniwa_quadrata" = "CHINA",
    "Bathynotus_kueichouensis" = "CHINA",
    "Jiulongshania_longispina" = "CHINA",
    "Lisania_paibiensis" = "CHINA",
    "Neodrepanura_premesnili" = "CHINA",
    "Palaeadotes_hunanensis" = "CHINA",
    "Xingrenaspis_xingrenensis" = "CHINA",
    
    # GONDWANA
    "Pagetia_ocellata" = "GONDWANA",
    "Agraulos_ceticephalus" = "GONDWANA",
    "Balcoracania_dailyi" = "GONDWANA",
    "Bigotina_bivallata" = "GONDWANA",
    "Choubertella_spinosa" = "GONDWANA",
    "Daguinaspis_ambroggii" = "GONDWANA",
    "Dolerolenus_zoppii" = "GONDWANA",
    "Eccaparadoxides_pusillus" = "GONDWANA",
    "Emuella_polymera" = "GONDWANA",
    "Estaingia_bilobata" = "GONDWANA",
    "Fallotaspis_bondoni" = "GONDWANA",
    "Hamatolenus_vincenti" = "GONDWANA",
    "Lemdadella_antarcticae" = "GONDWANA",
    "Megapharanaspis_nedini" = "GONDWANA",
    "Metadoxides_armatus" = "GONDWANA",
    "Protolenus_densigranulatus" = "GONDWANA",
    "Xystridura_templetonensis" = "GONDWANA",
    "Dinesus_ida" = "GONDWANA",
    "Fuchouia_fecunda" = "GONDWANA",
    "Sao_hirsuta" = "GONDWANA",
    "Iveria_iverensis" = "GONDWANA",
    "Hungaia_puelchana" = "GONDWANA",
    "Myopsolenites_altus" = "GONDWANA",
    "Papyriaspis_lanceola" = "GONDWANA",
    "Penarosa_netenta" = "GONDWANA",
    "Ptychoparia_striata" = "GONDWANA",
    "Tasmacephalus_platypus" = "GONDWANA",
    
    # LAURENTIA
    "Elliptocephala_asaphoides" = "LAURENTIA",
    "Montezumaspis_parallela" = "LAURENTIA",
    "Nephrolenellus_geniculatus" = "LAURENTIA",
    "Nevadia_weeksi" = "LAURENTIA",
    "Olenellus_gilberti" = "LAURENTIA",
    "Peachella_iddingsi" = "LAURENTIA",
    "Wanneria_walcottana" = "LAURENTIA",
    "Corynexochus_plumula" = "LAURENTIA",
    "Olenoides_serratus" = "LAURENTIA",
    "Oryctocephalus_indicus" = "LAURENTIA",
    "Coosella_kieri" = "LAURENTIA",
    "Elrathia_kingii" = "LAURENTIA",
    "Modocia_kohli" = "LAURENTIA",
    "Tricrepicephalus_texanus" = "LAURENTIA",
    "Asaphiscus_wheeleri" = "LAURENTIA",
    "Burnetiella_leechi" = "LAURENTIA",
    "Cedaria_minor" = "LAURENTIA",
    "Labiostria_westropi" = "LAURENTIA",
    "Loganellus_logani" = "LAURENTIA",
    "Orygmaspis_billingsi" = "LAURENTIA",
    "Pterocephalia_norfordi" = "LAURENTIA",
    "Bolaspidella_housensis" = "LAURENTIA",
    "Brachyaspidion_microps" = "LAURENTIA",
    "Catillicephala_shawi" = "LAURENTIA",
    "Cheilocephalus_brevilobus" = "LAURENTIA",
    "Dikelocephalus_minnesotensis" = "LAURENTIA",
    "Elvinia_roemeri" = "LAURENTIA",
    "Eurekia_rintintini" = "LAURENTIA",
    "Norwoodia_boninoi" = "LAURENTIA",
    "Plethopeltis_armatus" = "LAURENTIA"
)

# Define colors for palaeocontinents
palaeocontinent_colors <- c(
    "BALTICA" = "#4575B4",    # Deep blue
    "CHINA" = "#D73027",      # Deep red
    "GONDWANA" = "#91CF60",   # Green
    "LAURENTIA" = "#FEE090"   # Golden yellow
)

# Create mapping for specimens
specimens <- dimnames(pca_results$scores)[[1]]
specimen_palaeocontinents <- data.frame(
    specimen = specimens,
    palaeocontinent = NA_character_,
    stringsAsFactors = FALSE
)

# Map specimens to palaeocontinents
for(i in 1:nrow(specimen_palaeocontinents)) {
    spec_name <- specimens[i]
    for(taxon in names(taxa_to_palaeocontinent)) {
        if(grepl(taxon, spec_name, ignore.case = TRUE)) {
            specimen_palaeocontinents$palaeocontinent[i] <- taxa_to_palaeocontinent[[taxon]]
            break
        }
    }
}

# Function to draw gradient-edged hull
draw_gradient_hull <- function(x, y, continent) {
    continent_points <- specimen_palaeocontinents$palaeocontinent == continent
    if(sum(continent_points) > 2) {
        # Get hull points
        hpts <- chull(x[continent_points], y[continent_points])
        hpts <- c(hpts, hpts[1])
        
        # Create polygon with very light fill
        polygon(x[continent_points][hpts], y[continent_points][hpts],
               col = adjustcolor(palaeocontinent_colors[continent], alpha.f = 0.1),
               border = NA)
        
        # Draw border with gradient effect
        for(i in 1:(length(hpts)-1)) {
            x1 <- x[continent_points][hpts[i]]
            y1 <- y[continent_points][hpts[i]]
            x2 <- x[continent_points][hpts[i+1]]
            y2 <- y[continent_points][hpts[i+1]]
            
            # Draw multiple lines with decreasing opacity
            for(j in 1:3) {
                lines(c(x1, x2), c(y1, y2),
                      col = adjustcolor(palaeocontinent_colors[continent], 
                                      alpha.f = 0.8/j),
                      lwd = 4/j,
                      lty = j)
            }
        }
    }
}

# Save the mapping
save(specimen_palaeocontinents, palaeocontinent_colors, file = "specimen_palaeocontinent_mapping.RData")


# Create PCA plot
pdf("palaeocontinent_pca_morphospace.pdf", width = 8, height = 6)

# Calculate proportion of variance explained
prop_var <- pca_results$variance$proportion

# Set up empty plot
plot(NA, NA,
     xlim = range(pca_results$scores[,1]),
     ylim = range(pca_results$scores[,2]),
     xlab = paste0("PC1 (", round(prop_var[1]*100, 1), "%)"),
     ylab = paste0("PC2 (", round(prop_var[2]*100, 1), "%)"),
     main = "Palaeocontinent Distribution in Morphospace",
     asp = 1)

# Add zero lines
abline(h = 0, v = 0, lty = 2, col = "gray80")

# Draw gradient-edged hulls for each palaeocontinent
for(continent in names(palaeocontinent_colors)) {
    draw_gradient_hull(pca_results$scores[,1], pca_results$scores[,2], continent)
}

# Add points with buffer zones
for(continent in names(palaeocontinent_colors)) {
    continent_points <- specimen_palaeocontinents$palaeocontinent == continent
    
    # Add buffer zones (slightly larger points in lighter color)
    points(pca_results$scores[continent_points, 1], 
           pca_results$scores[continent_points, 2],
           pch = 21,
           bg = adjustcolor(palaeocontinent_colors[continent], alpha.f = 0.3),
           col = adjustcolor(palaeocontinent_colors[continent], alpha.f = 0.5),
           cex = 2.2)
    
    # Add main points
    points(pca_results$scores[continent_points, 1], 
           pca_results$scores[continent_points, 2],
           pch = 21,
           bg = palaeocontinent_colors[continent],
           col = "black",
           cex = 1.5)
}

# Add legend
legend("topright", 
       legend = names(palaeocontinent_colors),
       col = palaeocontinent_colors,
       pt.bg = palaeocontinent_colors,
       pch = 21,
       pt.cex = 1.5,
       cex = 0.8,
       bty = "n")

# Add sample size
mtext(paste0("n = ", nrow(pca_results$scores)), 
      side = 3, adj = 1, font = 3)

# Print summary
cat("\nPalaeocontinent Distribution Summary:\n")
cat("----------------------------------------\n")
print(table(specimen_palaeocontinents$palaeocontinent))

# Check for unmapped specimens
unmapped <- specimen_palaeocontinents[is.na(specimen_palaeocontinents$palaeocontinent), ]
if(nrow(unmapped) > 0) {
    cat("\nWarning: The following specimens couldn't be mapped to palaeocontinents:\n")
    print(unmapped)
}

dev.off()


# -------------------------------------------------------------------------------------------------
# 3.1.4 Thin-Plate Spline Deformation Grid
# -------------------------------------------------------------------------------------------------
# Set working directory
setwd("~/Desktop/dissertation_project_trilobite")

# Load required packages
if (!require("fields")) install.packages("fields")
library(geomorph)
library(fields)

# Load the pre-processed data
load("gpa_analysis_results.RData")
load("pca_analysis_results.RData")

# Function to create and transform grid
create_grid <- function(ref, target, n.grid=30) {
    rng.x <- range(ref[,1])
    rng.y <- range(ref[,2])
    pad.x <- diff(rng.x) * 0.2
    pad.y <- diff(rng.y) * 0.2
    
    x <- seq(rng.x[1] - pad.x, rng.x[2] + pad.x, length=n.grid)
    y <- seq(rng.y[1] - pad.y, rng.y[2] + pad.y, length=n.grid)
    
    xm <- matrix(rep(x, n.grid), n.grid, n.grid)
    ym <- matrix(rep(y, each=n.grid), n.grid, n.grid)
    
    tps.x <- Tps(ref, target[,1])
    new.x <- predict.Krig(tps.x, cbind(c(xm), c(ym)))
    
    tps.y <- Tps(ref, target[,2])
    new.y <- predict.Krig(tps.y, cbind(c(xm), c(ym)))
    
    list(
        x = matrix(new.x, n.grid, n.grid),
        y = matrix(new.y, n.grid, n.grid),
        orig.x = xm,
        orig.y = ym
    )
}

# Function to create deformation grid plots
create_deformation_plots <- function(pc_values, pc_number, pc_shape) {
    # Calculate all deformed shapes first to determine unified plot ranges
    all_shapes <- lapply(pc_values, function(pc_val) {
        deformed_shape <- mean_shape + pc_val * pc_shape
        deformed_shape[,2] <- -deformed_shape[,2]
        return(deformed_shape)
    })
    
    # Calculate unified ranges
    x_range <- range(sapply(all_shapes, function(s) range(s[,1])))
    y_range <- range(sapply(all_shapes, function(s) range(s[,2])))
    
    # Create plot layout
    layout <- c(3,2)
 
    # Set up the plotting device
    pdf(sprintf("pc%d_deformation_grids.pdf", pc_number), width = 12, height = 6)
    par(mfrow = layout, mar = c(1,1,2,1))
    
    # Create visualizations for each PC value
    for(pc_val in pc_values) {
        # Calculate deformed shape
        deformed_shape <- mean_shape + pc_val * pc_shape
        
        # Flip y coordinates
        deformed_shape[,2] <- -deformed_shape[,2]
        mean_shape_flipped <- mean_shape
        mean_shape_flipped[,2] <- -mean_shape_flipped[,2]
        
        # Set up plot with unified ranges
        plot(NULL, 
             xlim = x_range * 1.2,
             ylim = rev(y_range * 1.2),
             asp = 1,
             xlab = "",
             ylab = "",
             main = sprintf("PC%d = %g", pc_number, pc_val),
             axes = FALSE,
             frame.plot = FALSE)
        
        # Create and transform grid using flipped coordinates
        grid <- create_grid(mean_shape_flipped, deformed_shape)
        
        # Draw grid lines
        for(i in 1:nrow(grid$x)) {
            lines(grid$x[i,], grid$y[i,], col="forestgreen", lwd=0.5)
            lines(grid$x[,i], grid$y[,i], col="forestgreen", lwd=0.5)
        }
        # Draw curves as lines
        curves_start <- 17
        curves_per_structure <- c(95, 199, 173, 134)
        current_pos <- curves_start
        
        for(n_points in curves_per_structure) {
            curve_indices <- current_pos:(current_pos + n_points - 1)
            lines(deformed_shape[curve_indices,], col="black", lwd=2)
            current_pos <- current_pos + n_points
        }
        # Add landmarks visualisation (landmarks are points 1-16)
        landmarks_indices <- 1:16
        points(deformed_shape[landmarks_indices, 1], 
               deformed_shape[landmarks_indices, 2],
               col="red", 
               pch=19,  # Filled circles
               cex=0.6) # Size of points
    }
    
    dev.off()
    cat(sprintf("Created PC%d deformation grids with unified scaling for values:", pc_number), 
        paste(pc_values, collapse=", "), "\n")
    cat(sprintf("Output saved to 'pc%d_deformation_grids.pdf'\n", pc_number))
}

# Set up parameters for PC1 and PC2
pc1_values <- c(-0.2, 0, 0.2, 0.4)
pc2_values <- c(-0.3, -0.2, -0.1, 0, 0.1, 0.2)

# Get mean shape and PC coefficients
mean_shape <- mshape(gpan$coords)
pc1_shape <- matrix(pca_results$pca$rotation[,1], ncol=2, byrow=TRUE)
pc2_shape <- matrix(pca_results$pca$rotation[,2], ncol=2, byrow=TRUE)

# Create deformation plots for PC1 and PC2
create_deformation_plots(pc1_values, 1, pc1_shape)
create_deformation_plots(pc2_values, 2, pc2_shape)


# -------------------------------------------------------------------------------------------------
# 3.1.5 Total PC Vairances
# -------------------------------------------------------------------------------------------------
# ----------------------------------------
# Generating an excel sheet containing PCA analysis data of each individual specimen
# ----------------------------------------
# Set working directory
setwd("~/Desktop/dissertation_project_trilobite")

# Install and load required packages
install.packages("geomorph")
install.packages("writexl")
library(geomorph)
library(writexl)

# Load data
load("gpa_analysis_results.RData")
load("specimen_order_mapping.RData")

# Perform PCA analysis
pcan <- geomorph::gm.prcomp(gpan$coords)
variances <- pcan$sdev^2
prop_var <- variances/sum(variances)

# Extracting PCA scores
pca_scores <- pcan$x

# Create a data frame
pca_table <- data.frame(
    Species = rownames(pca_scores),
    Order = specimen_orders$order
)

# Adding PCA scores
pca_table <- cbind(pca_table, as.data.frame(pca_scores))

# Rename PCA columns
colnames(pca_table)[3:ncol(pca_table)] <- paste0("PC", 1:(ncol(pca_table)-2))

# Adding the Explanation of Variance attribute
attr(pca_table, "variance_explained") <- prop_var * 100

# Create a table of explained variation
variance_table <- data.frame(
    PC = paste0("PC", 1:length(prop_var)),
    Variance_Explained = prop_var * 100,
    Cumulative_Variance = cumsum(prop_var) * 100
)

# Save as Excel file
write_xlsx(list(
    "PCA_Scores" = pca_table,
    "Variance_Explained" = variance_table
), "pca_analysis_results.xlsx")

# ----------------------------------------
# Generate variance histogram for the first 10 PCs
# ----------------------------------------
# Load required packages
library(geomorph)
library(ggplot2)

# Load data 
load("pca_analysis_results.RData")

# Extract first 10 PCs and their variance
pc_data <- data.frame(
    PC = factor(paste0("PC", 1:10), levels = paste0("PC", 1:10)),
    Variance = pca_results$variance$proportion[1:10] * 100
)

# Create plot
pdf("pc_variance_histogram.pdf", width = 10, height = 6)

ggplot(pc_data, aes(x = PC, y = Variance)) +
    geom_bar(stat = "identity", fill = "steelblue", width = 0.7) +
    geom_text(aes(label = sprintf("%.1f%%", Variance)), 
              vjust = -0.5, 
              size = 3) +
    theme_minimal() +
    labs(title = "Variance Explained by First 10 Principal Components",
         x = "Principal Component",
         y = "Variance Explained (%)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5)) +
    ylim(0, max(pc_data$Variance) * 1.1) # Add 10% space for labels

# Print summary
cat("\nVariance explained by first 10 PCs:\n")
print(pc_data)
dev.off()


# -------------------------------------------------------------------------------------------------
# 3.1.6 SoV + Rate Plot (Based On Three Classification Schemes)
# -------------------------------------------------------------------------------------------------
# Set working directory
setwd("~/Desktop/dissertation_project_trilobite")

# Load required packages
library(ape)
library(phytools)
library(ggplot2)
library(gridExtra)
library(cowplot)

# ----------------------------------------
# Load previously saved data
# ----------------------------------------
# Load the rate tree data
load("phylogenetic_tree_rates_palaeocontinents_orders.RData")

# Load time mapping
load("specimen_time_mapping.RData")

# Load order mapping
load("specimen_order_mapping.RData")

# Load palaeocontinent mapping
load("specimen_palaeocontinent_mapping.RData")

# Load PCA results
load("pca_analysis_results.RData") 

# ----------------------------------------
# Define time intervals with their boundaries (from original code)
# ----------------------------------------
time_intervals <- list(
  "Stage_3" = c(521, 514),       # From Stage 3 start to end
  "Stage_4" = c(514, 509),       # From Stage 4 start to end
  "Wuliuan" = c(509, 504.5),     # From Wuliuan start to end
  "Drumian" = c(504.5, 500.5),   # From Drumian start to end
  "Guzhangian" = c(500.5, 497),  # From Guzhangian start to end
  "Furongian" = c(497, 485)      # From Paibian start to Stage 10 end (combined Furongian)
)

# ----------------------------------------
# Define group orders and colors
# ----------------------------------------
# Orders
orders <- c("EODISCIDA", "REDLICHIIDA", "CORYNEXOCHIDA", 
           "AULACOPLEURIDA", "OLENIDA", "UNCERTAIN")

order_colors <- c(
    "EODISCIDA" = "#E41A1C",      # Red
    "REDLICHIIDA" = "#4DAF4A",    # Green
    "CORYNEXOCHIDA" = "#377EB8",  # Blue
    "AULACOPLEURIDA" = "#FF7F00", # Orange
    "OLENIDA" = "#984EA3",        # Purple
    "UNCERTAIN" = "#999999"       # Gray
)

# Time periods
time_periods <- c("Stage_3", "Stage_4", "Wuliuan", 
                 "Drumian", "Guzhangian", "Furongian")

time_colors <- c(
    "Stage_3" = "#1B9E77",      # Dark teal
    "Stage_4" = "#D95F02",      # Orange
    "Wuliuan" = "#7570B3",      # Purple
    "Drumian" = "#E7298A",      # Pink
    "Guzhangian" = "#66A61E",   # Green
    "Furongian" = "#666666"     # Grey
)

# Palaeocontinents
palaeocontinents <- c("BALTICA", "CHINA", "GONDWANA", "LAURENTIA")

palaeocontinent_colors <- c(
    "BALTICA" = "#4575B4",    # Deep blue
    "CHINA" = "#D73027",      # Deep red
    "GONDWANA" = "#91CF60",   # Green
    "LAURENTIA" = "#FEE090"   # Golden yellow
)

# ----------------------------------------
# Extract tip rates from the rate_tree
# ----------------------------------------
tip_rates <- rate_tree$edge.length[rate_tree$edge[,2] <= length(rate_tree$tip.label)]
tip_indices <- rate_tree$edge[,2][rate_tree$edge[,2] <= length(rate_tree$tip.label)]
tip_names <- rate_tree$tip.label[tip_indices]

# Create data frame with tip names and rates
tip_data <- data.frame(
  specimen = tip_names,
  rate = tip_rates,
  stringsAsFactors = FALSE
)

# ----------------------------------------
# Map Time Periods, Orders, and Palaeocontinents to tip data
# ----------------------------------------
# Time Periods
tip_data$time_period <- NA
for (i in 1:nrow(tip_data)) {
  species_name <- tip_data$specimen[i]
  for (j in 1:nrow(specimen_times)) {
    if (specimen_times$specimen[j] == species_name) {
      tip_data$time_period[i] <- specimen_times$time_period[j]
      break
    }
  }
}

# Orders
tip_data$order <- NA
for (i in 1:nrow(tip_data)) {
  species_name <- tip_data$specimen[i]
  for (j in 1:nrow(specimen_orders)) {
    if (specimen_orders$specimen[j] == species_name) {
      tip_data$order[i] <- specimen_orders$order[j]
      break
    }
  }
}

# Palaeocontinents
tip_data$palaeocontinent <- NA
for (i in 1:nrow(tip_data)) {
  species_name <- tip_data$specimen[i]
  for (j in 1:nrow(specimen_palaeocontinents)) {
    if (specimen_palaeocontinents$specimen[j] == species_name) {
      tip_data$palaeocontinent[i] <- specimen_palaeocontinents$palaeocontinent[j]
      break
    }
  }
}

# ----------------------------------------
# Calculate mean rate for each group
# ----------------------------------------
# Time periods
time_means <- aggregate(rate ~ time_period, data = tip_data, FUN = mean)
time_sds <- aggregate(rate ~ time_period, data = tip_data, FUN = sd)
time_counts <- aggregate(rate ~ time_period, data = tip_data, FUN = length)

# Orders
order_means <- aggregate(rate ~ order, data = tip_data, FUN = mean)
order_sds <- aggregate(rate ~ order, data = tip_data, FUN = sd)
order_counts <- aggregate(rate ~ order, data = tip_data, FUN = length)

# Palaeocontinents
palaeocontinent_means <- aggregate(rate ~ palaeocontinent, data = tip_data, FUN = mean)
palaeocontinent_sds <- aggregate(rate ~ palaeocontinent, data = tip_data, FUN = sd)
palaeocontinent_counts <- aggregate(rate ~ palaeocontinent, data = tip_data, FUN = length)

# ----------------------------------------
# Create summary data frames
# ----------------------------------------
# Time periods
time_summary <- data.frame(
  time_period = time_means$time_period,
  mean_rate = time_means$rate,
  sd_rate = time_sds$rate,
  count = time_counts$rate,
  mid_point = sapply(time_means$time_period, function(tp) {
    interval <- time_intervals[[tp]]
    if(!is.null(interval)) {
      return(mean(interval))
    } else {
      return(NA)
    }
  }),
  stringsAsFactors = FALSE
)

# Order the time periods chronologically
time_summary$time_period <- factor(time_summary$time_period, 
                                  levels = c("Stage_3", "Stage_4", "Wuliuan", 
                                            "Drumian", "Guzhangian", "Furongian"))
time_summary <- time_summary[order(time_summary$time_period), ]

# Orders
order_summary <- data.frame(
  order = order_means$order,
  mean_rate = order_means$rate,
  sd_rate = order_sds$rate,
  count = order_counts$rate,
  stringsAsFactors = FALSE
)

# Palaeocontinents
palaeocontinent_summary <- data.frame(
  palaeocontinent = palaeocontinent_means$palaeocontinent,
  mean_rate = palaeocontinent_means$rate,
  sd_rate = palaeocontinent_sds$rate,
  count = palaeocontinent_counts$rate,
  stringsAsFactors = FALSE
)

# ----------------------------------------
# Calculate Sum of Variances (SoV) from PCA scores
# ----------------------------------------
# Function to calculate SoV for each specimen
calculate_specimen_sov <- function(data) {
  n_specimens <- nrow(data)
  sov <- numeric(n_specimens)
  
  for(i in 1:n_specimens) {
    specimen_data <- as.numeric(data[i,])
    sov[i] <- var(specimen_data, na.rm = TRUE)
  }
  
  return(sov)
}

# Calculate SoV from PCA scores
pca_scores <- as.matrix(pca_results$scores)
sov_by_specimen <- calculate_specimen_sov(pca_scores)

# Create data frame with SoV for each specimen
disparity_data <- data.frame(
  specimen = rownames(pca_scores),
  sov = sov_by_specimen
)

# Map time periods, orders, and palaeocontinents to disparity data
# Time Periods
disparity_data$time_period <- NA
for (i in 1:nrow(disparity_data)) {
  species_name <- disparity_data$specimen[i]
  for (j in 1:nrow(specimen_times)) {
    if (specimen_times$specimen[j] == species_name) {
      disparity_data$time_period[i] <- specimen_times$time_period[j]
      break
    }
  }
}

# Orders
disparity_data$order <- NA
for (i in 1:nrow(disparity_data)) {
  species_name <- disparity_data$specimen[i]
  for (j in 1:nrow(specimen_orders)) {
    if (specimen_orders$specimen[j] == species_name) {
      disparity_data$order[i] <- specimen_orders$order[j]
      break
    }
  }
}

# Palaeocontinents
disparity_data$palaeocontinent <- NA
for (i in 1:nrow(disparity_data)) {
  species_name <- disparity_data$specimen[i]
  for (j in 1:nrow(specimen_palaeocontinents)) {
    if (specimen_palaeocontinents$specimen[j] == species_name) {
      disparity_data$palaeocontinent[i] <- specimen_palaeocontinents$palaeocontinent[j]
      break
    }
  }
}

# Create factor variables with the correct order
disparity_data$time_period <- factor(disparity_data$time_period, 
                                   levels = c("Stage_3", "Stage_4", "Wuliuan", 
                                             "Drumian", "Guzhangian", "Furongian"))

disparity_data$order <- factor(disparity_data$order, 
                             levels = orders)

disparity_data$palaeocontinent <- factor(disparity_data$palaeocontinent, 
                                       levels = palaeocontinents)

# ----------------------------------------
# Calculate mean SoV for each group
# ----------------------------------------
# Time periods
time_sov_means <- aggregate(sov ~ time_period, data = disparity_data, FUN = mean, na.rm = TRUE)
time_sov_sds <- aggregate(sov ~ time_period, data = disparity_data, FUN = sd, na.rm = TRUE)
time_sov_counts <- aggregate(sov ~ time_period, data = disparity_data, FUN = function(x) sum(!is.na(x)))

# Orders
order_sov_means <- aggregate(sov ~ order, data = disparity_data, FUN = mean, na.rm = TRUE)
order_sov_sds <- aggregate(sov ~ order, data = disparity_data, FUN = sd, na.rm = TRUE)
order_sov_counts <- aggregate(sov ~ order, data = disparity_data, FUN = function(x) sum(!is.na(x)))

# Palaeocontinents
palaeocontinent_sov_means <- aggregate(sov ~ palaeocontinent, data = disparity_data, FUN = mean, na.rm = TRUE)
palaeocontinent_sov_sds <- aggregate(sov ~ palaeocontinent, data = disparity_data, FUN = sd, na.rm = TRUE)
palaeocontinent_sov_counts <- aggregate(sov ~ palaeocontinent, data = disparity_data, FUN = function(x) sum(!is.na(x)))

# ----------------------------------------
# Create visualisation for Time Periods
# ----------------------------------------
create_time_plot <- function() {
  # Set up the plotting parameters
  pdf("time_disparity_rate_combined.pdf", width = 10, height = 5)
  
  # Panel 1: SoV Boxplot with enhanced styling and secondary rate y-axis
  par(mar = c(4, 4.5, 3, 4.5))
  
  # Calculate y-axis limits for SoV
  sov_ylim <- c(min(disparity_data$sov, na.rm = TRUE) * 0.9,
               max(disparity_data$sov, na.rm = TRUE) * 1.2)
  
  # Empty plot for boxplot
  plot(NA, NA, 
       xlim = c(0.5, length(levels(disparity_data$time_period)) + 0.5),
       ylim = sov_ylim,
       xlab = "",
       ylab = "Sum of Variances (SoV)",
       main = "Trilobite Disparity and Evolutionary Rates Through Time",
       xaxt = "n")
  
  # Add x-axis labels
  axis(1, at = 1:length(levels(disparity_data$time_period)), 
       labels = levels(disparity_data$time_period), 
       las = 1)
  
  # Add grid
  grid(nx = NA, ny = NULL, lty = 2, col = "lightgray")
  
  # Draw enhanced SoV boxplots with time colors and notches
  boxplot(sov ~ time_period, data = disparity_data,
          add = TRUE,
          notch = TRUE,  # Add statistical notches
          outline = FALSE, 
          col = time_colors[levels(disparity_data$time_period)],
          border = "black",
          lwd = 0.5,
          xaxt = "n")
  
  # Add jittered points for better distribution visualisation
  for(i in 1:length(levels(disparity_data$time_period))) {
    time_pd <- levels(disparity_data$time_period)[i]
    subset_data <- disparity_data[disparity_data$time_period == time_pd, ]
    points(jitter(rep(i, nrow(subset_data)), amount = 0.2),
           subset_data$sov,
           pch = 20,
           col = adjustcolor("black", alpha.f = 0.3),
           cex = 0.8)
  }
  
  # Add sample sizes
  for(i in 1:nrow(time_sov_counts)) {
    time_pd <- time_sov_counts$time_period[i]
    pos <- which(levels(disparity_data$time_period) == time_pd)
    text(pos, sov_ylim[2] * 0.95, 
         paste0("n=", time_sov_counts$sov[i]), 
         cex = 0.8)
  }
  
  # Add right-side secondary y-axis
  par(new = TRUE)
  
  # Calculate rate line y-axis limits
  rate_ylim <- c(min(time_summary$mean_rate - time_summary$sd_rate) * 0.9,
                max(time_summary$mean_rate + time_summary$sd_rate) * 1.1)
  
  # Plot evolutionary rate line
  plot(NA, NA,
       xlim = c(0.5, length(levels(disparity_data$time_period)) + 0.5),
       ylim = rate_ylim,
       xlab = "",
       ylab = "",
       xaxt = "n",
       yaxt = "n")
  
  # Add right y-axis
  axis(4, col = "blue", col.axis = "blue")
  mtext("Evolutionary Rate", side = 4, col = "blue", line = 2.5)
  
  # Draw evolutionary rate line
  lines(1:nrow(time_summary), time_summary$mean_rate, 
        col = "blue", lwd = 2)
  points(1:nrow(time_summary), time_summary$mean_rate, 
         pch = 21, bg = "white", col = "blue", cex = 1.5)
  
  # Add error bars for rates
  for(i in 1:nrow(time_summary)) {
    arrows(i, time_summary$mean_rate[i] - time_summary$sd_rate[i],
           i, time_summary$mean_rate[i] + time_summary$sd_rate[i],
           length = 0.05, angle = 90, code = 3, col = "blue")
  }
  
  # Close the PDF device
  dev.off()
  
  # Return file path
  return("time_disparity_rate_combined.pdf")
}

# ----------------------------------------
# Create visualisation for Orders
# ----------------------------------------
create_order_plot <- function() {
  # Set up the plotting parameters
  pdf("order_disparity_rate_combined.pdf", width = 10, height = 5)
  
  # Panel 1: SoV Boxplot with enhanced styling and secondary rate y-axis
  par(mar = c(4, 4.5, 3, 4.5))
  
  # Calculate y-axis limits for SoV
  sov_ylim <- c(min(disparity_data$sov, na.rm = TRUE) * 0.9,
               max(disparity_data$sov, na.rm = TRUE) * 1.2)
  
  # Empty plot for boxplot
  plot(NA, NA, 
       xlim = c(0.5, length(levels(disparity_data$order)) + 0.5),
       ylim = sov_ylim,
       xlab = "",
       ylab = "Sum of Variances (SoV)",
       main = "Trilobite Disparity and Evolutionary Rates Across Orders",
       xaxt = "n")
  
  # Add x-axis labels
  axis(1, at = 1:length(levels(disparity_data$order)), 
       labels = levels(disparity_data$order), 
       las = 1,
       cex.axis = 0.7) 
  
  # Add grid
  grid(nx = NA, ny = NULL, lty = 2, col = "lightgray")
  
  # Draw enhanced SoV boxplots with order colors and notches
  boxplot(sov ~ order, data = disparity_data,
          add = TRUE,
          notch = TRUE,  # Add statistical notches
          outline = FALSE,
          col = order_colors[levels(disparity_data$order)],
          border = "black",
          lwd = 0.5,
          xaxt = "n")
  
  # Add jittered points for better distribution visualisation
  for(i in 1:length(levels(disparity_data$order))) {
    order_val <- levels(disparity_data$order)[i]
    subset_data <- disparity_data[disparity_data$order == order_val, ]
    points(jitter(rep(i, nrow(subset_data)), amount = 0.2),
           subset_data$sov,
           pch = 20,
           col = adjustcolor("black", alpha.f = 0.3),
           cex = 0.8)
  }
  
  # Add sample sizes
  for(i in 1:nrow(order_sov_counts)) {
    order_val <- order_sov_counts$order[i]
    pos <- which(levels(disparity_data$order) == order_val)
    text(pos, sov_ylim[2] * 0.95, 
         paste0("n=", order_sov_counts$sov[i]), 
         cex = 0.8)
  }
  
  # Add right-side secondary y-axis
  par(new = TRUE)
  
  # Match order_summary order with the factor levels from disparity_data
  order_summary$order_factor <- factor(order_summary$order, levels = levels(disparity_data$order))
  order_summary <- order_summary[order(order_summary$order_factor), ]
  
  # Calculate rate line y-axis limits
  rate_ylim <- c(min(order_summary$mean_rate - order_summary$sd_rate) * 0.9,
                max(order_summary$mean_rate + order_summary$sd_rate) * 1.1)
  
  # Plot evolutionary rate line
  plot(NA, NA,
       xlim = c(0.5, length(levels(disparity_data$order)) + 0.5),
       ylim = rate_ylim,
       xlab = "",
       ylab = "",
       xaxt = "n",
       yaxt = "n")
  
  # Add right y-axis
  axis(4, col = "blue", col.axis = "blue")
  mtext("Evolutionary Rate", side = 4, col = "blue", line = 2.5)
  
  # Draw evolutionary rate line
  lines(1:nrow(order_summary), order_summary$mean_rate, 
        col = "blue", lwd = 2)
  points(1:nrow(order_summary), order_summary$mean_rate, 
         pch = 21, bg = "white", col = "blue", cex = 1.5)
  
  # Add error bars for rates
  for(i in 1:nrow(order_summary)) {
    arrows(i, order_summary$mean_rate[i] - order_summary$sd_rate[i],
           i, order_summary$mean_rate[i] + order_summary$sd_rate[i],
           length = 0.05, angle = 90, code = 3, col = "blue")
  }
  
  # Close the PDF device
  dev.off()
  
  # Return file path
  return("order_disparity_rate_combined.pdf")
}

# ----------------------------------------
# Create visualisation for Palaeocontinents
# ----------------------------------------
create_palaeocontinent_plot <- function() {
  # Set up the plotting parameters
  pdf("palaeocontinent_disparity_rate_combined.pdf", width = 10, height = 5)
  
  # Panel 1: SoV Boxplot with enhanced styling and secondary rate y-axis
  par(mar = c(4, 4.5, 3, 4.5))
  
  # Calculate y-axis limits for SoV
  sov_ylim <- c(min(disparity_data$sov, na.rm = TRUE) * 0.9,
               max(disparity_data$sov, na.rm = TRUE) * 1.2)
  
  # Empty plot for boxplot
  plot(NA, NA, 
       xlim = c(0.5, length(levels(disparity_data$palaeocontinent)) + 0.5),
       ylim = sov_ylim,
       xlab = "",
       ylab = "Sum of Variances (SoV)",
       main = "Trilobite Disparity and Evolutionary Rates Across Palaeocontinents",
       xaxt = "n")
  
  # Add x-axis labels
  axis(1, at = 1:length(levels(disparity_data$palaeocontinent)), 
       labels = levels(disparity_data$palaeocontinent), 
       las = 1)
  
  # Add grid
  grid(nx = NA, ny = NULL, lty = 2, col = "lightgray")
  
  # Draw enhanced SoV boxplots with palaeocontinent colors and gradient effects
  # First draw the filled boxplot with alpha for transparency
  boxplot(sov ~ palaeocontinent, data = disparity_data,
          add = TRUE,
          notch = TRUE,
          outline = FALSE,
          col = adjustcolor(palaeocontinent_colors[levels(disparity_data$palaeocontinent)], alpha.f = 0.6),
          border = NA,
          xaxt = "n")
  
  # Then add multiple borders with decreasing opacity for gradient effect
  # Border 1 (thickest, most opaque)
  boxplot(sov ~ palaeocontinent, data = disparity_data,
          add = TRUE,
          notch = TRUE,
          outline = FALSE,
          col = NA,
          border = palaeocontinent_colors[levels(disparity_data$palaeocontinent)],
          lwd = 1.2,
          xaxt = "n")
  
  # Border 2 (medium thickness)
  boxplot(sov ~ palaeocontinent, data = disparity_data,
          add = TRUE,
          notch = TRUE,
          outline = FALSE,
          col = NA,
          border = adjustcolor(palaeocontinent_colors[levels(disparity_data$palaeocontinent)], alpha.f = 0.6),
          lwd = 0.8,
          xaxt = "n")
  
  # Border 3 (thinnest, least opaque)
  boxplot(sov ~ palaeocontinent, data = disparity_data,
          add = TRUE,
          notch = TRUE,
          outline = FALSE,
          col = NA,
          border = adjustcolor(palaeocontinent_colors[levels(disparity_data$palaeocontinent)], alpha.f = 0.4),
          lwd = 0.4,
          xaxt = "n")
  
  # Add jittered points for better distribution visualisation
  for(i in 1:length(levels(disparity_data$palaeocontinent))) {
    continent <- levels(disparity_data$palaeocontinent)[i]
    subset_data <- disparity_data[disparity_data$palaeocontinent == continent, ]
    points(jitter(rep(i, nrow(subset_data)), amount = 0.2),
           subset_data$sov,
           pch = 20,
           col = adjustcolor("black", alpha.f = 0.3),
           cex = 0.8)
  }
  
  # Add sample sizes
  for(i in 1:nrow(palaeocontinent_sov_counts)) {
    continent <- palaeocontinent_sov_counts$palaeocontinent[i]
    pos <- which(levels(disparity_data$palaeocontinent) == continent)
    text(pos, sov_ylim[2] * 0.95, 
         paste0("n=", palaeocontinent_sov_counts$sov[i]), 
         cex = 0.8)
  }
  
  # Add right-side secondary y-axis
  par(new = TRUE)
  
  # Match palaeocontinent_summary order with the factor levels from disparity_data
  palaeocontinent_summary$palaeocontinent_factor <- factor(palaeocontinent_summary$palaeocontinent, 
                                                         levels = levels(disparity_data$palaeocontinent))
  palaeocontinent_summary <- palaeocontinent_summary[order(palaeocontinent_summary$palaeocontinent_factor), ]
  
  # Calculate rate line y-axis limits
  rate_ylim <- c(min(palaeocontinent_summary$mean_rate - palaeocontinent_summary$sd_rate) * 0.9,
                max(palaeocontinent_summary$mean_rate + palaeocontinent_summary$sd_rate) * 1.1)
  
  # Plot evolutionary rate line
  plot(NA, NA,
       xlim = c(0.5, length(levels(disparity_data$palaeocontinent)) + 0.5),
       ylim = rate_ylim,
       xlab = "",
       ylab = "",
       xaxt = "n",
       yaxt = "n")
  
  # Add right y-axis
  axis(4, col = "blue", col.axis = "blue")
  mtext("Evolutionary Rate", side = 4, col = "blue", line = 2.5)
  
  # Draw evolutionary rate line
  lines(1:nrow(palaeocontinent_summary), palaeocontinent_summary$mean_rate, 
        col = "blue", lwd = 2)
  points(1:nrow(palaeocontinent_summary), palaeocontinent_summary$mean_rate, 
         pch = 21, bg = "white", col = "blue", cex = 1.5)
  
  # Add error bars for rates
  for(i in 1:nrow(palaeocontinent_summary)) {
    arrows(i, palaeocontinent_summary$mean_rate[i] - palaeocontinent_summary$sd_rate[i],
           i, palaeocontinent_summary$mean_rate[i] + palaeocontinent_summary$sd_rate[i],
           length = 0.05, angle = 90, code = 3, col = "blue")
  }
  
  # Close the PDF device
  dev.off()
  
  # Return file path
  return("palaeocontinent_disparity_rate_combined.pdf")
}

# ----------------------------------------
# Generate the plots
# ----------------------------------------
# Time period plots
time_plot_path <- create_time_plot()
cat("Time period plot saved to:", time_plot_path, "\n")

# Order plots
order_plot_path <- create_order_plot()
cat("Order plot saved to:", order_plot_path, "\n")

# Palaeocontinent plots
palaeocontinent_plot_path <- create_palaeocontinent_plot()
cat("Palaeocontinent plot saved to:", palaeocontinent_plot_path, "\n")

# Save the data for further analysis
save(time_summary, order_summary, palaeocontinent_summary, 
     disparity_data, time_sov_means, order_sov_means, palaeocontinent_sov_means,
     file = "combined_rate_disparity_data.RData")


# -------------------------------------------------------------------------------------------------
# 3.1.7 SoV + Rate Plot (Only Combined Palaeocontinents + Time)
# -------------------------------------------------------------------------------------------------
# Set working directory
setwd("~/Desktop/dissertation_project_trilobite")

# Load required packages
library(ape)
library(phytools)
library(ggplot2)
library(gridExtra)

# ----------------------------------------
# Load previously saved data
# ----------------------------------------
# Load the rate tree data
load("phylogenetic_tree_rates_palaeocontinents_orders.RData")

# Load time mapping
load("specimen_time_mapping.RData")

# Load order mapping
load("specimen_order_mapping.RData")

# Load palaeocontinent mapping
load("specimen_palaeocontinent_mapping.RData")

# Load PCA results
load("pca_analysis_results.RData") 

# ----------------------------------------
# Define time intervals with their boundaries
# ----------------------------------------
time_intervals <- list(
  "Stage_3" = c(521, 514),       # From Stage 3 start to end
  "Stage_4" = c(514, 509),       # From Stage 4 start to end
  "Wuliuan" = c(509, 504.5),     # From Wuliuan start to end
  "Drumian" = c(504.5, 500.5),   # From Drumian start to end
  "Guzhangian" = c(500.5, 497),  # From Guzhangian start to end
  "Furongian" = c(497, 485)      # From Paibian start to Stage 10 end (combined Furongian)
)

# ----------------------------------------
# Extract tip rates from the rate_tree
# ----------------------------------------
tip_rates <- rate_tree$edge.length[rate_tree$edge[,2] <= length(rate_tree$tip.label)]
tip_indices <- rate_tree$edge[,2][rate_tree$edge[,2] <= length(rate_tree$tip.label)]
tip_names <- rate_tree$tip.label[tip_indices]

# Create data frame with tip names and rates
tip_data <- data.frame(
  specimen = tip_names,
  rate = tip_rates,
  stringsAsFactors = FALSE
)

# ----------------------------------------
# Map Time Periods, Orders, and Palaeocontinents to tip data
# ----------------------------------------
# Time Periods
tip_data$time_period <- NA
for (i in 1:nrow(tip_data)) {
  species_name <- tip_data$specimen[i]
  for (j in 1:nrow(specimen_times)) {
    if (specimen_times$specimen[j] == species_name) {
      tip_data$time_period[i] <- specimen_times$time_period[j]
      break
    }
  }
}

# Orders
tip_data$order <- NA
for (i in 1:nrow(tip_data)) {
  species_name <- tip_data$specimen[i]
  for (j in 1:nrow(specimen_orders)) {
    if (specimen_orders$specimen[j] == species_name) {
      tip_data$order[i] <- specimen_orders$order[j]
      break
    }
  }
}

# Palaeocontinents
tip_data$palaeocontinent <- NA
for (i in 1:nrow(tip_data)) {
  species_name <- tip_data$specimen[i]
  for (j in 1:nrow(specimen_palaeocontinents)) {
    if (specimen_palaeocontinents$specimen[j] == species_name) {
      tip_data$palaeocontinent[i] <- specimen_palaeocontinents$palaeocontinent[j]
      break
    }
  }
}

# ----------------------------------------
# Calculate Sum of Variances (SoV) from PCA scores
# ----------------------------------------
# Function to calculate SoV for each specimen
calculate_specimen_sov <- function(data) {
  n_specimens <- nrow(data)
  sov <- numeric(n_specimens)
  
  for(i in 1:n_specimens) {
    specimen_data <- as.numeric(data[i,])
    sov[i] <- var(specimen_data, na.rm = TRUE)
  }
  
  return(sov)
}

# Calculate SoV from PCA scores
pca_scores <- as.matrix(pca_results$scores)
sov_by_specimen <- calculate_specimen_sov(pca_scores)

# Create data frame with SoV for each specimen
disparity_data <- data.frame(
  specimen = rownames(pca_scores),
  sov = sov_by_specimen
)

# Map time periods, orders, and palaeocontinents to disparity data
# Time Periods
disparity_data$time_period <- NA
for (i in 1:nrow(disparity_data)) {
  species_name <- disparity_data$specimen[i]
  for (j in 1:nrow(specimen_times)) {
    if (specimen_times$specimen[j] == species_name) {
      disparity_data$time_period[i] <- specimen_times$time_period[j]
      break
    }
  }
}

# Orders
disparity_data$order <- NA
for (i in 1:nrow(disparity_data)) {
  species_name <- disparity_data$specimen[i]
  for (j in 1:nrow(specimen_orders)) {
    if (specimen_orders$specimen[j] == species_name) {
      disparity_data$order[i] <- specimen_orders$order[j]
      break
    }
  }
}

# Palaeocontinents
disparity_data$palaeocontinent <- NA
for (i in 1:nrow(disparity_data)) {
  species_name <- disparity_data$specimen[i]
  for (j in 1:nrow(specimen_palaeocontinents)) {
    if (specimen_palaeocontinents$specimen[j] == species_name) {
      disparity_data$palaeocontinent[i] <- specimen_palaeocontinents$palaeocontinent[j]
      break
    }
  }
}

# Create factor variables with the correct order
disparity_data$time_period <- factor(disparity_data$time_period, 
                                   levels = c("Stage_3", "Stage_4", "Wuliuan", 
                                             "Drumian", "Guzhangian", "Furongian"))

# ----------------------------------------
# Define broader time periods
# ----------------------------------------
# Define broader time periods grouping
broader_time_periods <- list(
  "Series2" = c("Stage_3", "Stage_4"),
  "Miaolingian" = c("Wuliuan", "Drumian", "Guzhangian"),
  "Furongian" = c("Furongian")
)

# Add broader time period to disparity_data
disparity_data$broader_time_period <- NA
for (i in 1:nrow(disparity_data)) {
  if(!is.na(disparity_data$time_period[i])) {
    for(bp in names(broader_time_periods)) {
      if(disparity_data$time_period[i] %in% broader_time_periods[[bp]]) {
        disparity_data$broader_time_period[i] <- bp
        break
      }
    }
  }
}

# Add broader time period to tip_data
tip_data$broader_time_period <- NA
for (i in 1:nrow(tip_data)) {
  if(!is.na(tip_data$time_period[i])) {
    for(bp in names(broader_time_periods)) {
      if(tip_data$time_period[i] %in% broader_time_periods[[bp]]) {
        tip_data$broader_time_period[i] <- bp
        break
      }
    }
  }
}

# ----------------------------------------
# Calculate metrics for each palaeocontinent-time combination
# ----------------------------------------
# Define the palaeocontinents to analyze
palaeocontinents_to_analyze <- c("CHINA", "GONDWANA", "LAURENTIA")

# Create a function to calculate metrics for each palaeocontinent-time combination
calculate_metrics_by_continent_time <- function() {
  # Create a data frame to store the results
  results <- data.frame(
    palaeocontinent = character(),
    time_period = character(),
    disparity_mean = numeric(),
    disparity_sd = numeric(),
    disparity_n = numeric(),
    rate_mean = numeric(),
    rate_sd = numeric(),
    rate_n = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Calculate metrics for each palaeocontinent-time combination
  for(continent in palaeocontinents_to_analyze) {
    for(period in names(broader_time_periods)) {
      # Calculate disparity metrics
      disparity_subset <- disparity_data[disparity_data$palaeocontinent == continent & 
                                       disparity_data$broader_time_period == period, ]
      
      if(nrow(disparity_subset) > 0) {
        disparity_mean <- mean(disparity_subset$sov, na.rm = TRUE)
        disparity_sd <- sd(disparity_subset$sov, na.rm = TRUE)
        disparity_n <- nrow(disparity_subset)
      } else {
        disparity_mean <- NA
        disparity_sd <- NA
        disparity_n <- 0
      }
      
      # Calculate rate metrics
      rate_subset <- tip_data[tip_data$palaeocontinent == continent & 
                            tip_data$broader_time_period == period, ]
      
      if(nrow(rate_subset) > 0) {
        rate_mean <- mean(rate_subset$rate, na.rm = TRUE)
        rate_sd <- sd(rate_subset$rate, na.rm = TRUE)
        rate_n <- nrow(rate_subset)
      } else {
        rate_mean <- NA
        rate_sd <- NA
        rate_n <- 0
      }
      
      # Add to results
      results <- rbind(results, data.frame(
        palaeocontinent = continent,
        time_period = period,
        disparity_mean = disparity_mean,
        disparity_sd = disparity_sd,
        disparity_n = disparity_n,
        rate_mean = rate_mean,
        rate_sd = rate_sd,
        rate_n = rate_n,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  return(results)
}

# Calculate the metrics
continent_time_metrics <- calculate_metrics_by_continent_time()

# ----------------------------------------
# Create the complex combination plot with horizontal layout
# ----------------------------------------
continent_time_sov_rate <- function() {
  # Set plotting parameters
  pdf("continent_time_sov_rate_plot.pdf", width = 12, height = 10)
  
  # Create 3x1 layout (3 palaeocontinents in rows)
  layout(matrix(1:3, nrow = 3, ncol = 1))
  
  # Increase margins to accommodate labels
  par(oma = c(1, 1, 2, 1))
  
  # Calculate disparity and rate range including error bars
  disparity_min <- min(continent_time_metrics$disparity_mean - continent_time_metrics$disparity_sd, na.rm = TRUE)
  disparity_max <- max(continent_time_metrics$disparity_mean + continent_time_metrics$disparity_sd, na.rm = TRUE)
  disparity_range <- c(disparity_min * 0.9, disparity_max * 1.1)
  
  rate_min <- min(continent_time_metrics$rate_mean - continent_time_metrics$rate_sd, na.rm = TRUE)
  rate_max <- max(continent_time_metrics$rate_mean + continent_time_metrics$rate_sd, na.rm = TRUE)
  rate_range <- c(rate_min * 0.9, rate_max * 1.1)
  
  # Define time periods and positions
  time_periods <- c("Series2", "Miaolingian", "Furongian")
  time_positions <- c(1, 2, 3)
  
  # Define colors
  continent_colors <- c(
    "CHINA" = "#D73027",      # Deep red
    "GONDWANA" = "#91CF60",   # Green
    "LAURENTIA" = "#FEE090"   # Golden yellow
  )
  
  time_colors <- c(
    "Series2" = "#F8D568",     # Light orange
    "Miaolingian" = "#A2CD5A", # Light green
    "Furongian" = "#87CEEB"    # Light blue
  )
  
  # Create plots for each palaeocontinent
  for(i in 1:length(palaeocontinents_to_analyze)) {
    continent <- palaeocontinents_to_analyze[i]
    
    # Set margins
    if(i == 1) {
      # First row - include title
      par(mar = c(1, 5.5, 3, 4.5))
    } else if(i == length(palaeocontinents_to_analyze)) {
      # Last row - include x-axis labels
      par(mar = c(4, 5.5, 1, 4.5))
    } else {
      # Middle rows
      par(mar = c(1, 5.5, 1, 4.5))
    }
    
    # Create empty plot
    plot(NA, NA,
         xlim = c(0.7, 3.3),
         ylim = disparity_range,
         xlab = "",
         ylab = "Disparity (SoV)",
         main = ifelse(i == 1, "Trilobite Disparity and Evolutionary Rates Across Time and Palaeocontinents", ""),
         xaxt = "n",
         cex.lab = 1.3,
         cex.axis = 1.2,
         cex.main = 2)
    
    # Add time period labels on x-axis for the last row only
    if(i == length(palaeocontinents_to_analyze)) {
      axis(1, at = time_positions, labels = time_periods, cex.axis = 2.5, font = 2)
    }
    
    # Add grid lines
    grid(nx = NA, ny = NULL, lty = 2, col = "lightgray")
    
    # Add background color for each time period section
    for(j in 1:length(time_periods)) {
      rect(j - 0.5, par("usr")[3],
           j + 0.5, par("usr")[4],
           col = adjustcolor(time_colors[time_periods[j]], alpha.f = 0.3),
           border = NA)
    }
    
    # Add vertical lines between time periods
    for(j in 1:2) {
      abline(v = j + 0.5, col = "gray50", lty = 2)
    }
    
    # Get data for current palaeocontinent
    continent_data <- continent_time_metrics[continent_time_metrics$palaeocontinent == continent, ]
    
    # Add palaeocontinent label on left side
    mtext(continent, side = 2, line = 4, cex = 1.6, font = 2)
    
    # Fixed heights for median and mean labels
    median_disparity_height <- 0.0028   # Lower height for median labels
    mean_rate_height <- 0.0031          # Higher height for mean labels
    
    # Add disparity points, error bars, and labels for each time period
    for(j in 1:length(time_periods)) {
      period <- time_periods[j]
      period_data <- continent_data[continent_data$time_period == period, ]
      
      if(nrow(period_data) > 0 && !is.na(period_data$disparity_mean)) {
        # Add disparity point
        points(time_positions[j], period_data$disparity_mean,
               pch = 22, 
               bg = continent_colors[continent],
               col = "black",
               lwd = 2,
               cex = 2.5)
        
        # Add disparity error bar
        if(!is.na(period_data$disparity_sd)) {
          arrows(time_positions[j], period_data$disparity_mean - period_data$disparity_sd,
                 time_positions[j], period_data$disparity_mean + period_data$disparity_sd,
                 length = 0.1, angle = 90, code = 3, lwd = 2)
        }
        # Add sample size with larger font
        if(period_data$disparity_n > 0) {
          text(time_positions[j], period_data$disparity_mean,
               paste0("n=", period_data$disparity_n),
               pos = 4, offset = 1, cex = 1.6, font = 2)
        }
        # Add mean disparity label 
        text(time_positions[j], median_disparity_height,
             paste("Mean Disparity=", round(period_data$disparity_mean, 5)),
             cex = 1.5, font = 1, adj = 0.5, col = "black")
      }
    }
    
    # Add secondary y-axis for evolutionary rate
    par(new = TRUE)
    plot(NA, NA,
         xlim = c(0.7, 3.3),
         ylim = rate_range,
         xlab = "",
         ylab = "",
         xaxt = "n",
         yaxt = "n")
    
    # Add rate y-axis
    axis(4, col = "blue", col.axis = "blue", cex.axis = 1.2)
    mtext("Evolutionary Rate", side = 4, col = "blue", line = 2.8, cex = 1.3)
    
    # Convert the fixed heights to rate scale
    disparity_pos_median <- median_disparity_height
    rate_pos_median <- rate_range[1] + (disparity_pos_median - disparity_range[1]) * 
                (rate_range[2] - rate_range[1]) / (disparity_range[2] - disparity_range[1])
    
    disparity_pos_mean <- mean_rate_height
    rate_pos_mean <- rate_range[1] + (disparity_pos_mean - disparity_range[1]) * 
                (rate_range[2] - rate_range[1]) / (disparity_range[2] - disparity_range[1])
    
    # Add rate points, error bars, and labels
    has_previous_point <- FALSE
    prev_x <- NA
    prev_y <- NA
    
    for(j in 1:length(time_periods)) {
      period <- time_periods[j]
      period_data <- continent_data[continent_data$time_period == period, ]
      
      if(nrow(period_data) > 0 && !is.na(period_data$rate_mean)) {
        # Add rate point
        points(time_positions[j], period_data$rate_mean,
               pch = 21, 
               bg = "white",
               col = "blue",
               lwd = 2,
               cex = 2.5)
        
        # Add rate error bar
        if(!is.na(period_data$rate_sd)) {
          arrows(time_positions[j], period_data$rate_mean - period_data$rate_sd,
                 time_positions[j], period_data$rate_mean + period_data$rate_sd,
                 length = 0.1, angle = 90, code = 3, col = "blue", lwd = 2)
        }
        
        # Connect to previous point if exists
        if(has_previous_point) {
          lines(c(prev_x, time_positions[j]), c(prev_y, period_data$rate_mean), 
                col = "blue", lwd = 2.5)
        }
        
        # Add mean rate label at the specified height
        text(time_positions[j], rate_pos_mean,
             paste("Mean Rate =", round(period_data$rate_mean, 5)),
             cex = 1.5, font = 1, adj = 0.5, col = "blue")
        
        # Update previous point
        has_previous_point <- TRUE
        prev_x <- time_positions[j]
        prev_y <- period_data$rate_mean
      }
    }
  }
  
  # Add legend
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend("bottom", 
         legend = c("Disparity (SoV)", "Evolutionary Rate"),
         pch = c(22, 21),
         pt.bg = c("red", "white"),
         col = c("black", "blue"),
         pt.cex = 2,
         lwd = 2,
         cex = 1.2,
         horiz = TRUE,
         bty = "n")
  
  # Close PDF device
  dev.off()
  
  return("continent_time_sov_rate_plot.pdf")
}

# Generate the horizontal complex plot
continent_time_sov_rate_plot <- continent_time_sov_rate()
cat("Horizontal complex plot saved to:", continent_time_sov_rate_plot, "\n")

# Print a summary of the palaeocontinent-time metrics
cat("\nSummary of Palaeocontinent-Time Metrics:\n")
print(continent_time_metrics)

# Save the data for future reference
save(continent_time_metrics, file = "continent_time_metrics.RData")
cat("\nMetrics saved to: continent_time_metrics.RData\n")


# -------------------------------------------------------------------------------------------------
# 3.1.8 Statistical Analyses Of SoV + Rate
# -------------------------------------------------------------------------------------------------
# Set working directory
setwd("~/Desktop/dissertation_project_trilobite")

# Load required packages
if(!require(ggplot2)) install.packages("ggplot2")
if(!require(dplyr)) install.packages("dplyr")
if(!require(writexl)) install.packages("writexl")
if(!require(openxlsx)) install.packages("openxlsx")
if(!require(car)) install.packages("car")
if(!require(tidyr)) install.packages("tidyr")
if(!require(patchwork)) install.packages("patchwork")
if(!require(corrplot)) install.packages("corrplot")

library(ggplot2)
library(dplyr)
library(writexl)
library(openxlsx)
library(car)
library(tidyr)
library(patchwork)
library(corrplot)

# ----------------------------------------
# Load data
# ----------------------------------------
load("combined_rate_disparity_data.RData")
load("continent_time_metrics.RData")

cat("Data files loaded\n")
cat("List of objects:\n")
print(ls())

# -----------------------------------------------------------------------------
# Part 1: Statistical analysis of detailed categories
# -----------------------------------------------------------------------------
cat("\n=========== Part 1: Statistical analysis of detailed categories ===========\n")

# ----------------------------------------
# 1. ANOVA analysis: Testing the effect of different classifications on disparity
# ----------------------------------------
# Time Period
time_anova <- aov(sov ~ time_period, data = disparity_data)
time_anova_summary <- summary(time_anova)[[1]]
cat("\n1.1 ANOVA analysis of time period effect on disparity:\n")
print(time_anova_summary)

# Order
order_anova <- aov(sov ~ order, data = disparity_data)
order_anova_summary <- summary(order_anova)[[1]]
cat("\n1.2 ANOVA analysis of order effect on disparity:\n")
print(order_anova_summary)

# Palaeocontinent
continent_anova <- aov(sov ~ palaeocontinent, data = disparity_data)
continent_anova_summary <- summary(continent_anova)[[1]]
cat("\n1.3 ANOVA analysis of palaeocontinent effect on disparity:\n")
print(continent_anova_summary)

# Post-hoc tests
if(time_anova_summary[["Pr(>F)"]][1] < 0.05) {
  cat("\nTime period post-hoc test (Tukey HSD):\n")
  print(TukeyHSD(time_anova))
}

if(order_anova_summary[["Pr(>F)"]][1] < 0.05) {
  cat("\nOrder post-hoc test (Tukey HSD):\n")
  print(TukeyHSD(order_anova))
}

if(continent_anova_summary[["Pr(>F)"]][1] < 0.05) {
  cat("\nPalaeocontinent post-hoc test (Tukey HSD):\n")
  print(TukeyHSD(continent_anova))
}

# ----------------------------------------
# 2. Correlation analysis: Examining the relationship between disparity and evolutionary rate
# ----------------------------------------
# Create summary data frames for correlation analysis
# Time period summary
time_summary_corr <- data.frame(
  time_period = time_summary$time_period,
  disparity = time_sov_means$sov,
  rate = time_summary$mean_rate
)

# Order summary
order_summary_corr <- data.frame(
  order = order_summary$order,
  disparity = order_sov_means$sov,
  rate = order_summary$mean_rate
)

# Palaeocontinent summary
continent_summary_corr <- data.frame(
  palaeocontinent = palaeocontinent_summary$palaeocontinent,
  disparity = palaeocontinent_sov_means$sov,
  rate = palaeocontinent_summary$mean_rate
)

# Calculate correlations
cat("\n2.1 Correlation between time period disparity and evolutionary rate:\n")
time_pearson <- cor.test(time_summary_corr$disparity, time_summary_corr$rate, method = "pearson")
time_spearman <- cor.test(time_summary_corr$disparity, time_summary_corr$rate, method = "spearman")
print(time_pearson)
print(time_spearman)

cat("\n2.2 Correlation between order disparity and evolutionary rate:\n")
order_pearson <- cor.test(order_summary_corr$disparity, order_summary_corr$rate, method = "pearson")
order_spearman <- cor.test(order_summary_corr$disparity, order_summary_corr$rate, method = "spearman")
print(order_pearson)
print(order_spearman)

cat("\n2.3 Correlation between palaeocontinent disparity and evolutionary rate:\n")
continent_pearson <- cor.test(continent_summary_corr$disparity, continent_summary_corr$rate, method = "pearson")
continent_spearman <- cor.test(continent_summary_corr$disparity, continent_summary_corr$rate, method = "spearman")
print(continent_pearson)
print(continent_spearman)

# ----------------------------------------
# 3. Regression models: Exploring the predictive relationship between disparity and evolutionary rate
# ----------------------------------------
# Create a combined data frame, unifying disparity and rate sources
combined_data <- data.frame(
  category = c(as.character(time_summary_corr$time_period), 
               as.character(order_summary_corr$order),
               as.character(continent_summary_corr$palaeocontinent)),
  type = c(rep("Time", nrow(time_summary_corr)),
           rep("Order", nrow(order_summary_corr)),
           rep("Continent", nrow(continent_summary_corr))),
  disparity = c(time_summary_corr$disparity,
                order_summary_corr$disparity,
                continent_summary_corr$disparity),
  rate = c(time_summary_corr$rate,
           order_summary_corr$rate,
           continent_summary_corr$rate)
)

# Basic model
basic_model <- lm(rate ~ disparity, data = combined_data)
cat("\n3.1 Basic regression model (using only disparity):\n")
basic_model_summary <- summary(basic_model)
print(basic_model_summary)

# Model considering classification type
type_model <- lm(rate ~ disparity + type, data = combined_data)
cat("\n3.2 Regression model considering classification type:\n")
type_model_summary <- summary(type_model)
print(type_model_summary)

# Interaction model
interaction_model <- lm(rate ~ disparity * type, data = combined_data)
cat("\n3.3 Interaction model of disparity and classification type:\n")
interaction_model_summary <- summary(interaction_model)
print(interaction_model_summary)

# Model comparison
cat("\n3.4 Model comparison (ANOVA):\n")
models_comparison <- anova(basic_model, type_model, interaction_model)
print(models_comparison)

# -----------------------------------------------------------------------------
# Part 2: Broad time period and palaeocontinent interaction analysis
# -----------------------------------------------------------------------------
cat("\n=========== Part 2: Broad time period and palaeocontinent interaction analysis ===========\n")

# ----------------------------------------
# 1. Examine the structure of continent_time_metrics
# ----------------------------------------
cat("\n1. Structure of broad time-continent data:\n")
str(continent_time_metrics)

# ----------------------------------------
# 2. Two-way ANOVA - Analysing effects of palaeocontinent, time period, and their interaction on disparity and evolutionary rate
# ----------------------------------------
# Two-way ANOVA for disparity
disparity_two_way <- aov(disparity_mean ~ palaeocontinent * time_period, data = continent_time_metrics)
cat("\n2.1 Two-way ANOVA of palaeocontinent and time period effects on disparity:\n")
disparity_two_way_summary <- summary(disparity_two_way)
print(disparity_two_way_summary)

# Two-way ANOVA for evolutionary rate
rate_two_way <- aov(rate_mean ~ palaeocontinent * time_period, data = continent_time_metrics)
cat("\n2.2 Two-way ANOVA of palaeocontinent and time period effects on evolutionary rate:\n")
rate_two_way_summary <- summary(rate_two_way)
print(rate_two_way_summary)

# ----------------------------------------
# 3. Analysis of temporal trends within each palaeocontinent
# ----------------------------------------
cat("\n3. Temporal trends analysis within each palaeocontinent:\n")

# Create empty list to store results
continent_trends <- list()

# Analyze each palaeocontinent
for (continent in unique(continent_time_metrics$palaeocontinent)) {
  # Extract data for this palaeocontinent
  subset_data <- continent_time_metrics[continent_time_metrics$palaeocontinent == continent, ]
  
  # If enough data points for linear regression
  if (nrow(subset_data) >= 3) {
    # Temporal trend for disparity
    disparity_trend <- lm(disparity_mean ~ time_period, data = subset_data)
    disparity_trend_summary <- summary(disparity_trend)
    
    # Temporal trend for evolutionary rate
    rate_trend <- lm(rate_mean ~ time_period, data = subset_data)
    rate_trend_summary <- summary(rate_trend)
    
    # Store results
    continent_trends[[continent]] <- list(
      disparity_trend = disparity_trend_summary,
      rate_trend = rate_trend_summary
    )
    
    cat("\nPalaeocontinent:", continent, "\n")
    cat("Temporal trend for disparity:\n")
    print(disparity_trend_summary)
    cat("Temporal trend for evolutionary rate:\n")
    print(rate_trend_summary)
  } else {
    cat("\nPalaeocontinent:", continent, "- insufficient data points for regression analysis (n =", nrow(subset_data), ")\n")
  }
}

# ----------------------------------------
# 4. Analysis of palaeocontinent differences within each time period
# ----------------------------------------
cat("\n4. Analysis of palaeocontinent differences within each time period:\n")

# Create empty list to store results
time_differences <- list()

# Analyze each time period
for (period in unique(continent_time_metrics$time_period)) {
  # Extract data for this time period
  subset_data <- continent_time_metrics[continent_time_metrics$time_period == period, ]
  
  # If enough data points for ANOVA
  if (length(unique(subset_data$palaeocontinent)) >= 2) {
    # Effect of palaeocontinent on disparity
    disparity_diff <- aov(disparity_mean ~ palaeocontinent, data = subset_data)
    disparity_diff_summary <- summary(disparity_diff)
    
    # Effect of palaeocontinent on evolutionary rate
    rate_diff <- aov(rate_mean ~ palaeocontinent, data = subset_data)
    rate_diff_summary <- summary(rate_diff)
    
    # Store results
    time_differences[[period]] <- list(
      disparity_diff = disparity_diff_summary,
      rate_diff = rate_diff_summary
    )
    cat("\nTime period:", period, "\n")
    cat("Effect of palaeocontinent on disparity:\n")
    print(disparity_diff_summary)
    cat("Effect of palaeocontinent on evolutionary rate:\n")
    print(rate_diff_summary)
  } else {
    cat("\nTime period:", period, "- insufficient number of palaeocontinents for ANOVA (n =", 
        length(unique(subset_data$palaeocontinent)), ")\n")
  }
}

# ----------------------------------------
# 5. Relationship between disparity and evolutionary rate - Grouped by palaeocontinent and time period
# ----------------------------------------
cat("\n5. Relationship between disparity and evolutionary rate - grouped by palaeocontinent and time period:\n")

# Create a function to calculate group correlations
calculate_group_correlation <- function(data, group_var) {
  group_levels <- unique(data[[group_var]])
  result <- data.frame(
    Group = character(),
    Pearson_r = numeric(),
    Pearson_p = numeric(),
    Spearman_rho = numeric(),
    Spearman_p = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (level in group_levels) {
    subset_data <- data[data[[group_var]] == level, ]
    
    # If enough data points
    if (nrow(subset_data) >= 3) {
      pearson_test <- cor.test(subset_data$disparity_mean, subset_data$rate_mean, method = "pearson")
      spearman_test <- cor.test(subset_data$disparity_mean, subset_data$rate_mean, method = "spearman")
      
      result <- rbind(result, data.frame(
        Group = level,
        Pearson_r = pearson_test$estimate,
        Pearson_p = pearson_test$p.value,
        Spearman_rho = spearman_test$estimate,
        Spearman_p = spearman_test$p.value
      ))
    }
  }
  
  return(result)
}

# Correlations grouped by palaeocontinent
continent_correlations <- calculate_group_correlation(continent_time_metrics, "palaeocontinent")
cat("\n5.1 Disparity-rate correlations grouped by palaeocontinent:\n")
print(continent_correlations)

# Correlations grouped by time period
time_correlations <- calculate_group_correlation(continent_time_metrics, "time_period")
cat("\n5.2 Disparity-rate correlations grouped by time period:\n")
print(time_correlations)

# -----------------------------------------------------------------------------
# Visualisation Section
# -----------------------------------------------------------------------------
# ----------------------------------------
# 1. Prepare themes and colors
# ----------------------------------------
theme_paleontology <- function() {
  theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = "gray80")
    )
}

# Create theme with bold x-axis labels for boxplots
theme_paleontology_bold_x <- function() {
  theme_paleontology() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold")
    )
}

# Color configurations
time_colors <- c(
  "Stage_3" = "#1B9E77",
  "Stage_4" = "#D95F02",
  "Wuliuan" = "#7570B3",
  "Drumian" = "#E7298A",
  "Guzhangian" = "#66A61E",
  "Furongian" = "#666666"
)

order_colors <- c(
  "EODISCIDA" = "#E41A1C",
  "REDLICHIIDA" = "#4DAF4A",
  "CORYNEXOCHIDA" = "#377EB8",
  "AULACOPLEURIDA" = "#FF7F00",
  "OLENIDA" = "#984EA3", 
  "UNCERTAIN" = "#999999"
)

continent_colors <- c(
  "BALTICA" = "#4575B4",
  "CHINA" = "#D73027",
  "GONDWANA" = "#91CF60",
  "LAURENTIA" = "#FEE090"
)

# ----------------------------------------
# 2. Create correlation scatter plot function
# ----------------------------------------
create_correlation_plot <- function(x_var_data, y_var_data, labels, color_palette, title, legend_title, sample_sizes) {
  # Create correctly matched data frame
  plot_data <- data.frame(
    group = labels,
    disparity = x_var_data,
    rate = y_var_data,
    n = sample_sizes
  )
  
  # Calculate correlation coefficient
  corr_value <- cor(plot_data$disparity, plot_data$rate, method = "pearson")
  corr_test <- cor.test(plot_data$disparity, plot_data$rate, method = "pearson")
  p_value <- corr_test$p.value
  
  # Correlation coefficient label
  corr_label <- paste0("Pearson's r = ", round(corr_value, 3))
  if (p_value < 0.05) {
    if (p_value < 0.01) {
      corr_label <- paste0(corr_label, " **")
    } else {
      corr_label <- paste0(corr_label, " *")
    }
  } else {
    corr_label <- paste0(corr_label, " (ns)")
  }
  
  # Create scatter plot
  p <- ggplot(plot_data, aes(x = disparity, y = rate, color = group)) +
    geom_point(size = 3) +
    geom_smooth(method = "lm", se = TRUE, color = "blue", linetype = "dashed", formula = y ~ x) +
    scale_color_manual(values = color_palette, name = legend_title, labels = function(x) gsub("^a ", "", x)) +
    labs(
      title = title,
      subtitle = corr_label,
      x = "Disparity (SoV)",
      y = "Evolutionary Rate"
    ) +
    theme_paleontology()
  
  # Add custom label positions for specific points
  hjust_values <- rep(-0.2, nrow(plot_data))
  vjust_values <- rep(0.5, nrow(plot_data))
  
  # Set specific point label positions based on title
  if (title == "Time Periods: Disparity vs. Evolutionary Rate") {
    # Guzhangian and Stage_3 labels on the left
    guzhangian_index <- which(plot_data$group == "Guzhangian")
    if (length(guzhangian_index) > 0) {
      hjust_values[guzhangian_index] <- 1.2
    }
    stage3_index <- which(plot_data$group == "Stage_3")
    if (length(stage3_index) > 0) {
      hjust_values[stage3_index] <- 1.2
    }
  } else if (title == "Orders: Disparity vs. Evolutionary Rate") {
    # EODISCIDA label on the left
    eodiscida_index <- which(plot_data$group == "EODISCIDA")
    if (length(eodiscida_index) > 0) {
      hjust_values[eodiscida_index] <- 1.2
    }
  } else if (title == "Palaeocontinents: Disparity vs. Evolutionary Rate") {
    # CHINA label on the left
    china_index <- which(plot_data$group == "CHINA")
    if (length(china_index) > 0) {
      hjust_values[china_index] <- 1.2
    }
  }
  
  # Add labels using custom hjust and vjust
  for (i in 1:nrow(plot_data)) {
    p <- p + annotate("text", 
                      x = plot_data$disparity[i], 
                      y = plot_data$rate[i], 
                      label = plot_data$group[i],
                      hjust = hjust_values[i], 
                      vjust = vjust_values[i], 
                      color = color_palette[as.character(plot_data$group[i])],
                      size = 3)
  }
  
  return(p)
}

# ----------------------------------------
# 3. Create boxplot function
# ----------------------------------------
create_boxplot <- function(data, x_var, y_var, color_var, color_palette, title, legend_title) {
  # Calculate means and sample sizes for each group
  means_df <- aggregate(data[[y_var]], by = list(Category = data[[x_var]]), FUN = mean, na.rm = TRUE)
  names(means_df) <- c("Category", "Mean")
  
  counts_df <- aggregate(data[[y_var]], by = list(Category = data[[x_var]]), FUN = function(x) sum(!is.na(x)))
  names(counts_df) <- c("Category", "Count")
  
  # Fixed height for sample size annotation
  n_height <- 0.0035
  
  # Create the plot
  p <- ggplot(data, aes_string(x = x_var, y = y_var, fill = color_var)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA, notch = TRUE) +
    geom_jitter(width = 0.2, height = 0, alpha = 0.5, size = 1) +
    # Add mean points as bold red filled circles
    stat_summary(fun = mean, geom = "point", shape = 16, size = 3, 
                 color = "red", fill = "red") +
    scale_fill_manual(values = color_palette, name = legend_title, labels = function(x) gsub("^a ", "", x)) +
    labs(
      title = title,
      x = "",
      y = "Sum of Variances (SoV)"
    ) +
    theme_paleontology_bold_x()
  
  # Add sample size annotations at fixed height
  for (i in 1:nrow(counts_df)) {
    category_name <- as.character(counts_df$Category[i])
    category_index <- which(levels(data[[x_var]]) == category_name)
    
    # Add sample size annotation
    p <- p + annotate("text",
                     x = category_index,
                     y = n_height,
                     label = paste("n =", counts_df$Count[i]),
                     hjust = 0.5,
                     color = "black",
                     size = 3)
  }
  
  return(p)
}

# ----------------------------------------
# 4. Generate scatter plots and boxplots
# ----------------------------------------
# Get sample sizes for each group
time_counts <- aggregate(sov ~ time_period, data = disparity_data, FUN = function(x) sum(!is.na(x)))
order_counts <- aggregate(sov ~ order, data = disparity_data, FUN = function(x) sum(!is.na(x)))
continent_counts <- aggregate(sov ~ palaeocontinent, data = disparity_data, FUN = function(x) sum(!is.na(x)))

# Time Period
time_sov_means <- aggregate(sov ~ time_period, data = disparity_data, FUN = mean, na.rm = TRUE)
time_periods <- time_sov_means$time_period
time_sov_values <- time_sov_means$sov
time_rate_values <- time_summary$mean_rate[match(time_periods, time_summary$time_period)]
time_sample_sizes <- time_counts$sov[match(time_periods, time_counts$time_period)]

time_corr_plot <- create_correlation_plot(
  time_sov_values, 
  time_rate_values, 
  time_periods,
  time_colors[as.character(time_periods)],
  "Time Periods: Disparity vs. Evolutionary Rate",
  "Time Period",
  time_sample_sizes
)

# Order
order_sov_means <- aggregate(sov ~ order, data = disparity_data, FUN = mean, na.rm = TRUE)
orders <- order_sov_means$order
order_sov_values <- order_sov_means$sov
order_rate_values <- order_summary$mean_rate[match(orders, order_summary$order)]
order_sample_sizes <- order_counts$sov[match(orders, order_counts$order)]

order_corr_plot <- create_correlation_plot(
  order_sov_values, 
  order_rate_values, 
  orders,
  order_colors[as.character(orders)],
  "Orders: Disparity vs. Evolutionary Rate",
  "Order",
  order_sample_sizes
)

# Palaeocontinent
continent_sov_means <- aggregate(sov ~ palaeocontinent, data = disparity_data, FUN = mean, na.rm = TRUE)
continents <- continent_sov_means$palaeocontinent
continent_sov_values <- continent_sov_means$sov
continent_rate_values <- palaeocontinent_summary$mean_rate[match(continents, palaeocontinent_summary$palaeocontinent)]
continent_sample_sizes <- continent_counts$sov[match(continents, continent_counts$palaeocontinent)]

continent_corr_plot <- create_correlation_plot(
  continent_sov_values, 
  continent_rate_values, 
  continents,
  continent_colors[as.character(continents)],
  "Palaeocontinents: Disparity vs. Evolutionary Rate",
  "Palaeocontinent",
  continent_sample_sizes
)

# Create boxplots
time_box_plot <- create_boxplot(
  disparity_data, "time_period", "sov", "time_period", 
  time_colors, "Disparity Across Time Periods", "Time Period"
)

order_box_plot <- create_boxplot(
  disparity_data, "order", "sov", "order", 
  order_colors, "Disparity Across Orders", "Order"
)

continent_box_plot <- create_boxplot(
  disparity_data, "palaeocontinent", "sov", "palaeocontinent", 
  continent_colors, "Disparity Across Palaeocontinents", "Palaeocontinent"
)

# ----------------------------------------
# 5. Combine plots and save as PDF
# ----------------------------------------
# Combine all plots
combined_all_plot <- (time_corr_plot | order_corr_plot | continent_corr_plot) / 
                    (time_box_plot | order_box_plot | continent_box_plot)

# Combine all scatter plots
combined_corr_plot <- time_corr_plot | order_corr_plot | continent_corr_plot

# Combine all boxplots
combined_box_plot <- time_box_plot | order_box_plot | continent_box_plot

# Save as PDF
ggsave("all_boxplot_correlation_plots.pdf", combined_all_plot, width = 18, height = 12)
ggsave("combined_correlation_plots.pdf", combined_corr_plot, width = 18, height = 6)
ggsave("combined_boxplot_plots.pdf", combined_box_plot, width = 18, height = 6)

cat("\nPlots saved as three PDF files:\n")
cat("1. all_boxplot_correlation_plots.pdf - contains all plots\n")
cat("2. combined_correlation_plots.pdf - contains only correlation scatter plots\n")
cat("3. combined_boxplot_plots.pdf - contains only boxplots\n")

# -----------------------------------------------------------------------------
# Export all results to excel
# -----------------------------------------------------------------------------
cat("\n=========== Exporting results to excel ===========\n")

# ----------------------------------------
# 1. Prepare ANOVA results data frame
# ----------------------------------------
anova_results <- data.frame(
  Test = c("Time Period ANOVA", "Order ANOVA", "Palaeocontinent ANOVA",
           "Disparity Two-Way ANOVA (Palaeocontinent)", 
           "Disparity Two-Way ANOVA (Time Period)",
           "Disparity Two-Way ANOVA (Interaction)",
           "Rate Two-Way ANOVA (Palaeocontinent)", 
           "Rate Two-Way ANOVA (Time Period)",
           "Rate Two-Way ANOVA (Interaction)"),
  
  F_value = c(time_anova_summary$`F value`[1], 
              order_anova_summary$`F value`[1],
              continent_anova_summary$`F value`[1],
              disparity_two_way_summary[[1]]$`F value`[1],
              disparity_two_way_summary[[1]]$`F value`[2],
              disparity_two_way_summary[[1]]$`F value`[3],
              rate_two_way_summary[[1]]$`F value`[1],
              rate_two_way_summary[[1]]$`F value`[2],
              rate_two_way_summary[[1]]$`F value`[3]),
  
  p_value = c(time_anova_summary$`Pr(>F)`[1],
              order_anova_summary$`Pr(>F)`[1],
              continent_anova_summary$`Pr(>F)`[1],
              disparity_two_way_summary[[1]]$`Pr(>F)`[1],
              disparity_two_way_summary[[1]]$`Pr(>F)`[2],
              disparity_two_way_summary[[1]]$`Pr(>F)`[3],
              rate_two_way_summary[[1]]$`Pr(>F)`[1],
              rate_two_way_summary[[1]]$`Pr(>F)`[2],
              rate_two_way_summary[[1]]$`Pr(>F)`[3]),
  
  Significance = c(
    ifelse(time_anova_summary$`Pr(>F)`[1] < 0.001, "***", 
           ifelse(time_anova_summary$`Pr(>F)`[1] < 0.01, "**",
                  ifelse(time_anova_summary$`Pr(>F)`[1] < 0.05, "*", "ns"))),
    ifelse(order_anova_summary$`Pr(>F)`[1] < 0.001, "***", 
           ifelse(order_anova_summary$`Pr(>F)`[1] < 0.01, "**",
                  ifelse(order_anova_summary$`Pr(>F)`[1] < 0.05, "*", "ns"))),
    ifelse(continent_anova_summary$`Pr(>F)`[1] < 0.001, "***", 
           ifelse(continent_anova_summary$`Pr(>F)`[1] < 0.01, "**",
                  ifelse(continent_anova_summary$`Pr(>F)`[1] < 0.05, "*", "ns"))),
    ifelse(disparity_two_way_summary[[1]]$`Pr(>F)`[1] < 0.001, "***", 
           ifelse(disparity_two_way_summary[[1]]$`Pr(>F)`[1] < 0.01, "**",
                  ifelse(disparity_two_way_summary[[1]]$`Pr(>F)`[1] < 0.05, "*", "ns"))),
    ifelse(disparity_two_way_summary[[1]]$`Pr(>F)`[2] < 0.001, "***", 
           ifelse(disparity_two_way_summary[[1]]$`Pr(>F)`[2] < 0.01, "**",
                  ifelse(disparity_two_way_summary[[1]]$`Pr(>F)`[2] < 0.05, "*", "ns"))),
    ifelse(disparity_two_way_summary[[1]]$`Pr(>F)`[3] < 0.001, "***", 
           ifelse(disparity_two_way_summary[[1]]$`Pr(>F)`[3] < 0.01, "**",
                  ifelse(disparity_two_way_summary[[1]]$`Pr(>F)`[3] < 0.05, "*", "ns"))),
    ifelse(rate_two_way_summary[[1]]$`Pr(>F)`[1] < 0.001, "***", 
           ifelse(rate_two_way_summary[[1]]$`Pr(>F)`[1] < 0.01, "**",
                  ifelse(rate_two_way_summary[[1]]$`Pr(>F)`[1] < 0.05, "*", "ns"))),
    ifelse(rate_two_way_summary[[1]]$`Pr(>F)`[2] < 0.001, "***", 
           ifelse(rate_two_way_summary[[1]]$`Pr(>F)`[2] < 0.01, "**",
                  ifelse(rate_two_way_summary[[1]]$`Pr(>F)`[2] < 0.05, "*", "ns"))),
    ifelse(rate_two_way_summary[[1]]$`Pr(>F)`[3] < 0.001, "***", 
           ifelse(rate_two_way_summary[[1]]$`Pr(>F)`[3] < 0.01, "**",
                  ifelse(rate_two_way_summary[[1]]$`Pr(>F)`[3] < 0.05, "*", "ns")))
  ),
  
  Factor = c("Time Period", "Order", "Palaeocontinent",
             "Palaeocontinent", "Time Period", "Palaeocontinent:Time Period",
             "Palaeocontinent", "Time Period", "Palaeocontinent:Time Period"),
  
  Response = c("Disparity", "Disparity", "Disparity",
               "Disparity", "Disparity", "Disparity",
               "Rate", "Rate", "Rate")
)

# ----------------------------------------
# 2. Prepare correlation analysis results
# ----------------------------------------
correlation_results <- data.frame(
  Group = c("Time Period", "Time Period", "Order", "Order", "Palaeocontinent", "Palaeocontinent"),
  Test_type = c("Pearson", "Spearman", "Pearson", "Spearman", "Pearson", "Spearman"),
  Correlation = c(time_pearson$estimate, time_spearman$estimate,
                  order_pearson$estimate, order_spearman$estimate,
                  continent_pearson$estimate, continent_spearman$estimate),
  p_value = c(time_pearson$p.value, time_spearman$p.value,
              order_pearson$p.value, order_spearman$p.value,
              continent_pearson$p.value, continent_spearman$p.value),
  Significance = c(
    ifelse(time_pearson$p.value < 0.001, "***", 
           ifelse(time_pearson$p.value < 0.01, "**",
                  ifelse(time_pearson$p.value < 0.05, "*", "ns"))),
    ifelse(time_spearman$p.value < 0.001, "***", 
           ifelse(time_spearman$p.value < 0.01, "**",
                  ifelse(time_spearman$p.value < 0.05, "*", "ns"))),
    ifelse(order_pearson$p.value < 0.001, "***", 
           ifelse(order_pearson$p.value < 0.01, "**",
                  ifelse(order_pearson$p.value < 0.05, "*", "ns"))),
    ifelse(order_spearman$p.value < 0.001, "***", 
           ifelse(order_spearman$p.value < 0.01, "**",
                  ifelse(order_spearman$p.value < 0.05, "*", "ns"))),
    ifelse(continent_pearson$p.value < 0.001, "***", 
           ifelse(continent_pearson$p.value < 0.01, "**",
                  ifelse(continent_pearson$p.value < 0.05, "*", "ns"))),
    ifelse(continent_spearman$p.value < 0.001, "***", 
           ifelse(continent_spearman$p.value < 0.01, "**",
                  ifelse(continent_spearman$p.value < 0.05, "*", "ns")))
  )
)

# Add grouped correlation results
if(exists("continent_correlations") && exists("time_correlations")) {
  grouped_correlation_results <- rbind(
    # Add group identifier column
    cbind(Analysis_Type = "By Continent", continent_correlations),
    cbind(Analysis_Type = "By Time Period", time_correlations)
  )
}

# ----------------------------------------
# 3. Prepare regression model results
# ----------------------------------------
# Get coefficients and p-values - add error handling
get_model_coef <- function(model, coef_name) {
  if(coef_name %in% names(coef(model))) {
    return(coef(model)[coef_name])
  } else {
    return(NA)
  }
}

get_model_pvalue <- function(model, coef_name) {
  if(coef_name %in% rownames(summary(model)$coefficients)) {
    return(summary(model)$coefficients[coef_name, "Pr(>|t|)"])
  } else {
    return(NA)
  }
}

# Prepare regression model results - using basic model and type model
regression_models_results <- data.frame(
  Model = c("Basic Model", "Model with Type", "Interaction Model"),
  Intercept = c(coef(basic_model)[1], coef(type_model)[1], coef(interaction_model)[1]),
  Disparity_coefficient = c(coef(basic_model)[2], coef(type_model)[2], coef(interaction_model)[2]),
  Disparity_p_value = c(summary(basic_model)$coefficients[2,4],
                        summary(type_model)$coefficients[2,4],
                        summary(interaction_model)$coefficients[2,4]),
  R_squared = c(summary(basic_model)$r.squared,
               summary(type_model)$r.squared,
               summary(interaction_model)$r.squared),
  Adj_R_squared = c(summary(basic_model)$adj.r.squared,
                   summary(type_model)$adj.r.squared,
                   summary(interaction_model)$adj.r.squared),
  F_statistic = c(summary(basic_model)$fstatistic[1],
                 summary(type_model)$fstatistic[1],
                 summary(interaction_model)$fstatistic[1]),
  Model_p_value = c(
    pf(summary(basic_model)$fstatistic[1], 
       summary(basic_model)$fstatistic[2], 
       summary(basic_model)$fstatistic[3], lower.tail = FALSE),
    pf(summary(type_model)$fstatistic[1], 
       summary(type_model)$fstatistic[2], 
       summary(type_model)$fstatistic[3], lower.tail = FALSE),
    pf(summary(interaction_model)$fstatistic[1], 
       summary(interaction_model)$fstatistic[2], 
       summary(interaction_model)$fstatistic[3], lower.tail = FALSE)
  )
)

# ----------------------------------------
# 4. Prepare continent-time interaction data
# ----------------------------------------
# Data frame containing means for all continent-time combinations
continent_time_summary <- continent_time_metrics %>%
  select(palaeocontinent, time_period, disparity_mean, rate_mean) %>%
  arrange(palaeocontinent, time_period)

# ----------------------------------------
# 5. Prepare descriptive statistics data
# ----------------------------------------
# Calculate descriptive statistics for disparity by category
descriptive_stats <- data.frame(
  Category_Type = c(rep("Time Period", length(levels(disparity_data$time_period))),
                   rep("Order", length(levels(disparity_data$order))),
                   rep("Palaeocontinent", length(levels(disparity_data$palaeocontinent)))),
  
  Category = c(levels(disparity_data$time_period),
              levels(disparity_data$order),
              levels(disparity_data$palaeocontinent)),
  
  Count = c(tapply(disparity_data$sov, disparity_data$time_period, function(x) sum(!is.na(x))),
           tapply(disparity_data$sov, disparity_data$order, function(x) sum(!is.na(x))),
           tapply(disparity_data$sov, disparity_data$palaeocontinent, function(x) sum(!is.na(x)))),
  
  Mean_SoV = c(tapply(disparity_data$sov, disparity_data$time_period, mean, na.rm = TRUE),
              tapply(disparity_data$sov, disparity_data$order, mean, na.rm = TRUE),
              tapply(disparity_data$sov, disparity_data$palaeocontinent, mean, na.rm = TRUE)),
  
  Median_SoV = c(tapply(disparity_data$sov, disparity_data$time_period, median, na.rm = TRUE),
                tapply(disparity_data$sov, disparity_data$order, median, na.rm = TRUE),
                tapply(disparity_data$sov, disparity_data$palaeocontinent, median, na.rm = TRUE)),
  
  SD_SoV = c(tapply(disparity_data$sov, disparity_data$time_period, sd, na.rm = TRUE),
            tapply(disparity_data$sov, disparity_data$order, sd, na.rm = TRUE),
            tapply(disparity_data$sov, disparity_data$palaeocontinent, sd, na.rm = TRUE))
)

# ----------------------------------------
# 6. Create excel workbook and add worksheets
# ----------------------------------------
wb <- createWorkbook()

# ANOVA results worksheet
addWorksheet(wb, "ANOVA_Results")
writeData(wb, "ANOVA_Results", anova_results, startRow = 1, startCol = 1)
setColWidths(wb, "ANOVA_Results", cols = 1:6, widths = "auto")

# Correlation results worksheet
addWorksheet(wb, "Correlation_Results")
writeData(wb, "Correlation_Results", correlation_results, startRow = 1, startCol = 1)
if(exists("grouped_correlation_results")) {
  writeData(wb, "Correlation_Results", data.frame(Group = "Group-Specific Correlations"), startRow = nrow(correlation_results) + 3, startCol = 1)
  writeData(wb, "Correlation_Results", grouped_correlation_results, startRow = nrow(correlation_results) + 5, startCol = 1)
}
setColWidths(wb, "Correlation_Results", cols = 1:6, widths = "auto")

# Regression model results worksheet
addWorksheet(wb, "Regression_Models")
writeData(wb, "Regression_Models", regression_models_results, startRow = 1, startCol = 1)
setColWidths(wb, "Regression_Models", cols = 1:8, widths = "auto")

# Continent-time interaction data worksheet
addWorksheet(wb, "Continent_Time_Data")
writeData(wb, "Continent_Time_Data", continent_time_summary, startRow = 1, startCol = 1)
setColWidths(wb, "Continent_Time_Data", cols = 1:4, widths = "auto")

# Descriptive statistics worksheet
addWorksheet(wb, "Descriptive_Statistics")
writeData(wb, "Descriptive_Statistics", descriptive_stats, startRow = 1, startCol = 1)
setColWidths(wb, "Descriptive_Statistics", cols = 1:6, widths = "auto")

# Add worksheet for temporal trends within each palaeocontinent
if(exists("continent_trends")) {
  addWorksheet(wb, "Continent_Time_Trends")
  current_row <- 1

  for (continent in unique(continent_time_metrics$palaeocontinent)) {
    subset_data <- continent_time_metrics[continent_time_metrics$palaeocontinent == continent, ]
    
    if (nrow(subset_data) >= 3) {
      # Write heading
      writeData(wb, "Continent_Time_Trends", data.frame(Heading = paste0("Trends for ", continent)), 
                startRow = current_row, startCol = 1)
      current_row <- current_row + 1
      
      # Write disparity temporal trend results
      disparity_trend <- lm(disparity_mean ~ time_period, data = subset_data)
      disparity_summary <- summary(disparity_trend)
      
      writeData(wb, "Continent_Time_Trends", data.frame(Analysis = "Disparity Trend"), 
                startRow = current_row, startCol = 1)
      current_row <- current_row + 1
      
      # Create coefficient table
      coef_table <- as.data.frame(disparity_summary$coefficients)
      writeData(wb, "Continent_Time_Trends", coef_table, 
                startRow = current_row, startCol = 1)
      current_row <- current_row + nrow(coef_table) + 1
      
      writeData(wb, "Continent_Time_Trends", 
                data.frame(Stat = c("R-squared", "Adjusted R-squared", "F-statistic", "p-value"),
                          Value = c(disparity_summary$r.squared, 
                                   disparity_summary$adj.r.squared,
                                   disparity_summary$fstatistic[1],
                                   pf(disparity_summary$fstatistic[1], 
                                      disparity_summary$fstatistic[2], 
                                      disparity_summary$fstatistic[3], 
                                      lower.tail = FALSE))), 
                startRow = current_row, startCol = 1)
      current_row <- current_row + 5
      
      # Write evolutionary rate temporal trend results
      rate_trend <- lm(rate_mean ~ time_period, data = subset_data)
      rate_summary <- summary(rate_trend)
      
      writeData(wb, "Continent_Time_Trends", data.frame(Analysis = "Rate Trend"), 
                startRow = current_row, startCol = 1)
      current_row <- current_row + 1
      
      # Create coefficient table
      coef_table <- as.data.frame(rate_summary$coefficients)
      writeData(wb, "Continent_Time_Trends", coef_table, 
                startRow = current_row, startCol = 1)
      current_row <- current_row + nrow(coef_table) + 1
      
      writeData(wb, "Continent_Time_Trends", 
                data.frame(Stat = c("R-squared", "Adjusted R-squared", "F-statistic", "p-value"),
                          Value = c(rate_summary$r.squared, 
                                   rate_summary$adj.r.squared,
                                   rate_summary$fstatistic[1],
                                   pf(rate_summary$fstatistic[1], 
                                      rate_summary$fstatistic[2], 
                                      rate_summary$fstatistic[3], 
                                      lower.tail = FALSE))), 
                startRow = current_row, startCol = 1)
      current_row <- current_row + 7
    }
  }

  setColWidths(wb, "Continent_Time_Trends", cols = 1:5, widths = "auto")
}

# Save excel file
saveWorkbook(wb, "trilobite_statistical_results.xlsx", overwrite = TRUE)
cat("\nStatistical results exported to 'trilobite_statistical_results.xlsx'\n")

# ----------------------------------------
# Save R objects for future analysis
# ----------------------------------------
cat("\n=========== Saving R objects ===========\n")

# Collect all important analysis results to save
analysis_results <- list(
  # ANOVA results
  time_anova = time_anova,
  order_anova = order_anova,
  continent_anova = continent_anova,
  disparity_two_way = disparity_two_way,
  rate_two_way = rate_two_way,
  
  # Correlation results
  time_pearson = time_pearson,
  time_spearman = time_spearman,
  order_pearson = order_pearson,
  order_spearman = order_spearman,
  continent_pearson = continent_pearson,
  continent_spearman = continent_spearman,
  
  # Grouped correlations
  continent_correlations = if(exists("continent_correlations")) continent_correlations else NULL,
  time_correlations = if(exists("time_correlations")) time_correlations else NULL,
  
  # Regression models
  basic_model = basic_model,
  type_model = type_model,
  interaction_model = interaction_model,
  
  # Summary data
  anova_results = anova_results,
  correlation_results = correlation_results,
  regression_models_results = regression_models_results,
  
  # Continental trends
  continent_trends = if(exists("continent_trends")) continent_trends else NULL,
  time_differences = if(exists("time_differences")) time_differences else NULL
)

# Save the list of results
save(analysis_results, file = "disparity_rate_statistical_results.RData")
cat("Analysis results saved to 'disparity_rate_statistical_results.RData'\n")


# -------------------------------------------------------------------------------------------------
# 3.2 Phylogenetic Tree
# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------
# 3.2.1 Generate Morphological Tree (Using PCA Scores)
# -------------------------------------------------------------------------------------------------
# Set working directory
setwd("~/Desktop/dissertation_project_trilobite")

# Load required packages
library(ape)
library(phytools)

# Load PCA and tree data
load("pca_analysis_results.RData")
tree <- read.nexus("UCLN_CambrianTrilos.tre")
tree$tip.label[tree$tip.label == "Neocobboldia_chilinica"] <- "Neocobboldia_chinlinica"

# List of species to prune
species_to_prune <- c(
    "Bristolia_anteros", "Calodiscus_lobatus", "Serrodiscus_speciosus",
    "Anabaraspis_splendens", "Ichangia_ichangensis", "Tayanaspis_bulbosus",
    "Yinites_wanshougongensis", "Lonchopygella_megaspina", "Oligometopus_breviceps",
    "Tsinania_canens", "Zacanthopsis_palmeri", "Ajrikina_hunanensis",
    "Golasaphus_momedahensis", "Liostracina_tangwangzhaiensis",
    "Proceratopyge_truncata", "Notchpeakia_taylori", "Akoldinioidia_shanjiangensis", 
    "Asioptychaspis_subglobosa", "Conocoryphe_sulzeri", "Fenghuangella_coniforma",
    "Kuetsingocephalus_qiannanensis", "Mapania_striata", "Monkaspis_quadrata",
    "Onaraspis_rubra", "Parakoldinioidia_akerfeldti", "Prodamesella_punctata"
)

# Prune tree to match PCA dataset
pruned_tree <- drop.tip(tree, species_to_prune)

# Get PCA scores and match with tree tips
pca_scores <- pca_results$scores
matched_scores <- pca_scores[match(pruned_tree$tip.label, rownames(pca_scores)),]

# Calculate ancestral states at internal nodes
anc_states <- apply(matched_scores, 2, function(x) fastAnc(pruned_tree, x))

# Create morphological tree
morph_tree <- pruned_tree
for(i in 1:nrow(morph_tree$edge)) {
    node1 <- morph_tree$edge[i,1]  # Parent node
    node2 <- morph_tree$edge[i,2]  # Child node
    # Get states for both nodes
    if(node1 <= length(pruned_tree$tip.label)) {
        node1_state <- matched_scores[node1,]
    } else {
        node1_state <- anc_states[node1 - length(pruned_tree$tip.label),]
    }
    
    if(node2 <= length(pruned_tree$tip.label)) {
        node2_state <- matched_scores[node2,]
    } else {
        node2_state <- anc_states[node2 - length(pruned_tree$tip.label),]
    }
    
    # Calculate morphological change as Euclidean distance between states
    morph_tree$edge.length[i] <- sqrt(sum((node2_state - node1_state)^2))
}

# Visualise morphological change tree
pdf("morphological_tree.pdf", width=12, height=8)
cols <- colorRampPalette(c("blue","yellow","red"))(100)
plot(morph_tree, type="phylogram", cex=0.6)

# Add colored edge labels
br_cols <- cols[cut(morph_tree$edge.length, breaks=100)]
edgelabels(round(morph_tree$edge.length, 3), bg=br_cols, cex=0.5)

# Add legend
legend("topright",
       legend=round(seq(min(morph_tree$edge.length), 
                       max(morph_tree$edge.length), length.out=5), 3),
       fill=cols[c(1,25,50,75,100)],
       title="Morphological Change",
       cex=0.8)
dev.off()

# Save results
save(morph_tree, file="morphological_tree.RData")


# -------------------------------------------------------------------------------------------------
# 3.2.2 Phylognetic Tree With Evolutionary Rate
# -------------------------------------------------------------------------------------------------
# Set working directory
setwd("~/Desktop/dissertation_project_trilobite")

# Load required packages
library(ape)
library(phytools)

# Load data
load("morphological_tree.RData")
load("specimen_palaeocontinent_mapping.RData")  # Contains specimen_palaeocontinents and palaeocontinent_colors
load("specimen_order_mapping.RData")            # Contains specimen_orders and order_colors
tree <- read.nexus("UCLN_CambrianTrilos.tre")
tree$tip.label[tree$tip.label == "Neocobboldia_chilinica"] <- "Neocobboldia_chinlinica"

# List of species to prune
species_to_prune <- c(
  "Bristolia_anteros", "Calodiscus_lobatus", "Serrodiscus_speciosus",
  "Anabaraspis_splendens", "Ichangia_ichangensis", "Tayanaspis_bulbosus",
  "Yinites_wanshougongensis", "Lonchopygella_megaspina", "Oligometopus_breviceps",
  "Tsinania_canens", "Zacanthopsis_palmeri", "Ajrikina_hunanensis",
  "Golasaphus_momedahensis", "Liostracina_tangwangzhaiensis",
  "Proceratopyge_truncata", "Notchpeakia_taylori", "Akoldinioidia_shanjiangensis", 
  "Asioptychaspis_subglobosa", "Conocoryphe_sulzeri", "Fenghuangella_coniforma",
  "Kuetsingocephalus_qiannanensis", "Mapania_striata", "Monkaspis_quadrata",
  "Onaraspis_rubra", "Parakoldinioidia_akerfeldti", "Prodamesella_punctata"
)

# Get time tree with time intervals
time_tree <- drop.tip(tree, species_to_prune)

# Calculate evolutionary rates
rate_tree <- time_tree
rate_tree$edge.length <- morph_tree$edge.length / time_tree$edge.length

# Define vivid colors for palaeocontinents
palaeocontinent_colors <- c(
    "BALTICA" = "#0066FF",    # Bright blue
    "CHINA" = "#FF1E1E",      # Vivid red
    "GONDWANA" = "#00CC00",   # Bright green
    "LAURENTIA" = "#FFD700"   # Bright gold
)

# Define shapes for taxonomic orders
order_shapes <- c(
    "EODISCIDA" = 18,                 # Diamond plus
    "REDLICHIIDA (Olenellina)" = 17,  # Filled triangle
    "REDLICHIIDA (Redlichiina)" = 16, # Filled circle
    "CORYNEXOCHIDA" = 4,              # Cross
    "AULACOPLEURIDA" = 3,             # Plus
    "OLENIDA" = 15,                   # Filled square
    "UNCERTAIN" = 8                   # Asterisk
)

# Define geological stages with their boundaries
stages <- list(
  list(name="Stage 2", start=529, end=521),
  list(name="Stage 3", start=521, end=514),
  list(name="Stage 4", start=514, end=509),
  list(name="Wuliuan", start=509, end=504.5),
  list(name="Drumian", start=504.5, end=500.5),
  list(name="Guzhangian", start=500.5, end=497),
  list(name="Paibian", start=497, end=494),
  list(name="Jiangshanian", start=494, end=489.5),
  list(name="Stage 10", start=489.5, end=485)
)

# Function to convert geological time to x-coordinate
time_to_x <- function(time_ma) {
  plot_range * (530 - time_ma)/(530 - 480)
}

# Create visualisation
pdf("phylogenetic_tree_rates_palaeocontinents_orders.pdf", width=15, height=10)

# Set up plotting parameters
par(mar=c(5,1,2,8))  # Standard margins

# Create color scheme for rates
rate_vals <- rate_tree$edge.length
n_colors <- 100
cols <- hcl.colors(256, "Temps")
br_cols <- cols[cut(log(rate_tree$edge.length), 256)]

# Get time interval from tree
node_heights <- nodeHeights(time_tree)
max_height <- max(node_heights)

# Calculate time ranges
time_range <- 530 - 485             # Actual tree time range (530 to 485 Ma)
plot_range <- max_height * (50/45)  # Extend plot range to accommodate 480-485 Ma

# Plot tree
plot(time_tree, 
     type="phylogram",
     cex = 0.6,
     edge.width = 6,
     edge.col = br_cols,
     show.tip.label=TRUE,
     label.offset = 0.3,
     x.lim = c(0, plot_range),
     direction = "rightward")

# Create a subset for Redlichiida (Olenellina)
# Known Olenellina species
olenellina_species <- c(
  "Choubertella_spinosa", "Daguinaspis_ambroggii", "Elliptocephala_asaphoides",
  "Fallotaspis_bondoni", "Holmia_kjerulfi", "Judomia_tera",
  "Montezumaspis_parallela", "Nephrolenellus_geniculatus", "Nevadia_weeksi",
  "Olenellus_gilberti", "Peachella_iddingsi", "Wanneria_walcottana"
)

# Add colored shapes at tips for palaeocontinent (color) and order (shape)
last_plot_coords <- get("last_plot.phylo", envir = .PlotPhyloEnv)
tip_coords <- data.frame(
    x = last_plot_coords$xx[1:length(time_tree$tip.label)],
    y = last_plot_coords$yy[1:length(time_tree$tip.label)]
)

# Add colored shapes for each tip based on palaeocontinent and order
for(i in 1:length(time_tree$tip.label)) {
    species <- time_tree$tip.label[i]
    
    # Get palaeocontinent for color
    continent <- specimen_palaeocontinents$palaeocontinent[
        match(species, specimen_palaeocontinents$specimen)
    ]
    
    # Get order for shape
    taxon_order <- specimen_orders$order[specimen_orders$specimen == species]
    
    # Determine shape - special case for Redlichiida (Olenellina)
    shape <- if (!is.na(taxon_order)) {
        if (taxon_order == "REDLICHIIDA" && species %in% olenellina_species) {
            order_shapes["REDLICHIIDA (Olenellina)"]
        } else if (taxon_order == "REDLICHIIDA") {
            order_shapes["REDLICHIIDA (Redlichiina)"]
        } else {
            order_shapes[taxon_order]
        }
    }
    # If both are available, plot the point
    if(!is.na(continent) && !is.na(taxon_order)) {
        points(tip_coords$x[i], tip_coords$y[i],
               pch=shape,                              # Shape based on order with special case for Olenellina
               col=palaeocontinent_colors[continent],  # Color based on palaeocontinent
               cex=1.2,   
               lwd=1.2)   
    }
}

# Create time axis with geological stage boundaries
# Get all unique boundary times from stages
boundary_times <- unique(c(529, sapply(stages, function(s) s$end)))
boundary_times <- sort(boundary_times, decreasing = TRUE) # Sort from oldest to youngest

# Calculate positions for time axis
axis_positions <- sapply(boundary_times, time_to_x)

# Add time axis with geological boundaries
axis(1, 
     at = axis_positions,
     labels = boundary_times,
     las = 1)
title(xlab = "Time (Ma)", line = 2.5)

# Add stage names above the plot with adjusted position
for (stage in stages) {
  if (stage$start >= 480 && stage$end <= 530) {
    mid_time <- (stage$start + stage$end) / 2
    mid_pos <- time_to_x(mid_time)
    mtext(stage$name, side = 3, at = mid_pos, line = -44, cex = 0.8)
  }
}

# Get rate breaks and colors for legend
rate_breaks <- round(seq(min(rate_vals), max(rate_vals), length.out=5), 3)
legend_colors <- cols[seq(1, 256, length.out=length(rate_breaks)-1)]

# Add evolutionary rate legend
legend("topright", 
       legend=paste(rate_breaks[-length(rate_breaks)], "-", rate_breaks[-1]),
       fill=legend_colors,
       title="Evolutionary Rates",
       cex=0.8,
       bty="n",
       inset=c(0.3, 0.05))

# Add palaeocontinent legend
legend("topright",
       legend=names(palaeocontinent_colors),
       col=palaeocontinent_colors,
       pch=15, 
       pt.cex=1.2,
       title="Palaeocontinents",
       cex=0.8,
       bty="n",
       inset=c(0.18, 0.05))

# Add taxonomic order legend
legend("topright",
       legend=names(order_shapes),
       col="black", 
       pch=order_shapes,
       pt.cex=1.2,
       pt.lwd=1.2,
       title="Taxonomic Orders",
       cex=0.8,
       bty="n",
       inset=c(0, 0.05)) 

dev.off()

# Save results
save(rate_tree, file="phylogenetic_tree_rates_palaeocontinents_orders.RData")


# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 4 Variability Analyses
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# ----------------------------------------
# Set up environment
# ----------------------------------------
rm(list = ls())
library(geomorph)

# Set working directory
setwd("~/Desktop/variability_test")

# ----------------------------------------
# Process file names and import TriloMorph
# ----------------------------------------
# Get file list
files <- list.files(pattern = "\\.txt$", full.names = TRUE)
cat("Found", length(files), "text files\n")

# Process file names (Remove extra .txt extensions)
file_names <- basename(files)
clean_names <- gsub("\\.txt\\.txt$", ".txt", file_names)
clean_names <- gsub("\\.txt$", "", clean_names)

# Import TriloMorph functions
source("~/Downloads/TriloMorph_Functions.r")

# ----------------------------------------
# Set landmark template and define groups
# ----------------------------------------
# Set landmark template
nlms <- c(2, 16, 95, 199, 173, 134)

# Define variability test groups
variability_groups <- list(
    Intra_sample = paste0("Dinesus_ida_", 1:5),
    Intra_specific = paste0("Asaphiscus_wheeleri_", 1:5),
    Perturbed_landmark = paste0("Nevadia_weeksi_", 1:5)
)

# ----------------------------------------
# Process specimens
# ----------------------------------------
# Function to process specimens
process_new_specimens <- function(specimen_names) {
    # Read specimens
    new_specimens <- shapRead(specimen_names, subdir = ".")
    
    # Process specimens
    processed_specimens <- shapFix(new_specimens, nlms, lm.scale = TRUE)
    
    # Perform GPA
    gpa_results <- gpagen(processed_specimens, Proj = TRUE, PrinAxes = FALSE)
    
    return(gpa_results)
}

# Process each variability group
var_results <- lapply(variability_groups, process_new_specimens)

# ----------------------------------------
# Calculate Procrustes distances for box plots
# ----------------------------------------
# Function to calculate Procrustes distances
calc_proc_dist <- function(coords) {
    # Calculate mean shape
    mean_shape <- mshape(coords)
    
    # Calculate distances from mean
    n_specimens <- dim(coords)[3]
    pd <- numeric(n_specimens)
    
    for(i in 1:n_specimens) {
        pd[i] <- sqrt(sum((coords[,,i] - mean_shape)^2))
    }
    return(pd)
}

# Calculate distances for each group
proc_distances <- list()
for(group in names(var_results)) {
    proc_distances[[group]] <- calc_proc_dist(var_results[[group]]$coords)
}

# Create data frame for plotting
plot_data <- data.frame(
    Distance = unlist(proc_distances),
    Group = rep(names(proc_distances), sapply(proc_distances, length))
)

# ----------------------------------------
# Create box plot
# ----------------------------------------
# Set up color scheme
box_colors <- c(
    "Intra_sample" = "gold",
    "Intra_specific" = "salmon",
    "Perturbed_landmark" = "lightblue"
)

# Calculate mean, median and SD for each group
means <- tapply(plot_data$Distance, plot_data$Group, mean)
medians <- tapply(plot_data$Distance, plot_data$Group, median)
sds <- tapply(plot_data$Distance, plot_data$Group, sd)

# Format median and SD values for labels
median_labels <- sprintf("Median=%.4f", medians)
sd_labels <- sprintf("SD=%.4f", sds)

# Create box plot
pdf("variability_boxplot.pdf", width = 8, height = 5)
mar = c(5, 4, 4, 2)
boxplot(Distance ~ Group, data = plot_data,
        col = box_colors[levels(as.factor(plot_data$Group))],
        main = "Procrustes Distance Variation Among Groups",
        xlab = "Group",
        ylab = "Procrustes Distance")

# Add mean values to plot (as points)
points(1:length(means), means, pch = 19, col = "black")

# Add median and SD annotations
for(i in 1:length(medians)) {
    # Add median label
    text(i, 0.02, 
         median_labels[i], 
         cex = 0.8)
    
    # Add SD label below the median label
    text(i, 0.025, 
         sd_labels[i], 
         cex = 0.8)
}
dev.off()

# Save results
save(var_results, proc_distances, file = "variability_gpa_analysis_results.RData")

# ----------------------------------------
# Load original PCA results and project variability specimens
# ----------------------------------------
# Load original results
load("pca_analysis_results.RData")

# Find target species in original dataset
target_species <- c("Dinesus_ida", "Asaphiscus_wheeleri", "Nevadia_weeksi")
species_colors <- c(
    "Dinesus_ida" = "gold",
    "Asaphiscus_wheeleri" = "red",
    "Nevadia_weeksi" = "blue"
)

# Find indices of target species in original dataset
original_indices <- list()
for(species in target_species) {
    original_indices[[species]] <- grep(species, rownames(pca_results$scores))
}

# Correct projection function
project_onto_pca <- function(coords, pca_obj) {
    # Get dimensions
    n_landmarks <- dim(coords)[1]
    n_dims <- dim(coords)[2]
    n_specimens <- dim(coords)[3]
    
    # Create matrix to match original data structure
    specimens_matrix <- matrix(0, nrow = n_specimens, 
                             ncol = length(pca_obj$center))
    
    # Fill matrix with coordinates in correct order
    for(i in 1:n_specimens) {
        spec_coords <- coords[,,i]
        specimens_matrix[i,] <- as.vector(t(spec_coords))
    }
    
    # Center using original mean
    centered <- scale(specimens_matrix, center = pca_obj$center, scale = FALSE)
    
    # Project onto PC axes
    projected <- centered %*% pca_obj$rotation
    return(projected)
}

# Project variability specimens
projected_coords <- list()
for(group in names(var_results)) {
    projected_coords[[group]] <- project_onto_pca(
        var_results[[group]]$coords,
        pca_results$pca
    )
}

# ----------------------------------------
# Create PCA plot with variability specimens
# ----------------------------------------
# Calculate variance explained
prop_var <- pca_results$variance$proportion

# Create plot
pdf("pca_with_variability.pdf", width = 10, height = 8)

# Set up basic plot with original specimens (except target species)
non_target_indices <- setdiff(1:nrow(pca_results$scores), 
                            unlist(original_indices))

plot(pca_results$scores[non_target_indices,1], 
     pca_results$scores[non_target_indices,2],
     pch = 19,  # Solid circles for regular specimens
     col = "black",
     xlab = paste0("PC1 (", round(prop_var[1]*100, 1), "%)"),
     ylab = paste0("PC2 (", round(prop_var[2]*100, 1), "%)"),
     main = "PCA Morphospace with Variability Test Specimens",
     asp = 1)  # Keep aspect ratio 1:1

# Add zero lines
abline(h = 0, v = 0, lty = 2, col = "gray80")

# Add original target species with special colors
for(species in target_species) {
    points(pca_results$scores[original_indices[[species]],1],
           pca_results$scores[original_indices[[species]],2],
           pch = 19,  # Solid circles
           col = species_colors[species],
           cex = 1.5)  # Slightly larger
}

# Add variability specimens
group_colors <- c(
    "Intra_sample" = "gold",
    "Intra_specific" = "red",
    "Perturbed_landmark" = "blue"
)

for(group in names(projected_coords)) {
    points(projected_coords[[group]][,1], 
           projected_coords[[group]][,2],
           pch = 4,  # X symbol
           col = group_colors[group],
           lwd = 2,
           cex = 1.2)
}

# Add legend
legend("topright",
       legend = c("Original specimens",
                 "Dinesus_ida (original)",
                 "Asaphiscus_wheeleri (original)",
                 "Nevadia_weeksi (original)",
                 "Intra_sample_variability",
                 "Intra_specific_variability",
                 "Perturbed_landmark_variability"),
       pch = c(19, 19, 19, 19, 4, 4, 4),
       col = c("black", "gold", "red", "blue", "gold", "red", "blue"),
       bty = "n",
       cex = 0.8)

# Add sample size
mtext(paste0("n = ", nrow(pca_results$scores)), 
      side = 3, adj = 1, font = 3)
dev.off()

# Save results
save(projected_coords, file = "variability_pca_results.RData")

# Print coordinate summary for verification
cat("\nOriginal specimen coordinates:\n")
for(species in target_species) {
    cat("\n", species, ":\n")
    print(pca_results$scores[original_indices[[species]], 1:2])
}

cat("\nVariability specimen coordinates:\n")
for(group in names(projected_coords)) {
    cat("\n", group, ":\n")
    print(projected_coords[[group]][, 1:2])
}

# -----------------------------------------------------------------------------
# Complete statistical analysis for variability test
# -----------------------------------------------------------------------------
# Load required packages
library(geomorph)
library(vegan)
library(car)

# Load required data if not in environment
if(!exists("var_results")) load("variability_gpa_analysis_results.RData")
if(!exists("projected_coords")) load("variability_pca_results.RData")

# ----------------------------------------
# 1. Analysis of Procrustes distances
# ----------------------------------------
# Prepare data for Kruskal-Wallis test
kw_data <- data.frame(
    Distance = unlist(proc_distances),
    Group = rep(names(proc_distances), sapply(proc_distances, length))
)

# Perform Kruskal-Wallis test
kw_result <- kruskal.test(Distance ~ Group, data = kw_data)

# Print Kruskal-Wallis results
cat("\n=== 1. Kruskal-Wallis Test Results ===\n")
print(kw_result)

# Calculate and print summary statistics for each group
cat("\nProcrustes Distance Summary Statistics:\n")
group_stats <- do.call(rbind, lapply(names(proc_distances), function(group) {
    data.frame(
        Group = group,
        Mean = mean(proc_distances[[group]]),
        SD = sd(proc_distances[[group]]),
        Median = median(proc_distances[[group]]),
        Min = min(proc_distances[[group]]),
        Max = max(proc_distances[[group]])
    )
}))
print(group_stats)

# ----------------------------------------
# 2. Multivariate dispersion analysis
# ----------------------------------------
# Create a single matrix with all specimens
all_coords <- list()
group_labels <- character()

for(group in names(var_results)) {
    coords <- var_results[[group]]$coords
    n_specimens <- dim(coords)[3]
    
    for(i in 1:n_specimens) {
        all_coords[[length(all_coords) + 1]] <- c(coords[,,i])
        group_labels <- c(group_labels, group)
    }
}

# Convert to matrix and calculate distances
coord_matrix <- do.call(rbind, all_coords)
dist_matrix <- dist(coord_matrix)
group_factor <- factor(group_labels)

# Perform multivariate dispersion test
disp_test <- betadisper(dist_matrix, group_factor)
disp_anova <- anova(disp_test)

cat("\n=== 2. Multivariate Dispersion Test Results ===\n")
print(disp_anova)

# ----------------------------------------
# 3. Tukey Test for multiple comparisons
# ----------------------------------------
tukey_result <- TukeyHSD(disp_test)

cat("\n=== 3. Tukey Test Results ===\n")
print(tukey_result)

# Calculate and print group dispersions
group_dispersions <- data.frame(
    Group = levels(group_factor),
    Mean_Distance = tapply(disp_test$distances, group_factor, mean),
    SD_Distance = tapply(disp_test$distances, group_factor, sd)
)

cat("\nGroup Dispersion Statistics:\n")
print(group_dispersions)

# ----------------------------------------
# Save complete statistical results
# ----------------------------------------
statistical_results <- list(
    kruskal_wallis = list(
        test = kw_result,
        group_stats = group_stats
    ),
    multivar_dispersion = list(
        test = disp_anova,
        group_dispersions = group_dispersions
    ),
    tukey_test = tukey_result
)

save(statistical_results, file = "complete_statistical_results.RData")

# ----------------------------------------
# Create summary report
# ----------------------------------------
cat("\n=== COMPLETE STATISTICAL ANALYSIS SUMMARY ===\n")
cat("\n1. Overall Variability (Kruskal-Wallis):\n")
cat("Chi-squared =", kw_result$statistic,
    "\np-value =", kw_result$p.value)

cat("\n\n2. Group Differences (Multivariate Dispersion):\n")
cat("F =", disp_anova$F[1],
    "\np-value =", disp_anova$`Pr(>F)`[1])

cat("\n\n3. Pairwise Comparisons Summary:\n")
print(data.frame(
    Comparison = rownames(tukey_result$group),
    P_value = tukey_result$group[,"p adj"]
))


# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 5 Validation Analyses
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------------------
# 5.1 Validation Analyses Of Substitution Locus (Landmark 11-12)
# -------------------------------------------------------------------------------------------------
# Clear workspace
rm(list = ls())

# Intall and load required packages
install.packages("geomorph")
install.packages("ggplot2")
library(geomorph)
library(ggplot2)

# Set working directory
setwd("~/Desktop/validation_analyses_trilobite")

# -----------------------------------------------------------------------------
# GPA
# -----------------------------------------------------------------------------
# ----------------------------------------
# Step 1: Process file names
# ----------------------------------------
# Get file list
files <- list.files(".", pattern = "\\.txt$", full.names = TRUE)
cat("Found", length(files), "text files\n")

# Process file names
file_names <- basename(files)
clean_names <- gsub("\\.txt\\.txt$", ".txt", file_names)
clean_names <- gsub("\\.txt$", "", clean_names)

# ----------------------------------------
# Step 2: Set landmark template
# ----------------------------------------
nlms <- c(2,    # Dimensions
          16,   # Fixed landmarks
          95,   # Minimum points for glabella curve
          199,  # Minimum points for facial suture curve
          173,  # Minimum points for anterior margin curve
          134)  # Minimum points for posterior margin curve

# ----------------------------------------
# Step 3: Read and process data
# ----------------------------------------
# Import TriloMorph functions (Serra et al.2023)
source("~/Downloads/TriloMorph_Functions.r")

# Read all specimen data
all_specimens <- shapRead(clean_names, subdir = ".")

# Process and standardize data
processed_specimens <- shapFix(all_specimens, nlms, lm.scale = TRUE)

# ----------------------------------------
# Step 4: Define species groups and create subsets
# ----------------------------------------
# Define Olenellina species list
olenellina_species <- c(
    "Choubertella_spinosa",
    "Daguinaspis_ambroggii",
    "Elliptocephala_asaphoides",
    "Fallotaspis_bondoni",
    "Holmia_kjerulfi",
    "Judomia_tera",
    "Montezumaspis_parallela",
    "Nephrolenellus_geniculatus",
    "Nevadia_weeksi",
    "Olenellus_gilberti",
    "Peachella_iddingsi",
    "Wanneria_walcottana"
)

# Create non-Olenellina species list by excluding Olenellina species
all_specimen_names <- dimnames(processed_specimens)[[3]]
all_non_olenellina <- all_specimen_names[!all_specimen_names %in% olenellina_species]

# ----------------------------------------
# Step 5: Perform GPA
# ----------------------------------------
# Process GPA for all specimens
gpa_results <- gpagen(processed_specimens, 
                     Proj = TRUE,
                     PrinAxes = FALSE)

# Create and process non-Olenellina dataset
non_olenellina_indices <- which(dimnames(processed_specimens)[[3]] %in% all_non_olenellina)
non_olenellina_subset <- processed_specimens[,,non_olenellina_indices]

gpa_no_olenellina <- gpagen(non_olenellina_subset, 
                           Proj = TRUE,
                           PrinAxes = FALSE)

# ----------------------------------------
# Step 6: Save results
# ----------------------------------------
# Save complete dataset GPA results
save(gpa_results, file = "gpa_analysis_results.RData")

# Save non-Olenellina dataset GPA results
save(gpa_no_olenellina, 
     file = "gpa_analysis_results_no_olenellina.RData")

# ----------------------------------------
# Step 7: Print processing summary
# ----------------------------------------
cat("\nData Processing Summary:\n")
cat("----------------------------------------\n")
cat("Total number of specimens:", dim(processed_specimens)[3], "\n")
cat("Number of non-Olenellina specimens:", dim(non_olenellina_subset)[3], "\n")


# -----------------------------------------------------------------------------
# Validation Analyses Of Substitution Locus (Landmark 11 + Landmark 12)
# -----------------------------------------------------------------------------
load("gpa_analysis_results.RData")  # All specimens
load("gpa_analysis_results_no_olenellina.RData")  # Non-Olenellina specimens

# ----------------------------------------
# Step 1: Analysis function definition
# ----------------------------------------
analyze_all_landmarks <- function(gpa_data, group_name = "") {
    # Extract all fixed landmarks (1-16)
    n_landmarks <- 16
    landmark_coords <- gpa_data$coords[1:n_landmarks,,]
    n_specimens <- dim(landmark_coords)[3]
    
    # Calculate mean position (consensus) for each landmark
    landmark_consensus <- apply(landmark_coords, c(1,2), mean)
    
    # Initialise array for storing individual landmark distances
    landmark_proc_dist <- array(0, dim = c(n_specimens, n_landmarks))
    
    # Calculate procrustes distances for each landmark
    for(i in 1:n_specimens) {
        for(j in 1:n_landmarks) {
            landmark_proc_dist[i,j] <- sqrt(sum((landmark_coords[j,,i] - landmark_consensus[j,])^2))
        }
    }

    # Calculate summary statistics for each landmark
    landmark_stats <- list()
    for(i in 1:n_landmarks) {
        distances <- landmark_proc_dist[,i]
        landmark_stats[[i]] <- list(
            landmark = paste0("LM", i),
            mean = mean(distances),
            sd = sd(distances),
            cv = sd(distances)/mean(distances) * 100,
            min = min(distances),
            max = max(distances)
        )
    }
    
    results <- list(
        group = group_name,
        n_specimens = n_specimens,
        landmark_stats = landmark_stats
    )
    
    return(results)
}

# ----------------------------------------
# Step 2: Generate report function
# ----------------------------------------
generate_landmark_report <- function(results_all, results_no_olenellina) {
    cat("\n=== Complete Landmark Position Analysis Report ===\n")
    cat("\n1. Sample Sizes:\n")
    cat("All Specimens:", results_all$n_specimens, "\n")
    cat("Non-Olenellina:", results_no_olenellina$n_specimens, "\n")
    
    # Create comprehensive comparison dataframe
    comparison_df <- data.frame(
        Landmark = character(),
        Mean_All = numeric(),
        Mean_NonOlen = numeric(),
        SD_All = numeric(),
        SD_NonOlen = numeric(),
        CV_All = numeric(),
        CV_NonOlen = numeric(),
        CV_Difference = numeric(),
        stringsAsFactors = FALSE
    )
    
    for(i in 1:16) {
        stats_all <- results_all$landmark_stats[[i]]
        stats_no_olen <- results_no_olenellina$landmark_stats[[i]]
        
        comparison_df[i,] <- c(
            paste0("LM", i),
            stats_all$mean,
            stats_no_olen$mean,
            stats_all$sd,
            stats_no_olen$sd,
            stats_all$cv,
            stats_no_olen$cv,
            stats_all$cv - stats_no_olen$cv
        )
    }
    
    # Convert numeric columns
    numeric_cols <- 2:ncol(comparison_df)
    comparison_df[,numeric_cols] <- apply(comparison_df[,numeric_cols], 2, as.numeric)
    
    cat("\n2. Individual Landmark Analysis:\n")
    print(format(comparison_df, digits = 4))
    
    # Calculate statistics for non-substituted landmarks
    non_subst_indices <- c(1:10, 13:16)
    cv_non_subst_all <- comparison_df$CV_All[non_subst_indices]
    cv_non_subst_no_olen <- comparison_df$CV_NonOlen[non_subst_indices]
    
    cat("\n3. Analysis Summary:\n")

    # Performance of substituted landmarks
    cat("\na) Performance of Substituted Landmarks:\n")
    cat("LM11 CV: All =", round(comparison_df$CV_All[11], 2), 
        "%, Non-Olenellina =", round(comparison_df$CV_NonOlen[11], 2), "%\n")
    cat("LM12 CV: All =", round(comparison_df$CV_All[12], 2), 
        "%, Non-Olenellina =", round(comparison_df$CV_NonOlen[12], 2), "%\n")
    
    # Overall landmark variation
    cat("\nb) Overall Landmark Variation:\n")
    cat("Non-substituted landmarks CV range:\n")
    cat("All Specimens:", round(min(cv_non_subst_all), 2), "% to", 
        round(max(cv_non_subst_all), 2), "%\n")
    cat("Non-Olenellina:", round(min(cv_non_subst_no_olen), 2), "% to", 
        round(max(cv_non_subst_no_olen), 2), "%\n")
    
    # Mean CV for context
    cat("\nc) Mean CV for context:\n")
    cat("Mean CV (non-substituted landmarks) - All Specimens:", 
        round(mean(cv_non_subst_all), 2), "%\n")
    cat("Mean CV (non-substituted landmarks) - Non-Olenellina:", 
        round(mean(cv_non_subst_no_olen), 2), "%\n")
    
    # Save results
    write.csv(comparison_df, "landmark_complete_analysis.csv", row.names = FALSE)
}

# ----------------------------------------
# Step 3: Analysis execution
# ----------------------------------------
# Analyse both datasets
results_all <- analyze_all_landmarks(gpa_results, "All Specimens")
results_no_olenellina <- analyze_all_landmarks(gpa_no_olenellina, "Non-Olenellina")

# Generate comprehensive report
generate_landmark_report(results_all, results_no_olenellina)

# -------------------------------------------------------------------------------------------------
# 5.2 Validation Analyses Of Substitution Facial Suture (Connected by LM10 + LM12) 
# -------------------------------------------------------------------------------------------------
# Clear workspace
rm(list = ls())

# Intall and load required packages
install.packages("geomorph")
install.packages("ggplot2")
install.packages("gridExtra")
library(geomorph)
library(ggplot2)
library(gridExtra)

# Set working directory
setwd("~/Desktop/validation_analyses_trilobite")

# Load GPA results
load("gpa_analysis_results.RData")  # All specimens
load("gpa_analysis_results_no_olenellina.RData")  # Non-Olenellina specimens

# -----------------------------------------------------------------------------
# Calculate Substituted Facial Suture Procrustes Distance And Deviations
# -----------------------------------------------------------------------------
# ----------------------------------------
# Step 1: Function to calculate Procrustes Distances and deviations
# ----------------------------------------
calculate_deviations <- function(gpa_data, suture_range = 112:310) {
    # Calculate mean shape
    mean_shape <- apply(gpa_data$coords[suture_range,,], c(1,2), mean)
    
    # Calculate Procrustes distance from each specimen to mean shape
    n_specimens <- dim(gpa_data$coords)[3]
    proc_dist <- numeric(n_specimens)
    
    for(i in 1:n_specimens) {
        specimen_coords <- gpa_data$coords[suture_range,,i]
        proc_dist[i] <- sqrt(sum((specimen_coords - mean_shape)^2))
    }
    
    return(proc_dist)
}

# ----------------------------------------
# Step 2: Generate analysis report
# ----------------------------------------
generate_report <- function(gpa_data, name = "") {
    proc_dist <- calculate_deviations(gpa_data)
    
    report <- list(
        n_specimens = dim(gpa_data$coords)[3],
        mean_dist = mean(proc_dist),
        sd_dist = sd(proc_dist),
        cv = sd(proc_dist)/mean(proc_dist) * 100,
        min_dist = min(proc_dist),
        max_dist = max(proc_dist),
        range = max(proc_dist) - min(proc_dist)
    )
    
    cat("\n=== Facial Suture Analysis Report for", name, "===\n")
    cat("Number of specimens:", report$n_specimens, "\n")
    cat("Mean Procrustes distance:", round(report$mean_dist, 6), "\n")
    cat("Standard deviation:", round(report$sd_dist, 6), "\n")
    cat("Coefficient of variation (%):", round(report$cv, 2), "\n")
    cat("Range:", round(report$range, 6), "\n")
    cat("Min distance:", round(report$min_dist, 6), "\n")
    cat("Max distance:", round(report$max_dist, 6), "\n")
    
    return(report)
}

# Generate statistical reports
all_report <- generate_report(gpa_results, "All Specimens")
no_olen_report <- generate_report(gpa_no_olenellina, "Non-Olenellina Specimens")

# Save report as CSV
report_df <- data.frame(
    Metric = c("N_specimens", "Mean_distance", "SD", "CV", "Range", "Min", "Max"),
    All_Specimens = c(all_report$n_specimens, all_report$mean_dist, 
                     all_report$sd_dist, all_report$cv, all_report$range,
                     all_report$min_dist, all_report$max_dist),
    Non_Olenellina = c(no_olen_report$n_specimens, no_olen_report$mean_dist,
                      no_olen_report$sd_dist, no_olen_report$cv, no_olen_report$range,
                      no_olen_report$min_dist, no_olen_report$max_dist)
)

write.csv(report_df, "facial_suture_analysis_report.csv", row.names = FALSE)


# -----------------------------------------------------------------------------
# Visualisation of Facial Suture Procrustes Configuration and Deviations Histrogram
# -----------------------------------------------------------------------------
# ----------------------------------------
# Step 1: Facial suture visualisation function
# ----------------------------------------
plot_facial_sutures <- function(gpa_data, title = "") {
    n_specimens <- dim(gpa_data$coords)[3]
    
    # Get coordinates for LM10 and LM12
    lm10 <- data.frame(
        x = gpa_data$coords[10,1,1],
        y = gpa_data$coords[10,2,1],
        label = "LM10 (start)"
    )
    
    lm12 <- data.frame(
        x = gpa_data$coords[12,1,1],
        y = gpa_data$coords[12,2,1],
        label = "LM12 (end)"
    )
    
    # Define facial suture range
    suture_range <- 112:310
    
    # Store curve data
    all_curves <- data.frame()
    for(i in 1:n_specimens) {
        curve_coords <- gpa_data$coords[suture_range,,i]
        specimen_df <- data.frame(
            x = curve_coords[,1],
            y = curve_coords[,2],
            specimen = i,
            type = "Individual specimens"
        )
        all_curves <- rbind(all_curves, specimen_df)
    }

    # Calculate mean curve
    mean_curve <- data.frame(
        x = apply(gpa_data$coords[suture_range,1,], 1, mean),
        y = apply(gpa_data$coords[suture_range,2,], 1, mean),
        type = "Consensus shape"
    )

    # Create base plot with consistent size
    p <- ggplot() +
        # Add individual curves
        geom_path(data = all_curves, 
                 aes(x = x, y = y, group = specimen, color = type),
                 alpha = 0.3,
                 linewidth = 0.5) +
        # Add mean curve
        geom_path(data = mean_curve,
                 aes(x = x, y = y, color = "Consensus shape"),
                 linewidth = 1) +
        # Add landmark points
        geom_point(data = rbind(lm10, lm12),
                  aes(x = x, y = y),
                  color = "red", 
                  size = 3) +
        # Add landmark labels
        geom_text(data = rbind(lm10, lm12),
                 aes(x = x, y = y, label = label),
                 hjust = -0.1, 
                 vjust = -0.5,
                 size = 3) +
        # Set color scale
        scale_color_manual(values = c("Consensus shape" = "red",
                                    "Individual specimens" = "grey70")) +
        theme_minimal() +
        ggtitle(title) +
        coord_fixed(xlim = c(-0.06, 0.04), ylim = c(-0.08, 0.06)) +
        theme(plot.title = element_text(hjust = 0.5),
              legend.position.inside = c(0.85, 0.95),
              legend.title = element_blank(),    
              legend.background = element_blank(),  
              legend.key = element_blank(),      
              legend.text = element_text(size = 8))
    
    return(p)
}

# ----------------------------------------
# Step 2: Create a combined plot of facial suture configurations and deviations histrogram
# ----------------------------------------
plot_deviations <- function(proc_dist, title = "") {
    mean_dist <- mean(proc_dist)
    sd_val <- sd(proc_dist)
    cv_val <- sd_val/mean_dist * 100
    n_specimens <- length(proc_dist)
    
    # Create stats text
    stats_text <- sprintf("n = %d\nSD = %.3f\nCV = %.1f%%", 
                         n_specimens, sd_val, cv_val)
    
    # Create mean label
    mean_label <- sprintf("Mean = %.3f", mean_dist)
    
    # Calculate histogram breaks for proper y-axis scaling
    hist_data <- hist(proc_dist, plot = FALSE, breaks = 20)
    max_height <- max(hist_data$counts)
    
    p <- ggplot(data.frame(distance = proc_dist), aes(x = distance)) +
        geom_histogram(fill = "lightblue", color = "white", bins = 20) +
        geom_vline(xintercept = mean_dist, color = "red", linewidth = 1) +
        # Add mean label
        annotate("text", x = mean_dist, y = max_height,
                label = mean_label, hjust = -0.1) +
        # Add stats text in top right
        annotate("text", x = max(proc_dist), y = max_height,
                label = stats_text, hjust = 1, vjust = 1) +
        theme_minimal() +
        ggtitle(title) +
        xlab("Procrustes Distance") +
        ylab("Frequency") +
        theme(plot.title = element_text(hjust = 0.5))
    return(p)
}

# Function to create combined plot
create_combined_plot <- function(gpa_results, gpa_no_olenellina) {
    # Create individual plots
    p1 <- plot_facial_sutures(gpa_results, "All Orders Suture Configurations")
    p2 <- plot_facial_sutures(gpa_no_olenellina, "Without Olenellina Configurations")
    
    d1 <- plot_deviations(calculate_deviations(gpa_results), 
                         "All Orders Deviations")
    d2 <- plot_deviations(calculate_deviations(gpa_no_olenellina), 
                         "Without Olenellina Deviations")
    
    # Generate combined plot with equal sizes
    combined_plot <- grid.arrange(
        p1, p2,
        d1, d2,
        ncol = 2,
        heights = c(1.5, 1)
    )
    return(combined_plot)
}

# Save plots
ggsave("facial_suture_complete_analysis.pdf", 
       create_combined_plot(gpa_results, gpa_no_olenellina), 
       width = 12, height = 10)


# -------------------------------------------------------------------------------------------------
# 5.3 Validation Analyses Of Substitution Scheme For Suborder Olenellina In PCA 
# -------------------------------------------------------------------------------------------------
# ----------------------------------------
# 1. Set up environment
# ----------------------------------------
# Set working directory
setwd("~/Desktop/validation_analyses_trilobite")

# Load required packages
library(geomorph)
library(vegan)
library(ggplot2)

# ----------------------------------------
# 2. Data preparation
# ----------------------------------------
# Load specimen-order mapping information
load("specimen_order_mapping.RData")  # loads specimen_orders object

# Load complete dataset (81 specimens)
load("gpa_analysis_results.RData")  # loads as gpa_results
gpan_full <- gpa_results

# Load dataset without olenellina (69 specimens)
load("gpa_analysis_results_no_olenellina.RData")  # loads as gpa_no_olenellina
gpan_nonolen <- gpa_no_olenellina

# Define olenellina specimens
olenellina_species <- c(
    "Choubertella_spinosa",
    "Daguinaspis_ambroggii",
    "Elliptocephala_asaphoides",
    "Fallotaspis_bondoni",
    "Holmia_kjerulfi",
    "Judomia_tera",
    "Montezumaspis_parallela",
    "Nephrolenellus_geniculatus",
    "Nevadia_weeksi",
    "Olenellus_gilberti",
    "Peachella_iddingsi",
    "Wanneria_walcottana"
)

# ----------------------------------------
# 3. PCA analysis
# ----------------------------------------
# Perform PCA on complete dataset
pca_full <- gm.prcomp(gpan_full$coords)

# Extract indices for non-olenellina specimens
non_olen_indices <- !sapply(dimnames(gpan_full$coords)[[3]], function(x) {
    any(sapply(olenellina_species, function(y) grepl(y, x, ignore.case = TRUE)))
})

# Extract PCA scores for non-olenellina specimens from complete dataset
pca_scores_full_nonolen <- pca_full$x[non_olen_indices, ]

# Perform PCA on dataset without olenellina
pca_nonolen <- gm.prcomp(gpan_nonolen$coords)

# ----------------------------------------
# 4. Visualisation
# ----------------------------------------
# Calculate explained variance
var_full <- (pca_full$sdev^2)/sum(pca_full$sdev^2) * 100
var_nonolen <- (pca_nonolen$sdev^2)/sum(pca_nonolen$sdev^2) * 100

# Create a vector for olenellina indices in the full dataset
olenellina_indices <- sapply(dimnames(gpan_full$coords)[[3]], function(x) {
    any(sapply(olenellina_species, function(y) grepl(y, x, ignore.case = TRUE)))
})

# Set up plotting window with appropriate margins
dev.new(width = 8, height = 6)
par(mar = c(5, 5, 4, 2))  # Bottom, left, top, right margins

# Create main plot area for complete dataset
plot(NA, NA,
     xlim = range(pca_full$x[,1]),
     ylim = range(pca_full$x[,2]),
     xlab = paste0("PC1 (", round(var_full[1], 1), "%)"),
     ylab = paste0("PC2 (", round(var_full[2], 1), "%)"),
     main = "PCA Comparison",
     asp = 1)

# Add zero lines
abline(h = 0, v = 0, lty = 2, col = "gray")

# Add non-Olenellina points from complete dataset (black solid circles)
points(pca_full$x[!olenellina_indices,1], 
       pca_full$x[!olenellina_indices,2],
       pch = 19,  # Solid circle
       col = "black",
       cex = 1.2)

# Add Olenellina points from complete dataset (red solid circles)
points(pca_full$x[olenellina_indices,1], 
       pca_full$x[olenellina_indices,2],
       pch = 19,  # Solid circle
       col = "red",
       cex = 1.2)

# Add non-Olenellina subset points (gray open circles)
points(pca_nonolen$x[,1], 
       pca_nonolen$x[,2],
       pch = 1,   # Open circle
       col = "gray50",
       cex = 1.2)

# Add legend on the right side
legend("topright", 
       legend = c("Complete dataset (Non-Olenellina)", 
                 "Complete dataset (Olenellina)",
                 "Non-Olenellina subset"),
       col = c("black", "red", "gray50"),
       pch = c(19, 19, 1),
       cex = 0.7,
       bty = "n")

# Add sample size information at the top
mtext("Complete dataset: n = 81", 
      side = 3, line = 0, adj = 0, font = 3)
mtext("Non-Olenellina subset: n = 69", 
      side = 3, line = 1, adj = 0, font = 3, col = "gray50")

# Add Non-Olenellina PC variance information (gray text)
mtext(paste0("PC1 (", round(var_nonolen[1], 1), "%)"), 
      side = 1, line = 3, col = "gray50", adj = 0.75)  # Near right end of x-axis

mtext(paste0("PC2 (", round(var_nonolen[2], 1), "%)"), 
      side = 2, line = 3, col = "gray50", adj = 0.9)  # Near top of y-axis

# Save the plot
dev.copy2pdf(file = "combined_pca_comparison.pdf")

# ----------------------------------------
# 5. Statistical analyses
# ----------------------------------------
# ----------------------------------------
# 5.1 Compare PC variance explained
# ----------------------------------------
# Calculate variance explained
var_full <- (pca_full$sdev^2)/sum(pca_full$sdev^2) * 100
var_nonolen <- (pca_nonolen$sdev^2)/sum(pca_nonolen$sdev^2) * 100

var_comp <- data.frame(
    PC = paste0("PC", 1:4),
    Full_Dataset = var_full[1:4],
    Reduced_Dataset = var_nonolen[1:4]
)
print(var_comp)

# ----------------------------------------
# 5.2 Calculate PC correlations
# ----------------------------------------
pc1_correlation <- cor(pca_scores_full_nonolen[,1], pca_nonolen$x[,1])
pc2_correlation <- cor(pca_scores_full_nonolen[,2], pca_nonolen$x[,2])
print(paste("PC1 correlation between datasets:", round(pc1_correlation, 3)))
print(paste("PC2 correlation between datasets:", round(pc2_correlation, 3)))

# ----------------------------------------
# 5.3 Procrustes analysis
# ----------------------------------------
coords1 <- pca_scores_full_nonolen[,1:2]
coords2 <- pca_nonolen$x[,1:2]
proc_test <- protest(coords1, coords2)
print(proc_test)

# ----------------------------------------
# 5.4 Calculate centroid shifts
# ----------------------------------------
# Calculate centroids for each dataset
full_scores <- pca_scores_full_nonolen[,1:2]
nonolen_scores <- pca_nonolen$x[,1:2]
groups_full <- specimen_orders$order[non_olen_indices]

# Function to calculate centroids
calc_centroids <- function(scores, groups) {
    centroids <- aggregate(scores, list(Group = groups), mean)
    names(centroids)[2:3] <- c("X1", "X2")
    return(centroids)
}

centroids_full <- calc_centroids(full_scores, groups_full)
centroids_nonolen <- calc_centroids(nonolen_scores, groups_full)

# Ensure matching order of groups
centroids_nonolen <- centroids_nonolen[match(centroids_full$Group, centroids_nonolen$Group),]

# Calculate shifts
centroid_shifts <- data.frame(
    Order = centroids_full$Group,
    PC1_shift = centroids_full$X1 - centroids_nonolen$X1,
    PC2_shift = centroids_full$X2 - centroids_nonolen$X2
)

# Calculate total shift
centroid_shifts$Total_shift <- sqrt(centroid_shifts$PC1_shift^2 + centroid_shifts$PC2_shift^2)

# Sort by total shift
centroid_shifts <- centroid_shifts[order(centroid_shifts$Total_shift),]
print(centroid_shifts)

# ----------------------------------------
# 5.5 PERMANOVA test
# ----------------------------------------
# Create grouping variable
groups <- rep(c("Full", "Reduced"), each=nrow(coords1))

# Combine datasets
combined_scores <- rbind(coords1, coords2)

# Calculate distance matrix
combined_dist <- dist(combined_scores)

# Perform PERMANOVA test
perm_test <- adonis2(combined_dist ~ groups, permutations = 999)
print(perm_test)

