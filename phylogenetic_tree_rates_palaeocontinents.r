# -------------------------------------------------------------------------------------------------
# Generate Morphological Tree By Uslising PCA Scores
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
# Generate Time-Calibrated Phylogenetic Tree With Evolutionary Rates and Palaeocontinent Information
# -------------------------------------------------------------------------------------------------
# Set working directory
setwd("~/Desktop/dissertation_project_trilobite")

# Load required packages
library(ape)
library(phytools)

# Load data
load("morphological_tree.RData")
load("specimen_palaeocontinent_mapping.RData")  # Contains specimen_palaeocontinents and palaeocontinent_colors
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

# Create visualisation
pdf("phylogenetic_tree_rates_palaeocontinents.pdf", width=15, height=10)

# Set up plotting parameters
par(mar=c(5,1,2,8))  # Modified to accommodate time axis

# Create color scheme for rates
rate_vals <- rate_tree$edge.length
n_colors <- 100
cols <- hcl.colors(256, "Temps")
br_cols <- cols[cut(log(rate_tree$edge.length), 256)]

# Get time interval from tree
node_heights <- nodeHeights(time_tree)
max_height <- max(node_heights)

# Calculate time ranges
time_range <- 530 - 485  # Actual tree time range (530 to 485 Ma)
plot_range <- max_height * (50/45)  # Extend plot range to accommodate 480-485 Ma

# Plot tree with exact axis limit
plot(time_tree, 
     type="phylogram",
     cex = 0.6,
     edge.width = 6,
     edge.col = br_cols,
     show.tip.label=TRUE,
     label.offset = 0.3,
     x.lim = c(0, plot_range),  # Extended range for axis
     direction = "rightward")

# Add palaeocontinent dots at tips
last_plot_coords <- get("last_plot.phylo", envir = .PlotPhyloEnv)
tip_coords <- data.frame(
    x = last_plot_coords$xx[1:length(time_tree$tip.label)],
    y = last_plot_coords$yy[1:length(time_tree$tip.label)]
)

# Add colored dots for each tip based on palaeocontinent
for(i in 1:length(time_tree$tip.label)) {
    species <- time_tree$tip.label[i]
    continent <- specimen_palaeocontinents$palaeocontinent[
        match(species, specimen_palaeocontinents$specimen)
    ]
    
    if(!is.na(continent)) {
        points(tip_coords$x[i], tip_coords$y[i],
               pch=19,        # Solid circle
               col=palaeocontinent_colors[continent],
               cex=1.2)       # Dot size
    }
}

# Add vertical dashed lines at 521 Ma and 501 Ma
# Calculate x-positions for the lines
line_pos_521 <- plot_range * (530 - 521)/(530 - 480)
line_pos_501 <- plot_range * (530 - 501)/(530 - 480)

# Add vertical lines
abline(v = line_pos_521, col = "black", lty = 2)  # Dashed line for 521 Ma
abline(v = line_pos_501, col = "black", lty = 2)  # Dashed line for 501 Ma

# Add time axis (from 530 Ma to 480 Ma)
time_scale <- seq(480, 530, by=5)  # Full time scale including 480 Ma
# Calculate positions to match tree end at 485 Ma
axis_positions <- seq(plot_range, 0, length.out = length(time_scale))
axis(1, 
     at = axis_positions,
     labels = time_scale,
     las = 1)
title(xlab = "Time (Ma)", line = 2.5)

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
       inset=c(0.15, 0.05))

# Add palaeocontinent legend
legend("topright",
       legend=names(palaeocontinent_colors),
       col=palaeocontinent_colors,
       pch=19,
       pt.cex=1.2,
       title="Palaeocontinents",
       cex=0.8,
       bty="n",
       inset=c(0.05, 0.05))

dev.off()

# Save results
save(rate_tree, file="phylogenetic_tree_rates_palaeocontinents.RData")