# 0) Load required packages
library(ape)
library(dplyr)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)

# 1) Colour & shape scheme
pop_colors <- c(
  "inornatus E"         = "#DA3424",
  "inornatus W"         = "#F3A740",
  "amplus NW AZ"        = "#00A0DF",
  "amplus ColPlat"      = "#EA2E13",
  "amplus S AZ"         = "#AE95BF",
  "amplus Salt Gila"    = "#8F438F",
  "longimembris UpColR" = "#DA3424",
  "longimembris AZ hap5"= "#F4D690",
  "longimembris AZ hap6"= "#FAAA0A",
  "longimembris VirRiv" = "#DC9512",
  "longimembris GrBas"  = "#3CA0DB",
  "longimembris Mojave" = "#1A6AD2",
  "longimembris SW CA"  = "#84CAD9",
  "longimembris SE CA"  = "#C1E204",
  "longimembris SW AZ"  = "#FCEB09",
  "alticola"            = "#ED5C14",
  "flavescens"          = "#D3D3D3",
  "flavus"              = "#A9A9A9",
  "mollipilosus 1"      = "#CB731A",
  "mollipilosus 2"      = "#FCEB09"
)
pop_shapes <- c(
  "inornatus E"         = 15,
  "inornatus W"         = 16,
  "amplus NW AZ"        = 17,
  "amplus ColPlat"      = 18,
  "amplus S AZ"         = 19,
  "amplus Salt Gila"    = 20,
  "longimembris UpColR" = 21,
  "longimembris AZ hap5"= 22,
  "longimembris AZ hap6"= 23,
  "longimembris VirRiv" = 24,
  "longimembris GrBas"  = 25,
  "longimembris Mojave" = 7,
  "longimembris SW CA"  = 8,
  "longimembris SE CA"  = 6,
  "longimembris SW AZ"  = 14,
  "alticola"            = 9,
  "flavescens"          = 10,
  "flavus"              = 11,
  "mollipilosus 1"      = 13,
  "mollipilosus 2"      = 12
)

# 2) Read tree & metadata; clean pop keys
tree <- read.tree("Perognathus_RELEC_wASTRAL.newick.tree")
meta <- read.csv("map_data.csv", stringsAsFactors = FALSE)
meta$pop_key <- sub("^P_", "", meta$population)
meta$pop_key <- gsub("_", " ", meta$pop_key)

# 3) Compute tree‐depth for padding
tree_depth <- max(node.depth.edgelength(tree))

# 4) Download sf layers & define fixed map bounds
rivers   <- ne_download(scale = "large",
                        type = "rivers_lake_centerlines",
                        category = "physical",
                        returnclass = "sf")
world_sf <- ne_countries(scale = "medium", returnclass = "sf")
states_sf<- ne_states(returnclass = "sf")

fixed_xlim <- c(-125, -108)
fixed_ylim <- c(  27,   43)
bbox_sf    <- st_as_sfc(st_bbox(c(
  xmin = fixed_xlim[1],
  xmax = fixed_xlim[2],
  ymin = fixed_ylim[1],
  ymax = fixed_ylim[2]
), crs = 4326))

# 5) Repair & crop sf layers
states_sf    <- st_make_valid(states_sf)
world_clip   <- st_crop(world_sf,  bbox_sf)
states_clip  <- st_crop(states_sf, bbox_sf)
rivers_clip  <- st_crop(rivers,    bbox_sf)

# 6) Prepare tip keys
tip_raw <- gsub('^"|"$', '', tree$tip.label)
tip_key <- sub("^P_", "", tip_raw)
tip_key <- gsub("_", " ", tip_key)

# 7) Plot side‐by‐side
par(mfrow = c(1,2), oma = c(0,0,0,0))

# 7a) Phylogeny panel
par(mar = c(4,1,2,4))
plot(
  tree,
  direction      = "rightwards",
  tip.col        = pop_colors[tip_key],
  cex            = 0.8,
  font           = 3,
  show.tip.label = TRUE,
  x.lim          = c(0, tree_depth * 1.4)
)
tiplabels(
  pch    = pop_shapes[tip_key],
  col    = pop_colors[tip_key],
  cex    = 0.8,
  adj    = 1,
  offset = 0.4 * tree_depth
)

# 7b) Map panel with black background & white boundaries
par(mar = c(4,1,2,1))

# 1) draw black fill for the entire box
plot(bbox_sf, col = "black", border = "black", reset = FALSE)

# 2) overlay country borders in white
plot(st_geometry(world_clip), add = TRUE,
     col = NA, border = "white", lwd = 1.2)

# 3) overlay state borders in white
plot(st_geometry(states_clip), add = TRUE,
     col = NA, border = "white", lwd = 0.5)

# 4) rivers in steelblue
plot(st_geometry(rivers_clip), add = TRUE,
     col = "steelblue", lwd = 0.7)

# 5) individual samples
points(
  meta$lon, meta$lat,
  pch = pop_shapes[meta$pop_key],
  col = pop_colors[meta$pop_key],
  cex = 0.8
)

# Reset
par(mfrow = c(1,1))
