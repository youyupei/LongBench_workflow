# This script knit a Rmarkdown file with a custom output name: base_name_version(if specified)_date.html

# read from command line
input_file <- commandArgs(trailingOnly = TRUE)[1]

# get the base name of the input file
base_name <- tools::file_path_sans_ext(basename(input_file))

# get the version from the yaml front matter
version <- rmarkdown::yaml_front_matter(input_file)$version
if (is.null(version)) {
    output_name <- paste0(base_name, "_", Sys.Date(), ".html")
} else {
    output_name <- paste0(base_name, "_", version, "_", Sys.Date(), ".html")
}

# render the Rmarkdown file
message("Rendering ", input_file, " to ", output_name)
rmarkdown::render(input_file, output_file = output_name)