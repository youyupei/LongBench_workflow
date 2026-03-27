All Rmd file should follow the following:
1. All code need to be versioned and can be knit by versioned_knit.R
2. Put proper chunk names so the figure saved with the same name. 
3. Need a set param with cache_dir and fig.path, can leave it with Null if use versioned_knit.R to knit the file. It will be automatically set to  base_name + date + _cache and base_name + date + _figures
4. set
        knitr::opts_chunk$set(
        cache = !is.null(params$cache_dir),  # Enable caching only if cache_dir is provided
        fig.width = 15,                     # Optional: Set default figure width
        fig.height = 10,                    # Optional: Set default figure height
        fig.path = params$fig.path,         # Set figure path (NULL if not provided)
        dev = if (!is.null(params$fig.path)) c("png", "svg") else NULL  # This will save figures in both PNG and SVG formats if fig.path is provided, otherwise it won't save any figures
        )

5. All figure will need to directly plot  without saving to a file. 
6. Since the raw points for the figure will usually be required, use different chunk for generating the raw points and plotting the figure. The chunk for generating the raw points can be set with cache = TRUE, and the chunk for plotting the figure can be set with cache = FALSE.

## Output:
The versioned_knit.R will save the output in the same directory as the Rmd file with the name of base_name + date + .html. The cache and figure will be saved in the same directory with the name of base_name + date + _cache and base_name + date + _figures.


