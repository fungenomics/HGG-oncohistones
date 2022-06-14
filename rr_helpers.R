
# Initialization ----

#' Clean-up the rr repository of example files for initialization
#' 
#' A helper function for cleaning up a repository generated from the rr
#' template, by removing the example documents & their outputs
rr_initialize_cleanup <- function(dry_run = FALSE) {
    
    require(here)
    
    # For Yes/No prompts, treat a selection of 1 as TRUE, and 2 as FALSE
    s2l <- function(x) ifelse(x == 1, TRUE, FALSE)
    
    # Delete the example documents & HTMLs
    del_ex <- s2l(menu(c("Yes", "No"), title = paste0("Delete example Rmds and corresponding HTMLs? i.e.\n",
                                                      "  git rm analysis/01-first_step.{Rmd,html}\n",
                                                      "  git rm analysis/02-second_step.{Rmd,html}")))
    if (del_ex & !dry_run) {
        
        system(paste(paste0("git rm ", Sys.glob(here("analysis/01-first_step.*"))), collapse = ";"))
        system(paste(paste0("git rm ", Sys.glob(here("analysis/02-second_step.*"))), collapse = ";"))
        
    }
    
    # Delete the example output / figures
    del_out <- s2l(menu(c("Yes", "No"), title = paste0("Delete example outputs/figures? i.e.\n",
                                                       "  git rm output/01/mtcars.{tsv,desc}\n",
                                                       "  git rm figures/01/pressure-1.{png,pdf}\n",
                                                       "  git rm figures/02/figure2-*")))
    if (del_out & !dry_run) {
        
        system(paste(paste0("git rm ", Sys.glob(here("output/01/mtcars*"))), collapse = ";"))
        system(paste(paste0("git rm ", Sys.glob(here("figures/01/pressure-1*"))), collapse = ";"))
        system(paste(paste0("git rm ", Sys.glob(here("figures/02/figure2-*"))), collapse = ";"))
        
    }
    
    # Clear README file from rr repository
    del_readme <- s2l(menu(c("Yes", "No"), title = "Clear README.md?"))
    if (del_readme & !dry_run) {
        file.remove(here("README.md"))
        file.create(here("README.md"))
    }
    
    # Delete the images that are used for demo in the rr repository
    del_img <- s2l(menu(c("Yes", "No"), title = "Delete images used in rr repository README?"))
    # ... except for the lab logo, needed in the template
    img_paths <- setdiff(list.files(here("include/img"), full.names = TRUE),
                         here("include/img/kleinman_lab_logo.png"))
    if (del_img & !dry_run) {
        
        system(paste(paste0("git rm ", img_paths), collapse = ";"))
        
    }
    
}


#' Initialize a reproducible research repository
#' 
#' A function to setup a new repository from the rr template, by customizing
#' templates / headers with author info and project info. Interactively prompts
#' user for permission to delete files, and for the following information to
#' update templates:
#'  * Author name(s)
#'  * Contact email address
#'  * Short project name
#'  * Link to repository on GitHub (or other link)
#'
#' @param dry_run Logical, whether to actually perform the changes. Default: FALSE
#'
#' @return Nothing
rr_initialize <- function(dry_run = FALSE) {
    
    require(here)
    
    if (dry_run) message("Performing a dry-run of repository initialization...\n")
    
    rr_initialize_cleanup(dry_run = dry_run)
    
    # Replace the author name
    name <- readline("Author name(s): ") 
    # ....Construct sed command
    cmd_update_template_name <- paste0("sed -i '' '3s/Selin Jessa/", name, "/g' ", here("include/template.Rmd"))
    print(cmd_update_template_name)
    # ...Execute command
    if (!dry_run) system(cmd_update_template_name)
    
    # Replace author email
    email <- readline("Contact email: ")
    cmd_update_template_email <- paste0("sed -i '' '3s/selin.jessa@mail.mcgill.ca/", email, "/g' ", here("include/template.Rmd"))
    print(cmd_update_template_email)
    if (!dry_run) system(cmd_update_template_email)
    
    # Replace project name in header,
    # and change name of .Rproj file
    # and add to README
    project_name <- readline("Short project name (no spaces or slashes), e.g. matching GitHub repo name: ")
    cmd_update_proj_name <- paste0("sed -i '' 's/rr project/", project_name, " project/g' ", here("include/header.html"))
    print(cmd_update_proj_name)
    cmd_update_Rproj_name <- paste0("mv ", here("rr.Rproj"), project_name, ".Rproj")
    print(cmd_update_Rproj_name)
    if (!dry_run) {
        system(cmd_update_proj_name)
        system(cmd_update_Rproj_name)
        writeLines(project_name, here("README.md"))
    }
    
    # Replace github source link
    github_link <- readline("Link to GitHub repository: ")
    cmd_update_proj_link <- paste0("sed -i '' 's#https://github.com/sjessa/rr/#", github_link, "#g' ", here("include/header.html"))
    print(cmd_update_proj_link)
    if (!dry_run) system(cmd_update_proj_link)
    
    message("\nInitialization complete. To start an analysis, copy ",
            here("include/template.Rmd"),
            " to the analysis folder.")
    
}


# Helpers ----

#' The here() function always returns a full path from the root directory
#' This function returns a path from the project root for less clutter
#' and greater portability
#' 
#' @param path String, full path as generated by here::here("blah")
path_from_here <- function(path) {
    
    paste0(basename(here::here()), # Project root directory
           strsplit(path, basename(here::here()))[[1]][2]) # Chosen path relative to project root
    
}


# Savers / loaders ----

rr_save <- function() {
    
    
}


#' A wrapper function for writing a description with a TSV
#'
#' This funcion simply wraps readr::write_tsv, but also saves a user-provided
#' description alongside the tsv with the same filename but extension ".desc"
#'
#' @param df Data frame to write to tsv
#' @param path Path for output tsv file, as returned by here e.g. here("my/file.tsv")
#' @param desc String, brief description of file contents
#' @param verbose Logical, whether to print .desc file path to console
#'
#' @return Nothing
#'
#' @examples
#' mtcars %>% 
#'     rr_write_tsv(path = here("output/01/mtcars.tsv),
#'                  desc = "The mtcars dataset, verbatim")
rr_write_tsv <- function(df, path, desc, verbose = TRUE) {
    
    # Need readr to simplify table writing
    require(readr)
    readr::write_tsv(df, path)
    
    # Create the path for the description file, swapping .tsv extension to .desc
    desc_path <- gsub("tsv$", "desc", path)
    
    # Print the object description to the desc file
    cat(desc, file = desc_path, sep = "\n")
    
    # Output a message with path to desc file
    if (verbose) message("...writing description of ", basename(path), " to ", path_from_here(desc_path))
    
}


#' A wrapper function for saving source data along side a ggplot
#' 
#' This funcion simply wraps ggplot2::ggplot, but also saves the input data
#' alongside the figure, with the same filename but extension ".source_data.tsv".
#' This is extremely useful for being able to quickly extract the data needed to
#' regenerate the figure, sometimes also required for papers.
#' 
#' NOTE: Saving soure data only works if a ggplot is generated within a code chunk and the
#' document rendered by RMarkdown, otherwise, a warning is emitted and the function
#' proceeds with ggplot code. 
#'
#' @param df Data frame, input to ggplot2
#' @param plot_num Numeric, index of plot within R Markdown chunk, used to determine
#' the filname of the figure when the document is rendered
#' @param ... Additional parameters passed to ggplot2::ggplot, e.g. "aes(x = mpg, y = cyl)"
#'
#' @return A ggplot2 object, to which additional gg elements can be added with +,
#' same as ggplot2::ggplot
#' @export
#'
#' @examples
#' mtcars2 %>% 
#'     rr_ggplot(1, aes(x = disp, y = wt)) +
#'     geom_line() +
#'     theme_bw()
rr_ggplot <- function(df, plot_num, ...) {
    
    require(ggplot2)
    require(readr)
    
    if (!interactive()) {
        
        # TODO: Currently, it's not possible to not specify plot_num, because
        # it messes up the dots (...) which are passed to ggplot, so this if statement
        # is never evaluated, an error is thrown instead:
        
        # If the plot # is not provided
        if (missing(plot_num)) {
            
            plot_num <- 1
            # This is beacuse plots are named by their number in each chunk, but
            # that number cannot be accessed by this function
            warning("!! If more than one ggplot is generated in this chunk with rr_ggplot(),",
                    "only the source data for the first one will be saved.",
                    "Pass plot # explicitly to plot_num argument to correct this.")
            
        }
        
        # Get the figure path for the current chunk, without file extensions
        # https://github.com/yihui/knitr/issues/73#issuecomment-3514096
        fig_path <- knitr::fig_path(number = plot_num)
        
        # Make a path for the source data by appending a suffix to the figure path,
        # and write source data there as a TSV
        src_path <- paste0(fig_path, ".source_data.tsv")
        write_tsv(df, src_path)
        
        # Output a message with path to source data file
        message("...writing source data of ggplot to ", path_from_here(src_path))
        
    } else {
        
        warning("!! This function is being run in an interactive session ",
                "and the source data is NOT being saved. Render the document ",
                "to save source data.")
        
    }
    
    # Proceed with ggplot
    ggplot(data = df, ...)
    
}


rr_load <- function() {
    
    
}


#' A wrapper function for reading a TSV along with its metadata & description
#' 
#' This funcion simply wraps readr::read_tsv, but at the same time, prints some
#' information about the file to help with reproducibility & dependency tracking:
#'  * The description of the file, if one exists at the same filepath with ".desc" extension
#'  * The timestamp for when the file was last modified
#'  * The script that generated the file, under the assumption it was generated
#'    by a script within the analysis folder of this repository
#'    
#' NOTE: for files NOT produced in this repositoy, this function is not receommonded.
#' 
#' To produce a toggle button showing/hiding the output of this function in an R Markdown
#' HTML report, wrap the chunk in <div class="fold o"></div>
#' (in which case the outpout is hidden by default)
#'
#' @param path String, as returned by here::here("blah")
#' @param ... Additional parameters passed to readr::read_tsv()
#' 
#' @return A tibble, same as readr::read_tsv
#' 
#' @examples
#' mtcars %>% 
#'     rr_write_tsv(path = here("output/01/mtcars.tsv),
#'                  desc = "The mtcars dataset, verbatim")
#'                  
#' mtcars2 <- rr_read_tsv(path = here("output/01/mtcars.tsv))
rr_read_tsv <- function(path, ...) {
    
    require(readr)
    require(stringr)
    
    # Create the path for the description file, swapping .tsv extension to .desc
    desc_path <- gsub("tsv$", "desc", path)
    
    if(!file.exists(desc_path)) warning("!! No description file (.desc) found. ",
                                        "To automatically write a description file ",
                                        "when saving a tsv, use rr_write_tsv().")
    
    # Get the number of the analysis, e.g. "01"
    doc_idx <- stringr::str_extract(path, "(\\d)+")
    
    # Search analysis folder for .Rmd file matching doc_idx
    script <- list.files(here("analysis"), pattern = glob2rx(paste0(doc_idx, "*.Rmd")))
    
    # Get the timestamp for when the file contents were last modified
    timestamp <- file.info(path)$mtime
    
    # We use cat here because it returns to stdout, which will be picked
    # up by the code folding js script (hideOutput.js). Otherwise, there will
    # be one folding button per line of output
    # https://stackoverflow.com/questions/36699272/why-is-message-a-better-choice-than-print-in-r-for-writing-a-package
    cat(paste0(path_from_here(path), " info:\n",
               "...description : ",   ifelse(file.exists(desc_path), readLines(desc_path), "NOT SPECIFIED"),
               "\n...generated by: ", path_from_here(here("analysis", script)),
               "\n...last updated: ", timestamp))
    # e.g. output
    # /Users/selinjessa/Repos/rr/output/01/mtcars.tsv info:
    # ...description : The mtcars dataset, verbatim
    # ...generated by: rr/analysis/01-first_step.Rmd
    # ...last updated: 2020-05-23 22:31:53
    
    # Read the file and return as dataframe
    suppressMessages(readr::read_tsv(path, ...))
    
}
