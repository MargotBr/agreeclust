# Dataset binary_data_for_example
load("~/Documents/agreeclust/data/binary_data_for_example.rda")
binary_data_for_example <- pedagdataBin
binary_data_for_example <- thinkr::clean_names(binary_data_for_example)
rownames(binary_data_for_example) <- stringr::str_replace_all(tolower(rownames(binary_data_for_example)), "\\.", "_")
usethis::use_data(binary_data_for_example, overwrite = TRUE, version = 3, compress = "bzip2")

# Dataset continuous_data_for_example
load("~/Documents/agreeclust/data/continuous_data_for_example.rda")
continuous_data_for_example <- pedagdataCont
continuous_data_for_example <- thinkr::clean_names(continuous_data_for_example)
rownames(continuous_data_for_example) <- stringr::str_replace_all(tolower(rownames(continuous_data_for_example)), "\\.", "_")
usethis::use_data(continuous_data_for_example, overwrite = TRUE, version = 3, compress = "bzip2")
