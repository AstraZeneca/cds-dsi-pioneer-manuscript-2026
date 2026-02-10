grep_in_col <- function (dataframe = db, string_grep = "string")  {
    counts <- sapply(dataframe, function(col) sum(grepl(string_grep, 
        as.character(col), ignore.case = TRUE)))
    return(counts[counts > 0])
}


create_header <- function(name = "Spare", sep = "=", length_n = 60) {
    if (!hasArg(length_n)) {
        len_n <- nchar(name) + 8
    } else {
        len_n <- length_n
    }
    border <- paste(rep(sep, len_n), collapse = "")
    centered_title <- sprintf(
        "%s%s%s",
        paste(rep(" ", (len_n - nchar(name)) %/% 2), collapse = ""),
        name,
        paste(rep(" ", (len_n - nchar(name) + 1) %/% 2), collapse = "")
    )
    header <- paste0(
        "#", border, "#\n",
        "#", centered_title, "#\n",
        "#", border, "#\n"
    )
    cat(header)
}
