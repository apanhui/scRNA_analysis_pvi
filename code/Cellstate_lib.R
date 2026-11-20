WriteTable<-function(
  x,
  file,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  ...
) {
  write.table(
    x,
    file = file,
    quote = quote,
    sep = sep,
    row.names = row.names,
    col.names = col.names,
    ...
  )
}
Message<-function(..., type = c("Running", "Warning", "Stopped")) {
  type <- match.arg(arg = type)
  message(paste0("\n[", as.character(Sys.time()), "] [", type, "]: ", ...))
  if (type == "Stopped") {
    q(status = 1)
  }
}

Filter<-function (f, x)
{
    f <- match.fun(f)
    ind <- as.logical(unlist(lapply(x, f)))
    x[which(ind)]
}
validCharacters<-function(x) {
  return(
    nzchar(x = x, keepNA = FALSE) &
      !is.na(x = x) &
      is.atomic(x = x)
  )
}

FetchData2<-function(object, vars, ...) {
  new_names <- names(vars)
  df <- FetchData(object, vars = vars, ...)
  if (any(validCharacters(new_names))) {
    vars <- vars[validCharacters(new_names)]
    vars <- vars[vars %in% colnames(df)]
    df <- df %>%
      dplyr::rename(all_of(vars))
  }
  df <- df %>% mutate(Cells = rownames(.), .before = 1)
  return(df)
}

