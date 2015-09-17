mod.formula <- function(formula, meta = NULL){
  # to create new formula with center and scale from formula and meta information

  if (is.null(meta)){
    new_formula <- formula
    label_info <- NULL
  } else {
    # get old terms
    term <- attr(terms(formula), 'term.labels')

    # get the new formula with transformation
    factor_matrix <- as.data.frame(attr(terms(formula), 'factors'))[-1,]

    ## add transformation for numeric variables
    if (is.null(dim(factor_matrix))){
      .type <- meta$type[meta$name == term]
      .center <- meta$center[meta$name == term]
      .scale <- meta$scale[meta$name == term]
      .label <- meta$label[meta$name == term]
      .units <- meta$units[meta$name == term]

      new_term <- ifelse(.type == 'factor', term,
                         paste('I((', term, '-', ifelse(is.na(.center), 0, .center), ')', '/',
                               ifelse(is.na(.scale), 1, .scale), ')', sep = ''))
      label_info <- cbind(name = term,
                          label = .label,
                          order = 1,
                          type = .type,
                          units = ifelse(.type == 'factor', '', paste('[', .units, ']', sep = '')),
                          units_fit = ifelse(.type == 'factor', '',
                                             paste('[+', ifelse(is.na(.scale), 1, .scale), ' ', .units, ']', sep = '')))
      rownames(label_info) <- term
    } else {
      new_var <- sapply(rownames(factor_matrix), FUN = function(x){
        ifelse(meta$type[meta$name == x] == 'factor', x,
               paste('I((', x, '-', ifelse(is.na(meta$center[meta$name == x]), 0, meta$center[meta$name == x]), ')', '/',
                     ifelse(is.na(meta$scale[meta$name == x]), 1, meta$scale[meta$name == x]), ')', sep = ''))
      })
      new_term <- apply(factor_matrix, 2, FUN = function(x){paste(new_var[x == 1], collapse = ':')})

      label_info <- t(apply(factor_matrix, 2, FUN = function(x){
        .name <- rownames(factor_matrix)[x == 1]
        .label <- meta$label[meta$name %in% rownames(factor_matrix)[x == 1]]
        .type  <- meta$type[meta$name %in% rownames(factor_matrix)[x == 1]]
        .units <- meta$units[meta$name %in% rownames(factor_matrix)[x == 1]]
        .scale <- meta$scale[meta$name %in% rownames(factor_matrix)[x == 1]]

        c(name = paste(.name, collapse = ':'),
          label = paste(.label, collapse = ':'),
          order = sum(x),
          type = ifelse(sum(x) > 1 | 'factor' %in% .type, 'factor', .type),
          units = ifelse(sum(x) > 1 | 'factor' %in% .type, '', paste('[', .units, ']', sep = '')),
          units_fit = ifelse(sum(x) > 1 | 'factor' %in% .type, '', paste('[+', ifelse(is.na(.scale), 1, .scale), ' ', .units, ']', sep = '')))
      }))
    }

    # create new formula
    new_formula <- eval(parse(text = paste(formula[2], formula[1], paste(new_term, collapse = '+'))))
  }

  # add attributes
  attr(new_formula, 'label_info') <- label_info
  attr(new_formula, 'org_formula') <- formula

  # export new formula
  new_formula
}
