#'@importFrom methods setClassUnion
setClassUnion('characterOrnumericOrNULL', members = c('character', 'numeric', 'NULL'))
setClassUnion('numericOrNULL', members = c('numeric', 'NULL'))
setClassUnion('integerOrNULL', members = c('integer', 'NULL'))
