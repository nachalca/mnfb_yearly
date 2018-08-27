cleanFolder <- function(dir) {
  answer <- NA
  while(!(answer %in% c('y', 'n'))) {
    answer <- readline(paste("Clean ", dir ,"? (y/n)"))
  }
  if(answer == 'y') {
    rules <- c('.log', '.vrb', '.nav', '.snm', '.toc', '.nav',
               '-tikzDictionary', '.tex', '.synctex.gz')
    
    ll <- list.files(path=dir, pattern = paste0('\\',rules ,'$', collapse = '|'))
    file.remove( paste(dir, ll, sep='/') )
  }
}