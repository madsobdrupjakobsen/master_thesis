files <- dir("~/r/functions", pattern="/*", full.names=TRUE)
## If a subfolder for the current the working directory names "functions", then source also all files in that
if(length(dir("./functions")) > 0) files <- c(files,dir("./functions", pattern="/*", full.names=TRUE))

for(i in 1:length(files))
  {
    ## Only source it if it doens not end with a ~
    if(length(grep("~|#", files[i])) == 0)
      {
        source(files[i])
        print(paste("Sourced:",files[i]))
      }
  }
