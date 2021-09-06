##########################################################################################################
#### 1. Load packages, and set working directory
##########################################################################################################

    ### Load libraries

        library(Spectre)
        Spectre::package.check()    # Check that all required packages are installed
        Spectre::package.load()     # Load required packages

    ### Set PrimaryDirectory
        
        dirname(rstudioapi::getActiveDocumentContext()$path)
        setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
        getwd()
        PrimaryDirectory <- getwd()
        PrimaryDirectory
        
    ### Create output directory
        
        setwd(PrimaryDirectory)
        dir.create("Output_Spectre", showWarnings = FALSE)
        setwd("Output_Spectre")
        OutputDirectory <- getwd()
        setwd(PrimaryDirectory)

##########################################################################################################
#### 2. Import and prep data
##########################################################################################################

    ### Read in data
        
        setwd(PrimaryDirectory)
        setwd('data/')
        
        files <- list.files(getwd(), '.csv')
        files
        
        dat.list <- list()
        
        for(i in files){
            a <- gsub('.csv', '', i)
            dat.list[[a]] <- fread(i)
        }
        
        dat.list
    
##########################################################################################################
#### 3. Fixed Ly6C transformation
##########################################################################################################
        
        names(dat.list)
        
        dat.list[[1]] <- do.asinh(dat.list[[1]], 'Ly6C', cofactor = 5000)
        dat.list[[2]] <- do.asinh(dat.list[[2]], 'Ly6C', cofactor = 15)
        dat.list[[3]] <- do.asinh(dat.list[[3]], 'Ly6C', cofactor = 1000)
        
##########################################################################################################
#### 4. Transformations
##########################################################################################################

    for(i in names(dat.list)){
        # i <- 'aurora'
        
        setwd(OutputDirectory)
        dir.create(i)
        setwd(i)
        
        for(a in c(5, 15, 25, 50, 250, 500, 2500, 5000, 25000)){
            dat <- do.asinh(dat.list[[i]], 'B220', cofactor = a)
            make.colour.plot(dat, 'B220_asinh', 'Ly6C_asinh', filename = paste0(a, '.png'))
        }
    }
        
        