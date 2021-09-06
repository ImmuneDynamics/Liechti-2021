##########################################################################################################
#### 1. Load packages, and set working directory
##########################################################################################################

    ### Load libraries
    
        library(Spectre)
        Spectre::package.check()    # Check that all required packages are installed
        Spectre::package.load()     # Load required packages
        
    ### Set PrimaryDirectory
    
        dirname(rstudioapi::getActiveDocumentContext()$path)            # Finds the directory where this script is located
        setwd(dirname(rstudioapi::getActiveDocumentContext()$path))     # Sets the working directory to where the script is located
        getwd()
        PrimaryDirectory <- getwd()
        PrimaryDirectory
    
    ### Set 'input' directory
    
        setwd(PrimaryDirectory)
        setwd('data')
        InputDirectory <- getwd()
        InputDirectory
        setwd(PrimaryDirectory)
    
    ### Create output directory
        
        dir.create("Output_Spectre", showWarnings = FALSE)
        setwd("Output_Spectre")
        OutputDirectory <- getwd()
        setwd(PrimaryDirectory)

##########################################################################################################
#### 2. Import and prep data
##########################################################################################################
    
    ### Import data
    
        setwd(InputDirectory)
        cell.dat <- fread('data.csv')
        cell.dat
        
    ### Specify columns
                         
        as.matrix(names(cell.dat))
        cellular.cols <- names(cell.dat)[c(9:16)]
        cellular.cols
        
        as.matrix(names(cell.dat))
        cluster.cols <- names(cell.dat)[c(9:16)]
        cluster.cols

##########################################################################################################
#### xxx
##########################################################################################################
        
    ### Adjust clustering columns
        
        as.matrix(cluster.cols)
        
        cluster.cols <- cluster.cols[-4]
        as.matrix(cluster.cols)
        
    ### Arcsinh transformation
        
        cell.dat <- do.asinh(cell.dat, 'PECy5-5 CD3e', cofactor = 5000, append.cf = TRUE)
        cell.dat <- do.asinh(cell.dat, 'PECy5-5 CD3e', cofactor = 500, append.cf = TRUE)
        cell.dat <- do.asinh(cell.dat, 'PECy5-5 CD3e', cofactor = 50, append.cf = TRUE)
        cell.dat <- do.asinh(cell.dat, 'PECy5-5 CD3e', cofactor = 5, append.cf = TRUE)
        
        cell.dat <- run.fitsne(cell.dat, c(cluster.cols, 'PECy5-5 CD3e_asinh_cf5000'), fitsne.x.name = 'FItSNE_X_CD3_5000', fitsne.y.name = 'FItSNE_Y_CD3_5000')
        cell.dat <- run.fitsne(cell.dat, c(cluster.cols, 'PECy5-5 CD3e_asinh_cf500'), fitsne.x.name = 'FItSNE_X_CD3_500', fitsne.y.name = 'FItSNE_Y_CD3_500')
        cell.dat <- run.fitsne(cell.dat, c(cluster.cols, 'PECy5-5 CD3e_asinh_cf50'), fitsne.x.name = 'FItSNE_X_CD3_50', fitsne.y.name = 'FItSNE_Y_CD3_50')
        cell.dat <- run.fitsne(cell.dat, c(cluster.cols, 'PECy5-5 CD3e_asinh_cf5'), fitsne.x.name = 'FItSNE_X_CD3_5', fitsne.y.name = 'FItSNE_Y_CD3_5')
        cell.dat <- run.fitsne(cell.dat, c(cluster.cols, 'PECy5-5 CD3e'), fitsne.x.name = 'FItSNE_X_CD3_lin', fitsne.y.name = 'FItSNE_Y_CD3_lin')
        
    ### FItSNE plots
         
        setwd(OutputDirectory)
        dir.create("FItSNE plots")
        setwd("FItSNE plots")
        
        make.colour.plot(cell.dat, 'FItSNE_X_CD3_5000', 'FItSNE_Y_CD3_5000', 'Population')
        make.colour.plot(cell.dat, 'FItSNE_X_CD3_500', 'FItSNE_Y_CD3_500', 'Population')
        make.colour.plot(cell.dat, 'FItSNE_X_CD3_50', 'FItSNE_Y_CD3_50', 'Population')
        make.colour.plot(cell.dat, 'FItSNE_X_CD3_5', 'FItSNE_Y_CD3_5', 'Population')
        make.colour.plot(cell.dat, 'FItSNE_X_CD3_lin', 'FItSNE_Y_CD3_lin', 'Population')
        
        make.colour.plot(cell.dat, 'FItSNE_X_CD3_5000', 'FItSNE_Y_CD3_5000', 'PECy5-5 CD3e_asinh_cf5000')
        make.colour.plot(cell.dat, 'FItSNE_X_CD3_500', 'FItSNE_Y_CD3_500', 'PECy5-5 CD3e_asinh_cf500')
        make.colour.plot(cell.dat, 'FItSNE_X_CD3_50', 'FItSNE_Y_CD3_50', 'PECy5-5 CD3e_asinh_cf50')
        make.colour.plot(cell.dat, 'FItSNE_X_CD3_5', 'FItSNE_Y_CD3_5', 'PECy5-5 CD3e_asinh_cf5')
        
    ### Axis plots
        
        setwd(OutputDirectory)
        dir.create("Axis plots")
        setwd("Axis plots")
        
        make.colour.plot(cell.dat, 'PECy5-5 CD3e_asinh_cf5000', 'BV605 Ly6C_asinh', 'Population', plot.width = 7, plot.height = 6)
        make.colour.plot(cell.dat, 'PECy5-5 CD3e_asinh_cf500', 'BV605 Ly6C_asinh', 'Population', plot.width = 7, plot.height = 6)
        make.colour.plot(cell.dat, 'PECy5-5 CD3e_asinh_cf50', 'BV605 Ly6C_asinh', 'Population', plot.width = 7, plot.height = 6)
        make.colour.plot(cell.dat, 'PECy5-5 CD3e_asinh_cf5', 'BV605 Ly6C_asinh', 'Population', plot.width = 7, plot.height = 6)  
        
        make.colour.plot(cell.dat, 'PECy5-5 CD3e_asinh_cf5000', 'BV605 Ly6C_asinh', plot.width = 6, plot.height = 5)
        make.colour.plot(cell.dat, 'PECy5-5 CD3e_asinh_cf500', 'BV605 Ly6C_asinh', plot.width = 6, plot.height = 5)
        make.colour.plot(cell.dat, 'PECy5-5 CD3e_asinh_cf50', 'BV605 Ly6C_asinh', plot.width = 6, plot.height = 5)
        make.colour.plot(cell.dat, 'PECy5-5 CD3e_asinh_cf5', 'BV605 Ly6C_asinh', plot.width = 6, plot.height = 5)  
        