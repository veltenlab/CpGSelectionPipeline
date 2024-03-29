library(RnBeads)
theme_set(theme_bw())
xml.file <- "rnbOptions.xml"
arch <- new("ClusterArchitectureSGE")
arch <- setExecutable(arch,"R","/usr/bin/R")
arch <- setExecutable(arch,"Rscript","/usr/bin/Rscript")
rnb.cr <- new("RnBClusterRun",arch)
rnb.cr <- setModuleResourceRequirements(rnb.cr,c(h_vmem="120G",virtual_free="120G",h_rt="05:00:00"),"all")
rnb.cr <- setModuleResourceRequirements(rnb.cr,c(h_vmem="120G",virtual_free="120G",h_rt="05:00:00", queue='long-sl7'),"differential")
rnb.cr <- setModuleNumCores(rnb.cr,1L,"all")
rnb.cr <- setModuleNumCores(rnb.cr,1L,"exploratory")
run(rnb.cr, "rnB_Kulis", xml.file, queue="short-sl7")
