library(shiny)
library(igraph)
library(corrplot)

# To Do: read the labels from a file

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  
  output$pca_plot <- renderPlot({
    
    labels <- c("United Kingdom", "Spain", "Ireland", "Portugal", "Hungary", "Italy", "France", "Albania", 
                "Macedonia", "Belgium", "Switzerland", "Germany", "Croatia", "Bosnia", "Poland", "Czech Republic",
                "Serbia", "Netherlands", "Austria")
    
    if (input$isFinite == 'infinite') {
      upper_thr = 'Inf'
    } else {
      upper_thr = input$upper_thr
    }
    
    sims <- as.matrix(read.table(paste0("data/", input$datatype, "/popres_", input$lower_thr, "_", upper_thr, "_mean_sims.txt")))
    rownames(sims) <- labels
    colnames(sims) <- labels
    
    popres.pca <- prcomp(sims, center = TRUE, scale. = TRUE)
    plot(popres.pca$x[,1], popres.pca$x[,2], type = "n", xlab = "PC1", ylab = "PC2")
    text(popres.pca$x[,1], popres.pca$x[,2], labels = labels, cex=.7)
    
})
  
  output$sim_plot <- renderPlot({
    
    labels <- c("United Kingdom", "Spain", "Ireland", "Portugal", "Hungary", "Italy", "France", "Albania", 
                "Macedonia", "Belgium", "Switzerland", "Germany", "Croatia", "Bosnia", "Poland", "Czech Republic",
                "Serbia", "Netherlands", "Austria")
    
    if (input$isFinite == 'infinite') {
      upper_thr = 'Inf'
    } else {
      upper_thr = input$upper_thr
    }
    
    sims <- as.matrix(read.table(paste0("data/", input$datatype, "/popres_", input$lower_thr, "_", upper_thr, "_mean_sims.txt")))
    rownames(sims) <- labels
    colnames(sims) <- labels
    
    min2 <- min(sims[sims!=min(sims)])/2;
    corr_sims <- sims+min2
    is.na(corr_sims) <- do.call(cbind,lapply(corr_sims, is.infinite))
    corrplot(corr = log(corr_sims)+max(abs(log(corr_sims))), diag=TRUE, order= "hclust", is.corr = FALSE, addrect = 3, method = "circle",cl.pos = "n")

  })
  
  output$graph_plot <- renderPlot({
    
    labels <- c("United Kingdom", "Spain", "Ireland", "Portugal", "Hungary", "Italy", "France", "Albania", 
                "Macedonia", "Belgium", "Switzerland", "Germany", "Croatia", "Bosnia", "Poland", "Czech Republic",
                "Serbia", "Netherlands", "Austria")
    
    if (input$isFinite == 'infinite') {
      upper_thr = 'Inf'
    } else {
      upper_thr = input$upper_thr
    }
    
    sims <- as.matrix(read.table(paste0("data/", input$datatype, "/popres_", input$lower_thr, "_", upper_thr, "_mean_sims.txt")))
    rownames(sims) <- labels
    colnames(sims) <- labels
    
    net=graph.adjacency(sims,mode="undirected",weighted=TRUE,diag=FALSE)
    node.size <- 50*diag(sims)
    plot.igraph(net,vertex.label=V(net)$name,layout=layout.fruchterman.reingold, edge.color="black",edge.width=5*E(net)$weight,
                vertex.size=node.size)
  })
  
  output$maps_plot_migration <- renderImage({
    
    if (input$isFinite == 'infinite') {
      upper_thr = 'Inf'
    } else {
      upper_thr = input$upper_thr
    }
    file_name <- paste0("data/", input$datatype, "/plots_lowerBnd_", input$lower_thr, "_upperBnd_", upper_thr, "-mrates.png")
    
    list(src = file_name,
         contentType = 'image/png', width = 600, height = 480)
    
  }, deleteFile = FALSE)
    
  output$maps_plot_population <- renderImage({
    if (input$isFinite == 'infinite') {
      upper_thr = 'Inf'
    } else {
      upper_thr = input$upper_thr
    }
    file_name <- paste0("data/", input$datatype, "/plots_lowerBnd_", input$lower_thr, "_upperBnd_", upper_thr, "-popsizes.png")
    
    list(src = file_name,
         contentType = 'image/png', width = 600, height = 480)
    
  }, deleteFile = FALSE)
  
})