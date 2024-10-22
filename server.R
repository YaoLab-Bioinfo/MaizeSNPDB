
# options(shiny.maxRequestSize = 200*1024^2)

shinyServer(function(input, output, session) {
  
  # GBrowser
  observe({
    if (input$submit1>0) {
      if ( exists("GBrowser") ){
      } else {
        source("GBrowser.R")
      }
      
      if ( exists("anaReg") ){
      } else {
        source("anaReg.R")
      }
      isolate({
        myPos <- anaReg(input$regB)
        if ( exists("snpInfo") ){
        } else {
          source("snpInfo.R")
        }
        
        if ( exists("validReg") ){
        } else {
          source("validReg.R")
        }
        if (validReg(myPos)) {
          if (!is.null(myPos)) {
            snp.info <- snpInfo(chr=myPos$chr, start=myPos$start - input$GBUP, end=myPos$end + input$GBDOWN, 
                                accession = input$mychooserB$selected, mutType = input$GB_mut_group)
          } else {
            snp.info <- NULL
          }
          
          if (is.null(snp.info) || nrow(snp.info[[1]][[1]]) < 1) {
            sendSweetAlert(
              session = session,
              title = "Error input!", type = "error",
              text = "No SNPs are detected in the specified genomic region or the specified genomic region is too large!"
            )
          } else {
            GBplot <<- NULL
            output$gbrowser <- plotly::renderPlotly({
              GBplot <<- GBrowser(chr=myPos$chr, start=myPos$start - input$GBUP, 
                                  end=myPos$end + input$GBDOWN,
                                  accession = input$mychooserB$selected,
                                  mutType = input$GB_mut_group)
              GBplot[[2]]
            })
            
            ## Download PDF file of GBrowser
            output$downloadGB.pdf <- downloadHandler(
              filename <- function() { paste('GBrowser.pdf') },
              content <- function(file) {
                pdf(file, width = 900/72, height = 300/72)
                grid::grid.draw(GBplot[[1]])
                dev.off()
              }, contentType = 'application/pdf')
            
            # Download genotypes of seleceted SNPs
            output$downloadsnp.txt <- downloadHandler(
              filename = function() { "snp.geno.txt" },
              content = function(file) {
                write.table(snp.info[[1]][[1]], file, sep="\t", quote=F)
              })
            
            # Download information of SNPs
            output$downloadsnpInfo.txt <- downloadHandler(
              filename = function() { "snp.info.txt" },
              content = function(file) {
                write.table(snp.info[[2]], file, sep="\t", quote=F, row.names=F)
              })
          }
        } else {
          sendSweetAlert(
            session = session,
            title = "Error input!", type = "error",
            text = "Please input genomic region or gene model in appropriate format!"
          )
        }
        
      })
    } else {
      NULL
    }
  })
  
  observe({
    if (input$clearGB>0) {
      isolate({
        updateTextInput(session, "regB", value="")
      })
    } else {NULL}
  })
  
  observe({
    if (input$GBExam >0) {
      isolate({
        updateTextInput(session, "regB", value="chr1:29765419-29793053")
      })
    } else {NULL}
  })
  
  # LDheatmap
  observe({
    if (input$submit2>0) {
      if ( exists("anaReg") ){
      } else {
        source("anaReg.R")
      }
      if ( exists("ld.heatmap") ){
      } else {
        source("ld.heatmap.R")
      }
      
      if ( exists("fetchSnp") ){
      } else {
        source("fetchSnp.R")
      }
      isolate({
        ld.height <<- input$ldHeight
        ld.width <<- input$ldWidth
        myPos <- anaReg(input$regL)
        
        if (validReg(myPos)) {
          if (!is.null(myPos)) {
            snp.reg <- fetchSnp(chr=myPos$chr, start=myPos$start - input$ldUp * 1000, 
                                end=myPos$end + input$ldDown * 1000, accession = input$mychooserLD$selected,
                                mutType = input$ld_mut_group)[[1]]
          } else {
            snp.reg <- NULL
          }
          
          if (is.null(snp.reg) || nrow(snp.reg) < 5) {
            sendSweetAlert(
              session = session,
              title = "Error input!", type = "error",
              text = "Too few SNPs are detected in the specified genomic region or the specified genomic region is too large!"
            )
          } else {
            snp.pos <- as.numeric(unlist(strsplit(input$ldpos, split=",")))
            
            if (input$uploadLD == 1) {
              ld.snp.site <- NULL
            } else {
              if (!is.null(input$LD.snpsite)) {
                ld.snp.site <- readLines(input$LD.snpsite$datapath)
              } else {
                ld.snp.site <- NULL
              }
            }
            
            output$ldheatmap <- renderPlot({
              if (input$flip == "0") {
                ld.heatmap(chr=myPos$chr, start=myPos$start - input$ldUp * 1000, end=myPos$end + input$ldDown * 1000, text=c(FALSE, TRUE)[as.numeric(input$showText)+1],
                           snp.pos=snp.pos, gene=FALSE, flip=c(FALSE, TRUE)[as.numeric(input$flip)+1],
                           col=list(grey.colors(20), heat.colors(20))[[as.numeric(input$ldcol)]],
                           mutType = input$ld_mut_group, accession = input$mychooserLD$selected, 
                           snpSites = ld.snp.site)
              } else if (input$flip == "1") {
                if (input$LDshowGene) {
                  ld.heatmap(chr=myPos$chr, start=myPos$start - input$ldUp * 1000, end=myPos$end + input$ldDown * 1000, text=FALSE,
                             snp.pos=snp.pos, ld.y=input$ldY/100, ld.w=input$ldW/100, gene=TRUE, 
                             col=list(grey.colors(20), heat.colors(20))[[as.numeric(input$ldcol)]],
                             mutType = input$ld_mut_group, accession = input$mychooserLD$selected, 
                             snpSites = ld.snp.site)
                } else {
                  ld.heatmap(chr=myPos$chr, start=myPos$start - input$ldUp * 1000, end=myPos$end + input$ldDown * 1000, text=FALSE,
                             gene=FALSE, flip=c(FALSE, TRUE)[as.numeric(input$flip)+1],
                             col=list(grey.colors(20), heat.colors(20))[[as.numeric(input$ldcol)]],
                             mutType = input$ld_mut_group, accession = input$mychooserLD$selected, 
                             snpSites = ld.snp.site)
                }
              }
              
            }, height = ld.height, width = ld.width)
          }
        } else {
          sendSweetAlert(
            session = session,
            title = "Error input!", type = "error",
            text = "Please input genomic region or gene model in appropriate format!"
          )
        }
        
      })
    } else {
      NULL
    }
  })
  
  observe({
    if (input$clearLD>0) {
      isolate({
        updateTextInput(session, "regL", value="")
      })
    } else {NULL}
  })
  
  observe({
    if (input$LDExam >0) {
      isolate({
        updateTextInput(session, "regL", value="Zm00001d033673")
      })
    } else {NULL}
  })

	## Download PDF file of LDheatmap
	output$downloadLD.pdf <- downloadHandler(
	  filename <- function() { paste('LDheatmap.pdf') },
	  content <- function(file) {
	    myPos <- anaReg(input$regL)
	    
	    snp.reg <- fetchSnp(chr=myPos$chr, start=myPos$start - input$ldUp * 1000, 
	                        end=myPos$end + input$ldDown * 1000, accession = input$mychooserLD$selected,
	                        mutType = input$ld_mut_group)[[1]]
	    if (nrow(snp.reg) < 5) {
	      js_string <- 'alert("Too few SNPs in specified genomic region!");'
	      session$sendCustomMessage(type='jsCode', list(value = js_string))
	    } else {
	      pdf(file, width = input$ldWidth/72, height = input$ldHeight/72, onefile = FALSE)
	      
	      snp.pos <- as.numeric(unlist(strsplit(input$ldpos, split=",")))
	      
	      if (input$uploadLD == 1) {
	        ld.snp.site <- NULL
	      } else {
	        if (!is.null(input$LD.snpsite)) {
	          ld.snp.site <- readLines(input$LD.snpsite$datapath)
	        } else {
	          ld.snp.site <- NULL
	        }
	      }
	      
	      if (input$flip == "0") {
	        ld.heatmap(chr=myPos$chr, start=myPos$start - input$ldUp * 1000, end=myPos$end + input$ldDown * 1000, text=c(FALSE, TRUE)[as.numeric(input$showText)+1],
	                   snp.pos=snp.pos, gene=FALSE, flip=c(FALSE, TRUE)[as.numeric(input$flip)+1],
	                   col=list(grey.colors(20), heat.colors(20))[[as.numeric(input$ldcol)]],
	                   mutType = input$ld_mut_group, accession = input$mychooserLD$selected, 
	                   snpSites = ld.snp.site)
	      } else if (input$flip == "1") {
	        if (input$LDshowGene) {
	          ld.heatmap(chr=myPos$chr, start=myPos$start - input$ldUp * 1000, end=myPos$end + input$ldDown * 1000, text=FALSE,
	                     snp.pos=snp.pos, ld.y=input$ldY/100, ld.w=input$ldW/100, gene=TRUE, 
	                     col=list(grey.colors(20), heat.colors(20))[[as.numeric(input$ldcol)]],
	                     mutType = input$ld_mut_group, accession = input$mychooserLD$selected, 
	                     snpSites = ld.snp.site)
	        } else {
	          ld.heatmap(chr=myPos$chr, start=myPos$start - input$ldUp * 1000, end=myPos$end + input$ldDown * 1000, text=FALSE,
	                     gene=FALSE, flip=c(FALSE, TRUE)[as.numeric(input$flip)+1],
	                     col=list(grey.colors(20), heat.colors(20))[[as.numeric(input$ldcol)]],
	                     mutType = input$ld_mut_group, accession = input$mychooserLD$selected, 
	                     snpSites = ld.snp.site)
	        }
	      }
	      dev.off()
	    }
	    
	  }, contentType = 'application/pdf')
	
	## Download SVG file of LDheatmap
	output$downloadLD.svg <- downloadHandler(
	  filename <- function() { paste('LDheatmap.svg') },
	  content <- function(file) {
	    myPos <- anaReg(input$regL)
	    snp.reg <- fetchSnp(chr=myPos$chr, start=myPos$start - input$ldUp * 1000, 
	                        end=myPos$end + input$ldDown * 1000, accession = input$mychooserLD$selected,
	                        mutType = input$ld_mut_group)[[1]]
	    if (nrow(snp.reg) < 5) {
	      js_string <- 'alert("Too few SNPs in specified genomic region!");'
	      session$sendCustomMessage(type='jsCode', list(value = js_string))
	    } else {
	      svg(file, width = input$ldWidth/72, height = input$ldHeight/72)
	      
	      snp.pos <- as.numeric(unlist(strsplit(input$ldpos, split=",")))
	      
	      if (input$uploadLD == 1) {
	        ld.snp.site <- NULL
	      } else {
	        if (!is.null(input$LD.snpsite)) {
	          ld.snp.site <- readLines(input$LD.snpsite$datapath)
	        } else {
	          ld.snp.site <- NULL
	        }
	      }
	      
	      if (input$flip == "0") {
	        ld.heatmap(chr=myPos$chr, start=myPos$start - input$ldUp * 1000, end=myPos$end + input$ldDown * 1000, text=c(FALSE, TRUE)[as.numeric(input$showText)+1],
	                   snp.pos=snp.pos, gene=FALSE, flip=c(FALSE, TRUE)[as.numeric(input$flip)+1],
	                   col=list(grey.colors(20), heat.colors(20))[[as.numeric(input$ldcol)]],
	                   mutType = input$ld_mut_group, accession = input$mychooserLD$selected, 
	                   snpSites = ld.snp.site)
	      } else if (input$flip == "1") {
	        if (input$LDshowGene) {
	          ld.heatmap(chr=myPos$chr, start=myPos$start - input$ldUp * 1000, end=myPos$end + input$ldDown * 1000, text=FALSE,
	                     snp.pos=snp.pos, ld.y=input$ldY/100, ld.w=input$ldW/100, gene=TRUE, 
	                     col=list(grey.colors(20), heat.colors(20))[[as.numeric(input$ldcol)]],
	                     mutType = input$ld_mut_group, accession = input$mychooserLD$selected, 
	                     snpSites = ld.snp.site)
	        } else {
	          ld.heatmap(chr=myPos$chr, start=myPos$start - input$ldUp * 1000, end=myPos$end + input$ldDown * 1000, text=FALSE,
	                     gene=FALSE, flip=c(FALSE, TRUE)[as.numeric(input$flip)+1],
	                     col=list(grey.colors(20), heat.colors(20))[[as.numeric(input$ldcol)]],
	                     mutType = input$ld_mut_group, accession = input$mychooserLD$selected, 
	                     snpSites = ld.snp.site)
	        }
	      }
	      dev.off()
	    }
	    
	  }, contentType = 'image/svg')
	
	
	# Diversity
	observe({
	  if (input$submit4>0) {
	    if ( exists("anaReg") ){
	    } else {
	      source("anaReg.R")
	    }
	    
	    if ( exists("nucDiv") ){
	    } else {
	      source("nucDiv.R")
	    }
	    
	    if ( exists("fetchSnp") ){
	    } else {
	      source("fetchSnp.R")
	    }
	    
	    if ( exists("validReg") ){
	    } else {
	      source("validReg.R")
	    }
	    
	    isolate({
	      div.height <<- input$divHeight
	      div.width <<- input$divWidth
	      
	      myPos <- anaReg(input$regD)
	      
	      if (validReg(myPos)) {
	        div.up <- input$divUp * 1000
	        div.down <- input$divDown * 1000
	        div.group <- input$div_acc_group
	        div.step <- input$snpnumD
	        div.numerator <- input$nuc_numerator
	        div.denominator <- input$nuc_denominator
	        div.mut.group <- input$div_mut_group
	        
	        if (input$uploadDIV == 1) {
	          div.snp.site <- NULL
	        } else {
	          if (!is.null(input$DIV.snpsite)) {
	            div.snp.site <- readLines(input$DIV.snpsite$datapath)
	          } else {
	            div.snp.site <- NULL
	          }
	        }
	        
	        if (!is.null(myPos)) {
	          snp.reg <- fetchSnp(chr=myPos$chr, start=myPos$start - div.up, end=myPos$end + div.down,
	                              mutType=input$div_mut_group)[[1]]
	        } else {
	          snp.reg <- NULL
	        }
	        
	        if (is.null(snp.reg) || nrow(snp.reg) < 10) {
	          sendSweetAlert(
	            session = session,
	            title = "Error input!", type = "error",
	            text = "Too few SNPs are detected in the specified genomic region or the specified genomic region is too large!"
	          )
	        } else {
	          nuc.div.plot <<- NULL
	          output$diversity <- renderPlot({
	            nuc.div.plot <<- nucDiv(chr=myPos$chr, nuc.start=myPos$start - div.up, nuc.end=myPos$end + div.down, 
	                                    groups = div.group, step = div.step,
	                                    numerator = div.numerator, denominator = div.denominator, 
	                                    mutType = div.mut.group, snpSites = div.snp.site)
	            grid::grid.draw(grid.arrange(nuc.div.plot[[1]], nuc.div.plot[[2]], ncol=1, heights=c(2.3, 1)))
	          }, height = div.height, width = div.width)
	          
	          ## Download PDF file of Diversity
	          output$downloadDiv01 <- renderUI({
	            req(input$submit4, nuc.div.plot)
	            downloadButton("downloadDiv.pdf", "Download pdf-file")
	          })
	          
	          output$downloadDiv.pdf <- downloadHandler(
	            filename <- function() { paste('diversity.pdf') },
	            content <- function(file) {
	              pdf(file, width = input$divWidth/72, height = input$divHeight/72)
	              grid::grid.draw(grid.arrange(nuc.div.plot[[1]], nuc.div.plot[[2]], ncol=1, heights=c(2.3, 1)))
	              
	              dev.off()
	            }, contentType = 'application/pdf')
	          
	          ## Download SVG file of Diversity
	          output$downloadDiv02 <- renderUI({
	            req(input$submit4, nuc.div.plot)
	            downloadButton("downloadDiv.svg", "Download svg-file")
	          })
	          
	          output$downloadDiv.svg <- downloadHandler(
	            filename <- function() { paste('diversity.svg') },
	            content <- function(file) {
	              svg(file, width = input$divWidth/72, height = input$divHeight/72)
	              grid::grid.draw(grid.arrange(nuc.div.plot[[1]], nuc.div.plot[[2]], ncol=1, heights=c(2.3, 1)))
	              
	              dev.off()
	            }, contentType = 'image/svg')
	          
	          ## Download TXT file of diversity
	          output$downloadDiv03 <- renderUI({
	            req(input$submit4)
	            downloadButton("downloadDiv.txt", "Download TXT-file")
	          })
	          
	          output$downloadDiv.txt <- downloadHandler(
	            filename <- function() { paste('diversity.txt') },
	            content <- function(file) {
	              write.table(diVTxt, file, sep="\t", quote=F, row.names = F)
	            }, contentType = 'text/plain')
	        }
	      } else {
	        sendSweetAlert(
	          session = session,
	          title = "Error input!", type = "error",
	          text = "Please input genomic region or gene model in appropriate format!"
	        )
	      }
	    })
	  } else {
	    NULL
	  }
	})
	
	observe({
	  if (input$clearDiv>0) {
	    isolate({
	      updateTextInput(session, "regD", value="")
	    })
	  } else {NULL}
	})
	
	observe({
	  if (input$DivExam >0) {
	    isolate({
	      updateTextInput(session, "regD", value="Zm00001d033673")
	    })
	  } else {NULL}
	})
	
	
	# phylogenetics
	observe({
	  if (input$submit5>0) {
	    if ( exists("anaReg") ){
	    } else {
	      source("anaReg.R")
	    }
	    
	    if ( exists("phylo") ){
	    } else {
	      source("phylo.R")
	    }
	    isolate({
	      phy.height <<- input$phyHeight
	      phy.width <<- input$phyWidth
	      phy.up <- input$phyUp * 1000
	      phy.down <- input$phyDown * 1000
	      
	      myPos <- anaReg(input$regP)
	      
	      if (validReg(myPos)) {
	        phy.acc <- input$mychooserPhy$selected
	        phy.mut.group <- input$phy_mut_group
	        
	        if (input$uploadPHY == 1) {
	          phy.snp.site <- NULL
	        } else {
	          if (!is.null(input$PHY.snpsite)) {
	            phy.snp.site <- readLines(input$PHY.snpsite$datapath)
	          } else {
	            phy.snp.site <- NULL
	          }
	        }
	        
	        if (!is.null(myPos)) {
	          snp.reg <- fetchSnp(chr=myPos$chr, start=myPos$start - phy.up, end=myPos$end + phy.down,
	                              accession=phy.acc, mutType=phy.mut.group)[[1]]
	        } else {
	          snp.reg <- NULL
	        }
	        
	        if (is.null(snp.reg) || nrow(snp.reg) < 10) {
	          js_string <- 'alert("No SNPs are detected in the specified genomic region or the specified genomic region is too large!");'
	          session$sendCustomMessage(type='jsCode', list(value = js_string))
	        } else {
	          output$phylo <- renderPlot({
	            phylo(chr=myPos$chr, start=myPos$start-phy.up, end=myPos$end+phy.down,
	                  accession=phy.acc, mutType=phy.mut.group, snpSites = phy.snp.site)
	          }, height = phy.height, width = phy.width)
	          
	          ## Download PDF file of phylogenetics
	          output$downloadPhy01 <- renderUI({
	            req(input$submit5)
	            downloadButton("downloadPhylo.pdf", "Download pdf-file")
	          })
	          
	          output$downloadPhylo.pdf <- downloadHandler(
	            filename <- function() { paste('phylogenetics.pdf') },
	            content <- function(file) {
	              pdf(file, width = input$phyWidth/72, height = input$phyHeight/72)
	              print(figurecp)
	              dev.off()
	            }, contentType = 'application/pdf')
	          
	          ## Download NWK file of phylogenetics
	          output$downloadPhy02 <- renderUI({
	            req(input$submit5)
	            downloadButton("downloadPhylo.nwk", "Download Newick-file")
	          })
	          
	          output$downloadPhylo.nwk <- downloadHandler(
	            filename <- function() { paste('phylogenetics.nwk') },
	            content <- function(file) {
	              ape::write.tree(treNwk, file)
	            }, contentType = 'text/plain')
	        }
	      } else {
	        sendSweetAlert(
	          session = session,
	          title = "Error input!", type = "error",
	          text = "Please input genomic region or gene model in appropriate format!"
	        )
	      }

	    })
	  } else {
	    NULL
	  }
	})
	
	observe({
	  if (input$clearPhy>0) {
	    isolate({
	      updateTextInput(session, "regP", value="")
	    })
	  } else {NULL}
	})
	
	observe({
	  if (input$PhyExam >0) {
	    isolate({
	      updateTextInput(session, "regP", value="Zm00001d033673")
	    })
	  } else {NULL}
	})
	
	
	# accession information
	output$acc.info.txt <- downloadHandler(
	  filename = function() { "acc.info.txt" },
	  content = function(file) {
	    if (exists("acc.info")){
	    }else{
	      acc.info <- read.table("./data/all.acc.txt", head=T, as.is=T, sep="\t", quote="")
	    }
	    write.table(acc.info, file, sep = "\t", quote=FALSE, row.names = FALSE)
	}, contentType = 'text/plain')
	
	
	output$sel.acc.info.txt <- downloadHandler(
	  filename = function() { "sel.acc.info.txt" },
	  content = function(file) {
	    accession <- input$mychooserA$selected
	    accession <- gsub(",.+", "", accession)
	    accession <- sapply(accession, function(x){
	      if (x %in% c("Improved", "Landrace", "Parviglumis")) {
	        x.dat <- readLines(paste0("./data/", x, ".acc.txt"))
	        return(x.dat)
	      } else {
	        return(x)
	      }
	    })
	    accession <- unique(unlist(accession))
	    if (exists("acc.info")){
	    }else{
	      acc.info <- read.table("./data/all.acc.txt", head=T, as.is=T, sep="\t", quote="")
	    }
	    write.table(acc.info[acc.info$ID %in% accession, ], 
	                file, sep = "\t", quote=FALSE, row.names = FALSE)
	  }, contentType = 'text/plain')

	
	output$mytable1 = renderDataTable({
	  accession <- input$mychooserA$selected
	  accession <- gsub(",.+", "", accession)
	  accession <- sapply(accession, function(x){
	    if (x %in% c("Improved", "Landrace", "Parviglumis")) {
	      x.dat <- readLines(paste0("./data/", x, ".acc.txt"))
	      return(x.dat)
	    } else {
	      return(x)
	    }
	  })
	  accession <- unique(unlist(accession))
	  if (exists("acc.info")){
	  }else{
	    acc.info <- read.table("./data/all.acc.txt", head=T, as.is=T, sep="\t", quote="")
	  }
	  acc.info[acc.info$ID %in% accession, ]
	}, options = list(lengthMenu = c(5, 8, 10), pageLength = 10, searching = TRUE, autoWidth = FALSE), escape = FALSE
	)
	
	output$mytable3 = renderDataTable({
	  geneID <- read.table("maize.v3TOv4.geneIDhistory.txt", head=T, as.is=T)
	  
	  geneID
	}, options = list(lengthMenu = c(10, 20, 30, 40), pageLength = 15, searching = TRUE, autoWidth = FALSE), escape = FALSE
	)
	
	# Bulk download genotypes of seleceted SNPs
	observe({
	  if (input$submit6>0) {
	    if ( exists("anaReg") ){
	    } else {
	      source("anaReg.R")
	    }
	    
	    if ( exists("validReg") ){
	    } else {
	      source("validReg.R")
	    }
	    isolate({
	      myPos <- anaReg(input$regBB)
	      
	      if (validReg(myPos)) {
	        snp.info.down <<- NULL
	        
	        if (!is.null(myPos)) {
	          output$mytable2 = renderDataTable({
	            snp.info.down <<- snpInfo(chr=myPos$chr, start=myPos$start, end=myPos$end, 
	                                      accession = input$mychooserD$selected, mutType = input$down_mut_group)
	            snp.info.down[[2]]
	          }, options = list(lengthMenu = c(5, 8, 10), pageLength = 5, searching = TRUE, autoWidth = TRUE), escape = FALSE
	          )
	          
	          output$bulkdownloadsnp.txt <- downloadHandler(
	            filename = function() { "down.snp.geno.txt" },
	            content = function(file) {
	              write.table(snp.info.down[[1]][[1]], file, sep="\t", quote=F)
	            })
	          
	          # Bulk download information of SNPs
	          output$bulkdownloadsnpInfo.txt <- downloadHandler(
	            filename = function() { "down.snp.info.txt" },
	            content = function(file) {
	              write.table(snp.info.down[[2]], file, sep="\t", quote=F, row.names=F)
	            })
	          
	          # Bulk download gene annotation
	          output$bulkdownloadgene.txt <- downloadHandler(
	            filename = function() { "down.gene.info.txt" },
	            content = function(file) {
	              if (exists("gff")){
	              }else{
	                load("./data/gff.AGP.v4.RData")
	              }
	              gene.info <- gff[gff$chr==myPos$chr & gff$start>=myPos$start & gff$end<=myPos$end, ]
	              write.table(gene.info, file, sep="\t", quote=F, row.names=F)
	            })
	        } else {
	          NULL
	        }
	      } else {
	        sendSweetAlert(
	          session = session,
	          title = "Error input!", type = "error",
	          text = "Please input genomic region or gene model in appropriate format!"
	        )
	      }
	     
	    })
	  } else {
	    NULL
	  }
	})
	
	observe({
	  if (input$clearBDo>0) {
	    isolate({
	      updateTextInput(session, "regBB", value="")
	    })
	  } else {NULL}
	})
	
	observe({
	  if (input$BDoExam >0) {
	    isolate({
	      updateTextInput(session, "regBB", value="chr1:29611303-29639223")
	    })
	  } else {NULL}
	})

})


