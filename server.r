# shiny functionality
# retrieve input when changed
# add test dependent input (mu, sg, es)
# show specification
# provide output using functions of susanne
source("susanne.r")
server <- function(input, output) {

	getInput <- reactive({
		out <- list()
		validate(
			if(input$test=="t-test"){
				need(input$es, 'specify a cohen D')
			},
			if(input$test=="F-test"){
				need(input$mu, 'specify a vector of averages')
				need(input$sg, 'specify a standard deviation')
				need(input$sg > 0, 'standard deviations should be strictly positive')
			},
			need(input$ia, 'specify data points per stage')
		)
	   
		.ia <- as.numeric(unlist(gsub(" ","",unlist(strsplit(input$ia,";")))))
		validate(
			need(all(diff(c(1.99,.ia))>0),'data points should be increasing, values no less then 3')
		)
		# t1e
		if(input$test=="t-test"){
			out <- list(test=input$test,es=as.numeric(input$es),t1e=as.numeric(input$t1e),ia=.ia,nr=as.numeric(input$nrsim),type=input$type,side=as.numeric(input$side))
		}
		if(input$test=="F-test"){
			.mu <- as.numeric(unlist(gsub(" ","",unlist(strsplit(input$mu,";")))))
			validate(
				need(length(.mu)>1,'multiple group averages are required')
			)
			out <- list(test=input$test,mu=.mu,sg=as.numeric(input$sg),t1e=as.numeric(input$t1e),ia=.ia,nr=as.numeric(input$nrsim),type=input$type)
		}
		out
	})
	output$side <- renderUI({
		if(input$test=='t-test'){
			selectInput("side", "sided:", choices = c("1","2"))
		}
	})
	output$es <- renderUI({
		if(input$test=='t-test'){
			textInput('es','Cohens D (eg., 1)','1')
		}
	})
	output$mu <- renderUI({
		if(input$test=='F-test'){
			textInput('mu','group means (semicolon delimited, eg., 1;2;1)','1')
		}
	})
	output$sg <- renderUI({
		if(input$test=='F-test'){
			textInput('sg','st.dev.','1')
		}
	})

										   
						 
																																																
																												
	   
						  
	
	
	processInput <- reactive({
  
		myinput <- getInput()
		if(myinput$test=='t-test'){
			result <- simIntAtTest(myinput$es, myinput$ia, myinput$t1e, myinput$nr, sides = myinput$side, alphaSpendingType = myinput$type)
			tmp <- data.frame(do.call(rbind,lapply(result[c(1,2,5,6)],function(.x) as.vector(.x))))
			tmp <- rbind(tmp,NA)
			tmp <- rbind(tmp,NA)
														
			tmp <- cbind(c('cumulative power','alphas','cumulative alpha','target cumulative alpha','expected number','expected stop'),tmp)
			tmp <- cbind(tmp,NA)
														 
			tmp[5:6,ncol(tmp)] <- as.numeric(unlist(result[c(3,4)]))
			# tmp <- data.frame(lapply(tmp,as.character),stringsAsFactors=FALSE)
			# tmp[is.na(tmp)] <- ""
			# names(tmp) <- c("",paste0("test ",1:(ncol(tmp)-2)),"")
		}
		if(myinput$test=='F-test'){
			result <- simIntAnFtest(myinput$mu, myinput$sg, myinput$ia, myinput$t1e, myinput$nr)
			tmp <- data.frame(do.call(rbind,lapply(result[c(1,2,5,6)],function(.x) as.vector(.x))))
			tmp <- rbind(tmp,NA)
			tmp <- rbind(tmp,NA)
														
			tmp <- cbind(c('cumulative power','alphas','cumulative alpha','target cumulative alpha','expected number','expected stop'),tmp)
			tmp <- cbind(tmp,NA)
														 
			tmp[5:6,ncol(tmp)] <- as.numeric(unlist(result[c(3,4)]))
		}
		tmp <- data.frame(lapply(tmp,as.character),stringsAsFactors=FALSE)
		tmp[is.na(tmp)] <- ""
		names(tmp) <- c("",paste0("test ",1:(ncol(tmp)-2)),"")
		tmp
	})

	# add button
	simtF <- eventReactive(input$simulate, {
		out <- processInput()
		out
	}, ignoreNULL = TRUE)	
	

	output$setup <- renderTable({
		out <- simtF()
		out
	})
}
