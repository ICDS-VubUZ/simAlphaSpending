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
			},
			need(input$ia, 'specify interim analysis vector')
		)
		# t1e
		.ia <- as.numeric(unlist(gsub(" ","",unlist(strsplit(input$ia,";")))))
		if(input$test=="t-test"){
			out <- list(test=input$test,es=as.numeric(input$es),t1e=as.numeric(input$t1e),ia=.ia,nr=as.numeric(input$nrsim),type=input$type,side=as.numeric(input$side))
		}
		if(input$test=="F-test"){
			.mu <- as.numeric(unlist(gsub(" ","",unlist(strsplit(input$mu,";")))))
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

	# simtF <- eventReactive(input$simulate, {
		# myinput <- getInput()
		# if(myinput$test=="t-test") out <- simIntAtTest(myinput$es, myinput$ia, myinput$t1e, myinput$nr, sides = 1, alphaSpendingType = myinput$type)
		# if(myinput$test=="F-test") out <- simIntAnFtest(myinput$mu, myinput$sg, myinput$ia, myinput$t1e, myinput$nr)
		# out
	# }, ignoreNULL = FALSE)	
	

	output$setup <- renderTable({
		
		myinput <- getInput()
		if(myinput$test=='t-test'){
			result <- simIntAtTest(myinput$es, myinput$ia, myinput$t1e, myinput$nr, sides = myinput$side, alphaSpendingType = myinput$type)
			tmp <- data.frame(do.call(rbind,lapply(result[c(1,2,5,6)],function(.x) as.vector(.x))))
			tmp <- rbind(tmp,NA)
			tmp <- rbind(tmp,NA)
			# names(tmp) <- paste0("test ",1:length(input$nrobs))
			tmp <- cbind(c('cumulative power','alphas','cumulative alpha','target cumulative alpha','expected number','expected stop'),tmp)
			tmp <- cbind(tmp,NA)
			names(tmp) <- c("",paste0("test ",1:(ncol(tmp)-2)),"")
			tmp[5:6,ncol(tmp)] <- as.numeric(unlist(result[c(3,4)]))
		}
		if(myinput$test=='F-test'){
			result <- simIntAnFtest(myinput$mu, myinput$sg, myinput$ia, myinput$t1e, myinput$nr)
			tmp <- data.frame(do.call(rbind,lapply(result[c(1,2,5,6)],function(.x) as.vector(.x))))
			tmp <- rbind(tmp,NA)
			tmp <- rbind(tmp,NA)
			# names(tmp) <- paste0("test ",1:length(input$nrobs))
			tmp <- cbind(c('cumulative power','alphas','cumulative alpha','target cumulative alpha','expected number','expected stop'),tmp)
			tmp <- cbind(tmp,NA)
			names(tmp) <- c("",paste0("test ",1:(ncol(tmp)-2)),"")
			tmp[5:6,ncol(tmp)] <- as.numeric(unlist(result[c(3,4)]))
		}
		tmp
	})
}
