# shiny in and output
# input in sidebar, including action button
# output in main panel, including information and plot
ui <- fluidPage(

	titlePanel("Simulating Alpha Spending"),

	sidebarLayout(

		sidebarPanel(
			fluidRow(
				column(3,
					selectInput("test", "Test:", choices = c("t-test", "F-test"))
				),
				column(2,
					uiOutput('side'),
					uiOutput('sg')
				),
				column(7,
					uiOutput('es'),
					uiOutput('mu')
				)
			),
			fluidRow(
				column(3,
					numericInput("t1e", "type I error", min=0.005, max=.5, step=.005, value=.05)
				),
				column(9,
					textInput('ia', 'interim stages (semicolon delimited, eg., 0;3;9)', "")
				)
			),
			fluidRow(
				column(3,
					selectInput("type", "type:", choices = c("OBF", "Pocock","compromise"))
				),
				column(2,
					textInput('nrsim', '# sim (eg., 10000)', "10")
				),
				column(7#,
					# actionButton("simulate", "simulate the alpha spending function sample sizes")
				)
			),
			width=6
		),

		mainPanel(
			h4("set up used for simulation"),
			# actionButton("simulate", "SIMULATE the alpha spending function sample sizes"),
			# br(),
			# br(),
			tableOutput("setup"),
			width=6
		)

	)
)