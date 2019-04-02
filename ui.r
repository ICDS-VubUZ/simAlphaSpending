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
				column(3,
					uiOutput('side'),
					uiOutput('sg')
				),
				column(6,
					uiOutput('es'),
					uiOutput('mu')
				)
			),
			fluidRow(
				column(3,
					numericInput("t1e", "type I error", min=0.005, max=.5, step=.005, value=.05)
				),
				column(9,
					textInput('ia', 'data points per stage (semicolon delimited, eg., 3;9;17)', "")
				)
			),
			fluidRow(
				column(3,
					selectInput("type", "type:", choices = c("OBF", "Pocock","compromise"))
				),
				column(3,
					textInput('nrsim', '# sim (eg., 10000)', "10")
				),
				column(6#,
					# actionButton("simulate", "simulate the alpha spending function sample sizes")
				)
			),
			h5("Susanne Blotwijk: simulation algorithm"),
			h6("https://www.icds.be  (Wilfried Cools: shiny suit)"),
			width=6
		),

		mainPanel(
			# h4("set up used for simulation"),
			actionButton("simulate", "SIMULATE the alpha spending function sample sizes"),
			br(),
			br(),
			tableOutput("setup"),
			width=6
		)

	)
)