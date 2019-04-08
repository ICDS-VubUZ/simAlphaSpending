# shiny in and output
# input in sidebar, including action button
# output in main panel, including information and plot
ui <- fluidPage(

	titlePanel("Simulating Alpha Spending"),

	sidebarLayout(

		sidebarPanel(
			fluidRow(
				column(4,
					selectInput("test", "test:", choices = c("t-test", "F-test"))
				),
				column(4,
					uiOutput('side'),
					uiOutput('sg')
				),
				column(4,
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
				column(4,
					selectInput("type", "type:", choices = c("OBF", "Pocock","compromise"))
				),
				column(4,
					textInput('nrsim', '# sim (eg., 10000)', "10")
				),
				column(4#,
					# actionButton("simulate", "simulate the alpha spending function sample sizes")
				)
			),
			tags$a(href="manual.pdf", "open introduction in new window", target="new"),
			h5("Susanne Blotwijk: simulation algorithm"),
			h6("https://www.icds.be  (Wilfried Cools: shiny suit)"),
			hr(),
			h5("interpretation"),
			img(src='ToolOutput.png',width='100%'),
			br(),
			br(),
			h5("diagnostics"),
			img(src='MoreSimulationsNeeded.png',width='100%'),
			width=4
		),

		mainPanel(
			actionButton("simulate", "SIMULATE the alpha spending function sample sizes"),
			br(),
			br(),
			tableOutput("setup"),
			width=8
		)

	)
)