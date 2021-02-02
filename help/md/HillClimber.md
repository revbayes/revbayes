## name
HillClimber
## title
Hill-Climber analysis object
## description
The HillClimber analysis object keeps a model and the associated moves and monitors. The object is used to run Markov chain Monte Carlo (HillClimber) simulation on the model, using the provided moves, to obtain a sample of the posterior probability distribution. During the analysis, the monitors are responsible for sampling model parameters of interest.
## details
 The HillClimber analysis object produced by a call to this function keeps copies of the model and the associated moves and monitors. The HillClimber analysis object is used to run Markov chain Monte Carlo (HillClimber) simulation on the model, using the provided moves, to obtain a sample of the posterior probability distribution. During the analysis, the monitors are responsible for sampling model parameters of interest.
## authors
Sebastian Hoehna
## see_also
SimulatedAnnealing
## example
	# Create a simple model (unclamped)
	a ~ exponential(1)
	mymodel = model(a)
	
	# Create a move vector and a monitor vector
	moves[1] = mvScale(a, lambda=1.0, weight=1.0)
	monitors[1] = mnFile(a,"output/out.log")
	
	# Create an HillClimber object
	myHillClimberObject = HillClimber( mymodel, monitors, moves)
	
	# Run a short analysis
	myHillClimberObject.burnin( generations = 400, tuningInterval = 100)
	myHillClimberObject.run( generations = 400)
	
	# print the summary of the operators (now tuned)
	myHillClimberObject.operatorSummary()
	
## references
