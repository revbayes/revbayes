## name
model
## title
Create a model object

## description
Creates a model object that can be graphed or subjected to Bayesian inference.

## details
`model(x)` creates a model object by creating a copy of all elements and 
parameters that influence or are influenced by the likelihood of `x`.

Because `model` works with copies of objects, conducting an mcmc(mc) analysis
on a model object will not change the values of the objects in the RevBayes
workspace.

The model object can be modified to ignore specific data elements using the
method `ignoreData`.  Thus to run without the sequence data `phySeq` you
might specify:

   mymodel.ignoreData(phySeq)

Only clamped nodes can be ignored. To ignore all clamped nodes you can use
the method `ignoreAllData`:

   mymodel.ignoreAllData()

## authors
## see_also
## example
    # Create a simple model (unclamped)
    a ~ dnExponential(1)
    b ~ dnExponential(a)
    mymodel = model(b) # equivalent to model(a) or model(a, b)
    
    # Save a DOT visualization of the model to file
    mymodel.graph("mymodel.dot")
    
    # Create a move vector and a monitor vector
    moves = [ mvScale( a, lambda = 1.0, weight = 1.0 ) ]
    monitors = [ mnScreen(printgen = 10, a) ]
    
    # Create an mcmc object
    mymcmcObject = mcmc( mymodel, monitors, moves )
    
    # Print value of a
    print(a)
    
    # Run a short analysis
    mymcmcObject.run( generations = 100 )
    
    print(a) # Value is unchanged in the workspace - only the copy is modified
    
## references
