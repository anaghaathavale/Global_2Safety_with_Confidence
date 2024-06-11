'''
Tensorflow Example
====================

Top contributors (to current version):
  - Christopher Lazarus
  - Kyle Julian
  
This file is part of the Marabou project.
Copyright (c) 2017-2019 by the authors listed in the file AUTHORS
in the top-level source directory) and their institutional affiliations.
All rights reserved. See the file COPYING in the top-level source
directory for licensing information.
'''
import sys
sys.path.append('/home/ubuntu/Marabou')
from maraboupy import Marabou
import numpy as np

# %%
# This network has inputs x0, x1, and was trained to create outputs that approximate
# y0 = abs(x0) + abs(x1), y1 = x0^2 + x1^2
filename = "/home/ubuntu/frozen_models/frozen_graph_for_test.pb"
network = Marabou.read_tf(filename)

# %%
# Or, you can specify the operation names of the input and output operations.
# The default chooses the placeholder operations as input and the last operation as output
inputNames = ['x','y']
outputName = ['Identity','Identity_1']
network = Marabou.read_tf(filename = filename, inputNames = inputNames, outputNames = outputName, modelType='frozen')

# %%
# Get the input and output variable numbers; [0] since first dimension is batch size
inputVars = network.inputVars[0][0]
outputVars = network.outputVars[0][0]

# %%
# Set input bounds on both input variables
network.setLowerBound(network.inputVars[0][0][0], 0.0)
network.setUpperBound(network.inputVars[0][0][0], 1.0)
network.setLowerBound(network.inputVars[0][0][1], 0.0)
network.setUpperBound(network.inputVars[0][0][1], 1.0)
network.setLowerBound(network.inputVars[1][0][0], 0.0)
network.setUpperBound(network.inputVars[1][0][0], 1.0)
network.setLowerBound(network.inputVars[1][0][1], 0.0)
network.setUpperBound(network.inputVars[1][0][1], 1.0)

print ("output 1: ",network.outputVars[0][0][0])
print ("output 2: ",network.outputVars[1][0][0])
print ("input 1: ",  network.inputVars[0][0][0])
print ("input 2: ",network.inputVars[0][0][1])
print ("input 3: ",network.inputVars[1][0][0])
print ("input 4: ",network.inputVars[1][0][1])


# %%
# Set output bounds on the second output variable
network.setLowerBound(outputVars[0], 194.0)
network.setUpperBound(outputVars[0], 200.0)
network.setLowerBound(network.outputVars[1][0][0], 194.0)
network.setUpperBound(network.outputVars[1][0][0], 200.0)

#y0-y1 <= 0
#network.addInequality([network.outputVars[0][0][0],
#                               network.outputVars[1][0][0]],
#                              [1, -1], 0.1, isProperty=True)



network.addInequality([network.inputVars[0][0][0],
                               network.inputVars[1][0][0]],
                              [1, -1], 0.00000000001, isProperty=True)  

network.addInequality([network.inputVars[0][0][1],
                               network.inputVars[1][0][1]],
                              [1, -1], 0.00000000001, isProperty=True)  



# %%
# Call to C++ Marabou solver
#exitCode, vals, stats = network.solve("marabou.log")

#options = Marabou.createOptions(snc=True, verbosity=0, initialSplits=2);
exitCode, vals, stats = network.solve()
#print(exitCode)


#res, _, _ = network.solve("marabou.log")
#if res == 'sat':
#            # It is possible that y_correct <= y_i.
#            print("sat:not robust")
            
#else:
#           print("unsat:robust")        
