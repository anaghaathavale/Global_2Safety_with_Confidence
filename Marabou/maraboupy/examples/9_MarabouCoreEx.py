import sys
sys.path.append('/home/ubuntu/Marabou')
from maraboupy import MarabouCore
from maraboupy.Marabou import createOptions

inputQuery = MarabouCore.InputQuery()
inputQuery.setNumberOfVariables(4)

# %%
# Set the lower and upper bounds on variables

inputQuery.setLowerBound(0, 0)
inputQuery.setUpperBound(0, 1)
inputQuery.setLowerBound(1, 2)
inputQuery.setUpperBound(1, 3)
inputQuery.setLowerBound(2, 4)
inputQuery.setUpperBound(2, 5)
inputQuery.setLowerBound(3, -1.0)
inputQuery.setUpperBound(3, 10.0)

MarabouCore.addMaxConstraint(inputQuery, {0,1,2}, 3)

equation1 = MarabouCore.Equation()
equation1.addAddend(1, 3)
equation1.addAddend(-1, 0)
equation1.setScalar(0)
inputQuery.addEquation(equation1)

#input query 2:

inputQuery2 = MarabouCore.InputQuery()
inputQuery2.setNumberOfVariables(6)

inputQuery2.setLowerBound(0, 0)
inputQuery2.setUpperBound(0, 10)
inputQuery2.setLowerBound(1, 0)
inputQuery2.setUpperBound(1, 10)
inputQuery2.setLowerBound(2, 4)
inputQuery2.setUpperBound(2, 5)
inputQuery2.setLowerBound(3, 6.0)
inputQuery2.setUpperBound(3, 7.0)
inputQuery2.setLowerBound(4, 12.0)
inputQuery2.setUpperBound(4, 15.0)
inputQuery2.setLowerBound(5, -1.0)
inputQuery2.setUpperBound(5, 50.0)

MarabouCore.addMaxConstraint(inputQuery2, {0,1,2,3,4}, 5)
MarabouCore.addReluConstraint(inputQuery2, 0, 1)

equation2 = MarabouCore.Equation()
equation2.addAddend(1, 5)
equation2.addAddend(-1, 4)
equation2.setScalar(0)
inputQuery2.addEquation(equation2)


options = createOptions()
exitCode, vars, stats = MarabouCore.solve(inputQuery, options, "")
print(exitCode)
if exitCode == "sat":
    print(vars)

exitCode2, vars2, stats2 = MarabouCore.solve(inputQuery2, options, "")
print(exitCode2)
if exitCode2 == "sat":
    print(vars2)

