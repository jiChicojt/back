from flask import Flask
from flask_cors import CORS
import json
from flask import request, jsonify, Response
from sympy import *
from sympy.parsing.sympy_parser import standard_transformations, implicit_multiplication_application, convert_xor

app = Flask(__name__)
app.config["DEBUG"] = True
CORS(app)

from routes import solver

if __name__ == "__main__":
    app.run(debug = True,threaded = True)
    
@app.route('/message', methods = ['GET'])
def message():
    return {"message": "All good"}

@app.route('/solver', methods = ['POST'])
def solveEquqtion():
    data = json.loads(request.data.decode())
    equation = data['equation']
    hom = equation.split(sep="=")[0]
    part = equation.split(sep="=")[1]

    splitEquation = list(filter(None, hom.split(sep="'")))

    solution = solveEquation(splitEquation, part)
    steps = getSteps(splitEquation)

    return {"solution": solution, "steps": steps}

def getSteps(splitEquation):
    number = len(splitEquation)
    x, m = symbols('x m')
    y = [x**m]
    transformations = (standard_transformations + (implicit_multiplication_application,) + (convert_xor,))
    substitutionEq = 0
    auxEq = 0

    for num in range(0, number):
        y.append(y[0].diff(x, number - (num + 1)))
        splitEquation[num] = parse_expr(splitEquation[num].replace('y', str(y[0].diff(x, number - (num + 1)))),
                                        transformations=transformations)
        substitutionEq += splitEquation[num]
        auxEq += splitEquation[num] / x**m

    for num in range(0, len(y)):
        y[num] = transformXpr(y[num])

    del y[0]

    expandedAuxEq = expand(auxEq)
    solutionsAuxEq = solveset(Eq(expandedAuxEq, 0), m)
    solutionsAuxEq = str(solutionsAuxEq).replace('{', '').replace('}', '').split(',')

    return {"derivatives": y, "substitution": transformXpr(substitutionEq), "auxEq": transformXpr(auxEq),
            "xAuxEq": transformXpr(expandedAuxEq), "solutions": transformXpr(solutionsAuxEq)}

def solveEquation(splitEquation, part):
    x = symbols('x')
    f = Function('f')
    transformations = (standard_transformations + (implicit_multiplication_application,) + (convert_xor,))
    number = len(splitEquation)
    homEq = 0
    partEq = parse_expr(part, transformations=transformations)

    for num in range(0, number):
        print(splitEquation[num])
        homEq += parse_expr(splitEquation[num].replace('y', '1'), transformations=transformations)*f(x).diff(x, number - (num + 1))
        print(homEq)

    diffEq = Eq(homEq, partEq)
    sol = dsolve(diffEq, f(x))
    strSol = 'f(x) = ' + str(sol)[9:len(str(sol))-1].replace('**', '^').replace('*', '')

    return strSol

def transformXpr(xpr):
    return str(xpr).replace('**', '^').replace('*', '')