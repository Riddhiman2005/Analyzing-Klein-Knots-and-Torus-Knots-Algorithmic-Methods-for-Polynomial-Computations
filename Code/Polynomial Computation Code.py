


# importing necessary functions and symbols from the SymPy library. 
# Importing the symbol "A" and "t" for algebraic calculations.



from sympy import *

A = symbols("A")
t = symbols("t")

# TLComponent class to store the order of the generator and provides methods to manipulate and access the generator's properties.

class TLComponent:
    # representation of a single generator of TL_n algebra.
    
    def __init__(self, i, n):
        # initiate, in TL_n, the ith generator.
        self.i = i
        self.order = list(range(1, n+1))
        self.order[i - 1 ] = i+1
        self.order[i] = i
    
    def __getitem__(self, key):
        temp = self.order[key]
        return temp

    def sum(self):
        return sum(self.order)

    def getFirstNonzero(self):
        # Finds first non-traveled stretch on this generator.

        for i in range(len(self.order)):
            if self.order[i] != 0:
                return i + 1

        return -1

class TL:
    # representation of an entire tangle.

    def __init__(self, order, n):
        # initiate in TL_n, given an order of generators.
        self.n = n
        self.originalOrder = order
        self.order = []
        self.count = 0
        for i in order:
            self.order.append(TLComponent(i, n))




class Klein:
    # Klein representation in code.
  # The Klein class represents a Klein knot and stores the parameters m and n of the knot
  # Provides methods to calculate the braid representation of the Klein knot.

    def __init__(self, m, n):
        # Initiate a Klein knot K(m, n). Also calculates the writhe of the knot.
        self.m = m
        self.n = n
        self.writhe = m * (n - 1)
        for i in range(1, n):
            self.writhe -= i

    def returnBraidRep(self):
        # Calculate the braid representation of this Klein knot K(m, n).
        returnList = []
        for i in range(self.m):
            for j in range(1, self.n):
                returnList.append(j)

        for i in range(0, self.n - 1):
            for j in range(self.n-1, i, -1):
                returnList.append(-j)

        return removeInverses(returnList)


class PolyPart:
    # Class for implementing non-commutative polynomial multiplication. Useful for keeping the TLn algebra generators in correct order.
  # PolyPart class represents a part of a polynomial. 
  # Provides methods for non-commutative polynomial multiplication and handles the order of generators of the TLn algebra

    def __init__(self, generator):
        # Initiate the polynomial representation.
        # It is possible to initiate with a list, in which case the polynomial is simply defined using that list,
        # or it is also possible to initiate with a generator of a TLn algebra, which the code then calculates.
        if type(generator) == list:
            self.poly = generator
            return

        if(generator > 0):
            self.poly = ["A", "A**-1 x e"+ str(int(generator))]
        else:
            self.poly = ["A**-1", "A x e" + str(int(-generator))]

    def __mul__(self, other):
        # Special multiplication method to preserve order of generators of TLn algebra.
        product = []
        len1 = len(self.poly)
        len2 = len(other.poly)
        for i in range(len1):
            for j in range(len2):
                if("e" in self.poly[i] or "e" in other.poly[j]):
                    # if there are any generators..
                    selfreplacement = self.poly[i].split(" x ")
                    otherreplacement = other.poly[j].split(" x ")
                    if(len(selfreplacement) > 1):
                        vars1 = selfreplacement[1]
                    else:
                        vars1 = ""
                    if(len(otherreplacement) > 1):
                        vars2 = otherreplacement[1]
                    else:
                        vars2 = ""
                    replacement = selfreplacement[0] + "*" +  otherreplacement[0]
                    variables = vars1 + vars2
                    product.append(replacement + " x " + variables)
                else:
                    # no generators! Simple, just multiply.
                    product.append(self.poly[i] + "*" + other.poly[j])
        return PolyPart(product)


def removeInverses(inputList):
    # Simplify braid polynomials to make sure that an element and its inverse cancel each other out if adjacent.
  
    changed = True
    while(changed):
        length = len(inputList)
        changed = False
        for i in range(length - 1):
            if inputList[i] == -1*inputList[i + 1]:
                changed = True
                inputList.pop(i)
                inputList.pop(i)
                break
    return inputList

def turnListIntoTLN(inputList, n):
    # Takes a braid polynomial (represented as a list) and represents as polynomials using the rho homomorphism.
    polyList = []
    for i in range(0, len(inputList)):
        polyList.append(PolyPart(inputList[i]))

    if(polyList != []):
        polyProduct = polyList[0]
        for i in range(1, len(polyList)):
            polyProduct = polyProduct*polyList[i]
    else: 
        polyProduct = PolyPart(["1"])

    return polyProduct

def evaluateTLNPolynomial(poly, n):
    # Evaluate the polynomial product. Each term in poly.poly is meant to be added together, but before we do so, 
    # we must first evaluate each of them.
    # If there are any generators in a term, the tangle they represent is calculated, and the trace is found. 
    # The remaining coefficients of As are multiplied together.
    print("Evaluating: ", end="")
    for i in range(0, len(poly.poly)):
        print(str(i) + " ", end="")
        term = poly.poly[i]
        if(" x " in term):
            # this means it has a generator. Must take care.
            term = term.split(" x ")
            coeffs = eval(term[0])
            variableOrders = term[1].split('e')[1:]
            for j in range(0, len(variableOrders)):
                variableOrders[j] = int(variableOrders[j])
            TLObject = TL(variableOrders, n)
            trace = TLObject.countComponents()
            variableBracket = (-A**-2 - A**2)**(trace - 1)
            poly.poly[i] = variableBracket*coeffs
        else:
            # no generator! But still, it's now the identity element.
            coeffs = eval(term)
            variableOrders = []
            TLObject = TL([], n)
            trace = TLObject.countComponents()
            variableBracket = (-A**-2 - A**2)**(trace - 1)
            poly.poly[i] = variableBracket*coeffs

    print()

    bracket = 0

    print("Adding: ", end="")
    for i in range(0, len(poly.poly)):
        print(str(i) + " ", end="")
        bracket = bracket + poly.poly[i]

    # Bracket polynomial!
    print()

    return bracket

def main():
    print("Which Klein link's polynomial would you like to compute?")
    m = int(input("Enter m now: "))
    n = int(input("Enter n now: "))

    print("Computing Klein representation...")

    klein = Klein(m, n)

    if n == 0:
        print("This Klein link is simply " + str(m) + " unknots.")
        print("Evaluating polynomial...")
        bracket = (-A**2 - A**-2) ** (m - 1)
        klein.writhe = 0
    else:
        print("Computing Klein representation...")
        klein = Klein(m, n)
        print("This knot has writhe " + str(klein.writhe))
        print("Finding braid polynomial...")
        braid = klein.returnBraidRep()
        print("Turning braid into a product of polynomials...")
        polyProduct = turnListIntoTLN(braid, n)
        print("Evaluating each polynomial and multiplying...")
        bracket = evaluateTLNPolynomial(polyProduct, n)

    print("Normalizing...")
    # normalize the bracket polynomial according to the writhe of the Klein knot.
    normalizedFinalBracket = simplify((-A**3)**(-klein.writhe) * bracket)
    print("The normalized bracket polynomial of this Klein knot is " + str(normalizedFinalBracket))
    jones = expand(normalizedFinalBracket.subs(A, t**(-1/4)))
    print("The Jones polynomial for this Klein knot is " + str(jones))
    print("==================== Success! ======================")
    return normalizedFinalBracket, jones


while True:
    main()
