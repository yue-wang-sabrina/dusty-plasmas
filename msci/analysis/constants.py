import numpy
import math
from scipy import special

# Some constants
radd = 4.5 * 10 ** (-6)  # radius of dust particle
e = 1.60217662 * 10 ** (-19)  # electron charge
e0 = 8.85418782 * 10 ** (-12)  # permittivity free space
Te = 46000  # Electrons usually a few EV. Non-thermal plasma, E field couples with electron more efficiently, higher temperature
Ti = 310  # Ion thermal temperature, room temperature
kb = 1.38064852 * 10 ** (-23)  # boltzmans constant
ne0 = 1. * 10 ** (15)  # initial electron density
ni0 = ne0  # Quasi-neutral initial conditions
md = 1000. * (4 / 3) * math.pi * radd ** 3  # Mass of dust particle
mi = 39.948 * 1.66053904 * 10 ** (-27)  # argon mass ion
me = 9.10938356 * 10 ** (-31)  # electron mass
cs = numpy.sqrt(kb * Te / mi)  # Bohm speed = sound speed
lambdade = (kb * Te * e0 / (ne0 * e ** 2)) ** 0.5  # Electron debye rad
lambdadi = (kb * Ti * e0 / (ni0 * e ** 2)) ** 0.5  # Ion debye rad
lambdaD = 1. / (1. / (lambdadi ** 2) + 1. / (lambdade ** 2)) ** 0.5  # dusty debye radius
boxr = 523 * lambdaD  # cylinder radius
boxz = 0.001  # cylinder height
g = -9.8  # gravitational acceleration
dt = 0.0001  # Much longer than charging time which is order nanoseconds
sheathd = 10 * lambdaD
electrodeV = abs((kb * Te / (2 * e)) * (numpy.log(2 * math.pi * me / mi)))  # potential at electrode
wallV = electrodeV  # cylindrical sides of wall same potential
radinfluence = 10 * lambdaD
dipolea = boxr / 100.
mu0 = 4 * math.pi * 10 ** (-7)  # Permeability free space
gamma = 5000.  # damping gamma
mu0 = 4 * math.pi * 10 ** (-7)  # Permeaility free space
Bmom = ((2 * math.pi * (0.00038) ** 3) * 0.014 / mu0) * numpy.array(
    [0, 0, 1])  # Nm/T #At 1cm away I want the B to be 0.014T
magBmom = numpy.sqrt(Bmom[0] ** 2 + Bmom[1] ** 2 + Bmom[2] ** 2)
Bmomhat = numpy.array(Bmom) / magBmom
dipolepos = [0, 0, -0.0005]


def OLMsol():  # Solve for dust grain surface potential and dust grain charge
    x0 = 3.3 * Te / e
    xnew = x0

    def f(x):
        return math.sqrt(me * Ti / (mi * Te)) * (1. - (Te / Ti) * x) - numpy.exp(x)  # x is the normalised potential

    def fprime(x):
        return math.sqrt(me * Ti / (mi * Te)) * (-Te / Ti) - numpy.exp(x)

    x0 = 0.  # Value approximated as number given in M.Lampe paper
    xnew = x0

    for _ in numpy.arange(1000):  # number of iterations
        xnew += - f(xnew) / fprime(xnew)

    surfpot = xnew  # normalised potential
    chargedust = (surfpot * radd * 4 * math.pi * e0) * (Te * kb / e) / e
    return surfpot * (Te * kb / e), chargedust  # This is in Volts and Coloumbs, i.e. SI units


# phi at surface of dust and the corresponding dust charge
phia, Zd = OLMsol()


def K(x):
    return x * numpy.arctan(x) + (numpy.sqrt(math.pi / 2) - 1) * (x ** 2 / (1 + x ** 2)) - numpy.sqrt(
        math.pi / 2) * numpy.log(1 + x ** 2)


def chandra(x):
    return (special.erf(x) - 2 * x * numpy.exp(-x ** 2)) / (2 * x ** 2 * numpy.sqrt(math.pi))