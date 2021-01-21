
def Cpcal(x1):
    x2 = 1 -x1
    cp1 = 144
    cp2 = 154
    Cp = x1* cp1 + x2 * cp2
    return Cp

def kwhToJ(kwh):
    J = kwh * 1000 * 3600
    return J



def EnergyBalace(x1, kwhboiler, T2, T1, grams):
    hvap = 39840
    moles = grams/60.1
    Cp = Cpcal(x1)
    J = kwhToJ(kwhboiler)
    time = moles * Cp * (T2 - T1)/ J * 3600
    vapperhour = J/hvap
    energyreleased = -J
    return time, vapperhour, energyreleased

# time, vaper, energyreleased = EnergyBalace(.7, 3.29, 83, 25, 1589)
# print("Time to boil:", round(time), "seconds\nVapor produced per hour:", round(vaper),  "moles\nDuty needed for condenser: ", round(energyreleased),
#       "Joules per hour")

print(Cpcal(.3))

