alpha=0.7836
Tr=296.7
Ts=61.97 
def fx(a,b,c):
    return alpha*(-b*Tr-c*Tr*Tr + 2*c*Tr*Ts)

def gx(a,b,c):
    return alpha*(b-2*c*Ts)

def hx(a,b,c):
    return c*alpha


params = {
    """
    From eq. S2, Ex in form of a+bT+cT^2
    """
    'H':{
        'a': -22.657,
        'b': 0.15379,
        'c': -0.00025597
    },
    'A':{
        'a': -23.364,
        'b': 0.15876,
        'c': -0.00026696
    },
    'O':{
        'a': 2.1607,
        'b': -0.015064,
        'c': 0.000026
    },
    'P':{
        'a': 10.475,
        'b': 0.071482,
        'c': 0.0001201
    },
    'C':{
        'a': 8.5997,
        'b': -0.057676,
        'c': 0.000093317
    }
}
for key, val in params.items():
    # print(key, val.values())
    a,b,c = val.values()
    print(f"lambda({key}) = lambda_HPS + {fx(a,b,c):.8f} + {gx(a,b,c):.8f}*T + {hx(a,b,c):.8f}*T^2")
