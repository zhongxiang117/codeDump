
def rng(seed=8):
    """
    Method: Linear Congruent Method

    ri := (MULTIPLIER * ri-1 + INCREMENT) mod MODULUS
    r0 supplied by user
    
    assume: MULTIPLIER = 4,  INCREMENT = 5,  MODULUS = 9
    """
    return (seed*4 + 5) % 9



def rng_packmol(seed=129128787):
    """beteen (0,1)"""
    def mult(p,q):
        p1 = int(p/10000)
        p0 = p % 10000
        q1 = int(q/10000)
        q0 = q % 10000
        num = (p0*q1+p1*q0) % 10000
        num = (num*10000+p0*q0) % 100000000
        return num
    s = (mult(seed,3141581)+1) % 100000000
    return s/100000000.0




