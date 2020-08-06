import random

def flashsort_classical(t,m=None):
    """
    Flash Sort Method -- classical
        K(Ai) = 1 + int((m-1)*(Ai-Amin)/(Amax-Amin))
    t : input array
    m : number of buckets
    """
    tmin = min(t)
    tmax = max(t)
    if tmax - tmin <= 0: return t

    # determine number of buckets
    if m is None:
        tmp = int(0.1*len(t))
        m = tmp if tmp != 0 else len(t)

    dt = (m-1.0) / (tmax-tmin)
    k = [int(dt*(i-tmin)) for i in t]

    # buckets
    L = [0 for i in range(m)]
    for i in k: L[i] += 1
    # accumulating
    i = 1
    while i < len(L):
        L[i] += L[i-1]
        i += 1

    # key-idea: make tmax always at the end
    nmax = t.index(tmax)
    tmp = t[0]
    t[0] = tmax
    t[nmax] = tmp

    cnt = 0
    # accumulating list index
    ndx = 0
    # real position index
    j = 0
    while cnt < len(t):
        while j > L[ndx]-1:
            j += 1
            ndx = int((t[j]-tmin)*dt)
        
        # if j == L[ndx], which means value is placed at its position, else do swap
        flash = t[j]
        while j != L[ndx]:
            ndx = int((flash-tmin)*dt)
            # care: python index starting from zero
            tmp = t[L[ndx]-1]
            t[L[ndx]-1] = flash
            flash = tmp
            # update accumulating buckets
            L[ndx] = L[ndx] - 1
            # update counting
            cnt += 1

    # now use straight insertion method resort t
    i = len(t) - 2
    while i>=0:
        if t[i+1] <= t[i]:
            j = i
            v = t[j]
            while t[j+1] < v:
                t[j] = t[j+1]
                j += 1
            t[j] = v
        i -= 1

    return t



def flashsort_custom(t,m=None):
    """
    Flash Sort Method -- optimalize
        K(Ai) = 1 + int((m-1)*(Ai-Amin)/(Amax-Amin))
    t : input array
    m : number of buckets
    """
    def func_bubbleSort(a):
        """sort from small to big"""
        i = 0
        while i < len(a)-1:
            j = i + 1
            while j < len(a):
                if a[j] < a[i]:
                    tmp = a[i]
                    a[i] = a[j]
                    a[j] = tmp
                j += 1
            i += 1
        return a

    tmin = min(t)
    tmax = max(t)
    if tmax - tmin <= 0: return t

    # determine number of buckets
    if m is None:
        tmp = int(0.1*len(t))
        m = tmp if tmp != 0 else len(t)

    # 2D buckets
    L = [[] for i in range(m)]

    dt = (m-1.0) / (tmax-tmin)
    for i in t:
        ndx = int(dt*(i-tmin))
        L[ndx].append(i)

    tsort = []
    for i in L:
        if len(i) != 0:
            tsort += func_bubbleSort(i)

    return tsort



if __name__ == '__main__':
    for i in range(100):
        tlist = [random.randrange(1,100) for i in range(5000)]

        #print('input-- ',tlist)
        t = flashsort_custom(tlist,300)
        #t = flashsort_classical(tlist)

        sort = sorted(tlist)
        for i,j in enumerate(sort):
            if t[i] != j:
                print('Error')
                exit()

    print('DONE')



