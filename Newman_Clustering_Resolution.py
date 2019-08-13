Ki = np.sum(A,0)                            #in-degree
Ko = np.sum(A,1)                            #out-degree
m = np.sum(Ki)                           	#number of edges
b = A - gamma*np.transpose(Ko*Ki)/m;
B = b + np.transpose(b)                              #directed modularity matrix
Ci = np.ones((len(A),1))                        #community indices
cn = 1                                   #number of communities
U = np.array([1,0])                                #array of unexamined communites

ind = np.arange(N)
Bg = B
Ng = N

while U[0]:
    w, v = np.linalg.eig(Bg)
    #i1 = list(w).index(max(np.real(w)))         #maximal positive (real part of) eigenvalue of Bg
    i1 = np.where(w == max(np.real(w)))[0]
    #print(i1)
    v1 = v[:,i1]
    
    S = np.ones((Ng,1))
    S[v1<0] = -1
    q=np.transpose(S)*Bg*S
    
    if (q > 1e-10):
        qmax = q
        np.fill_diagonal(Bg,0.)
        indg = np.ones((Ng,1))
        Sit = S
        while np.count_nonzero(~np.isnan(indg)) > 0:
            Qit = qmax - 4*np.multiply(Sit,(Bg*Sit))
            qmax = np.nanmax(np.multiply(Qit,indg))
            imax = np.where(np.multiply(Qit,indg) == np.nanmax(np.multiply(Qit,indg)))[0]
            #imax = list(np.multiply(Qit,indg)).index(np.nanmax(np.multiply(Qit,indg)))
            Sit[imax] = -Sit[imax]
            indg[imax] = float('nan')
            
            if (qmax > q):
                q = qmax              
                S = Sit
                
        if (np.abs(sum(S)) == Ng):
            del U[0]
            
        else:
            cn = cn+1
            #Ci[list(S).index(1)] = U[0]         #split old U(1) into new U(1) and into cn
            Ci[np.where(S == 1)[0]] = U[0]
            Ci[np.where(S == -1)[0]] = cn
            #Ci(list(S).index(-1)) = cn
            U = np.hstack((cn,U))  
    else:
        U[0] = 0
        
    ind = np.where(Ci == U[0])[0]                #indices of unexamined community U(1)
    bg = B[:,ind][ind]
    Bg = bg - np.diag(np.matrix(sum(bg)))                #modularity matrix for U(1)
    Ng = len(ind)
    
s = np.tile(Ci,N)                   #compute modularity
Q = np.multiply(np.logical_not(s - np.transpose(s)),B/(2*m))
Q = sum(Q)
        
            
        
        
    