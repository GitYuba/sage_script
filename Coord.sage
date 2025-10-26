def coord_bary(rep_aff, pt):
    # Définir les points
    A = vector(rep_aff[0])
    B = vector(rep_aff[1])
    C = vector(rep_aff [2])
    P = vector(pt)
    var('a,b,c')

    # Définir les equations
    eq0 = a * A[0] + b * B[0] + c * C[0] == P[0] 
    eq1 = a * A[1] + b * B[1] + c * C[1] == P[1]
    # Calculer les coordonnées barycentriques de P par rapport à A, B, C
    s = solve([eq0, eq1, a + b + c == 1],a, b, c)

    # Get the right-hand side of solutions
    return show("P(-",abs(s[0][0].rhs()),",", s[0][1].rhs(),",", s[0][2].rhs(),")")

#coord_bary([[0,0],[1,0],[0,1]],[1,2])
