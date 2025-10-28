
###
## Importing SAGE files "pythonically" from github
## load(path_to_sage_file):
## load("https://raw.githubusercontent.com/medbednaoui/sage_script/refs/heads/main/Geometries/coord.sage")

###
##
def coord_bary(rep_aff, pt):
    """
    Calculer les coord bary du point pt dans le repère affine rep_aff
    
    Args:
    rep_aff = [[0,0],[1,0],[2,1]]
    pt = [1,2]
    
    Returns:
    Donne les coord bary absolues du point pt dans le repère affine rep_aff : [a,b,c] (list)
    """
    
    # Définir les points
    A = vector(rep_aff[0])
    B = vector(rep_aff[1])
    C = vector(rep_aff[2])
    P = vector(pt)
    var('a,b,c')
    # Justifier le fait que rep_aff est un repère affine
    d = det(matrix([B-A,C-A]))
    if d == 0:
        print("Les points entrés sont alignés !!!")
        print("Entrer des points non alignés.")
        return
        
    # Définir les equations
    eq0 = a * A[0] + b * B[0] + c * C[0] == P[0] 
    eq1 = a * A[1] + b * B[1] + c * C[1] == P[1]
    # Calculer les coordonnées barycentriques de P par rapport à A, B, C
    s = solve([eq0, eq1, a + b + c == 1],a, b, c)

    # Get the right-hand side of solutions
    return [s[0][0].rhs(), s[0][1].rhs(), s[0][2].rhs()]
    
###
##
def base_adaptee(rep_proj):
    """
    Trouver la base adaptée au repère projectif.
    
    Args:
    rep_proj = [[1,2,1],[1,-1,2],[0,1,3],[1,1,3]] (
    A = vector([1,2,1]) •••
    Où les quatre points forment un repère projectif du plan projectif
    
    Returns:
    La base adaptée (e1,e2,e3,e4) (list of vectors) au repère projectif (A,B,C,D)
    """
    
    # Définir les points
    A = vector(rep_proj[0])
    B = vector(rep_proj[1])
    C = vector(rep_proj[2])
    D = vector(rep_proj[3])
    # rep_proj doit être un repère projectif
    d_A = det(matrix([B,C,D]))
    d_B = det(matrix([A,C,D]))
    d_C = det(matrix([A,B,D]))
    d_D = det(matrix([A,B,C]))
    if d_A*d_B*d_C*d_D == 0:
        print("Les points entrés ne forment pas un repère projectif !!!")
        print("Entrer des points qui forment un repère projectif.")
        quit()    
    # Définir les variables 
    var('a,b,c')
    # Expression des équations
    eq0 = a * A[0] + b * B[0] + c * C[0] == D[0] 
    eq1 = a * A[1] + b * B[1] + c * C[1] == D[1]
    eq2 = a * A[2] + b * B[2] + c * C[2] == D[2]
    # Résolution du système
    s = solve([eq0, eq1, eq2], a, b, c)
    # Déduire la base adaptée (e1,e2,e3,e4) au repère projectif (A,B,C,D) (Get the right-hand side of solutions)
    e1 = s[0][0].rhs() * A
    e2 = s[0][1].rhs() * B
    e3 = s[0][2].rhs() * C
    e4 = e1 + e2 + e3 
    return [e1,e2,e3,e4]
