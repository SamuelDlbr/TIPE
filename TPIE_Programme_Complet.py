import math as ma
import matplotlib.pyplot as plt
import numpy as np
from numpy import diff
from numpy.linalg import norm, det
import circle_fit as cf

def Bezier(L,x):
    
    "retourne le polynôme de Bézier quadratique associé aux points de la"
    "liste L (de longueur 3) en x (pour x entre 0 et 1)"
    
    Polynome_de_bezier=np.array([0.0,0.0])
    
    for j in range(3): 
        Polynome_Bernstein_en_x=2/(ma.factorial(j)*ma.factorial(2-j))*(x**j)*(1-x)**(2-j)
        Polynome_de_bezier+=Polynome_Bernstein_en_x*L[j]
        
    return(Polynome_de_bezier)
    
def Raccordement_polynomes_de_Bezier(L,dx):
    
    "découpe une longue liste de points en petites liste de 3 points et"
    "rassemble les listes retournées par ListeBézier en une liste"
    
    Liste_des_positions=[]
    
    for i in range(int((len(L)-1)/2)):
        Liste_des_positions+=Bezier_sous_forme_de_liste(L[(i)*2:(i+1)*2+1],dx)
        
    return(np.array(Liste_des_positions))
    
def Bezier_sous_forme_de_liste(L,dx):
    
    "retourne une liste de points à intervalles réguliers appartenant au"
    "polynôme de Bézier quadratique associé à la liste L (de longueur 3)"
    
    Liste_des_positions=[]
    nombre_de_points=int(2/dx)
    
    for i in range(nombre_de_points):
        Liste_des_positions.append(Bezier(L,dx*i/2))
        
    return(Liste_des_positions)
    
def Calcul_du_mouvement_angulaire(angle_initial,vitesse_angulaire_initiale,Liste_points_courbe_Bezier,g,rayon,coef_friction,masse,vitesse):
    
    "résout l'équation du mouvement avec la méthode d'Euler"
    
    Vitesse_angulaire=[vitesse_angulaire_initiale]
    Angle=[angle_initial]
    Forces_d_inertie=Calcul_forces_inertie(Liste_points_courbe_Bezier,vitesse)
    Nombre_de_points_calcules=int(len(Forces_d_inertie))
    
    for i in range (1,Nombre_de_points_calcules):
        
        "le pas correspond à l'intervalle de temps entre les points de calcul"
        "de la méthode d'Euler. Il correspond au temps mit, à vitesse"
        "constante, pour relier deux points consécutifs de la liste"
        "des points de la suite de courbes de Bézier générées auparavant"
        
        pas=norm(Liste_points_courbe_Bezier[i+1]-Liste_points_courbe_Bezier[i])/vitesse
        Angle_precedent=Angle[-1]
        Vitesse_angulaire_precedente=Vitesse_angulaire[-1]
        Vitesse_angulaire.append(Vitesse_angulaire_precedente+pas*(Forces_d_inertie[i]/rayon*np.cos(Angle_precedent)
        -g/rayon*np.sin(Angle_precedent)-coef_friction/masse*Vitesse_angulaire_precedente))
        Angle.append(Angle_precedent+Vitesse_angulaire_precedente*pas)
        
    return(Angle)

def Calcul_forces_inertie(Liste_points_courbe_Bezier,vitesse):
    
    "déterminer les forces d'inertie le long du parcours (à vitesse constante)"
    
    Forces_inertie=[]
    for i in range(len(Liste_points_courbe_Bezier)-2):
        Points = Liste_points_courbe_Bezier[i:i+3]
        rayon_de_courbure = cf.least_squares_circle(Points)[2]
        Forces_inertie.append(vitesse**2/rayon_de_courbure)
        
    return (Forces_inertie)

def Affichage_trajectoire(angle_initial,vitesse_angulaire_initiale,Liste_points_courbe_Bezier,g,rayon,coef_friction,masse,vitesse):
    
    "affiche la trajectoire de la masse vue de haut, son angle le long du"
    "parcours étant calculé par Calcul_du_mouvement_angulaire"
    
    Angle=Calcul_du_mouvement_angulaire(angle_initial,vitesse_angulaire_initiale,Liste_points_courbe_Bezier,g,rayon,coef_friction,masse,vitesse)
    
    "on veut un vecteur unitaire orthonormal à la courbe pour pouvoir afficher"
    "la trajectoire de la bille vue de haut. On calcule le vecteur"
    "orthonormal à deux points consécutifs de la courbe"
    
    delta = diff(Liste_points_courbe_Bezier, axis = 0)
    norme = norm(delta, axis = 1)
    
    "en mettant en commentaire les lignes 145, 206 et 207 on affiche le graphe"
    "de l'angle le long de la trajectoire"
    "en mettant en commentaire la lignes 146 on affiche la trajectoire vue de"
    "haut ainsi que les points pointés"
    
    xcoord = [Liste_points_courbe_Bezier[x][0]+rayon*np.sin(-Angle[x])*delta[x,1]/norme[x] for x in range(2,len(Liste_points_courbe_Bezier)-3)]
    ycoord = [Liste_points_courbe_Bezier[x][1]+rayon*np.sin(Angle[x])*delta[x,0]/norme[x] for x in range(2,len(Liste_points_courbe_Bezier)-3)]
    
    plt.plot(xcoord,ycoord, c='r')
    "plt.plot(Angle)"
    
def Points_milieu(Liste):
    
    "ajoute à une liste les points au milieu des points présents afin de"
    "garantir des raccordements fluides aux courbes de Bézier"
    
    Liste_et_points_milieux=[]
    for i in range(len(Liste)-1):
        Liste_et_points_milieux.append([(Liste[i][0]+Liste[i+1][0])/2,(Liste[i][1]+Liste[i+1][1])/2])
        Liste_et_points_milieux.append(Liste[i+1])
    return(np.array(Liste_et_points_milieux[:-1]))
    
def General(pas,angle_initial,vitesse_angulaire_initiale,Bord_Droit,g,Bord_Gauche,coef_friction,masse,vitesse):
    
    "la fonction générale regroupe les résultats et affiche les données"
    
    CENTRE = (Bord_Droit+Bord_Gauche)/2
    "CENTRE est la liste des points au centre du tube"
    
    dist = norm(Bord_Droit-Bord_Gauche, axis = 1)
    rayon = np.mean(dist)/2
    "r est le rayon approximatif du tube"
    
    "scatter affiche sous forme de points les bords du tube"
    plt.scatter(Bord_Gauche[:,0],Bord_Gauche[:,1])
    plt.scatter(Bord_Droit[:,0],Bord_Droit[:,1])
    
    Points_centraux=Points_milieu(CENTRE)
    
    "pointsmilieu adapte la liste pour l'interpolation de Béziers"
    
    Liste_points_courbe_Bezier=Raccordement_polynomes_de_Bezier(Points_centraux,pas)
    
    "Raccordement_polynomes_de_Bezier renvoie une suite de polynomes de"
    "Bézier et trajgrapheliste calcule et affiche la trajectoire sur cette courbe"
    
    Affichage_trajectoire(angle_initial,vitesse_angulaire_initiale,Liste_points_courbe_Bezier,g,rayon,coef_friction,masse,vitesse)

if __name__ == '__main__':
    
    angle_initial=0
    vitesse_angulaire_initiale=0
    g=9.81
    pas=0.01
    coef_friction=0.016
    
    "Bord_Droit et Bord_Gauche sont les listes des points des bords du tube pointés avec Tracker"
    "Trajectoire_Pointée est la liste des points de la trajectoire de la bilee dans le tube"
    
    Bord_Droit = np.load('Bord_Droit.npy')
    Bord_Gauche = np.load('Bord_Gauche.npy')
    Trajectoire_Pointée = np.load('Trajectoire_Pointée.npy')

    plt.scatter(Trajectoire_Pointée[:,0],Trajectoire_Pointée[:,1], c='k', s=1)
    
    "la double boucle lance le progamme plusieurs  fois pour différentes"
    "masses de passagers et coefficients de friction"

    masse = 1
    vitesse = 8
    General(pas,angle_initial,vitesse_angulaire_initiale,Bord_Droit,g,Bord_Gauche,coef_friction,masse/1000,vitesse/5)
