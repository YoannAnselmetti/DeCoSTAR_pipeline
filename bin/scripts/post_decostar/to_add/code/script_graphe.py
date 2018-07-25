# composantes connexes d'un graphe
def composantes(graphe):
    comp=[]
    liste_des_sommets = graphe.keys()
    marques = {}
    premier_non_marque=0
    while premier_non_marque < len(liste_des_sommets):
        premier_sommet = liste_des_sommets[premier_non_marque]
        comp.append([premier_sommet])
        pile=[premier_sommet]
        marques[premier_sommet] = 'M'
        while pile:
            voisins = graphe[pile[0]]
            for v in voisins:
                if not marques.has_key(v):
                    pile.append(v)
                    comp[-1].append(v)
                    marques[v] = 'M'
            del pile[0]
        while (premier_non_marque < len(liste_des_sommets) and
               marques.has_key(liste_des_sommets[premier_non_marque])):
            premier_non_marque=premier_non_marque+1
return comp