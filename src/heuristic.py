#import networkx as nx

def find_cover(G, vertex_cover_set):
    biclique_cover = []
    
    for u,v in G.edges:
        G.edges[u,v]["visited"] = False
    
    #visited_edges = [False for u,v in G.edges]
    
    for vertex in vertex_cover_set:
        new_list = []
        for v in G.neighbors(vertex):
            if vertex < v and not G.edges[vertex,v]["visited"]:
                new_list.append((vertex,v))
                G.edges[vertex,v]["visited"] = True
            elif v < vertex and not G.edges[v,vertex]["visited"]: 
                new_list.append((v,vertex))
                G.edges[v,vertex]["visited"] = True
                
        biclique_cover.append(new_list)

    return biclique_cover    
                
    
    