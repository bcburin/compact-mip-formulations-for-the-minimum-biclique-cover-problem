

def find_cover(g, vertex_cover_set):
    biclique_cover = []
    for u, v in g.edges:
        g.edges[u, v]["visited"] = False
    for vertex in vertex_cover_set:
        new_list = []
        for v in g.neighbors(vertex):
            if vertex < v and not g.edges[vertex, v]["visited"]:
                new_list.append((vertex, v))
                g.edges[vertex, v]["visited"] = True
            elif v < vertex and not g.edges[v, vertex]["visited"]:
                new_list.append((v, vertex))
                g.edges[v, vertex]["visited"] = True
                
        biclique_cover.append(new_list)
    return biclique_cover
