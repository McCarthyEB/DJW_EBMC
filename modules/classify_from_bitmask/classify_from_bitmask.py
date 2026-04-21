import numpy as np


def build_neighbour_list(atoms, top_layer, cutoff=3.5):
    neighbours = {}
    for i_local, i_global in enumerate(top_layer):
        neigh = set()

        for j_local, j_global in enumerate(top_layer):
            if i_local == j_local:
                continue

            d = atoms.get_distance(i_global, j_global, mic=True)

            if d < cutoff:
                neigh.add(j_local)

        neighbours[i_local] = neigh

    return neighbours


def classify_from_bitmask(bitmask, neighbours):
    #  Classify vacancy pattern based on connectivity to other vacancies
    vac = [i for i, b in enumerate(bitmask) if b == 1]

    if len(vac) == 0:
        return "full_surface"

    # Build graph
    G = {i: set() for i in vac}

    for i in vac:
        for j in neighbours[i]:
            if j in vac:
                G[i].add(j)

    # Connected components
    seen = set()
    clusters = []

    for i in vac:
        if i not in seen:
            stack = [i]
            comp = set()

            while stack:
                x = stack.pop()
                if x in seen:
                    continue

                seen.add(x)
                comp.add(x)
                stack.extend(G[x] - seen)

            clusters.append(comp)

    # Classification
    labels = []

    for c in clusters:
        n = len(c)

        if n == 1:
            labels.append("isolated")
        elif n == 2:
            labels.append("pair")
        elif n == 3:
            #sort out triangles or chain by how many neighbours each has!
            degrees = []
            for node in c:
                deg = len(G[node].intersection(c))
                degrees.append(deg)

            degrees.sort()

            if degrees == [2, 2, 2]:
                labels.append("triangle")
            elif degrees == [1, 1, 2]:
                labels.append("chain3")
            else:
                labels.append("cluster3")
        else:
            labels.append(f"cluster{n}")

    return "+".join(sorted(labels))


def classify_structure(bitmask, atoms, top_layer, cutoff=3.0):
    neighbours = build_neighbour_list(atoms, top_layer)
    return classify_from_bitmask(bitmask, neighbours)
