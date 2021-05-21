package mcmc.modules.graph

import graph.components.Edge

fun edges(vertices: Pathway) =
    (0..vertices.size - 2).map { i -> Edge(vertices[i], vertices[i + 1]) }
