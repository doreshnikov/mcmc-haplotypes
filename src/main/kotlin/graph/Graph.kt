package graph

import graph.components.WeightedEdge
import graph.components.AlignedVertex
import kotlinx.collections.immutable.ImmutableList
import kotlinx.collections.immutable.toImmutableList

abstract class Graph(
    vertices: List<AlignedVertex>,
    edges: List<List<WeightedEdge>>
): ImmutableList<AlignedVertex> by vertices.toImmutableList() {

    val edges = edges.map { it.toImmutableList() }.toImmutableList()
    val topologicalSort: List<Int>

    init {
        check(edges.size == vertices.size) { "Adjacency list should have same size as vertices" }

        val visited = MutableList(vertices.size) { false }
        val topologicalSort = mutableListOf<Int>()
        fun dfs(v: Int) {
            visited[v] = true
            edges[v].forEach { u ->
                if (!visited[u.target]) {
                    dfs(u.target)
                }
            }
            topologicalSort.add(v)
        }

        vertices.indices.forEach { s ->
            if (!visited[s]) {
                dfs(s)
            }
        }
        topologicalSort.reverse()
        this.topologicalSort = topologicalSort.toImmutableList()
    }

}