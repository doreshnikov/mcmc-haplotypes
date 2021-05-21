package graph

import graph.components.WeightedEdge
import graph.components.AlignedVertex
import kotlinx.collections.immutable.ImmutableList
import kotlinx.collections.immutable.toImmutableList
import utils.Randseed

abstract class Graph(
    vertices: List<AlignedVertex>,
    val edges: List<List<WeightedEdge>>
) : ImmutableList<AlignedVertex> by vertices.toImmutableList() {

    val reversedEdges = MutableList<MutableList<WeightedEdge>>(vertices.size) { mutableListOf() }
        .also { re ->
            edges.forEach { neighbors ->
                neighbors.forEach { e ->
                    re[e.target].add(WeightedEdge(e.target, e.source, e.weight))
                }
            }
        }
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

    fun precalcPathCounts() {
        fun dfsForward(v: Int) {
            if (get(v).pathCountRight != 0L) return
            if (edges[v].isEmpty()) {
                get(v).pathCountRight = 1
            } else {
                for (e in edges[v]) {
                    dfsForward(e.target)
                    get(v).pathCountRight += get(e.target).pathCountRight
                }
                get(v).rightSelector = Randseed.INSTANCE.scoreIntegerEnum(
                    edges[v].map { get(it.target).pathCountRight.toDouble() }
                )
            }
        }

        fun dfsBackward(v: Int) {
            if (get(v).pathCountLeft != 0L) return
            if (reversedEdges[v].isEmpty()) {
                get(v).pathCountLeft = 1
            } else {
                for (e in reversedEdges[v]) {
                    dfsBackward(e.target)
                    get(v).pathCountLeft += get(e.target).pathCountLeft
                }
                get(v).leftSelector = Randseed.INSTANCE.scoreIntegerEnum(
                    reversedEdges[v].map { get(it.target).pathCountLeft.toDouble() }
                )
            }
        }

        forEach {
            it.pathCountLeft = 0L
            it.pathCountRight = 0L
        }
        indices.forEach {
            dfsForward(it)
            dfsBackward(it)
        }
    }

}