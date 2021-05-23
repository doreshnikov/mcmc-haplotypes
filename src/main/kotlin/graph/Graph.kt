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
        forEach {
            it.pathCountLeft = 0.0
            it.pathCountRight = 0.0
        }

        val groups = indices.groupBy { get(it).position }
        val maxPos = groups.keys.maxOrNull()!!
        for (i in maxPos downTo 0) {
            val group = groups.getValue(i)
            var total = 0.0
            for (v in group) {
                var opts = edges[v].map { e ->
                    get(e.target).pathCountRight
                }
                if (opts.isEmpty()) {
                    opts = listOf(1.0)
                }
                get(v).pathCountRight = opts.sum()
                get(v).rightSelector = Randseed.INSTANCE.scoreIntegerEnum(opts)
                total += get(v).pathCountRight
            }
            for (v in group) {
                get(v).pathCountRight /= total
            }
        }
        for (i in 0..maxPos) {
            val group = groups.getValue(i)
            var total = 0.0
            for (v in group) {
                var opts = reversedEdges[v].map { e ->
                    get(e.target).pathCountLeft
                }
                if (opts.isEmpty()) {
                    opts = listOf(1.0)
                }
                get(v).pathCountLeft = opts.sum()
                get(v).leftSelector = Randseed.INSTANCE.scoreIntegerEnum(opts)
                total += get(v).pathCountLeft
            }
            for (v in group) {
                get(v).pathCountLeft /= total
            }
        }
    }

}