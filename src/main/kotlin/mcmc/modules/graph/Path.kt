package mcmc.modules.graph

import graph.components.Edge
import kotlinx.collections.immutable.ImmutableList
import kotlinx.collections.immutable.toImmutableList

class Path(
    vertices: List<Int>,
    var weight: Double
) : ImmutableList<Int> by vertices.toImmutableList() {

    val edges get() = (0 until size - 1).map(this::edge)

    fun edge(i: Int): Edge {
        return Edge(get(i), (i + 1))
    }

    infix fun withWeight(weight: Double): Path {
        return Path(this, weight)
    }

    operator fun unaryMinus(): Path {
        return withWeight(-weight)
    }

}