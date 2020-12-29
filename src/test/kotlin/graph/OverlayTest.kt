package graph

import ManualCheckTest
import mcmc.modules.graph.PathsOverlay
import org.jetbrains.numkt.core.ExperimentalNumkt

@ExperimentalNumkt
class OverlayTest(
    private val graph: Graph,
    private val overlay: PathsOverlay
) : ManualCheckTest("Overlay test", {
    println("Graph:")
    graph.indices.forEach { v ->
        println(" - v$v (${graph[v].data}) -> " +
                graph.edges[v].joinToString(", ") { it.target.toString() })
    }
    println("Overlay log likelihood = ${overlay.logLikelihood()}")
})