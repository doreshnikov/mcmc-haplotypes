package graph.components

data class WeightedEdge(
    val source: Int,
    val target: Int,
    val weight: Double
) {

    val locator get() = Edge(source, target)

}

typealias AdjacencyList = List<List<WeightedEdge>>