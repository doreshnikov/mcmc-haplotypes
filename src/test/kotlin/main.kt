import aligner.SnipsOnlyAligner
import graph.AlignmentTest
import aligner.Reference
import graph.AlignedDeBruijnGraph
import org.jetbrains.numkt.core.ExperimentalNumkt

fun printGraph(graph: AlignedDeBruijnGraph) {
    println("Graph: ")
    graph.indices.forEach { v ->
        println("$v: ${graph[v]}")
        graph.edges[v].forEach { e ->
            println("\t-> $e")
        }
    }
}

@ExperimentalNumkt
fun main() {

    val testReference: Reference = "GATTACA"
    val testReads = listOf(
        "GCT" to 1, "GAT" to 4,
        "CTT" to 1, "ATT" to 3,
        "TTA" to 2,
        "TAT" to 1, "TAC" to 4,
        "ATA" to 1, "ACA" to 3
    )
    val allTestReads = testReads.flatMap { entry -> MutableList(entry.second) { entry.first } }

    AlignmentTest(
        testReference,
        allTestReads.distinct()
    )

    val aligner = SnipsOnlyAligner(testReference)
    val testGraph = AlignedDeBruijnGraph.build(allTestReads.map { aligner.align(it) }, k = 2)

    ManualCheckTest.executeAll()

    printGraph(testGraph)
    printGraph(testGraph.Normalizer().normalize())


}