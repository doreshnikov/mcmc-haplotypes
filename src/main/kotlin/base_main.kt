import aligner.SnipsOnlyAligner
import datagen.MutationSimulator
import datagen.ReadSimulator
import datagen.generateReference
import graph.AlignedDeBruijnGraph
import mcmc.Engine
import mcmc.modules.graph.DistributionConfig
import mcmc.modules.graph.PathsOverlay
import scoring.PrimitiveScorer
import utils.Randseed
import java.io.PrintWriter

fun main() {

    val reference = generateReference(30)
    val haplotypes = MutationSimulator(reference, 1.2).generate(4)

    val weightedHaplotypes = haplotypes.mapIndexed { i, it -> it to (i + 1).toDouble() / 10 }
    PrintWriter("src/main/resources/_tmp/_reference").use { pw ->
        weightedHaplotypes.forEach {
            pw.println("${it.first} ${it.second}")
        }
    }

    val reads = ReadSimulator(
        weightedHaplotypes,
        500,
        Randseed.INSTANCE.uniformInteger(12, 15),
        0.01
    ).simulateAll()

    val alignments = SnipsOnlyAligner(reference).run {
        reads.map { align(it) }
    }

    val graph = AlignedDeBruijnGraph.build(alignments, 10)
    println("Paths: ${graph.validate()}")
    val cutGraph = graph.cutErrorTails()
    println("Paths after cutting: ${cutGraph.validate()}")
    val normalizedGraph = cutGraph.Normalizer().normalize()
    normalizedGraph.precalcPathCounts()

    val model = PathsOverlay(normalizedGraph, DistributionConfig(haplotypes.size.toDouble()), 0.001, 10)
    val engine = Engine(normalizedGraph, model)
    engine.simulate(10000)

    val results = listOf(model.extractResult(), engine.bestResult!!)
    results.forEachIndexed { idx, result ->
        PrintWriter("src/main/resources/_tmp/my${idx}").use { pw ->
            result.forEach {
                pw.println("${it.first} ${it.second}")
            }
        }

        println("==============")
        val score = PrimitiveScorer().score(weightedHaplotypes, result)
        for (i in weightedHaplotypes.indices) {
            println(score[i].joinToString(" ") { it.toString().padStart(6) } +
                    " -> " + weightedHaplotypes[i].second)
        }
        println(result.joinToString(" ") { String.format("%.4f", it.second) })
    }

}