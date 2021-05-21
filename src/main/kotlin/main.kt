import aligner.SnipsOnlyAligner
import datagen.MutationSimulator
import datagen.ReadSimulator
import datagen.generateReference
import graph.AlignedDeBruijnGraph
import mcmc.Engine
import mcmc.modules.graph.DistributionConfig
import mcmc.modules.graph.PathsOverlay
import utils.Randseed

fun main() {

    val reference = generateReference(30)
    val haplotypes = MutationSimulator(reference, 1.2).generate(4)
    haplotypes.forEach { println(it) }

    val reads = ReadSimulator(
        haplotypes.mapIndexed { i, it -> it to (i + 1).toDouble() / 10 },
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

    val model = PathsOverlay(normalizedGraph, DistributionConfig(haplotypes.size.toDouble()), 0.001)
    val engine = Engine(normalizedGraph, model)
    engine.simulate(10000)

    val result = model.collectResult()
    println(result)

}