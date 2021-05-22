import aligner.SnipsOnlyAligner
import aligner.alignment.PrimitiveAlignment
import graph.AlignedDeBruijnGraph
import kotlinx.coroutines.async
import kotlinx.coroutines.awaitAll
import kotlinx.coroutines.runBlocking
import mcmc.Engine
import mcmc.modules.graph.DistributionConfig
import mcmc.modules.graph.PathsOverlay
import java.io.File
import java.io.PrintWriter
import java.util.concurrent.atomic.AtomicInteger
import java.util.concurrent.locks.ReentrantLock
import kotlin.streams.toList
import kotlin.system.measureTimeMillis

fun readFasta(fileName: String): List<String> {
    File(fileName).bufferedReader().useLines { lines ->
        return@readFasta lines.toList()
            .filter { it.isNotEmpty() && it[0] != '>' }
            .map { it.toUpperCase() }
    }
}

fun main(args: Array<String>) = runBlocking {

    val dir = args[0]
    val referenceFile = "$dir/_reference"
    val readsFile = "$dir/_reads"

    val reference = readFasta(referenceFile).joinToString("")
    val alignments = mutableListOf<PrimitiveAlignment>()
    val k = 10

    if (!File("$dir/_aligned").exists()) {
        val reads = File(readsFile).bufferedReader().lines().toList()
        val lock = ReentrantLock()
        println("Aligning ${reads.size} reads...")

        val time = measureTimeMillis {
            val cnt = AtomicInteger(1)
            val aligner = SnipsOnlyAligner(reference)
            print("[=> ")
            val jobs = reads.map { read ->
                async {
                    val size = synchronized(lock) {
                        alignments.add(aligner.align(read))
                        alignments.size
                    }
                    if (size > reads.size / 20 * cnt.get()) {
                        cnt.incrementAndGet()
                        print("-")
                    }
                }
            }
            jobs.awaitAll()
            println(" ]")
        }

        println("Total time: ${time / 1000}s, average ${time / reads.size}ms/read")
        val sortedAlignments = alignments.sortedBy { it.start }
        if (sortedAlignments.first().start > 0) {
            val start = sortedAlignments.first().start
            for (i in sortedAlignments.indices) {
                if (sortedAlignments[i].start > start) {
                    break
                }
                alignments.add(
                    PrimitiveAlignment(
                        0, reference.substring(0 until start) +
                                sortedAlignments[i].read.substring(0 until k - 1)
                    )
                )
            }
        }
        if (sortedAlignments.last().end < reference.length) {
            val end = sortedAlignments.last().end
            for (i in sortedAlignments.indices.reversed()) {
                if (sortedAlignments[i].end < end) {
                    break
                }
                val readSize = sortedAlignments[i].read.length
                alignments.add(
                    PrimitiveAlignment(
                        end - k, sortedAlignments[i].read.substring(readSize - k + 1 until readSize) +
                                reference.substring(end until reference.length)
                    )
                )
            }
        }

        File("$dir/_aligned").printWriter().use { pw ->
            alignments.forEach {
                pw.println("${it.read} ${it.start}")
            }
        }
    } else {
        println("Using preprocessed alignments from _aligned")
        File("$dir/_aligned").bufferedReader().use {
            it.lines().forEach { line ->
                val (read, start) = line.split(" ")
                alignments.add(PrimitiveAlignment(start.toInt(), read))
            }
        }
    }

    val lambda = 10.0
    val graph = AlignedDeBruijnGraph.build(alignments, 10)
    println("Paths: ${graph.validate()}")
    val cutGraph = graph.cutErrorTails()
    println("Paths after cutting: ${cutGraph.validate()}")
    val normalizedGraph = cutGraph.Normalizer().normalize()
    normalizedGraph.precalcPathCounts()

    val model = PathsOverlay(normalizedGraph, DistributionConfig(lambda), 0.001)
    val engine = Engine(normalizedGraph, model)
    engine.simulate(10000)

    val results = listOf(model.extractResult(), engine.bestResult!!)
    results.forEachIndexed { idx, result ->
        PrintWriter("src/main/resources/_tmp/my${idx}").use { pw ->
            result.forEach {
                pw.println("${it.first} ${it.second}")
            }
        }
    }

}