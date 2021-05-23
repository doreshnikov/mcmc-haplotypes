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

fun main(args: Array<String>) {

    val dir = args[0]
    val referenceFile = "$dir/_reference"
    val readsFile = "$dir/_reads"

    val reference = readFasta(referenceFile).joinToString("")
    val alignments = mutableListOf<PrimitiveAlignment>()

    if (!File("$dir/_aligned").exists()) {
        val reads = File(readsFile).bufferedReader().lines().toList()
        val lock = ReentrantLock()
        println("Aligning ${reads.size} reads...")

        val time = runBlocking {
            measureTimeMillis {
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
        }

        println("Total time: ${time / 1000}s, average ${time / reads.size}ms/read")
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
                val s = start.toInt()
                alignments.add(PrimitiveAlignment(s, read, read == reference.subSequence(s, s + read.length)))
            }
        }
    }

    val k = 6

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
                            sortedAlignments[i].read.substring(0 until k - 1), false
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
                            reference.substring(end until reference.length), false
                )
            )
        }
    }

    val referenceFilter = sortedAlignments.filter { it.full }
    var failed = false
    for (i in 1 until referenceFilter.size) {
        if (referenceFilter[i].start > referenceFilter[i - 1].end - k + 1) {
            fun repr(al: PrimitiveAlignment): String {
                return "${al.start}..${al.end}"
            }
            println("Fail from ${repr(referenceFilter[i - 1])} to ${repr(referenceFilter[i])}")
            failed = true
        }
    }
    if (failed) return

    val lambda = 200.0
    val graph = AlignedDeBruijnGraph.build(alignments, k)
    println("Paths (mod MAX_LONG): ${graph.validate()}")
    val cutGraph = graph.cutErrorTails()
    println("Paths after cutting (mod MAX_LONG): ${cutGraph.validate()}")
    val normalizedGraph = cutGraph.Normalizer().normalize()
    normalizedGraph.precalcPathCounts()

    val model = PathsOverlay(normalizedGraph, DistributionConfig(lambda), 0.001)
    val engine = Engine(normalizedGraph, model)
    val traceBest = true
    engine.simulate(
        1, timeLimit = 20 * 60 * 1000, criteria = "tl",
        verboseLevel = -1000, traceBest = traceBest
    )

    val results = listOf(model.extractResult(), engine.bestResult!!)
    val suf = listOf("last", "opt")
    results.forEachIndexed { idx, result ->
        if (!traceBest && idx > 0) return@forEachIndexed
        File("${dir}/my_${suf[idx]}").printWriter().use { pw ->
            result.forEach {
                pw.println("${it.first} ${it.second}")
            }
        }
    }

}