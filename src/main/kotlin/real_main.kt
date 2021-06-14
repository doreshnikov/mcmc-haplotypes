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
import kotlin.math.sqrt
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
                if (s + read.length <= reference.length) {
                    alignments.add(
                        PrimitiveAlignment(
                            s, read,
                            read == reference.subSequence(s, s + read.length)
                        )
                    )
                }
            }
        }
    }

    val k = 7
    val sortedAlignments = alignments.sortedBy { it.start }
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
    println("Reference check: ${if (failed) "empty or failed" else "passed"}")
//    if (failed) return

    val lambda = 3.0
    val graph = AlignedDeBruijnGraph.build(alignments, k)
    val oldPaths = graph.validate()
    val cutGraph = graph.cutErrorTails()
    println("Paths after cut: $oldPaths+${cutGraph.validate() - oldPaths}")
    println("Vertices count: ${cutGraph.size}")
    val normalizedGraph = cutGraph.Normalizer().normalize()
//    normalizedGraph.precalcPathCounts()

    var results: List<List<Pair<String, Double>>>
    val traceBest = true
    val model = PathsOverlay(normalizedGraph, DistributionConfig(lambda), 0.001, 3 * (k + 1))
//    val model = PathsOverlay(normalizedGraph, DistributionConfig(lambda), 0.001, 4)
    val engine = Engine(normalizedGraph, model)

    if (oldPaths == 1L) {
        results = listOf(listOf(Pair(reference, 1.0)))
    } else {
        try {
            engine.simulate(
                1, timeLimit = 6 * 60 * 60 * 1000, criteria = "tl",
                verboseLevel = -100, traceBest = traceBest, temp = { 1.0 }
            )
        } catch (e: Exception) {
            println(e.message)
            e.printStackTrace()
        } finally {
            model.log.close()
        }

        results = listOf(model.extractResult(), engine.bestResult!!)
    }
    val suf = listOf("last", "opt")
    results.forEachIndexed { idx, result ->
        if (!traceBest && idx > 0) return@forEachIndexed
        var name = "${dir}/my_${suf[idx]}"
        if (model.SHIFT_DELETE != 0.0 || model.SHIFT_NEW != 0.0) {
            name += "_scale"
        }
        File(name).printWriter().use { pw ->
            result.forEach {
                pw.println("${it.first} ${it.second}")
            }
        }
    }

    File("${dir}/_loglhistory").printWriter().use { pw ->
        engine.logLHistory.forEach { pw.println(it) }
    }
    println("Acceptance rate: ${engine.accepted.toDouble() / (engine.iter + 1)}")

}