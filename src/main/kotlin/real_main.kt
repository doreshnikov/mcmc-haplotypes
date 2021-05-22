import aligner.SnipsOnlyAligner
import graph.AlignedDeBruijnGraph
import mcmc.Engine
import mcmc.modules.graph.DistributionConfig
import mcmc.modules.graph.PathsOverlay
import scoring.PrimitiveScorer
import java.io.File
import java.io.PrintWriter

fun readFasta(fileName: String): List<String> {
    File(fileName).bufferedReader().useLines { lines ->
        return@readFasta lines.toList()
            .filter { it.isNotEmpty() && it[0] != '>' && it[0] != '@' }
            .map { it.toUpperCase() }
    }
}

fun main(args: Array<String>) {

    val dir = args[0]
    val referenceFile = "$dir/_reference"
    val readsFile = "$dir/_origin"

    val reference = readFasta(referenceFile).joinToString("")
    val reads = readFasta(readsFile)

    val alignments = SnipsOnlyAligner(reference).run {
        reads.map { align(it) }
    }

}