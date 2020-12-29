package graph

import ManualCheckTest
import aligner.SnipsOnlyAligner

class AlignmentTest(
    private val reference: String,
    private val reads: List<String>
) : ManualCheckTest("Alignment test", {
    println("Reference: $reference")
    val aligner = SnipsOnlyAligner(reference)
    reads.forEach { read ->
        println("$read -> ${aligner.align(read)}")
    }
})