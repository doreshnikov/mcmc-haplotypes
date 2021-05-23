package aligner

import aligner.alignment.PrimitiveAlignment
import aligner.util.MutationPolicy
import utils.matrixOf

import kotlin.math.min

class SnipsOnlyAligner(reference: Reference) : ReferenceAligner<PrimitiveAlignment>(
    reference,
    MutationPolicy.SnipsOnly
) {

    private fun alignScore(read: String): Pair<PrimitiveAlignment, Int> {
        val (n, m) = Pair(reference.length, read.length)
        val scoring = matrixOf(n + 1, m + 1) { _, _ -> 0 }

        for (i in 1..n) {
            for (j in 1..min(m, i)) {
                scoring[i][j] = scoring[i - 1][j - 1] - substitutionPenalty(reference[i - 1], read[j - 1])
            }
        }

        var end = 0
        for (i in 1..n) {
            if (scoring[i][m] >= scoring[end][m]) {
                end = i
            }
        }
        val full = scoring[end][m] == 5 * read.length
        return PrimitiveAlignment(end - read.length, read, full) to scoring[end][m]
    }

    private fun complementRead(read: String) = read.map { c ->
        when (c) {
            'A' -> 'T'
            'C' -> 'G'
            'G' -> 'C'
            'T' -> 'A'
            else -> c
        }
    }.joinToString("")

    override fun align(read: String): PrimitiveAlignment {
        val complement = complementRead(read)
        val reads = listOf(read, complement, read.reversed(), complement.reversed())
        val alignments = reads.map { alignScore(it) }

        val res = alignments.maxByOrNull { it.second }!!
//        println(res.second)
        return res.first
    }

}