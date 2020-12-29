@file:Suppress("NON_EXHAUSTIVE_WHEN")

package aligner

import aligner.util.Mutation
import aligner.util.MutationPolicy
import utils.matrixOf

abstract class ReferenceAligner<A : Any>(
    val reference: Reference,
    val mutationPolicy: MutationPolicy
) {

    protected companion object {
        const val gapPenalty = 2

        fun idx(c: Char) = when (c) {
            'A' -> 0
            'C' -> 1
            'G' -> 2
            'T' -> 3
            else -> error("Unknown nucleotide")
        }

        val substitutionMatrix = arrayOf(
            arrayOf(-3, 3, 3, 3),
            arrayOf(3, -3, 3, 3),
            arrayOf(3, 3, -3, 3),
            arrayOf(3, 3, 3, -3)
        )

        fun substitutionPenalty(c1: Char, c2: Char) =
            substitutionMatrix[idx(c1)][idx(c2)]
    }

    abstract fun align(read: String): A

}