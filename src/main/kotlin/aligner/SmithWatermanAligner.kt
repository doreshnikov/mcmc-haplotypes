@file:Suppress("NON_EXHAUSTIVE_WHEN")

package aligner

import aligner.alignment.FullAlignment
import aligner.util.Mutation
import aligner.util.MutationPolicy
import utils.matrixOf

@Deprecated("Not to be used until FullAlignment actually makes any sense")
class SmithWatermanAligner(reference: Reference) : ReferenceAligner<FullAlignment>(
    reference,
    MutationPolicy.NoLimits
) {

    override fun align(read: String): FullAlignment {
        val (n, m) = Pair(reference.length, read.length)
        val scoring = matrixOf(n + 1, m + 1) { _, _ -> 0 }
        val backtrack = matrixOf(n + 1, m + 1) { _, _ -> Mutation.NONE }

        fun evaluate(i: Int, j: Int, step: Mutation) = when (step) {
            Mutation.NONE -> 0
            Mutation.SNP -> scoring[i - 1][j - 1] - substitutionPenalty(reference[i - 1], read[j - 1])
            Mutation.IN -> scoring[i][j - 1] - gapPenalty
            Mutation.DEL -> scoring[i - 1][j] - gapPenalty
        }

        for (i in 1..n) {
            for (j in 1..m) {
                for (step in mutationPolicy.mutations) {
                    val newValue = evaluate(i, j, step)
                    if (newValue >= scoring[i][j]) {
                        scoring[i][j] = newValue
                        backtrack[i][j] = step
                    }
                }
            }
        }

        var (i, j) = Pair(0, m)
        for (_i in 1..n) {
            if (scoring[_i][m] > scoring[i][m]) {
                i = _i
            }
        }
        val mutations = mutableListOf<Mutation>()
        while (scoring[i][j] > 0) {
            mutations.add(backtrack[i][j])
            when (backtrack[i][j]) {
                Mutation.SNP -> {
                    i--; j--
                }
                Mutation.IN -> j--
                Mutation.DEL -> i--
            }
        }
        return FullAlignment(-1, -1)
    }

}