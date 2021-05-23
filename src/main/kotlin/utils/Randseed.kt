package utils

import org.apache.commons.math3.distribution.*
import org.apache.commons.math3.random.RandomGenerator
import org.apache.commons.math3.random.Well19937c
import org.apache.commons.math3.util.Pair

class Randseed private constructor(private var seed: Int) {

    companion object {
        val INSTANCE = Randseed(42)
    }

    private val _rng = Well19937c(seed)
    private val rng: RandomGenerator
        get() = _rng.also {
            seed++
            _rng.setSeed(seed)
        }

    fun randint() =
        rng.nextInt()

    fun normal(mean: Double, sd: Double) =
        NormalDistribution(rng, mean, sd)

    fun uniformReal(lower: Double, upper: Double) =
        UniformRealDistribution(rng, lower, upper)

    fun uniformInteger(lower: Int, upper: Int) =
        UniformIntegerDistribution(rng, lower, upper)

    fun poisson(lambda: Double) =
        PoissonDistribution(
            rng, lambda,
            PoissonDistribution.DEFAULT_EPSILON, PoissonDistribution.DEFAULT_MAX_ITERATIONS
        )

    fun <T> probEnum(data: List<T>, probs: List<Double>) =
        EnumeratedDistribution(rng, data.zip(probs) { d, p -> Pair(d, p) })

    fun <T> probEnum(data: List<kotlin.Pair<T, Double>>) =
        EnumeratedDistribution(rng, data.map { (d, p) -> Pair(d, p) })

    fun <T> scoreEnum(data: List<T>, scores: List<Double>): EnumeratedDistribution<T> {
        val totalScore = scores.sum()
        check(totalScore != 0.0) { "At least one option should be available" }
        return EnumeratedDistribution(
            data.zip(scores) { i, s -> Pair(i, s / totalScore) }
        )
    }

    fun probIntegerEnum(probs: List<Double>) =
        EnumeratedIntegerDistribution(rng, IntArray(probs.size) { it }, probs.toDoubleArray())

    fun scoreIntegerEnum(scores: List<Double>): EnumeratedIntegerDistribution {
        val totalScore = scores.sum()
        check(totalScore != 0.0) {
            "At least one option should be available"
        }
        return probIntegerEnum(scores.map { s -> s / totalScore })
    }

}

