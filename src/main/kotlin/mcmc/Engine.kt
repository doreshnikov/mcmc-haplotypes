package mcmc

import org.apache.commons.math3.distribution.UniformRealDistribution
import kotlin.math.ln

open class Engine<E : Any>(
    val entity: E,
    model: Model<E, *>
) {

    private val distribution = UniformRealDistribution(0.0, 1.0)
    var model: Model<E, *> = model
        private set

    fun simulate(iterations: Int) {
        repeat(iterations) {
            val candidate = model.proposeCandidate()
            if (candidate.logLikelihoodDelta() <= ln(distribution.sample())) {
                candidate.accept()
            }
            println("Iteration $it: log likelihood ${model.logLikelihood()}")
        }
    }

}