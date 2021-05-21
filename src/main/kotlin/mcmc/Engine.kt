package mcmc

import utils.Randseed
import kotlin.math.exp
import kotlin.math.ln

open class Engine<E : Any>(
    val entity: E,
    model: Model<E, *>
) {

    private val distribution = Randseed.INSTANCE.uniformReal(0.0, 1.0)
    var model: Model<E, *> = model
        private set

    fun simulate(iterations: Int) {
        repeat(iterations) {
            val candidate = model.proposeCandidate()
            val logAcceptanceRate = candidate.logLikelihoodDelta + candidate.logJumpDensity
            if (ln(distribution.sample()) < logAcceptanceRate) {
                candidate.accept()
                println("Iteration $it: log likelihood ${model.logLikelihood()}")
            } else {
                println("Iteration $it: candidate rejected with AR of ${exp(logAcceptanceRate)}")
            }
        }
    }

}