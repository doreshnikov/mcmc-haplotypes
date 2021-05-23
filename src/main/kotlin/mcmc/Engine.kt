package mcmc

import mcmc.modules.graph.PathsOverlay
import utils.Randseed
import kotlin.math.ln

open class Engine<E : Any, R : Any>(
    val entity: E,
    model: Model<E, *, R>
) {

    private val distribution = Randseed.INSTANCE.uniformReal(0.0, 1.0)
    var model: Model<E, *, R> = model
        private set

    var bestResult: R? = null
    var bestLogLikelihood: Double = Double.NEGATIVE_INFINITY

    fun simulate(
        iterations: Int, timeLimit: Int? = null, criteria: String = "iter",
        verboseLevel: Int = 0, traceBest: Boolean = false
    ) {
        val startTime = System.currentTimeMillis()
        var iterLimit = iterations
        var iter = 0
        while (iter < iterLimit) {
            val candidate = model.proposeCandidate()
            val logLD = candidate.logLikelihoodDelta
            val logJD = candidate.logJumpDensity
            val logAcceptanceRate = logLD + logJD
            val accept = ln(distribution.sample()) < logAcceptanceRate
            if (accept) {
                candidate.accept()
                if (traceBest) {
                    val logL = model.logLikelihood()
                    if (logL > bestLogLikelihood) {
                        bestLogLikelihood = logL
                        bestResult = model.extractResult()
                    }
                }
            }
            val timePassed = System.currentTimeMillis() - startTime

            if (verboseLevel != 0) {
                val report = """Iteration $iter:
 - logL : ${model.logLikelihood()}
 - paths: ${(model as PathsOverlay).paths.size}
 - time : ${timePassed / 1000}s
 - logLD: $logLD
 - logJD: $logJD
 - acc  : $accept ${candidate.javaClass.canonicalName.split(".").last()}"""
                when (verboseLevel) {
                    1 -> if (accept) println(report)
                    2 -> println(report)
                    else -> if ((iter + 1) % -verboseLevel == 0) {
                        println(
                            """Iteration ${(iter + 1) / -verboseLevel} x ${-verboseLevel}:
 - logL : ${model.logLikelihood()}
 - paths: ${(model as PathsOverlay).paths.size}"""
                        )
                    }
                }
            }

            iter++
            when (criteria) {
                "tl" -> {
                    if (timePassed >= timeLimit!!) return
                    if (iter == iterLimit) iterLimit *= 2
                }
                "both" -> if (iter == iterLimit) {
                    if (timePassed >= timeLimit!!) return
                    iterLimit++
                }
                "any" -> {
                    if (timePassed >= timeLimit!!) return
                }
            }
        }
    }

}