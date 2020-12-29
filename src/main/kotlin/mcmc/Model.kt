package mcmc

abstract class Model<E : Any, D : Model<E, D>.Delta>(
    val entity: E
) {

    abstract inner class Delta {
        abstract fun logLikelihoodDelta(): Double
        abstract fun accept()
    }

    abstract fun logLikelihood(): Double
    abstract fun proposeCandidate(): Delta

}