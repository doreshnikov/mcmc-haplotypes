package graph.components

import org.apache.commons.math3.distribution.EnumeratedIntegerDistribution

data class AlignedVertex(
    val position: Int,
    val data: String
) {

    var pathCountRight: Double = 0.0
    var pathCountLeft: Double = 0.0

    lateinit var rightSelector: EnumeratedIntegerDistribution
    lateinit var leftSelector: EnumeratedIntegerDistribution

}