package datagen

import org.apache.commons.math3.distribution.IntegerDistribution
import utils.Randseed
import kotlin.math.min

class ReadSimulator(
    sequences: List<Pair<String, Double>>,
    val readNumber: Int,
    val readLengthDistribution: IntegerDistribution,
    val readErrorRate: Double
) {

    private val sequenceDistribution = Randseed.INSTANCE.probEnum(sequences)

    private val readStartDistribution = Randseed.INSTANCE.uniformInteger(0, sequences[0].first.length - 1)
    private val errorDistribution = Randseed.INSTANCE.uniformReal(0.0, 1.0)

    fun simulate(): String {
        val reference = sequenceDistribution.sample()
        val start = readStartDistribution.sample()
        val length = min(readLengthDistribution.sample(), reference.length - start)
        val read = reference.substring(start, start + length).toCharArray()
        read.indices.forEach { i ->
            if (errorDistribution.sample() < readErrorRate) {
                read[i] = generateNucleotide()
            }
        }
        return read.concatToString()
    }

    fun simulateAll(): List<String> {
        return MutableList(readNumber) { simulate() }
    }

}