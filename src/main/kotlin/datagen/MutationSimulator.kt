package datagen

import aligner.Reference
import utils.Randseed

class MutationSimulator(
    val reference: Reference,
    val mutationRate: Double
) {

    private val mutationNumberGenerator = Randseed.INSTANCE.poisson(mutationRate)
    private val indexGenerator = Randseed.INSTANCE.uniformInteger(0, reference.length - 1)

    fun generate(number: Int): List<Reference> {
        val result = mutableListOf(reference)
        repeat(number) {
            val origin = result[Randseed.INSTANCE.uniformInteger(0, result.size - 1).sample()]
            val item = origin.toCharArray().also { item ->
                repeat(mutationNumberGenerator.sample()) {
                    var newNucleotide = generateNucleotide()
                    val index = indexGenerator.sample()
                    while (newNucleotide == item[index]) {
                        newNucleotide = generateNucleotide()
                    }
                    item[index] = newNucleotide
                }
            }.concatToString()
            result.add(item)
        }
        return result
    }

}