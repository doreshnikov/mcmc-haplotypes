package datagen

import utils.Randseed

val nucleotides = listOf('A', 'C', 'G', 'T')
val nucleotidesDistribution = Randseed.INSTANCE.uniformInteger(0, nucleotides.size - 1)

fun generateNucleotide(): Char {
    return nucleotides[nucleotidesDistribution.sample()]
}

fun generateReference(length: Int): String {
    return CharArray(length) { generateNucleotide() }
        .concatToString()
}