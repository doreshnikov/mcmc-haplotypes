package scoring

interface Scorer<T> {

    fun score(
        expect: List<Pair<String, Double>>,
        actual: List<Pair<String, Double>>
    ): T

}