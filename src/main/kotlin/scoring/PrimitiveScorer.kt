package scoring

class PrimitiveScorer : Scorer<List<List<Int>>> {

    private fun distance(s: String, t: String): Int {
        return s.length - s.mapIndexed { i, c ->
            if (c == t[i]) 1 else 0
        }.sum()
    }

    override fun score(
        expect: List<Pair<String, Double>>,
        actual: List<Pair<String, Double>>
    ): List<List<Int>> {
        val res = MutableList(expect.size) { MutableList(actual.size) { 0 } }
        for (i in expect.indices) {
            for (j in actual.indices) {
                res[i][j] = distance(expect[i].first, actual[j].first)
            }
        }
        return res
    }

}