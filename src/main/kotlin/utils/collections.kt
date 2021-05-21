package utils

inline fun <reified T> matrixOf(n: Int, m: Int, init: (Int, Int) -> T): MutableList<MutableList<T>> {
    return MutableList(n) { i -> MutableList(m) { j -> init(i, j) } }
}