package graph.components

open class Edge(
    val source: Int,
    val target: Int
) {

    operator fun component1(): Int = source
    operator fun component2(): Int = target

    override fun hashCode(): Int {
        return source.hashCode() xor (target.hashCode() shl 7)
    }

    override fun equals(other: Any?): Boolean {
        if (this === other) return true
        if (javaClass != other?.javaClass) return false

        other as Edge

        if (source != other.source) return false
        if (target != other.target) return false

        return true
    }

}
