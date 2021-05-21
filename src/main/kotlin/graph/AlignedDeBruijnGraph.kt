package graph

import aligner.alignment.PrimitiveAlignment
import graph.components.AdjacencyList
import graph.components.WeightedEdge
import graph.components.AlignedVertex
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation
import kotlin.math.pow

class AlignedDeBruijnGraph private constructor(
    vertices: List<AlignedVertex>,
    edges: AdjacencyList,
    val originLength: Int,
    val k: Int
) : Graph(vertices, edges) {

    companion object {

        private class VertexFactory {
            private val factory = hashMapOf<AlignedVertex, Int>()
            private val vertices = mutableListOf<AlignedVertex>()

            operator fun get(position: Int, data: String): Int {
                val vertex = AlignedVertex(position, data)
                return factory.getOrPut(vertex) {
                    vertices.add(vertex)
                    vertices.size - 1
                }
            }

            operator fun invoke() = vertices
        }

        private class EdgeFactory {
            private val factory = hashMapOf<Int, HashMap<Int, Int>>()

            operator fun plusAssign(edge: Pair<Int, Int>) {
                factory
                    .getOrPut(edge.first) { hashMapOf() }
                    .compute(edge.second) { _, v -> (v ?: 0) + 1 }
            }

            operator fun invoke(totalVertices: Int): AdjacencyList {
                val edges = MutableList(totalVertices) { mutableListOf<WeightedEdge>() }
                for (v in 0 until totalVertices) {
                    for (u in factory.getOrElse(v) { emptyMap() }) {
                        edges[v].add(WeightedEdge(v, u.key, u.value.toDouble()))
                    }
                }
                return edges
            }
        }

        fun build(alignments: List<PrimitiveAlignment>, k: Int): AlignedDeBruijnGraph {
            val vertexFactory = VertexFactory()
            val edgeFactory = EdgeFactory()

            for (alignment in alignments) {
                val verticesIds = mutableListOf<Int>()
                for (end in k..alignment.read.length) {
                    val start = end - k
                    verticesIds.add(
                        vertexFactory[
                                start + alignment.start,
                                alignment.read.substring(start until end)
                        ]
                    )
                }
                for (v in 1 until verticesIds.size) {
                    edgeFactory += verticesIds[v - 1] to verticesIds[v]
                }
            }

            val vertices = vertexFactory()
            return AlignedDeBruijnGraph(
                vertices, edgeFactory(vertices.size),
                alignments.maxOf { it.end }, k
            )
        }

    }

    fun isRightmost(v: Int): Boolean {
        return get(v).position + k == originLength
    }

    val sources = indices.filter { get(it).position == 0 }
    val targets = indices.filter { isRightmost(it) }

    @Suppress("DuplicatedCode")
    fun cutErrorTails(): AlignedDeBruijnGraph {
        val failIds = ArrayDeque<Int>()
        val failed = MutableList(size) { false }
        val rm = MutableList(size) { 0 }
        val reversedRm = MutableList(size) { 0 }

        for (id in indices) {
            if (get(id).position > 0 && get(id).position + k < originLength) {
                if (edges[id].isEmpty() || reversedEdges[id].isEmpty()) {
                    failIds.addLast(id)
                    failed[id] = true
                }
            }
        }

        while (failIds.isNotEmpty()) {
            val id = failIds.removeFirst()
            for (edge in edges[id]) {
                reversedRm[edge.target]++
                if (!failed[edge.target] && reversedEdges[edge.target].size == reversedRm[edge.target]) {
                    failIds.addLast(edge.target)
                    failed[edge.target] = true
                }
            }
            for (edge in reversedEdges[id]) {
                rm[edge.target]++
                if (!failed[edge.target] && edges[edge.target].size == rm[edge.target]) {
                    failIds.addLast(edge.target)
                    failed[edge.target] = true
                }
            }
        }

        val idsMapper = HashMap<Int, Int>()
        for (id in indices) {
            if (!failed[id]) {
                idsMapper[id] = idsMapper.size
            }
        }
        val newEdges = MutableList(idsMapper.size) { mutableListOf<WeightedEdge>() }
        for (source in indices) {
            if (!failed[source]) {
                for (edge in edges[source]) {
                    if (!failed[edge.target]) {
                        val s = idsMapper[source]!!
                        val t = idsMapper[edge.target]!!
                        newEdges[s].add(WeightedEdge(s, t, edge.weight))
                    }
                }
            }
        }
        return AlignedDeBruijnGraph(
            this.filterIndexed { i, v -> !failed[i] }, newEdges,
            originLength, k
        )
    }

    fun validate(): Int {
        val visited = MutableList(size) { false }
        val ok = MutableList(size) { 0 }

        fun dfs(v: Int) {
            visited[v] = true
            if (get(v).position + k == originLength) {
                ok[v] = 1
            }

            for (edge in edges[v]) {
                if (!visited[edge.target]) {
                    dfs(edge.target)
                }
                ok[v] += ok[edge.target]
            }
        }

        var answer = 0
        for (v in indices) {
            if (!visited[v]) {
                dfs(v)
                if (get(v).position == 0) {
                    answer += ok[v]
                }
            }
        }
        return answer
    }

    inner class Normalizer {
        private var pathError = 0.0
        private var bucketError = 0.0

        val error get() = pathError + bucketError

        private val inDegree = MutableList(size) { 0 }

        init {
            for (v in indices) {
                for (e in edges[v]) {
                    inDegree[e.target]++
                }
            }
        }

        @Deprecated("There is no point in doing 2 additional normalizations to artificially reduce error")
        private fun horizontalNormalization(localEdges: AdjacencyList): AdjacencyList {
            val paths = mutableListOf<MutableList<WeightedEdge>>()
            val visited = MutableList(size) { false }
            fun dfs(edge: WeightedEdge, pathId: Int) {
                val v = edge.target
                if (edge.source > -1) {
                    paths[pathId].add(edge)
                }
                if (visited[v]) {
                    return
                }
                visited[v] = true
                val keepPath = inDegree[v] == 1 && localEdges[v].size == 1
                localEdges[v].forEach { u ->
                    val nextPathId = if (keepPath) pathId else paths.size
                    if (nextPathId >= paths.size) {
                        paths.add(mutableListOf())
                    }
                    dfs(u, nextPathId)
                }
            }

            val sources = indices.filter { inDegree[it] == 0 }
            sources.forEach {
                paths.add(mutableListOf())
                dfs(WeightedEdge(-1, it, 0.0), paths.size - 1)
            }

            val pathAverages = paths.map { it.sumByDouble { e -> e.weight } / it.size }
            pathError = 0.0
            val newEdges = MutableList(size) { mutableListOf<WeightedEdge>() }
            paths.forEachIndexed { i, path ->
                path.forEach { edge ->
                    newEdges[edge.source].add(edge.copy(weight = pathAverages[i]))
                    pathError += (pathAverages[i] - edge.weight).pow(2)
                }
            }
            return newEdges
        }

        private fun verticalNormalization(localEdges: AdjacencyList): AdjacencyList {
            val buckets = mutableListOf<MutableList<WeightedEdge>>()
            val vertexDepth = MutableList(size) { -1 }
            fun dfs(edge: WeightedEdge, depth: Int) {
                val v = edge.target
                if (depth > -1) {
                    while (depth >= buckets.size) {
                        buckets.add(mutableListOf())
                    }
                    buckets[depth].add(edge)
                }
                if (vertexDepth[v] > -1) {
                    return
                }
                vertexDepth[v] = depth + 1
                localEdges[v].forEach { u ->
                    dfs(u, depth + 1)
                }
            }

            val sources = indices.filter { inDegree[it] == 0 }
            sources.forEach {
                dfs(WeightedEdge(-1, it, 0.0), -1)
            }

            val bucketTotals = buckets.map { it.sumByDouble { e -> e.weight } }
            bucketError = StandardDeviation().evaluate(bucketTotals.toDoubleArray())
            val newEdges = MutableList(size) { mutableListOf<WeightedEdge>() }
            buckets.forEach { bucket ->
                bucket.forEach { edge ->
                    val newWeight = edge.weight / bucketTotals[vertexDepth[edge.source]]
                    newEdges[edge.source].add(edge.copy(weight = newWeight))
                }
            }
            return newEdges
        }

        fun normalize(): AlignedDeBruijnGraph {
//            val pathNormalizedEdges = horizontalNormalization(edges)
            val pathNormalizedEdges = edges
            val bucketNormalizedEdges = verticalNormalization(pathNormalizedEdges)
            return AlignedDeBruijnGraph(
                this@AlignedDeBruijnGraph, bucketNormalizedEdges,
                originLength, k
            )
        }
    }

}