package aligner.util

sealed class MutationPolicy {

    abstract val mutations: List<Mutation>

    object SnipsOnly : MutationPolicy() {
        override val mutations = listOf(Mutation.NONE, Mutation.SNP)
    }

    object NoLimits : MutationPolicy() {
        override val mutations = Mutation.values().toList()
    }

}
