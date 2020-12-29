package aligner.alignment

data class PrimitiveAlignment(
    val start: Int,
    val read: String
) {

    val end get() = start + read.length
    val indices get() = start until end

}
