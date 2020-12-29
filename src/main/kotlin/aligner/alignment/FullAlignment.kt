package aligner.alignment

@Deprecated("Not to be used until implemented with IN/DELs check")
data class FullAlignment(
    val start: Int,
    val end: Int
) {

    val range = start until end

    operator fun compareTo(alignment: FullAlignment): Int {
        return if (start != alignment.start)
            start.compareTo(alignment.start)
        else
            end.compareTo(alignment.end)
    }

}