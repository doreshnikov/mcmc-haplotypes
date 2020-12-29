open class ManualCheckTest(val name: String, val action: () -> Any) {

    init {
        @Suppress("LeakingThis")
        allTests.add(this)
    }

    fun execute() {
        println("[${name}]")
        try {
            val result = action()
            if (result !is Unit) {
                println("Result: $result")
            }
        } catch (e: Exception) {
            println("Error: $e")
            e.printStackTrace()
        }
        println("-".repeat(20))
    }

    companion object {
        private val allTests = mutableListOf<ManualCheckTest>()

        fun executeAll() {
            allTests.forEach(ManualCheckTest::execute)
        }
    }

}