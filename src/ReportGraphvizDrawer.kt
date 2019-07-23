package edu.nus

import java.io.File
import kotlin.random.Random

typealias KeyValue = Pair<String, String>
class ReportGraphvizDrawer(private val file: File = File("C:/Users/tao/Desktop/report")) {
    private val spaceChars = """\s+""".toRegex()
    data class Segment(val identifier:String, val forward:Boolean){
        constructor(identifier:String): this(identifier.substring(1), identifier[0]=='+')
        companion object{
            private val SegmentPattern = Regex("""([^(]+)\((\d+),(\d+)\)\((\d+),(\d+)\)""") // 1:contig, 2,3,4,5:interval
            val pool = HashMap<String,Int>()
            fun getNodeID(identifier:String) = pool.getOrPut(identifier) {pool.size}
            fun beautifyID(identifier:String): String {
                val groups = SegmentPattern.find(identifier)!!.groups.map { it!!.value }
                return groups[1] + if (groups[2]==groups[3] && groups[4]==groups[5]) "[${groups[5]}]" else "[${groups[2]},${groups[3]}][${groups[4]},${groups[5]}]"
            }
        }
        fun port() = "node${getNodeID(identifier)}:${if (forward) "w" else "e"}"
    }
    class Graph(name:String = ""){
        private var dot = "graph $name {\n"
        private fun attrToDOT(kv: Array<out KeyValue>) = kv.joinToString(separator = " ") { "[${it.first}=${it.second}]" }
        fun setNodes(vararg kv: KeyValue) = addNode("node", *kv)
        fun setEdges(vararg kv: KeyValue) = addNode("edge", *kv)
        fun setGraph(vararg kv: KeyValue) = addNode("graph", *kv)
        fun addNode(name:String, vararg kv: KeyValue) { dot += "\t$name ${attrToDOT(kv)}\n" }
        fun addEdge(source:String, target: String, vararg kv: KeyValue) { dot += "\t$source -- $target ${attrToDOT(kv)}\n" }
        fun toDOT() = "$dot}"
    }
    operator fun String.minus(other: String) = Pair(this,other)
    private fun randomColorString() = "\"#%2x%2x%2X\"".format(100+Random.nextInt(100),100+Random.nextInt(100),100+Random.nextInt(100))
    fun draw(){
        val graph = Graph()
        graph.setNodes("shape"-"box")
        file.forEachLine{ line ->
            if (!line.startsWith("+") && !line.startsWith("-")) return@forEachLine
            val tokens = line.split(spaceChars)
            val ref = Segment(tokens[0])
            val qry = Segment(tokens[1])
            val s = "%.2f".format(tokens[2].toDouble())
            val randomColor = randomColorString()
            graph.addEdge(ref.port() , qry.port(), "label"-s, "color"-randomColor, "fontcolor"-randomColor)
        }
        for ((node, id) in Segment.pool) graph.addNode("node$id","label"-"\"${Segment.beautifyID(node)}\"")
        println(graph.toDOT())
//        Runtime.getRuntime().exec("dot -Tpng scaffold.dot -o scaffold.png").waitFor()
    }
}