package edu.nus

import me.tongfei.progressbar.ProgressBar
import htsjdk.samtools.reference.*
import java.io.File

class Edge(val overhang:Int, val overhangOverlap:Int, val s:Double, val target:SiteNode)
class SiteNode(val forward: Int, val contigName: String, private val contigSeq: String,
               val pos: Int, var flex: Int = 0, var partner:SiteNode?=null):ArrayList<Edge>() {
    var link:Edge? = null
    fun extpos() = pos - forward * flex
    private fun leftEx() = if (forward==1) extpos() else partner!!.extpos()
    private fun rightEx() = if (forward==-1) extpos() else partner!!.extpos()
    fun prune(len: Int) { flex -= len }
    fun seq(): String {
        val rawseq = contigSeq.substring(leftEx() - 1, rightEx())
        return if (forward == 1) rawseq else rawseq.reversed().map {
            Util.complementary.getOrDefault(it, 'N')
        }.joinToString(separator = "")
    }
    fun len() = Math.abs(partner!!.extpos()-extpos())+1
    fun getName(): String {
        val l = if (leftEx()==1) "S" else leftEx().toString()
        val r = if (rightEx()==contigSeq.length) "E" else rightEx().toString()
        return if (forward==1) "$contigName(${l}_$r)" else "$contigName(${r}_$l)"
    }
    fun lossSeq(): String = if (forward==1) contigSeq.substring(pos-1, extpos()) else contigSeq.substring(extpos()-1, pos)
    override fun hashCode() = "$contigName($pos)".hashCode()
    override fun equals(other: Any?) = hashCode()==other.hashCode()
}
class NodePool(private val contigs: MutableMap<String, String>) :HashMap<String, SiteNode>(){
    private val SegmentPattern = Regex("""([+-])([^(]+)\((\d+),(\d+)\)\((\d+),(\d+)\)""") // 1:end, 2:contig, 3,4,5,6:interval
    fun getBySegment(segment:String): SiteNode {
        val groups = SegmentPattern.find(segment)!!.groups.map { it!!.value }
        val name1 = "+${groups[2]}(${groups[4]})"
        val name2 = "-${groups[2]}(${groups[5]})"
        val name = if (groups[1]=="+") name1 else name2
        if (this.containsKey(name)) return this[name]!!
        val node1 = SiteNode(1, groups[2], contigs[groups[2]]!!,
            groups[4].toInt(), groups[4].toInt()-groups[3].toInt())
        val node2 = SiteNode(-1, groups[2], contigs[groups[2]]!!,
            groups[5].toInt(), groups[6].toInt()-groups[5].toInt())
        node1.partner = node2
        node2.partner = node1
        this[name1] = node1
        this[name2] = node2
        return if (groups[1]=="+") node1 else node2
    }
    fun getNode(forward: Int, contig: String, pos: Int, flex: Int): SiteNode {
        val name = if (forward==1) "+$contig($pos)" else "-$contig($pos)"
        if (!this.containsKey(name))
            this[name] = SiteNode(forward, contig, contigs[contig]!!, pos, flex)
        return this[name]!!
    }
    fun remove(node:SiteNode) = this.remove("${if (node.forward==1) '+' else '-'}${node.contigName}(${node.pos})")
}
class LinkageGroupDock: HashMap<SiteNode,SiteNode>() {
    var LGCounter = 0
    private fun find(node: SiteNode): SiteNode{
        var s = this[node]!!
        if (s==node) return s
        s = this.find(s)
        this[node] = s
        return s
    }
    private fun update(node: SiteNode, edge: Edge){
        node.link = edge
        edge.target.link = edge.target.find { it.target==node }
        val A = find(edge.target)
        val B = find(node)
        if (A==B) return
        this[B] = A
        LGCounter -= 1
    }
    fun place(node: SiteNode) {
        val partner = node.partner!!
        this[node] = node
        this[partner] = node
        LGCounter += 1
        val left = node.filter { this.contains(it.target) } .maxBy { it.s }
        val right = partner.filter { this.contains(it.target) } .maxBy { it.s }
//        println("=== ${node.getName()}\n" +
//                "left=${left?.target?.getName()} left.link=${left?.target?.link}\n" +
//                "right=${right?.target?.getName()} right.link=${right?.target?.link}\n")
        if (left!=null && right==null && left.target.link==null) update(node, left)
        if (left==null && right!=null && right.target.link==null) update(partner, right)
        if (left!=null && right!=null){
            if (left.target.link?.target==right.target && left.target.link!!.overhangOverlap==0 // locate between 2 segments
                || left.target.link==null && right.target.link==null && find(left.target)!=find(right.target)) { // link to 2 lg
                    update(node, left)
                    update(partner, right)
                if (left.target.link?.target==right.target) LGCounter += 1
            }
        }
    }
}
class Scaffolder(draftFile:File, private val reportFile:File, outDir:File) {
    private val javaRuntime = Runtime.getRuntime()
    private val blastTmpDir = File(outDir, "blast_tmp")
    private val tempRef = File(blastTmpDir, "ref.fasta")
    private val tempQry = File(blastTmpDir, "qry.fasta")
    private val lossSeq by lazy{ java.io.PrintWriter(File(blastTmpDir, "LossSeq.fasta"))}
    private val scaffoldInfo by lazy{ java.io.PrintWriter(File(blastTmpDir, "scaffold_info.log"))}
    private val outFasta = FastaWriter(File(outDir, "CAST.fasta"))
    private val contigs = loadFasta(draftFile) //Draft Genome
    private val nodePool = NodePool(contigs) //global dict of SiteNode for no duplicate and identity

    private fun BLASTSeqPair(ref:String, qry:String): List<List<Int>> {
        tempRef.writeText(">reference\n$ref\n")
        tempQry.writeText(">query\n$qry\n")
        javaRuntime.exec("makeblastdb -in ${tempRef.absolutePath} -dbtype nucl").waitFor()
        val cmd = javaRuntime.exec(arrayOf("blastn", "-num_threads", Runtime.getRuntime().availableProcessors().toString(), "-db", tempRef.absolutePath, "-query", tempQry.absolutePath, "-outfmt", "6 sstart send qstart qend"))
        return cmd.inputStream.bufferedReader().readLines().map{ line -> line.split('\t').map{it.toInt() } }
        // ref XXXXXXXXXXX
        //            XXXXXXXXX query
    }
    private fun anchorByBLAST(ref: SiteNode, qry: SiteNode): List<Int> {
        val hsps = BLASTSeqPair(ref.partner!!.seq(),qry.seq())
        val valid = hsps.filter {hsp -> hsp[0]<hsp[1] && hsp[2]<hsp[3]}
        if (valid.isEmpty()) return listOf(0,0,0,0) // no valid alignment
        var best = valid.minBy{ Math.max(ref.len()-ref.flex-it[1],0) + Math.max(it[2]-1-qry.flex,0)}!!
        val loss = Math.max(ref.len()-ref.flex-best[1],0) + Math.max(best[2]-1-qry.flex,0)
        scaffoldInfo.write("[BLAST] ${ref.getName()} ${qry.getName()} loss=$loss"
                + best.joinToString(",", " (",") ")
                + if (loss<=best[0]-best[1]+1) "PASS\n" else " FAIL\n")
        best = listOf(ref.len()-best[0], ref.len()-best[1], best[2]-1, best[3]-1)
        return if (loss<=best[0]-best[1]+1) best else listOf(0,0,0,0)
        // best=(ref_overlap+overhang, ref_overhang, qry_overhang, qry_overhang+overlap)
    }

    private fun mergeAndWrite() {
        val visit = mutableSetOf<SiteNode>()
        fun traverse(start: SiteNode){
            val contig = ContigBuffer()
            var node = start
            while (true) {
                val partner = node.partner!!
                visit.add(node)
                visit.add(partner)
                if (partner.link==null) break
                val link = partner.link!!
                val next = link.target
                val nextlink = next.link!!
                val overlapSep = next.seq().substring(nextlink.overhang, link.overhangOverlap)
                val overlap = nextlink.overhangOverlap-nextlink.overhang
                partner.prune(link.overhang)
                next.prune(nextlink.overhangOverlap)
                val leftloss = Math.max(0,-partner.flex)
                val rightloss = Math.max(0,-next.flex)
                contig.append(node.getName(), node.seq())
                if (leftloss>0) lossSeq.write(">${node.getName()}_${next.getName()}_LeftLoss$leftloss\n${node.lossSeq()}\n")
                if (rightloss>0) lossSeq.write(">${node.getName()}_${next.getName()}_RightLoss$rightloss\n${next.lossSeq()}\n")
                if (overlap>0) lossSeq.write(">${node.getName()}_${next.getName()}_Overlap$overlap\n$overlapSep\n")
                else contig.append("", "NNNNNNNNNN".repeat(5))
                node = next
            }
            contig.append(node.getName(), node.seq())
            outFasta.writeRecord(contig.name.toString(), contig.buffer.toString())
        }
        for ( node in nodePool.values) if (node.link==null && !visit.contains(node)) traverse(node) // for chain
        for ( node in nodePool.values) if (!visit.contains(node)) println("[ERROR]") // traverse(node) for circle
        nodePool.values.forEach{s -> contigs.remove(s.contigName) }
        contigs.forEach{(k,v) -> outFasta.writeRecord(k,v) }
    }
    private fun placement() {
        val lgDock = LinkageGroupDock()
        val siteNodes = nodePool.values.filter { it.forward==1 }.sortedByDescending { it.len() }
        for (node in siteNodes)
            lgDock.place(node)
        println("${siteNodes.size} segments are merged into ${lgDock.LGCounter} scaffolds")
    }
    private fun rescueDebris() { //to double check\
        val contigGroup = nodePool.values.groupBy{it.contigName}
        contigGroup.forEach { (contig, siteList) ->
            val tiglen = contigs[contig]!!.length
            val sites = siteList.sortedBy { it.pos }
            var lastNode: SiteNode = nodePool.getNode(1, contig, 1, 0)
            for (node in sites) {
                if (node.isEmpty() && node.pos!=1 && node.pos != tiglen) { nodePool.remove(node); continue}
                if (node.forward==1) { //left end
                    if (lastNode == node) continue
                    val newNode = nodePool.getNode(-1, contig, node.extpos(), node.flex)
                    lastNode.partner = newNode
                    newNode.partner = lastNode
                    lastNode = node
                } else { //right end
                    lastNode.partner = node
                    node.partner = lastNode
                    lastNode = if (node.pos == tiglen) node else nodePool.getNode(1, contig, node.extpos(), node.flex)
                }
            }
            if (lastNode.pos!=tiglen) {
                val newNode = nodePool.getNode(-1, contig, tiglen, 0)
                lastNode.partner = newNode
                newNode.partner = lastNode
            }
        }
        println("[Info] ${contigGroup.size} contigs are split into ${nodePool.size/2} segments.")
    }
    fun correctAndScaffoldFasta() {
        if (!blastTmpDir.exists()) blastTmpDir.mkdirs()
        loadEdge(reportFile)
        rescueDebris()
        placement()
        mergeAndWrite()
        outFasta.close()
        lossSeq.close()
        scaffoldInfo.close()
    }
    private fun loadEdge(reportFile:File) {
        val spaceChars = """\s+""".toRegex()
        val links = reportFile.readLines().filter { it.startsWith("+") || it.startsWith("-") }
        val pb = ProgressBar("BLAST overlap",links.size.toLong())
        links.forEach { line ->
            val tokens = line.split(spaceChars)
            val ref = nodePool.getBySegment(tokens[0])
            val qry = nodePool.getBySegment(tokens[1])
            val s = tokens[2].toDouble()
            val overlap = anchorByBLAST(ref,qry)
            ref.add(Edge(overlap[1],overlap[0],s,qry))
            qry.add(Edge(overlap[2],overlap[3],s,ref))
            pb.step()
        }
        pb.close()
        println("[Info] ${links.size} segment pairs loaded")
    }
    private fun loadFasta(draftFile:File): MutableMap<String, String> {
        val contigs = mutableMapOf<String,String>()
        val fasta = FastaSequenceFile(draftFile,false)
        var contig: ReferenceSequence? = fasta.nextSequence()
        while (contig!=null){
            contigs[contig.name] = contig.baseString.toUpperCase()
            contig = fasta.nextSequence()
        }
        return contigs
    }
}
