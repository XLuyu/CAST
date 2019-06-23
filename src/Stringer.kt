package edu.nus

import me.tongfei.progressbar.ProgressBar
import htsjdk.samtools.reference.*
import java.io.File

class FastaWriter(file: File){
    private val writer = file.bufferedWriter()
    fun writeRecord(name:String, seq:String) {
        writer.write(">$name\n")
        for ( line in seq.chunked(80)) writer.write(line+"\n")
    }
    fun close() { writer.close()}
}
class ContigBuffer {
    var name = StringBuilder()
    val buffer = StringBuilder()
    fun append(id:String, seq:String){
        name.append(id)
        buffer.append(seq)
    }
}

//class SiteNode(val contigName: String, val contigSeq: String, var leftEx: Int, var left: Int, var right: Int, var rightEx: Int) {
//    fun seq() = contigSeq.substring(leftEx-1,rightEx)
//    fun len() = rightEx-leftEx+1
//    var joint = arrayOf<Pair<Link,Int>?>(null,null)
//    operator fun get(idx:Int) = joint[idx]
//    fun updateJoint(idx:Int, linkJoint:Pair<Link,Int>?) {
//        joint[idx] = linkJoint
//        linkJoint?.let { (link, end) ->
//            link.joint[end] = OrientedSegment(this,idx!=end)
//        }
//    }
//}
//class OrientedSegment(val siteNode:SiteNode, private val forward:Boolean){
//    fun len() = siteNode.len()
//    fun seq() = if (forward) siteNode.seq() else siteNode.seq().reversed() .map { Util.complementary.getOrDefault(it,'N') } .joinToString(separator="")
//    fun leftLossSeq() = if (forward) siteNode.contigSeq.substring(siteNode.left-1,siteNode.leftEx-1) else siteNode.contigSeq.substring(siteNode.rightEx,siteNode.right)
//    fun rightLossSeq() = if (forward) siteNode.contigSeq.substring(siteNode.rightEx,siteNode.right) else siteNode.contigSeq.substring(siteNode.left-1,siteNode.leftEx-1)
//    fun leftExLen() = if (forward) siteNode.left-siteNode.leftEx else siteNode.rightEx-siteNode.right
//    fun rightExLen() = if (!forward) siteNode.left-siteNode.leftEx else siteNode.rightEx-siteNode.right
//    fun pruneL(loss:Int) = if (!forward) siteNode.rightEx -= loss else siteNode.leftEx += loss
//    fun pruneR(loss:Int) = if (forward) siteNode.rightEx -= loss else siteNode.leftEx += loss
//    fun bestL() = if (forward) siteNode[0] else siteNode[1]
//    fun bestR() = if (forward) siteNode[1] else siteNode[0]
//    fun setBestL(v:Pair<Link,Int>?) = if (forward) siteNode.joint[0]=v else siteNode.joint[1]=v
//    fun setBestR(v:Pair<Link,Int>?) = if (forward) siteNode.joint[1]=v else siteNode.joint[0]=v
//    fun getName(): String {
//        val l = if (siteNode.leftEx==1) "S" else siteNode.leftEx.toString()
//        val r = if (siteNode.rightEx==siteNode.contigSeq.length) "E" else siteNode.rightEx.toString()
//        return if (forward) "${siteNode.contigName}(${l}_$r)" else "${siteNode.contigName}(${r}_$l)"
//    }
//    fun reverse() = OrientedSegment(siteNode, !forward)
//}
//class Link(A: OrientedSegment, B: OrientedSegment, val s:Double){
//    val joint = arrayOf(A,B)
//    operator fun get(idx:Int) = joint[idx]
//    fun asOnRight(idx:Int) = if (idx==0) joint[0].reverse() else joint[1]
//    var overlap:List<Int>? = null
//}
//class SegmentPool(private val contigs: MutableMap<String, String>) :HashMap<String,SiteNode>(){
//    private val IDpattern = Regex("""([^(]+)\((\d+),(\d+)\)\((\d+),(\d+)\)""")
//    fun getByName(name:String): OrientedSegment {
//        val siteNode = this.getOrPut(name.substring(1), {
//            val groups = IDpattern.find(name.substring(1))!!.groups.map { it!!.value }
//            SiteNode(groups[1], contigs[groups[1]]!!, groups[2].toInt(), groups[3].toInt(), groups[4].toInt(), groups[5].toInt())
//        })
//        return OrientedSegment(siteNode, name[0]=='+')
//    }
//    fun createWithoutPut(contigName: String, leftEx: Int, left: Int, right: Int, rightEx: Int) =
//        SiteNode(contigName, contigs[contigName]!!, leftEx, left, right, rightEx)
//
//}
//class Stringer(draftFile:File, reportFile:File, outDir:File) {
//    private val javaRuntime = Runtime.getRuntime()
//    private val blastTmpDir = File(outDir, "blast_tmp")
//    private val tempRef = File(blastTmpDir, "ref.fasta")
//    private val tempQry = File(blastTmpDir, "qry.fasta")
//    private val lossSeq by lazy{ java.io.PrintWriter(File(blastTmpDir, "LossSeq.fasta"))}
//    private val scaffoldInfo by lazy{ java.io.PrintWriter(File(blastTmpDir, "scaffold_info.log"))}
//    private val outFasta = FastaWriter(File(outDir, "CAST.fasta"))
//    private val contigs = loadFasta(draftFile) //Draft Genome
//    private val segmentPool = SegmentPool(contigs) //global dict of SiteNode for no duplicate and identity
//    private var links = loadEdge(reportFile)  //Edge
//
//    private fun BLASTSeqPair(ref:String, qry:String): List<List<Int>> {
//        tempRef.writeText(">reference\n$ref\n")
//        tempQry.writeText(">query\n$qry\n")
//        javaRuntime.exec("makeblastdb -in ${tempRef.absolutePath} -dbtype nucl").waitFor()
//        val cmd = javaRuntime.exec(arrayOf("blastn", "-num_threads", Runtime.getRuntime().availableProcessors().toString(), "-db", tempRef.absolutePath, "-query", tempQry.absolutePath, "-outfmt", "6 sstart send qstart qend"))
//        return cmd.inputStream.bufferedReader().readLines().map{ line -> line.split('\t').map{it.toInt() } }
//        // ref XXXXXXXXXXX
//        //            XXXXXXXXX query
//    }
//    private fun anchorByBLAST(ref:OrientedSegment, qry:OrientedSegment): List<Int>? {
//        val hsps = BLASTSeqPair(ref.seq(),qry.seq())
//        val valid = hsps.filter {hsp -> hsp[0]<hsp[1] && hsp[2]<hsp[3]}
//        if (valid.isEmpty()) return null // no valid alignment
//        var best = valid.minBy{ Math.max(ref.len()-ref.rightExLen()-it[1],0) + Math.max(it[2]-1-qry.leftExLen(),0)}!!
//        val loss = Math.max(ref.len()-ref.rightExLen()-best[1],0) + Math.max(best[2]-1-qry.leftExLen(),0)
//        scaffoldInfo.write("[BLAST] ${ref.getName()} ${qry.getName()} $loss" + best.joinToString(",", " (",") "))
//        scaffoldInfo.write(if (loss<=Math.min(ref.len(),qry.len())*0.05) "PASS\n" else " FAIL\n")
//        best = listOf(ref.len()-best[0]+1, ref.len()-best[1]+1, best[2], best[3])
//        return if (loss<=Math.min(ref.len(),qry.len())*0.05) best else null
//        // best=(dist[sStart,length], dist[sEnd,length], qStart, qEnd)
//    }
//    private fun linkFilter(): List<SiteNode> {
//        val pb = ProgressBar("BLAST overlap",links.size.toLong())
//        for ( link in links) {
//            val ref = link[0]
//            val qry = link[1]
//            link.overlap = anchorByBLAST(ref,qry)
//            pb.step()
//            if (link.overlap!=null) {
//                if (ref.bestR()==null || ref.bestR()!!.first.s<link.s) ref.setBestR(Pair(link,0))
//                if (qry.bestL()==null || qry.bestL()!!.first.s<link.s) qry.setBestL(Pair(link,1))
//            }
//        }
//        pb.close()
//        println("[Info] ${links.count{it.overlap!=null}} links found overlap by BLAST")
//        for ( link in links) if (link.overlap!=null) {
//            val ref = link[0]
//            val qry = link[1]
//            scaffoldInfo.write("[Mutual] ${ref.getName()} ${qry.getName()} ref:${ref.bestR()==Pair(link,0)} qry:${qry.bestL()==Pair(link,1)}\n")
//            if (ref.bestR()==Pair(link,0) && qry.bestL()!=Pair(link,1)) ref.setBestR(null)
//            if (ref.bestR()!=Pair(link,0) && qry.bestL()==Pair(link,1)) qry.setBestL(null)
//            if (ref.bestR()!=Pair(link,0) || qry.bestL()!=Pair(link,1)) link.overlap = null
//        }
//        println("[Info] ${links.count{it.overlap!=null}} links agrees on both ends")
//        val segmentsByContig = segmentPool.values.groupBy{it.contigName}.flatMap { (contig, segmentList) ->
//            val segments = segmentList.sortedBy { it.leftEx }
//            var last = Pair(1,1)
//            var lastJoint:Pair<Link,Int>? = null
//            val result = mutableListOf<SiteNode>()
//            for (siteNode in segments) {
//                if (siteNode[0]!=null) {
//                    if (last != Pair(siteNode.leftEx, siteNode.left)) {
//                        result.add(segmentPool.createWithoutPut(contig, last.first, last.second, siteNode.leftEx, siteNode.left))
//                        result.last().updateJoint(0,lastJoint)
//                        result.last().joint[1] = null
//                    }
//                    last = Pair(siteNode.leftEx, siteNode.left)
//                    lastJoint = siteNode[0]
//                }
//                if (siteNode[1]!=null){
//                    result.add(segmentPool.createWithoutPut(contig, last.first, last.second, siteNode.right, siteNode.rightEx))
//                    result.last().updateJoint(0,lastJoint)
//                    result.last().updateJoint(1,siteNode[1])
//                    last = Pair(siteNode.right,siteNode.rightEx)
//                    lastJoint = null
//                }
//            }
//            val tiglen = contigs[contig]!!.length
//            if (last.first!=tiglen) {
//                result.add(segmentPool.createWithoutPut(contig, last.first, last.second, tiglen, tiglen))
//                result.last().updateJoint(0, lastJoint)
//                result.last().joint[1] = null
//            }
//            result
//        }
//        println("[Info] ${segmentsByContig.map{ it.contigName}.toSet().size} contigs are split into ${segmentsByContig.size} segments.")
//        return segmentsByContig
//    }
//    private fun mergeAndWrite(segments: Iterable<SiteNode>) {
//        val visit = mutableSetOf<SiteNode>()
//        fun traverse(os:OrientedSegment, contig:ContigBuffer){
//            if (visit.contains(os.siteNode)) return
//            visit.add(os.siteNode)
//            if (os.bestR()!=null){
//                val (link,idx) = os.bestR()!!
//                os.pruneR(link.overlap!![1+idx]-1)
////        link(idx^1,1).pruneL(link.overlap(3-idx*3))
//                val theOtherEnd = link.asOnRight(idx xor 1)
//                val seq = theOtherEnd.seq()
//                theOtherEnd.pruneL(link.overlap!![2-idx])
//                val rightloss = Math.max(0,-os.rightExLen())
//                val leftloss = Math.max(0,-theOtherEnd.leftExLen())
//                if (rightloss>20) lossSeq.write(">${os.getName()}_${theOtherEnd.getName()}_RightLoss$rightloss\n${os.rightLossSeq()}\n")
//                if (leftloss>20) lossSeq.write(">${os.getName()}_${theOtherEnd.getName()}_LeftLoss$leftloss\n${theOtherEnd.leftLossSeq()}\n")
//                lossSeq.write(">${os.getName()}_${theOtherEnd.getName()}_Overlap\n${seq.substring(link.overlap!![2-idx], link.overlap!![3-idx*3])}\n")
//                theOtherEnd.pruneL(link.overlap!![3-idx*3]-link.overlap!![2-idx])
//                contig.append(os.getName(),os.seq())
//                traverse(theOtherEnd, contig)
//            } else contig.append(os.getName(),os.seq())
//        }
//        for ( siteNode in segments) if (siteNode.joint.contains(null) && !visit.contains(siteNode)) { // for chain
//            val contig = ContigBuffer()
//            traverse(OrientedSegment(siteNode,siteNode[0]==null), contig)
//            outFasta.writeRecord(contig.name.toString(), contig.buffer.toString())
//        }
//        for ( siteNode in segments) if (!visit.contains(siteNode)) { // for circle
//            val contig = ContigBuffer()
//            traverse(OrientedSegment(siteNode,true), contig)
//            outFasta.writeRecord(contig.name.toString(), contig.buffer.toString())
//        }
//        segments.forEach{s -> contigs.remove(s.contigName) }
//        contigs.forEach{(k,v) -> outFasta.writeRecord(k,v) }
//    }
//
//    fun correctAndScaffoldFasta() {
//        if (!blastTmpDir.exists()) blastTmpDir.mkdirs()
//        val filteredSegments = linkFilter()
//        mergeAndWrite(filteredSegments)
//        outFasta.close()
//        lossSeq.close()
//        scaffoldInfo.close()
//    }
//    private fun loadEdge(reportFile:File): MutableList<Link> {
//        val spaceChars = """\s+""".toRegex()
//        val link = mutableListOf<Link>()
//        reportFile.forEachLine { line ->
//            if (line.startsWith("+") || line.startsWith("-")){
//                val tokens = line.split(spaceChars)
//                link.add(Link(segmentPool.getByName(tokens[0]).reverse(), segmentPool.getByName(tokens[1]),tokens[2].toDouble()))
//            }
//        }
//        return link
//    }
//    private fun loadFasta(draftFile:File): MutableMap<String, String> {
//        val contigs = mutableMapOf<String,String>()
//        val fasta = FastaSequenceFile(draftFile,false)
//        var contig: ReferenceSequence? = fasta.nextSequence()
//        while (contig!=null){
//            contigs[contig.name] = contig.baseString.toUpperCase()
//            contig = fasta.nextSequence()
//        }
//        return contigs
//    }
//}
