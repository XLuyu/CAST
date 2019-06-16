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
class Contig {
    var name = StringBuilder()
    val buffer = StringBuilder()
    fun append(id:String, seq:String){
        name.append(id)
        buffer.append(seq)
    }
}
class Stringer(private val draftFile:File, private val reportFile:File, outDir:File) {
    private val javaRuntime = Runtime.getRuntime()
    private val IDpattern = Regex("""([^(]+)\((\d+),(\d+)\)\((\d+),(\d+)\)""")
    private val blastTmpDir = File(outDir, "blast_tmp")
    private val tempRef = File(blastTmpDir, "ref.fasta")
    private val tempQry = File(blastTmpDir, "qry.fasta")
    private val lossSeq by lazy{ java.io.PrintWriter(File(blastTmpDir, "LossSeq.fasta"))}
    private val scaffoldInfo by lazy{ java.io.PrintWriter(File(blastTmpDir, "scaffold_info.log"))}
    private val outFasta = FastaWriter(File(outDir, "CAST.fasta"))
    private val contigs = mutableMapOf<String,String>() //Draft Genome
    private var links = mutableListOf<Link>()  //Edge
    private val segmentPool = mutableMapOf<String,Segment>() //global dict of Segment for no duplicate and identity
    fun MutableMap<String,Segment>.getByName(name:String): OrientedSegment {
        val segment = this.getOrPut(name.substring(1), { Segment(name.substring(1)) })
        return OrientedSegment(segment, name[0]=='+')
    }
    inner class Segment(identifier:String) {
        private val groups = IDpattern.find(identifier)!!.groups.map { it!!.value }
        var contig = groups[1]
        var leftEx = groups[2].toInt()
        var left = groups[3].toInt()
        var right = groups[4].toInt()
        var rightEx = groups[5].toInt()
        fun seq() = contigs[contig]!!.substring(leftEx-1,rightEx)
        fun len() = rightEx-leftEx+1
        var joint = arrayOf<Pair<Link,Int>?>(null,null)
        operator fun get(idx:Int) = joint[idx]
        init{ println(identifier)}
        fun updateJoint(idx:Int, linkJoint:Pair<Link,Int>?) {
            joint[idx] = linkJoint
            linkJoint?.let { (link, end) ->
                link.joint[end] = OrientedSegment(this,idx!=end)
            }
        }
    }
    inner class OrientedSegment(val segment:Segment, private val forward:Boolean){
        fun len() = segment.len()
        fun seq() = if (forward) segment.seq() else segment.seq().reversed() .map { Util.complementary.getOrDefault(it,'N') } .joinToString(separator="")
        fun leftLossSeq() = if (forward) contigs[segment.contig]!!.substring(segment.left-1,segment.leftEx-1) else contigs[segment.contig]!!.substring(segment.rightEx,segment.right)
        fun rightLossSeq() = if (forward) contigs[segment.contig]!!.substring(segment.rightEx,segment.right) else contigs[segment.contig]!!.substring(segment.left-1,segment.leftEx-1)
        fun leftExLen() = if (forward) segment.left-segment.leftEx else segment.rightEx-segment.right
        fun rightExLen() = if (!forward) segment.left-segment.leftEx else segment.rightEx-segment.right
        fun pruneL(loss:Int) = if (!forward) segment.rightEx -= loss else segment.leftEx += loss
        fun pruneR(loss:Int) = if (forward) segment.rightEx -= loss else segment.leftEx += loss
        fun bestL() = if (forward) segment[0] else segment[1]
        fun bestR() = if (forward) segment[1] else segment[0]
        fun setBestL(v:Pair<Link,Int>?) = if (forward) segment.joint[0]=v else segment.joint[1]=v
        fun setBestR(v:Pair<Link,Int>?) = if (forward) segment.joint[1]=v else segment.joint[0]=v
        fun getName(): String {
            val l = if (segment.leftEx==1) "S" else segment.leftEx.toString()
            val r = if (segment.rightEx==contigs[segment.contig]!!.length) "E" else segment.rightEx.toString()
            return if (forward) "${segment.contig}(${l}_$r)" else "${segment.contig}(${r}_$l)"
        }
        fun reverse() = OrientedSegment(segment, !forward)
    }
    inner class Link(A:String, B:String, val s:Double){
        val joint = arrayOf(segmentPool.getByName(A).reverse(), segmentPool.getByName(B))
        operator fun get(idx:Int) = joint[idx]
        fun asOnRight(idx:Int) = if (idx==0) joint[0].reverse() else joint[1]
        var overlap:List<Int>? = null
    }

    private fun BLASTSeqPair(ref:String, qry:String): List<List<Int>> {
        tempRef.writeText(">reference\n$ref\n")
        tempQry.writeText(">query\n$qry\n")
        javaRuntime.exec("makeblastdb -in ${tempRef.absolutePath} -dbtype nucl").waitFor()
        val cmd = javaRuntime.exec(arrayOf("blastn", "-num_threads", Runtime.getRuntime().availableProcessors().toString(), "-db", tempRef.absolutePath, "-query", tempQry.absolutePath, "-outfmt", "6 sstart send qstart qend"))
        return cmd.inputStream.bufferedReader().readLines().map{ line -> line.split('\t').map{it.toInt() } }
        // ref XXXXXXXXXXX
        //            XXXXXXXXX query
    }
    private fun anchorByBLAST(ref:OrientedSegment, qry:OrientedSegment): List<Int>? {
        val hsps = BLASTSeqPair(ref.seq(),qry.seq())
        val valid = hsps.filter {hsp -> hsp[0]<hsp[1] && hsp[2]<hsp[3]}
        if (valid.isEmpty()) return null // no valid alignment
        var best = valid.minBy{ Math.max(ref.len()-ref.rightExLen()-it[1],0) + Math.max(it[2]-1-qry.leftExLen(),0)}!!
        val loss = Math.max(ref.len()-ref.rightExLen()-best[1],0) + Math.max(best[2]-1-qry.leftExLen(),0)
        scaffoldInfo.write("[BLAST] ${ref.getName()} ${qry.getName()} $loss" + best.joinToString(",", " (",") "))
        scaffoldInfo.write(if (loss<=Math.min(ref.len(),qry.len())*0.05) "PASS\n" else " FAIL\n")
        best = listOf(ref.len()-best[0]+1, ref.len()-best[1]+1, best[2], best[3])
        return if (loss<=Math.min(ref.len(),qry.len())*0.05) best else null
        // best=(dist[sStart,length], dist[sEnd,length], qStart, qEnd)
    }
    private fun linkFilter(): List<Segment> {
        val pb = ProgressBar("BLAST overlap",links.size.toLong())
        for ( link in links) {
            val ref = link[0]
            val qry = link[1]
            link.overlap = anchorByBLAST(ref,qry)
            pb.step()
            if (link.overlap!=null) {
                if (ref.bestR()==null || ref.bestR()!!.first.s<link.s) ref.setBestR(Pair(link,0))
                if (qry.bestL()==null || qry.bestL()!!.first.s<link.s) qry.setBestL(Pair(link,1))
            }
        }
        pb.close()
        println("[Info] ${links.count{it.overlap!=null}} links found overlap by BLAST")
        for ( link in links) if (link.overlap!=null) {
            val ref = link[0]
            val qry = link[1]
            scaffoldInfo.write("[Mutual] ${ref.getName()} ${qry.getName()} ref:${ref.bestR()==Pair(link,0)} qry:${qry.bestL()==Pair(link,1)}\n")
            if (ref.bestR()==Pair(link,0) && qry.bestL()!=Pair(link,1)) ref.setBestR(null)
            if (ref.bestR()!=Pair(link,0) && qry.bestL()==Pair(link,1)) qry.setBestL(null)
            if (ref.bestR()!=Pair(link,0) || qry.bestL()!=Pair(link,1)) link.overlap = null
        }
        println("[Info] ${links.count{it.overlap!=null}} links agrees on both ends")
        val segmentsByContig = segmentPool.values.groupBy{it.contig}.flatMap { (contig, segmentList) ->
            val segments = segmentList.sortedBy { it.leftEx }
            var last = Pair(1,1)
            var lastJoint:Pair<Link,Int>? = null
            val result = mutableListOf<Segment>()
            for (segment in segments) {
                if (segment[0]!=null) {
                    if (last != Pair(segment.leftEx, segment.left)) {
                        result.add(Segment("$contig(${last.first},${last.second})(${segment.leftEx},${segment.left})"))
                        result.last().updateJoint(0,lastJoint)
                        result.last().joint[1] = null
                    }
                    last = Pair(segment.leftEx, segment.left)
                    lastJoint = segment[0]
                }
                if (segment[1]!=null){
                    result.add(Segment("$contig(${last.first},${last.second})(${segment.right},${segment.rightEx})"))
                    result.last().updateJoint(0,lastJoint)
                    result.last().updateJoint(1,segment[1])
                    last = Pair(segment.right,segment.rightEx)
                    lastJoint = null
                }
            }
            val tiglen = contigs[contig]!!.length
            if (last.first!=tiglen) {
                result.add(Segment("$contig(${last.first},${last.second})($tiglen,$tiglen)"))
                result.last().updateJoint(0, lastJoint)
                result.last().joint[1] = null
            }
            result
        }
        println("[Info] ${segmentsByContig.map{ it.contig}.toSet().size} contigs are split into ${segmentsByContig.size} segments.")
        return segmentsByContig
    }
    private fun mergeAndWrite(segments: Iterable<Segment>) {
        val visit = mutableSetOf<Segment>()
        fun traverse(os:OrientedSegment, contig:Contig){
            if (visit.contains(os.segment)) return
            visit.add(os.segment)
            if (os.bestR()!=null){
                val (link,idx) = os.bestR()!!
                os.pruneR(link.overlap!![1+idx]-1)
//        link(idx^1,1).pruneL(link.overlap(3-idx*3))
                val theOtherEnd = link.asOnRight(idx xor 1)
                val seq = theOtherEnd.seq()
                theOtherEnd.pruneL(link.overlap!![2-idx])
                val rightloss = Math.max(0,-os.rightExLen())
                val leftloss = Math.max(0,-theOtherEnd.leftExLen())
                if (rightloss>20) lossSeq.write(">${os.getName()}_${theOtherEnd.getName()}_RightLoss$rightloss\n${os.rightLossSeq()}\n")
                if (leftloss>20) lossSeq.write(">${os.getName()}_${theOtherEnd.getName()}_LeftLoss$leftloss\n${theOtherEnd.leftLossSeq()}\n")
                lossSeq.write(">${os.getName()}_${theOtherEnd.getName()}_Overlap\n${seq.substring(link.overlap!![2-idx], link.overlap!![3-idx*3])}\n")
                theOtherEnd.pruneL(link.overlap!![3-idx*3]-link.overlap!![2-idx])
                contig.append(os.getName(),os.seq())
                traverse(theOtherEnd, contig)
            } else contig.append(os.getName(),os.seq())
        }
        for ( segment in segments) if (segment.joint.contains(null) && !visit.contains(segment)) { // for chain
            val contig = Contig()
            traverse(OrientedSegment(segment,segment[0]==null), contig)
            outFasta.writeRecord(contig.name.toString(), contig.buffer.toString())
        }
        for ( segment in segments) if (!visit.contains(segment)) { // for circle
            val contig = Contig()
            traverse(OrientedSegment(segment,true), contig)
            outFasta.writeRecord(contig.name.toString(), contig.buffer.toString())
        }
        segments.forEach{s -> contigs.remove(s.contig) }
        contigs.forEach{(k,v) -> outFasta.writeRecord(k,v) }
    }

    fun correctAndScaffoldFasta() {
        if (!blastTmpDir.exists()) blastTmpDir.mkdirs()
        loadFastaAndEdge()
        val filteredSegments = linkFilter()
        mergeAndWrite(filteredSegments)
        outFasta.close()
        lossSeq.close()
        scaffoldInfo.close()
    }
    private fun loadFastaAndEdge() {
        val fasta = FastaSequenceFile(draftFile,false)
        var contig: ReferenceSequence? = fasta.nextSequence()
        while (contig!=null){
            contigs[contig.name] = contig.baseString.toUpperCase()
            contig = fasta.nextSequence()
        }
        reportFile.forEachLine { line ->
            if (line.startsWith("+") || line.startsWith("-")){
                val tokens = line.split("""\s+""".toRegex())
                links.add(Link(tokens[0],tokens[1],tokens[2].toDouble()))
            }
        }
    }
}
