package edu.nus

import me.tongfei.progressbar.ProgressBar
import htsjdk.samtools.reference.*
import java.io.File

class Scaffolder(private val draftFile:File, private val reportFile:File, outDir:File) {
    private val javaRuntime = Runtime.getRuntime()
    private val IDpattern = Regex("""([^(]+)\((\d+),(\d+)\)\((\d+),(\d+)\)""")
    private val blastTmpDir = File(outDir, "blast_tmp")
    private val tempRef = File(blastTmpDir, "ref.fasta")
    private val tempQry = File(blastTmpDir, "qry.fasta")
    private val lossSeq by lazy{ java.io.PrintWriter(File(blastTmpDir, "LossSeq.fasta"))}
    private val scaffoldInfo by lazy{ java.io.PrintWriter(File(blastTmpDir, "scaffold_info.log"))}
    private val outFasta = FastaWriter(File(outDir, "CAST.fasta"))
    private val contigs = loadFasta(draftFile) //Draft Genome
    private val segmentPool = SegmentPool(contigs) //global dict of Segment for no duplicate and identity
    private var links = loadEdge(reportFile)  //Edge

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

    private fun mergeAndWrite(segments: Iterable<Segment>) {
        val visit = mutableSetOf<Segment>()
        fun traverse(os:OrientedSegment, contig:ContigBuffer){
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
            val contig = ContigBuffer()
            traverse(OrientedSegment(segment,segment[0]==null), contig)
            outFasta.writeRecord(contig.name.toString(), contig.buffer.toString())
        }
        for ( segment in segments) if (!visit.contains(segment)) { // for circle
            val contig = ContigBuffer()
            traverse(OrientedSegment(segment,true), contig)
            outFasta.writeRecord(contig.name.toString(), contig.buffer.toString())
        }
        segments.forEach{s -> contigs.remove(s.contigName) }
        contigs.forEach{(k,v) -> outFasta.writeRecord(k,v) }
    }
    private fun placement(){

    }
    fun correctAndScaffoldFasta() {
        if (!blastTmpDir.exists()) blastTmpDir.mkdirs()
        outFasta.close()
        lossSeq.close()
        scaffoldInfo.close()
    }
    private fun loadEdge(reportFile:File): MutableList<Link> {
        val spaceChars = """\s+""".toRegex()
        val links = mutableListOf<Link>()
        reportFile.forEachLine { line ->
            if (line.startsWith("+") || line.startsWith("-")){
                val tokens = line.split(spaceChars)
                links.add(Link(segmentPool.getByName(tokens[0]).reverse(), segmentPool.getByName(tokens[1]),tokens[2].toDouble()))
            }
        }
        val pb = ProgressBar("BLAST overlap",links.size.toLong())
        for (link in links) {
            val ref = link[0]
            val qry = link[1]
            link.overlap = anchorByBLAST(ref,qry)
            pb.step()
        }
        pb.close()
        println("[Info] ${links.count{it.overlap!=null}} links found overlap by BLAST")
        return links
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
